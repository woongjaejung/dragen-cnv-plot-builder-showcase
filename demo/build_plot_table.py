#!/usr/bin/env python3
"""Join simplified DRAGEN-shaped folders into one long plotting table.

The input files are the synthetic stand-ins produced by
``generate_synthetic_run.py``.  They intentionally resemble the shape of
DRAGEN output without attempting to be byte-level replicas of those files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import time
from dataclasses import dataclass
from pathlib import Path


TABLE_COLUMNS = [
    "sample_id",
    "result_dir",
    "contig",
    "start",
    "end",
    "target_name",
    "annotation_target_name",
    "gene",
    "gene_strand",
    "exon_number",
    "tn_log2",
    "normalized_ratio",
    "improper_pairs",
    "segment_start",
    "segment_end",
    "segment_probes",
    "segment_mean",
    "call_type",
    "call_copy_number",
    "call_filter",
    "call_qual",
    "call_start",
    "call_end",
    "call_count",
    "call_pass_count",
]


@dataclass(frozen=True)
class Segment:
    contig: str
    start: int
    end: int
    probes: int
    mean: str


@dataclass(frozen=True)
class Call:
    call_id: str
    contig: str
    call_type: str
    start: int
    end: int
    copy_number: int
    filter_tag: str
    qual: str

    @property
    def is_pass(self) -> bool:
        return self.filter_tag == "PASS"


def read_delimited(path: Path, *, compressed: bool = False) -> list[dict[str, str]]:
    """Read a small TSV while accepting a leading comment-style header."""

    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    rows = [row for row in rows if row and any(cell.strip() for cell in row)]
    if not rows:
        return []
    header = rows[0]
    header[0] = header[0].lstrip("#")
    return [dict(zip(header, row)) for row in rows[1:]]


def find_result_dir(sample_dir: Path) -> Path:
    """Resolve a sample root or accept the result directory directly."""

    if list(sample_dir.glob("*.dragen.tn.tsv.gz")):
        return sample_dir
    sample_id = sample_dir.name
    candidate = sample_dir / "vcf" / "dragen" / sample_id / "cnv_exome_varcaller"
    if not candidate.is_dir():
        raise FileNotFoundError(f"No cnv_exome_varcaller directory under {sample_dir}")
    return candidate


def required_target_file(result_dir: Path, sample_id: str) -> Path:
    exact = result_dir / f"{sample_id}.dragen.tn.tsv.gz"
    if exact.exists():
        return exact
    matches = sorted(result_dir.glob("*.dragen.tn.tsv.gz"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one *.dragen.tn.tsv.gz in {result_dir}")
    return matches[0]


def optional_path(result_dir: Path, suffix: str) -> Path | None:
    matches = sorted(result_dir.glob(f"*.{suffix}"))
    return matches[0] if matches else None


def read_annotation(path: Path) -> dict[tuple[str, int, int], dict[str, str]]:
    annotations: dict[tuple[str, int, int], dict[str, str]] = {}
    for row in read_delimited(path):
        key = (row["chrom"], int(row["start"]), int(row["end"]))
        annotations[key] = row
    return annotations


def read_segments(path: Path | None) -> list[Segment]:
    if path is None:
        return []
    return [
        Segment(
            row["contig"],
            int(row["start"]),
            int(row["end"]),
            int(row["num_probes"]),
            row["segment_mean"],
        )
        for row in read_delimited(path)
    ]


def parse_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for field in value.split(";"):
        if "=" in field:
            key, attr_value = field.split("=", 1)
            attributes[key] = attr_value
    return attributes


def read_calls(path: Path | None) -> list[Call]:
    if path is None:
        return []
    calls: list[Call] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attributes = parse_attributes(fields[8])
            calls.append(
                Call(
                    attributes.get("ID", "synthetic-call"),
                    fields[0],
                    fields[2],
                    int(fields[3]),
                    int(fields[4]),
                    int(attributes["CN"]),
                    attributes["Filter"],
                    attributes["QUAL"],
                )
            )
    return calls


def overlaps(start: int, end: int, other_start: int, other_end: int) -> bool:
    return start <= other_end and end >= other_start


def choose_reported_call(calls: list[Call]) -> Call | None:
    """Select PASS first, then stable coordinates, regardless of input order."""

    if not calls:
        return None
    return sorted(
        calls,
        key=lambda call: (0 if call.is_pass else 1, call.start, call.end, call.call_id),
    )[0]


def find_segment(target: dict[str, str], segments: list[Segment]) -> Segment | None:
    candidates = [
        segment
        for segment in segments
        if segment.contig == target["contig"]
        and overlaps(int(target["start"]), int(target["end"]), segment.start, segment.end)
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda segment: (min(int(target["end"]), segment.end) - max(int(target["start"]), segment.start), -segment.start))


def row_for_target(
    sample_id: str,
    result_dir: Path,
    target: dict[str, str],
    annotation: dict[str, str],
    segment: Segment | None,
    target_calls: list[Call],
) -> dict[str, str]:
    reported = choose_reported_call(target_calls)
    tn_log2 = float(target["tn_log2"])
    row = {
        "sample_id": sample_id,
        "result_dir": str(result_dir),
        "contig": target["contig"],
        "start": target["start"],
        "end": target["end"],
        "target_name": target["target_name"],
        "annotation_target_name": annotation.get("name", ""),
        "gene": annotation.get("GeneSymbol", ""),
        "gene_strand": annotation.get("strand", ""),
        "exon_number": annotation.get("exonNumber", ""),
        "tn_log2": target["tn_log2"],
        "normalized_ratio": f"{2 ** tn_log2:.6f}",
        "improper_pairs": target["improper_pairs"],
        "segment_start": str(segment.start) if segment else "",
        "segment_end": str(segment.end) if segment else "",
        "segment_probes": str(segment.probes) if segment else "",
        "segment_mean": segment.mean if segment else "",
        "call_type": reported.call_type if reported else "",
        "call_copy_number": str(reported.copy_number) if reported else "",
        "call_filter": reported.filter_tag if reported else "",
        "call_qual": reported.qual if reported else "",
        "call_start": str(reported.start) if reported else "",
        "call_end": str(reported.end) if reported else "",
        "call_count": str(len(target_calls)),
        "call_pass_count": str(sum(call.is_pass for call in target_calls)),
    }
    return row


def build_rows(
    sample_dirs: list[Path],
    annotation_path: Path,
) -> tuple[list[dict[str, str]], int, int, float]:
    started = time.monotonic()
    annotations = read_annotation(annotation_path)
    annotation_cache: dict[tuple[str, int, int], dict[str, str]] = {}
    cache_hits = 0
    cache_misses = 0
    rows: list[dict[str, str]] = []

    for sample_dir in sample_dirs:
        result_dir = find_result_dir(sample_dir)
        sample_id = sample_dir.name
        target_rows = read_delimited(required_target_file(result_dir, sample_id), compressed=True)
        segments = read_segments(optional_path(result_dir, "dragen.seg"))
        calls = read_calls(optional_path(result_dir, "dragen.cnv.gff3"))
        for target in target_rows:
            key = (target["contig"], int(target["start"]), int(target["end"]))
            if key in annotation_cache:
                annotation = annotation_cache[key]
                cache_hits += 1
            else:
                annotation = annotations.get(key, {})
                annotation_cache[key] = annotation
                cache_misses += 1
            if not annotation:
                raise KeyError(f"Target {key} is absent from {annotation_path}")
            target_calls = [
                call
                for call in calls
                if call.contig == target["contig"]
                and overlaps(int(target["start"]), int(target["end"]), call.start, call.end)
            ]
            rows.append(
                row_for_target(
                    sample_id,
                    result_dir,
                    target,
                    annotation,
                    find_segment(target, segments),
                    target_calls,
                )
            )
    return rows, cache_hits, cache_misses, time.monotonic() - started


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_dir", nargs="+", type=Path, help="sample roots or result directories")
    parser.add_argument("--annotation", required=True, type=Path, help="synthetic panel annotation TSV")
    parser.add_argument("-o", "--output", required=True, type=Path, help="long CSV output path")
    parser.add_argument("--slim", action="store_true", help="omit the repeated result_dir column")
    parser.add_argument("--dry-run", action="store_true", help="validate and summarize without writing CSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows, cache_hits, cache_misses, elapsed = build_rows(args.sample_dir, args.annotation)
    fieldnames = [column for column in TABLE_COLUMNS if not (args.slim and column == "result_dir")]
    print(
        f"Annotation cache: {cache_hits} hits, {cache_misses} unique targets resolved; "
        f"elapsed {elapsed:.4f}s"
    )
    print(f"Rows: {len(rows)}; columns: {len(fieldnames)}")
    if args.dry_run:
        print(f"Dry run: would write {args.output}")
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row[key] for key in fieldnames} for row in rows)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
