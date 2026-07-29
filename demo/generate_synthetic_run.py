#!/usr/bin/env python3
"""Generate a small synthetic run shaped like DRAGEN CNV exome output.

The files written here are simplified stand-ins for DRAGEN's real outputs, not
byte-level replicas.  They contain no patient data and are deliberately small
enough to inspect by hand.  The seeded generator plants a handful of events so
the later table and plot steps have useful examples to show.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import random
from dataclasses import dataclass
from pathlib import Path


DEFAULT_SEED = 20260729
DEFAULT_OUT = Path("demo/_work/synthetic_run")


@dataclass(frozen=True)
class GeneSpec:
    name: str
    contig: str
    strand: str
    exon_count: int


@dataclass(frozen=True)
class Target:
    contig: str
    start: int
    end: int
    name: str
    gene: str
    exon_number: int
    strand: str


@dataclass(frozen=True)
class PlantedEvent:
    gene: str
    sample: str
    start_exon: int
    end_exon: int
    ratio: float
    call_type: str
    copy_number: int
    filter_tag: str
    qual: float


GENES = (
    GeneSpec("BRCA1", "chr17", "-", 23),
    GeneSpec("BRCA2", "chr13", "+", 27),
    GeneSpec("TP53", "chr17", "-", 11),
    GeneSpec("MLH1", "chr3", "+", 19),
    GeneSpec("MSH2", "chr2", "+", 16),
    GeneSpec("PTEN", "chr10", "+", 9),
    GeneSpec("EGFR", "chr7", "+", 28),
    GeneSpec("DMD", "chrX", "-", 20),
    GeneSpec("CFTR", "chr7", "+", 27),
    GeneSpec("STK11", "chr19", "+", 9),
    GeneSpec("PMS2", "chr7", "-", 15),
    GeneSpec("CDKN2A", "chr9", "-", 3),
)

GENE_BY_NAME = {gene.name: gene for gene in GENES}

SAMPLE_IDS = tuple(
    [f"DEMO-RUN01-S{i:02d}" for i in range(1, 13)]
    + ["DEMO-RUN01-NTC-A", "DEMO-RUN01-NTC-B"]
)

EVENTS = (
    PlantedEvent("BRCA1", "DEMO-RUN01-S03", 1, 23, 0.50, "DEL", 1, "PASS", 48.0),
    PlantedEvent("EGFR", "DEMO-RUN01-S07", 18, 21, 1.50, "DUP", 3, "PASS", 43.0),
    PlantedEvent("PTEN", "DEMO-RUN01-S05", 6, 6, 0.45, "DEL", 1, "cnvLength", 7.5),
    PlantedEvent("DMD", "DEMO-RUN01-S09", 8, 13, 0.50, "DEL", 1, "PASS", 39.0),
    PlantedEvent("TP53", "DEMO-RUN01-S02", 4, 6, 0.55, "DEL", 1, "cnvQual", 8.2),
)


def event_for(gene: str, sample: str) -> PlantedEvent | None:
    """Return the patient-specific event, if this gene/sample has one."""

    for event in EVENTS:
        if event.gene == gene and event.sample == sample:
            return event
    return None


def target_ratio(
    gene: GeneSpec,
    exon_number: int,
    sample: str,
    rng: random.Random,
) -> float:
    """Make one target ratio, adding event plateaus after baseline noise."""

    is_ntc = "NTC" in sample
    if is_ntc:
        # NTCs are low and relatively scattered compared with a diploid ratio.
        ratio = max(0.02, rng.gauss(0.10, 0.04))
        event = None
    else:
        ratio = rng.gauss(1.00, 0.06)
        event = event_for(gene.name, sample)

    if event and event.start_exon <= exon_number <= event.end_exon:
        ratio = max(0.02, rng.gauss(event.ratio, event.ratio * 0.025))

    # The CFTR artifact is shared across the batch.  For an NTC it is applied
    # multiplicatively so the control remains near zero while the same target
    # footprint is still present in every sample folder.
    if gene.name == "CFTR" and 10 <= exon_number <= 11:
        ratio = max(0.02, ratio * 0.65)

    return ratio


def make_targets(rng: random.Random) -> list[Target]:
    """Lay out 1-4 capture targets per exon with deterministic realistic gaps."""

    cursor_by_contig = {contig: 1_000_000 for contig in {g.contig for g in GENES}}
    targets: list[Target] = []
    for gene in GENES:
        cursor = cursor_by_contig[gene.contig] + rng.randint(2_000, 6_000)
        for genomic_exon_index in range(1, gene.exon_count + 1):
            exon_number = (
                genomic_exon_index
                if gene.strand == "+"
                else gene.exon_count - genomic_exon_index + 1
            )
            cursor += rng.randint(450, 1_500)
            target_count = rng.randint(1, 4)
            for within_exon in range(1, target_count + 1):
                length = rng.randint(120, 200)
                start = cursor
                end = start + length - 1
                targets.append(
                    Target(
                        gene.contig,
                        start,
                        end,
                        f"{gene.name}_exon_{exon_number:02d}_target_{within_exon}",
                        gene.name,
                        exon_number,
                        gene.strand,
                    )
                )
                cursor = end + rng.randint(90, 260)
            cursor += rng.randint(500, 1_400)
        cursor_by_contig[gene.contig] = cursor + rng.randint(2_000, 6_000)
    return targets


def targets_by_gene(targets: list[Target]) -> dict[str, list[Target]]:
    grouped: dict[str, list[Target]] = {gene.name: [] for gene in GENES}
    for target in targets:
        grouped[target.gene].append(target)
    return grouped


def affected_spans(gene: GeneSpec, sample: str) -> list[tuple[int, int]]:
    """Return event spans used to split segment files into readable plateaus."""

    spans: list[tuple[int, int]] = []
    event = event_for(gene.name, sample)
    if event and "NTC" not in sample:
        spans.append((event.start_exon, event.end_exon))
    if gene.name == "CFTR":
        spans.append((10, 11))
    return sorted(spans)


def make_segments(
    gene: GeneSpec,
    targets: list[Target],
    sample: str,
    ratios: dict[str, float],
) -> list[dict[str, str]]:
    """Create one segment per contiguous normal or planted exon span."""

    spans = affected_spans(gene, sample)
    boundaries = {1, gene.exon_count + 1}
    for start, end in spans:
        boundaries.add(start)
        boundaries.add(end + 1)
    sorted_boundaries = sorted(boundaries)
    segments: list[dict[str, str]] = []

    for left, right_exclusive in zip(sorted_boundaries, sorted_boundaries[1:]):
        members = [
            target
            for target in targets
            if left <= target.exon_number < right_exclusive
        ]
        if not members:
            continue
        mean = sum(ratios[target.name] for target in members) / len(members)
        segments.append(
            {
                "contig": gene.contig,
                "start": str(members[0].start),
                "end": str(members[-1].end),
                "num_probes": str(len(members)),
                "segment_mean": f"{mean:.4f}",
            }
        )
    return segments


def write_panel(path: Path, targets: list[Target]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "name", "GeneSymbol", "exonNumber", "strand"])
        for target in targets:
            writer.writerow(
                [
                    target.contig,
                    target.start,
                    target.end,
                    target.name,
                    target.gene,
                    target.exon_number,
                    target.strand,
                ]
            )


def write_target_file(path: Path, targets: list[Target], ratios: dict[str, float], rng: random.Random) -> None:
    text = io.StringIO()
    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
    writer.writerow(["#contig", "start", "end", "target_name", "tn_log2", "improper_pairs"])
    for target in targets:
        ratio = ratios[target.name]
        writer.writerow(
            [
                target.contig,
                target.start,
                target.end,
                target.name,
                f"{math.log2(ratio):.6f}",
                rng.randint(0, 3),
            ]
        )
    path.write_bytes(gzip.compress(text.getvalue().encode("utf-8"), mtime=0))


def write_segment_file(path: Path, segments: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "start", "end", "num_probes", "segment_mean"])
        for segment in segments:
            writer.writerow([segment[key] for key in ("contig", "start", "end", "num_probes", "segment_mean")])


def call_coordinates(gene: GeneSpec, targets: list[Target], event: PlantedEvent) -> tuple[int, int]:
    members = [
        target
        for target in targets
        if event.start_exon <= target.exon_number <= event.end_exon
    ]
    return members[0].start, members[-1].end


def write_call_file(path: Path, targets_by_gene_map: dict[str, list[Target]], sample: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")
        handle.write("# Synthetic calls are illustrative and not clinical results.\n")
        for event in EVENTS:
            if event.sample != sample:
                continue
            gene = GENE_BY_NAME[event.gene]
            start, end = call_coordinates(gene, targets_by_gene_map[event.gene], event)
            attributes = (
                f"ID={sample}_{gene.name}_{event.call_type};"
                f"CN={event.copy_number};Filter={event.filter_tag};QUAL={event.qual:.1f}"
            )
            row = [gene.contig, "synthetic", event.call_type, start, end, ".", gene.strand, ".", attributes]
            handle.write("\t".join(str(value) for value in row) + "\n")


def generate_run(out_dir: Path, seed: int) -> tuple[int, int]:
    """Write all synthetic inputs and return the sample and target counts."""

    rng = random.Random(seed)
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = out_dir / "synthetic_panel_targets.tsv"
    targets = make_targets(rng)
    write_panel(panel_path, targets)
    by_gene = targets_by_gene(targets)

    for sample in SAMPLE_IDS:
        result_dir = out_dir / sample / "vcf" / "dragen" / sample / "cnv_exome_varcaller"
        result_dir.mkdir(parents=True, exist_ok=True)
        ratios: dict[str, float] = {}
        segments: list[dict[str, str]] = []
        for gene in GENES:
            gene_targets = by_gene[gene.name]
            gene_ratios = {
                target.name: target_ratio(gene, target.exon_number, sample, rng)
                for target in gene_targets
            }
            ratios.update(gene_ratios)
            segments.extend(make_segments(gene, gene_targets, sample, gene_ratios))
        write_target_file(result_dir / f"{sample}.dragen.tn.tsv.gz", targets, ratios, rng)
        write_segment_file(result_dir / f"{sample}.dragen.seg", segments)
        write_call_file(result_dir / f"{sample}.dragen.cnv.gff3", by_gene, sample)

    return len(SAMPLE_IDS), len(targets)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="synthetic run directory")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="seed for deterministic random data")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_count, target_count = generate_run(args.out, args.seed)
    print(f"Generated {sample_count} synthetic samples and {target_count} capture targets in {args.out}")
    print(f"Seed: {args.seed}")


if __name__ == "__main__":
    main()
