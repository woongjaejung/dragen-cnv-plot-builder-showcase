#!/usr/bin/env python3
"""Turn the long synthetic table into a lazy, interactive Plotly page.

This HTML builder consumes the reduced table emitted by
``build_plot_table.py``.  The JSON is split into one script element per gene
so a browser parses only the group a reviewer opens; it is not a byte-level
replica of any private production format.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
from collections import OrderedDict
from pathlib import Path
from typing import Any


V8_STRING_LIMIT = 536_870_888
DEFAULT_PLOTLY_SRC = "../vendor/plotly-2.35.2.min.js"
DEFAULT_PLOTLY_PATH = Path(__file__).resolve().parents[1] / "docs/vendor/plotly-2.35.2.min.js"


def as_int(value: str) -> int:
    return int(value) if value else 0


def as_float(value: str) -> float | None:
    return float(value) if value else None


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def overlaps(start: int, end: int, other_start: int, other_end: int) -> bool:
    return start <= other_end and end >= other_start


def empty_group(gene: str, row: dict[str, str]) -> dict[str, Any]:
    return {
        "gene": gene,
        "strand": row["gene_strand"],
        "contig": row["contig"],
        "targets": {
            "start": [],
            "end": [],
            "target_name": [],
            "exon_number": [],
            "x": [],
        },
        "samples": [],
        "calls": [],
        "_target_map": {},
        "_sample_map": OrderedDict(),
        "_call_map": OrderedDict(),
    }


def build_payloads(rows: list[dict[str, str]], ntc_regex: str) -> list[dict[str, Any]]:
    """Hoist target metadata and build parallel arrays for each sample."""

    pattern = re.compile(ntc_regex)
    groups: OrderedDict[str, dict[str, Any]] = OrderedDict()
    for row in rows:
        gene = row["gene"]
        group = groups.setdefault(gene, empty_group(gene, row))
        target_key = (row["contig"], as_int(row["start"]), as_int(row["end"]))
        group["_target_map"].setdefault(
            target_key,
            {
                "start": as_int(row["start"]),
                "end": as_int(row["end"]),
                "target_name": row["annotation_target_name"] or row["target_name"],
                "exon_number": as_int(row["exon_number"]),
            },
        )
        sample_map = group["_sample_map"]
        if row["sample_id"] not in sample_map:
            sample_map[row["sample_id"]] = {
                "sample_id": row["sample_id"],
                "ntc": bool(pattern.search(row["sample_id"])),
                "values": {},
            }
        sample_map[row["sample_id"]]["values"][target_key] = {
            "tn_log2": float(row["tn_log2"]),
            "segment_mean": as_float(row["segment_mean"]),
        }

        if row["call_type"]:
            call_key = (
                row["sample_id"],
                row["call_type"],
                as_int(row["call_copy_number"]),
                row["call_filter"],
                row["call_qual"],
                as_int(row["call_start"]),
                as_int(row["call_end"]),
            )
            group["_call_map"].setdefault(
                call_key,
                {
                    "sample_id": row["sample_id"],
                    "type": row["call_type"],
                    "copy_number": as_int(row["call_copy_number"]),
                    "filter": row["call_filter"],
                    "qual": as_float(row["call_qual"]),
                    "start": as_int(row["call_start"]),
                    "end": as_int(row["call_end"]),
                },
            )

    payloads: list[dict[str, Any]] = []
    for group in groups.values():
        target_keys = sorted(group["_target_map"], key=lambda key: (key[1], key[2]))
        for index, key in enumerate(target_keys):
            target = group["_target_map"][key]
            for field in ("start", "end", "target_name", "exon_number"):
                group["targets"][field].append(target[field])
            group["targets"]["x"].append(index)

        for sample in group["_sample_map"].values():
            group["samples"].append(
                {
                    "sample_id": sample["sample_id"],
                    "ntc": sample["ntc"],
                    "tn_log2": [sample["values"][key]["tn_log2"] for key in target_keys],
                    "segment_mean": [sample["values"][key]["segment_mean"] for key in target_keys],
                }
            )

        for call in group["_call_map"].values():
            call["exons"] = sorted(
                {
                    group["_target_map"][key]["exon_number"]
                    for key in target_keys
                    if overlaps(call["start"], call["end"], key[1], key[2])
                }
            )
            group["calls"].append(call)

        for private_key in ("_target_map", "_sample_map", "_call_map"):
            del group[private_key]
        payloads.append(group)
    return payloads


def compact_json(value: Any) -> str:
    """Serialize without whitespace; stable insertion order keeps output deterministic."""

    return json.dumps(value, ensure_ascii=True, separators=(",", ":"), allow_nan=False)


def html_document(
    payloads: list[dict[str, Any]],
    title: str,
    default_group: str,
    ntc_pattern: str,
    runtime: str,
    plotly_script: str,
) -> tuple[str, int]:
    group_scripts: list[str] = []
    largest = 0
    for index, payload in enumerate(payloads):
        serialized = compact_json(payload).replace("</script", "<\\/script")
        largest = max(largest, len(serialized))
        group_scripts.append(
            f'<script type="application/json" id="cnvgrp-{index}">{serialized}</script>'
        )
    config = compact_json(
        {
            "defaultGroup": default_group,
            "ntcPattern": ntc_pattern,
            "groupCount": len(payloads),
            "groupIds": {payload["gene"]: f"cnvgrp-{index}" for index, payload in enumerate(payloads)},
        }
    )
    title_text = html.escape(title)
    source_script = plotly_script
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title_text}</title>
<style>
:root {{
  color-scheme: light;
  --ink: #17212b;
  --muted: #61707e;
  --line: #d9e1e7;
  --panel: #ffffff;
  --wash: #f4f7f9;
  --accent: #0f6b78;
  --accent-soft: #e2f1f3;
  --radius: 10px;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--wash); color: var(--ink); font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; font-size: 14px; }}
main {{ max-width: 1540px; margin: 0 auto; padding: 26px 34px 46px; }}
header {{ display: flex; align-items: flex-end; justify-content: space-between; gap: 24px; margin-bottom: 18px; }}
h1 {{ margin: 0; font-size: clamp(25px, 2.2vw, 34px); letter-spacing: -0.035em; line-height: 1.08; font-weight: 700; }}
.kicker {{ margin: 0 0 8px; color: var(--accent); font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px; letter-spacing: .13em; text-transform: uppercase; }}
.notice {{ margin: 9px 0 0; color: var(--muted); }}
.notice code {{ color: var(--ink); background: #e8eef1; border-radius: 4px; padding: 2px 5px; }}
.summary {{ color: var(--muted); font-size: 12px; text-align: right; white-space: nowrap; }}
.toolbar, .plot-panel, .calls-panel {{ background: var(--panel); border: 1px solid var(--line); border-radius: var(--radius); }}
.toolbar {{ display: flex; align-items: center; flex-wrap: wrap; gap: 14px 18px; padding: 13px 16px; margin-bottom: 16px; }}
.control {{ display: inline-flex; align-items: center; gap: 7px; color: var(--muted); }}
.control label {{ font-size: 12px; font-weight: 600; }}
input[type="search"], select {{ border: 1px solid #bcc9d1; border-radius: 7px; color: var(--ink); background: #fff; padding: 8px 10px; font: inherit; min-height: 35px; }}
input[type="search"] {{ width: 180px; }}
select {{ min-width: 165px; }}
input:focus, select:focus {{ outline: 3px solid var(--accent-soft); border-color: var(--accent); }}
.checks {{ display: inline-flex; align-items: center; flex-wrap: wrap; gap: 12px; margin-left: auto; }}
.check {{ display: inline-flex; align-items: center; gap: 6px; color: var(--muted); font-size: 12px; cursor: pointer; white-space: nowrap; }}
.check input {{ accent-color: var(--accent); }}
.plot-panel {{ padding: 8px 12px 0; min-height: 520px; }}
#plot {{ width: 100%; height: 520px; }}
.status {{ padding: 0 16px 12px; color: var(--muted); font-size: 12px; min-height: 17px; }}
.calls-panel {{ margin-top: 16px; padding: 19px 20px 20px; }}
.calls-heading {{ display: flex; justify-content: space-between; align-items: baseline; gap: 12px; margin-bottom: 12px; }}
h2 {{ margin: 0; font-size: 17px; letter-spacing: -.02em; }}
.strand-chip {{ display: inline-flex; align-items: center; gap: 6px; border-radius: 999px; background: var(--accent-soft); color: var(--accent); padding: 4px 9px; font-size: 11px; font-weight: 700; }}
.strand-chip .arrow {{ font-size: 14px; line-height: 1; }}
.table-wrap {{ overflow-x: auto; }}
table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
th {{ color: var(--muted); font-weight: 700; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--line); padding: 8px 9px; }}
td {{ border-bottom: 1px solid #edf1f3; padding: 9px; vertical-align: middle; white-space: nowrap; }}
tr:last-child td {{ border-bottom: 0; }}
.sample-cell {{ display: inline-flex; align-items: center; gap: 7px; }}
.sample-swatch {{ display: inline-block; width: 9px; height: 9px; border-radius: 50%; }}
.filter-tag {{ display: inline-block; border-radius: 4px; padding: 3px 6px; background: #edf2f4; color: var(--ink); font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px; }}
.filter-tag.pass {{ background: #e4f3e7; color: #24643b; }}
.empty {{ color: var(--muted); padding: 18px 9px; }}
.definitions {{ margin-top: 17px; border-top: 1px solid var(--line); padding-top: 12px; color: var(--muted); font-size: 12px; }}
.definitions p {{ margin: 4px 0; }}
.definitions strong {{ color: var(--ink); font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px; }}
.error {{ color: #a53d3d; font-weight: 600; }}
@media (max-width: 760px) {{
  main {{ padding: 18px 14px 30px; }}
  header {{ display: block; }}
  .summary {{ margin-top: 10px; text-align: left; }}
  .checks {{ width: 100%; margin-left: 0; }}
  #plot {{ height: 460px; }}
  .plot-panel {{ min-height: 460px; }}
}}
</style>
</head>
<body>
<main>
  <header>
    <div>
      <p class="kicker">CNV review aid / synthetic run</p>
      <h1>{title_text}</h1>
      <p class="notice">Synthetic data only. Coordinates and calls are illustrative, not clinical results.</p>
    </div>
    <div class="summary">{len(payloads)} genes / 14 samples / lazy per-gene payloads</div>
  </header>
  <section class="toolbar" aria-label="Plot controls">
    <div class="control"><label for="gene-search">Gene</label><input id="gene-search" type="search" list="gene-list" placeholder="Search genes" autocomplete="off"><datalist id="gene-list">{''.join(f'<option value="{html.escape(payload["gene"], quote=True)}"></option>' for payload in payloads)}</datalist></div>
    <div class="control"><label for="metric">Y metric</label><select id="metric"><option value="normalized_ratio">normalized_ratio</option><option value="tn_log2">tn_log2</option><option value="segment_mean">segment_mean</option></select></div>
    <div class="checks">
      <label class="check"><input id="toggle-reference" type="checkbox" checked> reference lines</label>
      <label class="check"><input id="toggle-bands" type="checkbox" checked> PASS bands</label>
      <label class="check"><input id="toggle-strand" type="checkbox" checked> strand reversal</label>
      <label class="check"><input id="toggle-ntc" type="checkbox" checked> NTC styling</label>
    </div>
  </section>
  <section class="plot-panel" aria-label="Exon-level copy-number plot">
    <div id="plot" aria-live="polite"></div>
    <div id="status" class="status"></div>
  </section>
  <section class="calls-panel" aria-labelledby="calls-heading">
    <div class="calls-heading"><h2 id="calls-heading">Calls in <span id="current-gene">{html.escape(default_group)}</span></h2><span id="strand-chip" class="strand-chip"></span></div>
    <div id="calls-table" class="table-wrap"></div>
    <div id="definitions" class="definitions"></div>
  </section>
</main>
{''.join(group_scripts)}
<script id="cnv-config" type="application/json">{config}</script>
{source_script}
<script>{runtime}</script>
</body>
</html>
"""
    return document, largest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path, help="long CSV from build_plot_table.py")
    parser.add_argument("-o", "--output", required=True, type=Path, help="HTML output path")
    parser.add_argument("--title", default="DRAGEN CNV Plot Builder", help="page title")
    parser.add_argument("--default-group", default="", help="gene selected on first load")
    parser.add_argument("--ntc-pattern", default=r"(?i)\bntc\b", help="regex used to identify NTC samples")
    parser.add_argument("--plotly-js-path", type=Path, help="local Plotly file to inline")
    parser.add_argument("--plotly-src", default=DEFAULT_PLOTLY_SRC, help="Plotly script URL or relative path")
    parser.add_argument("--inline-plotly", action="store_true", help="embed the local Plotly JavaScript")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_rows(args.table)
    if not rows:
        raise ValueError(f"No rows found in {args.table}")
    payloads = build_payloads(rows, args.ntc_pattern)
    if not payloads:
        raise ValueError("No gene groups found in the table")
    default_group = args.default_group or payloads[0]["gene"]
    if default_group not in {payload["gene"] for payload in payloads}:
        raise ValueError(f"Default group {default_group!r} is absent from the table")

    runtime_path = Path(__file__).with_name("plot_runtime.js")
    runtime = runtime_path.read_text(encoding="utf-8")
    if args.inline_plotly:
        plotly_path = args.plotly_js_path or DEFAULT_PLOTLY_PATH
        plotly_script = f"<script>{plotly_path.read_text(encoding='utf-8')}</script>"
    else:
        plotly_script = ""
        if args.plotly_js_path:
            raise ValueError("--plotly-js-path requires --inline-plotly; use --plotly-src for a script path")

    if args.inline_plotly:
        document, largest = html_document(
            payloads, args.title, default_group, args.ntc_pattern, runtime, plotly_script
        )
    else:
        script_tag = f'<script src="{html.escape(args.plotly_src, quote=True)}"></script>'
        document, largest = html_document(
            payloads, args.title, default_group, args.ntc_pattern, runtime, script_tag
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(document, encoding="utf-8", newline="\n")
    percentage = (largest / V8_STRING_LIMIT) * 100
    print(f"Largest embedded JSON string: {largest} chars / {V8_STRING_LIMIT} limit ({percentage:.6f}% used)")
    print(f"Wrote {args.output} ({len(document.encode('utf-8'))} bytes)")


if __name__ == "__main__":
    main()
