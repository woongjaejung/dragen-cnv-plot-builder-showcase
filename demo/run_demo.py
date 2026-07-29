#!/usr/bin/env python3
"""Run the complete deterministic synthetic CNV demo pipeline."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEMO_DIR = REPO_ROOT / "demo"
WORK_DIR = DEMO_DIR / "_work"
SYNTHETIC_DIR = WORK_DIR / "synthetic_run"
TABLE_PATH = WORK_DIR / "cnv_plot_table.csv"
HTML_PATH = REPO_ROOT / "docs/demo/index.html"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="keep the generated sample folders as well as the long table",
    )
    return parser.parse_args()


def run(command: list[str]) -> None:
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def main() -> None:
    args = parse_args()
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    run([sys.executable, str(DEMO_DIR / "generate_synthetic_run.py"), "--out", str(SYNTHETIC_DIR)])

    sample_dirs = sorted(
        path for path in SYNTHETIC_DIR.iterdir() if path.is_dir()
    )
    run(
        [
            sys.executable,
            str(DEMO_DIR / "build_plot_table.py"),
            *[str(path) for path in sample_dirs],
            "--annotation",
            str(SYNTHETIC_DIR / "synthetic_panel_targets.tsv"),
            "-o",
            str(TABLE_PATH),
            "--slim",
        ]
    )
    run(
        [
            sys.executable,
            str(DEMO_DIR / "build_plot_html.py"),
            str(TABLE_PATH),
            "-o",
            str(HTML_PATH),
            "--title",
            "DRAGEN CNV Plot Builder \u2013 synthetic demo",
            "--default-group",
            "BRCA1",
            "--plotly-src",
            "../vendor/plotly-2.35.2.min.js",
        ]
    )
    if not args.keep_work:
        shutil.rmtree(SYNTHETIC_DIR)
        print(f"Removed generated sample folders; kept {TABLE_PATH} for inspection")
    print(f"Demo page: {HTML_PATH}")


if __name__ == "__main__":
    main()
