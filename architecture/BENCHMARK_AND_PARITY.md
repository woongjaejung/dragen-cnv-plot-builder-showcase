# Golden-Master Parity & Benchmark Report (`v3` Python vs `v4` Rust)

> **Disclaimer:** Independent, unofficial open-source benchmark report — not affiliated with or endorsed by Illumina, Inc. All measurements were collected on de-identified/synthetic exome structures with zero patient identifiers.

## 1. Differential Verification Methodology

Inspired by behavior-first migration harnesses (`rust-anything`), a differential verification script parses `<script id="plot-meta">` and every `<script id="g0">...<script id="g7767">` payload from both the Python reference (`v3`) and the Rust engine (`v4`), asserting exact equality across:

- Sample ordering, NTC regex classification, group list, and skipped unannotated row counts
- Per-gene target metadata (`targets`), sample indices (`s`), exon ranks (`x`), call indices (`c`), and dominant strand (`strand`)
- Exon tick midpoints (`tv`), tick labels (`tt`), and exon separator lines (`bp`)
- Sorted `PASS` and `FILTERED` call tables (`calls`), exon range summaries, and target counts
- Floating-point arrays `v` (`tn_log2`) and `m` (`segment_mean`) at `0.0` tolerance (`0 ULP` difference)

### Parity Results

| Dataset | Mode | Gene Groups | Plotted Points | Overlapping Calls | Float Arrays Checked | ULP Diffs | Verdict |
|---|---|---:|---:|---:|---:|---:|---|
| **4-sample exome batch** | Raw `.gtf.gz` | 7,768 | 653,911 | 941 | 15,536 | **0** | **PASS** |
| **4-sample exome batch** | `GCV4` cache | 7,768 | 653,911 | 941 | 15,536 | **0** | **PASS** |
| **30-sample synthetic load batch** | `GCV4` cache | 7,768 | 4,908,720 | 7,620 | 15,536 | **0** | **PASS** |

> **Defect caught during verification:** Initial testing with single-member `flate2::read::GzDecoder` halted at the first 64 KB BGZF block boundary of `*.dragen.tn.tsv.gz`, silently yielding only 40 genes on `chr1p36`. The golden-master check caught the group-count mismatch immediately; switching to `MultiGzDecoder` restored all `7,768` gene groups with `0 ULP` differences.

---

## 2. Wall-Clock Runtime & Resource Measurements

Measured on Apple Silicon (`arm64`, local NVMe SSD) via `/usr/bin/time -l`:

| Workload | Pipeline | Step 1 | Step 2 | **Total Wall Time** | **Speedup** | Peak RSS | Temp Disk CSV |
|---|---|---:|---:|---:|---:|---:|---:|
| **4 samples** (`653,911` pts) | `v3` Python 3.14 | 48.37 s | 5.74 s | **54.11 s** | 1.0x | 1,489.2 MB | Written & deleted |
| **4 samples** (`653,911` pts) | `v4` Rust (`--no-cache`) | — | — | **3.36 s** | **16.1x** | **651.9 MB** | **0 MB** |
| **4 samples** (`653,911` pts) | **`v4` Rust (`GCV4` cache)** | — | — | **1.56 s** | **34.7x** | **642.1 MB** | **0 MB** |
| **30 samples** (`4,908,720` pts) | `v3` Python 3.14 | 128.11 s | 31.96 s | **160.07 s** | 1.0x | 1,835.4 MB | **948.9 MB** |
| **30 samples** (`4,908,720` pts) | `v4` Rust (`--no-cache`) | — | — | **10.66 s** | **15.0x** | 2,329.6 MB | **0 MB** |
| **30 samples** (`4,908,720` pts) | **`v4` Rust (`GCV4` cache)** | — | — | **6.59 s** | **24.3x** | 2,378.6 MB | **0 MB** |
