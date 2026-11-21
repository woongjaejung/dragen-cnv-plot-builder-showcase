# PRD — Exon-Level CNV Plot Builder (v4 · Pure Rust + Tauri v2)

> **Disclaimer:** Independent, unofficial open-source project — not affiliated with, endorsed by, or sponsored by Illumina, Inc. Not an in-vitro diagnostic (IVD) device.

## 1. Problem Statement

Segment-level CNV callers emit a copy number, a quality score (`QUAL`), and a filter tag (`PASS`, `cnvQual`, `dinucQual`, `cnvLikelihoodRatio`). Determining whether an exon-level call is a **true biological CNV** or a **panel capture/GC artifact** requires five concurrent signals that raw tabular outputs split across three files per sample:

1. **Per-exon target coverage** rather than segment averages alone.
2. **Batch-wide overlay** across all samples sequenced on the same capture panel.
3. **Expected copy-state baselines** (`1.0` diploid, `0.5` het del, `1.5` het dup).
4. **No-template control (NTC) behavior** at identical coordinates.
5. **Both `PASS` and `FILTERED` calls** alongside filter rationales and `5' → 3'` transcript orientation for negative-strand genes.

While `v2`–`v3` solved browser rendering scale (eliminating the 20-sample V8 string limit), the Python + Electron implementation still took **~160 seconds** at 30 samples, wrote and re-scanned nearly **1 GB** of temporary CSV data per run, and distributed as a **146 MB** (`417 MB` unpacked) bundle.

---

## 2. v4 Product Goals & Target KPIs

| Metric | Legacy `v3` (Python + Electron) | `v4` Target (Pure Rust + Tauri v2) | Achieved |
|---|---:|---:|---:|
| **4-sample batch build time** | 54.11 s | < 5 s | **1.56 s (34.7x)** |
| **30-sample batch build time** | 160.07 s | < 10 s | **6.59 s (24.3x)** |
| **Temporary CSV disk I/O** | 948.9 MB | 0 MB | **0 MB** |
| **Portable Windows `.exe` size** | 146 MB (417 MB unpacked) | < 30 MB single file | **23 MB (0 MB unpack)** |
| **Output numerical parity** | Baseline | 100% identical | **0 ULP diffs (`4.91M` pts)** |
| **Offline capability** | Full offline | Full offline | **100% offline** |

---

## 3. Target Users

1. **Genomics Reviewers & Variant Scientists** operating on locked-down Windows 10/11 workstations without Python or internet access.
2. **Bioinformatics Engineers** running headless batch pipelines via the standalone Rust CLI (`cnvplot-cli`).

---

## 4. Core Functional Scope

- **3-Step Guided Native Desktop Workflow:** Recursive folder discovery (up to depth 4, native OS folder drag-and-drop), KPI summary cards (`Samples`, `SEG`, `CNV GFF3`, `NTC`), live NTC regex preview, 4-stage pipeline progress bar, and one-click offline HTML launch.
- **Multi-Format Annotation Support:** Embedded pre-indexed hg38 GENCODE v38 (`GCV4` binary index) plus external `.gtf(.gz)`, `.gff3(.gz)`, `.bed(.gz)`, `.tsv`, `.csv`, and `.xlsx` panel annotations.
- **Self-Contained Interactive Report:** Inlined Plotly 2.35.2, lazy per-gene `<script>` parsing, metric switching (`Normalized Ratio`, `TN Log2`, `Segment Mean`), merged `PASS` bands, negative-strand `5' → 3'` reversal, and `PASS` + `FILTERED` call tables.
