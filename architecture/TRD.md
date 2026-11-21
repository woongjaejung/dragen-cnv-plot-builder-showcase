# TRD — Technical Architecture Specification (v4 · Pure Rust + Tauri v2)

> **Disclaimer:** Independent, unofficial open-source project — not affiliated with, endorsed by, or sponsored by Illumina, Inc.

## 1. Architectural Shift (`v3` → `v4`)

| Layer | Legacy `v3` | `v4` Native Architecture | Engineering Rationale |
|---|---|---|---|
| **Desktop Shell** | Electron 41 (bundled Chromium) | **Tauri v2** (OS-native WebView2 / WKWebView) | Drops ~120 MB Chromium payload |
| **Compute Engine** | Bundled CPython (`23 MB`) spawning two scripts | **Single-process Rust** (`rayon` + `flate2` `zlib-rs`) | Multi-core execution, zero subprocess overhead |
| **Inter-step Data** | `1–1.6 GB` disk CSV written once & read twice | **In-memory columnar structs** | Eliminates 3x disk serialization passes |
| **Default Annotation** | `46.5 MB` `.gtf.gz` text (~3.2M lines parsed/run) | **Lossless `GCV4` binary index (`13.7 MB`)** | Loads `1,499,012` exons in ~`0.1 s` |
| **Distribution** | Self-extracting archive (`417 MB` `%TEMP%` unpack) | **Single PE32+ `.exe` (`23 MB`)** via `include_bytes!` | Instant cold start, zero `%TEMP%` extraction |

---

## 2. Lossless `GCV4` Binary Exon Index

`gencode.v38.annotation.gtf.gz` contains ~3.2M lines, of which `1,499,012` are `exon` features. Only nine ranking attributes affect target assignment. `GCV4` interns repeated gene/transcript strings into a compact `StringPool` (`u32` IDs) and packs each exon into a 40-byte in-memory struct (`OverlapExon`):

- **Bitpacked flags (`u8`):** `MANE_Select`, `Ensembl_canonical`, `appris_principal_1`, `basic`, `protein_coding`, and strand (`+`/`-`).
- **9-tier deterministic comparator:**
  `(-overlap_len, MANE_Select, Ensembl_canonical, appris_principal_1, basic, protein_coding, level, exon_number_int, transcript_name, file_order)` — mathematically isomorphic to the Python reference `annotation_rank` tuple.
- **Custom fallback:** Supplying any external `.gtf(.gz)`, `.gff3(.gz)`, `.bed(.gz)`, `.tsv`, `.csv`, or `.xlsx` file bypasses `GCV4` and streams through the corresponding Rust parser (`calamine`, `csv`, `MultiGzDecoder`).

---

## 3. Deterministic Multi-Core Parallelism (`rayon`)

To parallelize across all CPU cores without changing a single output byte:

1. **Sample parse map:** `samples.par_iter().map(process_sample_files).collect::<Result<Vec<_>>>()` decompresses concatenated 64 KB BGZF blocks (`flate2::read::MultiGzDecoder` backed by SIMD `zlib-rs`) and matches `.seg` / `.cnv.gff3` intervals via 100 kb genomic bins while preserving exact input sample order.
2. **Batch-unique target annotation:** Unique `(contig, start, end)` panel coordinates (~190k keys) are deduplicated once and annotated in parallel (`into_par_iter()`), maintaining `O(targets)` scaling regardless of batch size.
3. **IEEE-754 `round4` parity:** Stack-buffered Dragon4 formatting (`write!("{:.10}", raw)` → `write!("{:.4}", step1)`) mirrors Python's two-stage `:.10f` + `round(x, 4)` round-half-to-even behavior with `0 ULP` divergence.
4. **Parallel group JSON serialization:** All ~7,768 gene groups serialize concurrently via `group_names.par_iter()` before streaming sequentially into `<script id="g{idx}" type="application/json">` tags.

---

## 4. Cross-Compilation Pipeline (macOS → Windows x64)

Windows `x86_64-pc-windows-msvc` binaries and user-mode NSIS installers are cross-compiled directly from Apple Silicon macOS using `cargo-xwin` + LLVM (`clang-cl`, `lld-link`, `llvm-rc`) + `makensis`, with `--remap-path-prefix` stripping host filesystem paths from `.rdata` panic locations.
