# Synthetic CNV Plot Demo (Pure Rust)

> **Unofficial showcase.** Independent open-source project — not affiliated with or endorsed by Illumina, Inc. All coordinates, ratios, and calls generated here are deterministic synthetic stand-ins containing zero patient data or proprietary panel designs.

This directory contains a self-contained, reproducible **pure-Rust** synthetic pipeline (`cnvplot-showcase-demo`) mirroring the v4 single-pass multi-core architecture:

1. **`src/synth.rs`** — Deterministically synthesizes 14 DRAGEN-shaped sample folders (`12` samples + `2` NTC controls across `12` genes) from a portable `SplitMix64` PRNG (`seed = 20260729`).
2. **`src/engine.rs`** — Parallel multi-core decompression (`rayon` + `flate2::read::MultiGzDecoder`), 100 kb binned interval overlap matching, single-pass in-memory gene grouping (zero intermediate disk CSV required), and lazy `<script type="application/json">` HTML emission.
3. **`src/plot_runtime.js`** — Browser-side interactivity embedded into the generated page.

## Run from the repository root

```bash
cargo run --release --manifest-path demo/Cargo.toml
```

Useful flags:

```bash
# Retain the synthetic sample folders and also emit the optional long CSV table:
cargo run --release --manifest-path demo/Cargo.toml -- --keep-work --emit-csv

# Build a single-file offline HTML with Plotly inlined directly:
cargo run --release --manifest-path demo/Cargo.toml -- --inline-plotly -o /tmp/cnv-demo-offline.html
```

## Planted synthetic teaching examples

| Gene | Sample | Planted pattern | Review takeaway |
|---|---|---|---|
| **BRCA1** | `DEMO-RUN01-S03` | Whole-gene ratio ~0.50, `PASS` `DEL` (`CN=1`) | Clean exon-wide shift in one sample only |
| **EGFR** | `DEMO-RUN01-S07` | Exons 18–21 ~1.50, `PASS` `DUP` (`CN=3`) | Focal gain flanked by flat diploid exons |
| **PTEN** | `DEMO-RUN01-S05` | Single exon 6 ~0.45, filtered `cnvLength` | Noisy single target — why filtered calls remain visible |
| **CFTR** | *all samples* | Exons 10–11 dip ~0.65 across the entire batch | Panel capture/GC artifact, immediately obvious via batch overlay |
| **DMD** | `DEMO-RUN01-S09` | Exons 8–13 ~0.50 `PASS` `DEL` (`-` strand) | `5' -> 3'` strand reversal reads in transcript order |
| **TP53** | `DEMO-RUN01-S02` | Exons 4–6 ~0.55, filtered `cnvQual` | Borderline rejected event surfaced for human review |
