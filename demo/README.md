# Synthetic CNV plot demo

This directory contains a small, reproducible pipeline for the public showcase.
It generates synthetic data, joins the data into one long table, and writes an
interactive Plotly page. The generated data contains no patient data, clinical
accession IDs, credentials, or real panel coordinates.

The files produced here are simplified stand-ins for DRAGEN's real outputs,
not byte-level replicas. The gene symbols are real symbols, but the genomic
coordinates and every measurement are illustrative and synthetic.

## Run the demo

From the repository root:

```bash
python3 demo/run_demo.py
```

This writes:

```text
demo/_work/synthetic_run/       generated sample folders and panel TSV
demo/_work/cnv_plot_table.csv   the step-one long table
docs/demo/index.html            the published demo page
```

The default run removes the generated sample folders after the page is built,
but leaves the long table for inspection and for testing the HTML builder. Use
`--keep-work` to retain every generated input file:

```bash
python3 demo/run_demo.py --keep-work
```

The default seed is `20260729`. The generator can be run independently with a
different output directory and seed:

```bash
python3 demo/generate_synthetic_run.py --out demo/_work/synthetic_run --seed 20260729
```

Step one accepts one or more sample directories:

```bash
python3 demo/build_plot_table.py demo/_work/synthetic_run/DEMO-RUN01-S01 \
  demo/_work/synthetic_run/DEMO-RUN01-S02 \
  --annotation demo/_work/synthetic_run/synthetic_panel_targets.tsv \
  -o demo/_work/cnv_plot_table.csv --slim
```

`--dry-run` parses and joins the inputs, prints the annotation-cache summary,
and does not write a CSV. `--slim` removes the repeated `result_dir` column.

Step two can reference the vendored Plotly file or inline it:

```bash
python3 demo/build_plot_html.py demo/_work/cnv_plot_table.csv \
  -o docs/demo/index.html --plotly-src ../vendor/plotly-2.35.2.min.js

python3 demo/build_plot_html.py demo/_work/cnv_plot_table.csv \
  -o /tmp/cnv-demo-offline.html --inline-plotly
```

The published page uses the relative vendor path. That is a local repository
resource, not a CDN request. `--inline-plotly` reads
`docs/vendor/plotly-2.35.2.min.js` and embeds the same Plotly build in the HTML
for a single-file offline example.

## Synthetic input formats

The generated run has this shape:

```text
demo/_work/synthetic_run/
  DEMO-RUN01-S01/
    vcf/dragen/DEMO-RUN01-S01/cnv_exome_varcaller/
      DEMO-RUN01-S01.dragen.tn.tsv.gz
      DEMO-RUN01-S01.dragen.seg
      DEMO-RUN01-S01.dragen.cnv.gff3
  synthetic_panel_targets.tsv
```

The same layout is repeated for all twelve patient-like synthetic samples and
the two NTC samples. The `*.dragen.tn.tsv.gz` file is required. It is a gzipped
TSV whose header starts with `#` and whose columns are:

```text
contig  start  end  target_name  tn_log2  improper_pairs
```

`tn_log2` is the tangent-normalized log2 measurement. Step one derives
`normalized_ratio` as `2 ** tn_log2`.

The `*.dragen.seg` file is a TSV with a header and these columns:

```text
contig  start  end  num_probes  segment_mean
```

The `*.dragen.cnv.gff3` file is a small GFF3-shaped call file. Column 3 is
`DEL` or `DUP`. Its column 9 attributes include all of the following keys:

```text
ID=<id>;CN=<int>;Filter=<tag>;QUAL=<float>
```

The synthetic filter tags are `PASS`, `cnvQual`, and `cnvLength`. The demo
runtime explains only the tags present for the selected gene.

The shared `synthetic_panel_targets.tsv` annotation file is a TSV with this
header and columns:

```text
chrom  start  end  name  GeneSymbol  exonNumber  strand
```

Each target belongs to one of twelve genes and one exon. Each exon has one to
four capture targets between roughly 120 and 200 bases long with gaps between
targets. Coordinates ascend within each contig and do not overlap between
genes. They are illustrative coordinates rather than a real panel design.

## Table and browser data layout

Step one emits these columns in this exact order unless `--slim` removes
`result_dir`:

```text
sample_id, result_dir, contig, start, end, target_name,
annotation_target_name, gene, gene_strand, exon_number, tn_log2,
normalized_ratio, improper_pairs, segment_start, segment_end,
segment_probes, segment_mean, call_type, call_copy_number, call_filter,
call_qual, call_start, call_end, call_count, call_pass_count
```

Annotation lookup is cached by `(contig, start, end)` across samples. A target
can overlap more than one call: the count fields include every overlap, while
the reported call fields choose a `PASS` call before a filtered call.

Step two puts one compact JSON payload in one `<script>` element per gene. The
payload hoists target metadata once, stores sample points in parallel arrays,
and omits `normalized_ratio`; the browser derives it from `tn_log2`. The page
parses a gene payload on first view rather than parsing the whole batch at
startup.

## Planted examples

The seeded batch contains these synthetic teaching examples:

| Gene | Sample | Synthetic pattern |
| --- | --- | --- |
| BRCA1 | S03 | Whole-gene ratio near 0.50 with a PASS DEL, CN 1 |
| EGFR | S07 | Exons 18-21 near 1.50 with a PASS DUP, CN 3 |
| PTEN | S05 | Exon 6 near 0.45 with a `cnvLength` call |
| CFTR | all samples | Exons 10-11 near 0.65 in non-NTC samples, with no call |
| DMD | S09 | Exons 8-13 near 0.50 with a PASS DEL on a negative strand gene |
| TP53 | S02 | Exons 4-6 near 0.55 with a low-QUAL `cnvQual` call |

All other gene/sample combinations are diploid noise. NTC samples are near a
small positive ratio floor with larger relative scatter. The CFTR footprint is
applied multiplicatively to NTCs so they remain near zero while retaining the
same shared target footprint.
