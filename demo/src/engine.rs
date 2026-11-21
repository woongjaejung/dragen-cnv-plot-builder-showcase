//! Single-pass multi-core Rust showcase engine (`rayon` + `MultiGzDecoder` + 100kb bins).
//!
//! Demonstrates the v4 architecture on synthetic DRAGEN folders:
//! - Target-keyed `O(targets)` annotation resolution
//! - Parallel per-sample BGZF/gzip decompression and 100kb binned interval join
//! - Single-pass in-memory per-gene aggregation (no intermediate multi-hundred-MB CSV required)
//! - Per-gene lazy `<script type="application/json">` HTML emission to stay far below V8's
//!   `String::kMaxLength` (`536,870,888` characters).

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use regex::Regex;
use rustc_hash::FxHashMap;
use serde::Serialize;

pub const V8_STRING_LIMIT: usize = 536_870_888;
const BIN_SIZE: i64 = 100_000;
static PLOT_RUNTIME_JS: &str = include_str!("plot_runtime.js");

#[derive(Clone, Debug)]
struct PanelAnnotation {
    gene: String,
    exon_number: i64,
    target_name: String,
    strand: String,
}

#[derive(Clone, Debug)]
struct SegInterval {
    start: i64,
    end: i64,
    probes: i64,
    mean: f64,
}

#[derive(Clone, Debug)]
struct CallInterval {
    start: i64,
    end: i64,
    call_type: String,
    copy_number: i64,
    filter: String,
    qual: f64,
}

#[derive(Clone, Debug)]
struct TargetRow {
    contig: String,
    start: i64,
    end: i64,
    target_name: String,
    tn_log2: f64,
    improper_pairs: i64,
    seg: Option<SegInterval>,
    best_call: Option<CallInterval>,
    call_count: usize,
    call_pass_count: usize,
}

struct SampleData {
    sample_id: String,
    result_dir: PathBuf,
    rows: Vec<TargetRow>,
}

#[derive(Serialize)]
struct TargetColumns {
    start: Vec<i64>,
    end: Vec<i64>,
    target_name: Vec<String>,
    exon_number: Vec<i64>,
    x: Vec<usize>,
}

#[derive(Serialize)]
struct SampleTracePayload {
    sample_id: String,
    ntc: bool,
    tn_log2: Vec<f64>,
    segment_mean: Vec<Option<f64>>,
}

#[derive(Serialize)]
struct CallPayload {
    sample_id: String,
    #[serde(rename = "type")]
    call_type: String,
    copy_number: i64,
    filter: String,
    qual: f64,
    start: i64,
    end: i64,
    exons: Vec<i64>,
}

#[derive(Serialize)]
struct GeneGroupPayload {
    gene: String,
    strand: String,
    contig: String,
    targets: TargetColumns,
    samples: Vec<SampleTracePayload>,
    calls: Vec<CallPayload>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PageConfig {
    default_group: String,
    ntc_pattern: String,
    group_count: usize,
    group_ids: BTreeMap<String, String>,
}

fn load_panel_tsv(path: &Path) -> Result<FxHashMap<(String, i64, i64), PanelAnnotation>> {
    let reader = BufReader::new(File::open(path)?);
    let mut out = FxHashMap::default();
    for (idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if idx == 0 || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 7 {
            continue;
        }
        let key = (
            cols[0].to_string(),
            cols[1].parse::<i64>()?,
            cols[2].parse::<i64>()?,
        );
        out.insert(
            key,
            PanelAnnotation {
                target_name: cols[3].to_string(),
                gene: cols[4].to_string(),
                exon_number: cols[5].parse::<i64>()?,
                strand: cols[6].to_string(),
            },
        );
    }
    Ok(out)
}

fn load_segments(path: &Path) -> Result<FxHashMap<String, FxHashMap<i64, Vec<SegInterval>>>> {
    let mut bins: FxHashMap<String, FxHashMap<i64, Vec<SegInterval>>> = FxHashMap::default();
    if !path.is_file() {
        return Ok(bins);
    }
    let reader = BufReader::new(File::open(path)?);
    for (idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if idx == 0 || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        let seg = SegInterval {
            start: cols[1].parse()?,
            end: cols[2].parse()?,
            probes: cols[3].parse()?,
            mean: cols[4].parse()?,
        };
        let contig_bins = bins.entry(cols[0].to_string()).or_default();
        for b in seg.start.div_euclid(BIN_SIZE)..=seg.end.div_euclid(BIN_SIZE) {
            contig_bins.entry(b).or_default().push(seg.clone());
        }
    }
    Ok(bins)
}

fn load_calls(path: &Path) -> Result<FxHashMap<String, FxHashMap<i64, Vec<CallInterval>>>> {
    let mut bins: FxHashMap<String, FxHashMap<i64, Vec<CallInterval>>> = FxHashMap::default();
    if !path.is_file() {
        return Ok(bins);
    }
    let reader = BufReader::new(File::open(path)?);
    for line_res in reader.lines() {
        let line = line_res?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let mut cn = 2i64;
        let mut filter = String::new();
        let mut qual = 0.0f64;
        for item in cols[8].split(';') {
            if let Some((k, v)) = item.split_once('=') {
                match k {
                    "CN" | "CopyNumber" => cn = v.parse().unwrap_or(2),
                    "Filter" => filter = v.to_string(),
                    "QUAL" | "Qual" => qual = v.parse().unwrap_or(0.0),
                    _ => {}
                }
            }
        }
        let call = CallInterval {
            start: cols[3].parse()?,
            end: cols[4].parse()?,
            call_type: cols[2].to_string(),
            copy_number: cn,
            filter,
            qual,
        };
        let contig_bins = bins.entry(cols[0].to_string()).or_default();
        for b in call.start.div_euclid(BIN_SIZE)..=call.end.div_euclid(BIN_SIZE) {
            contig_bins.entry(b).or_default().push(call.clone());
        }
    }
    Ok(bins)
}

fn parse_one_sample(sample_dir: &Path) -> Result<SampleData> {
    let sample_id = sample_dir
        .file_name()
        .and_then(|n| n.to_str())
        .context("Invalid sample directory")?
        .to_string();
    let caller_dir = sample_dir
        .join("vcf")
        .join("dragen")
        .join(&sample_id)
        .join("cnv_exome_varcaller");

    let tn_path = caller_dir.join(format!("{sample_id}.dragen.tn.tsv.gz"));
    let seg_path = caller_dir.join(format!("{sample_id}.dragen.seg"));
    let gff_path = caller_dir.join(format!("{sample_id}.dragen.cnv.gff3"));

    let seg_bins = load_segments(&seg_path)?;
    let call_bins = load_calls(&gff_path)?;

    let reader = BufReader::new(MultiGzDecoder::new(File::open(&tn_path)?));
    let mut rows = Vec::new();
    for line_res in reader.lines() {
        let line = line_res?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
            continue;
        }
        let contig = cols[0].to_string();
        let start: i64 = cols[1].parse()?;
        let end: i64 = cols[2].parse()?;
        let bin = start.div_euclid(BIN_SIZE);

        let seg = seg_bins
            .get(&contig)
            .and_then(|cb| cb.get(&bin))
            .and_then(|list| list.iter().find(|s| s.start <= end && start <= s.end).cloned());

        let mut matches: Vec<&CallInterval> = call_bins
            .get(&contig)
            .and_then(|cb| cb.get(&bin))
            .map(|list| list.iter().filter(|c| c.start <= end && start <= c.end).collect())
            .unwrap_or_default();

        let call_count = matches.len();
        let call_pass_count = matches.iter().filter(|c| c.filter == "PASS").count();
        matches.sort_by_key(|c| (if c.filter == "PASS" { 0u8 } else { 1u8 }, c.start));
        let best_call = matches.first().map(|c| (*c).clone());

        rows.push(TargetRow {
            contig,
            start,
            end,
            target_name: cols[3].to_string(),
            tn_log2: cols[4].parse()?,
            improper_pairs: cols[5].parse()?,
            seg,
            best_call,
            call_count,
            call_pass_count,
        });
    }

    Ok(SampleData {
        sample_id,
        result_dir: sample_dir.to_path_buf(),
        rows,
    })
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

pub fn build_showcase_html(
    sample_dirs: &[PathBuf],
    annotation_tsv: &Path,
    out_html: &Path,
    emit_csv: Option<&Path>,
    title: &str,
    default_group: &str,
    ntc_pattern: &str,
    plotly_src: &str,
    inline_plotly_path: Option<&Path>,
) -> Result<(usize, usize)> {
    if sample_dirs.is_empty() {
        bail!("No sample directories provided");
    }

    let panel = load_panel_tsv(annotation_tsv)?;
    let ntc_re = Regex::new(ntc_pattern)?;

    // Parallel decompression + interval join across all samples via Rayon
    let parsed_samples: Vec<SampleData> = sample_dirs
        .par_iter()
        .map(|d| parse_one_sample(d))
        .collect::<Result<Vec<_>>>()?;

    if let Some(csv_path) = emit_csv {
        if let Some(parent) = csv_path.parent() {
            fs::create_dir_all(parent)?;
        }
        let mut w = BufWriter::new(File::create(csv_path)?);
        writeln!(
            w,
            "sample_id,contig,start,end,target_name,annotation_target_name,gene,gene_strand,exon_number,tn_log2,normalized_ratio,improper_pairs,segment_start,segment_end,segment_probes,segment_mean,call_type,call_copy_number,call_filter,call_qual,call_start,call_end,call_count,call_pass_count"
        )?;
        for s in &parsed_samples {
            for r in &s.rows {
                let Some(ann) = panel.get(&(r.contig.clone(), r.start, r.end)) else {
                    continue;
                };
                let (ss, se, sp, sm) = match &r.seg {
                    Some(seg) => (
                        seg.start.to_string(),
                        seg.end.to_string(),
                        seg.probes.to_string(),
                        format!("{:.4}", seg.mean),
                    ),
                    None => (String::new(), String::new(), String::new(), String::new()),
                };
                let (ct, ccn, cf, cq, cs, ce) = match &r.best_call {
                    Some(c) => (
                        c.call_type.clone(),
                        c.copy_number.to_string(),
                        c.filter.clone(),
                        format!("{:.1}", c.qual),
                        c.start.to_string(),
                        c.end.to_string(),
                    ),
                    None => (
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                    ),
                };
                writeln!(
                    w,
                    "{},{},{},{},{},{},{},{},{},{:.6},{:.6},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    s.sample_id,
                    r.contig,
                    r.start,
                    r.end,
                    r.target_name,
                    ann.target_name,
                    ann.gene,
                    ann.strand,
                    ann.exon_number,
                    r.tn_log2,
                    2.0_f64.powf(r.tn_log2),
                    r.improper_pairs,
                    ss,
                    se,
                    sp,
                    sm,
                    ct,
                    ccn,
                    cf,
                    cq,
                    cs,
                    ce,
                    r.call_count,
                    r.call_pass_count
                )?;
            }
        }
        w.flush()?;
    }

    // Group in-memory by gene in stable gene-appearance order
    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_seen: BTreeSet<String> = BTreeSet::new();
    for r in &parsed_samples[0].rows {
        if let Some(ann) = panel.get(&(r.contig.clone(), r.start, r.end)) {
            if gene_seen.insert(ann.gene.clone()) {
                gene_order.push(ann.gene.clone());
            }
        }
    }

    let mut payloads: Vec<GeneGroupPayload> = Vec::with_capacity(gene_order.len());
    for gene in &gene_order {
        let mut target_keys: Vec<(String, i64, i64)> = parsed_samples[0]
            .rows
            .iter()
            .filter_map(|r| {
                let ann = panel.get(&(r.contig.clone(), r.start, r.end))?;
                if &ann.gene == gene {
                    Some((r.contig.clone(), r.start, r.end))
                } else {
                    None
                }
            })
            .collect();
        target_keys.sort_by_key(|k| (k.1, k.2));
        target_keys.dedup();

        let first_ann = &panel[&target_keys[0]];
        let contig = target_keys[0].0.clone();
        let strand = first_ann.strand.clone();

        let mut t_cols = TargetColumns {
            start: Vec::with_capacity(target_keys.len()),
            end: Vec::with_capacity(target_keys.len()),
            target_name: Vec::with_capacity(target_keys.len()),
            exon_number: Vec::with_capacity(target_keys.len()),
            x: Vec::with_capacity(target_keys.len()),
        };
        for (idx, k) in target_keys.iter().enumerate() {
            let ann = &panel[k];
            t_cols.start.push(k.1);
            t_cols.end.push(k.2);
            t_cols.target_name.push(ann.target_name.clone());
            t_cols.exon_number.push(ann.exon_number);
            t_cols.x.push(idx);
        }

        let mut sample_traces = Vec::with_capacity(parsed_samples.len());
        let mut calls_out: Vec<CallPayload> = Vec::new();
        let mut seen_calls: BTreeSet<(String, String, i64, String, i64, i64)> = BTreeSet::new();

        for s in &parsed_samples {
            let mut row_by_key: FxHashMap<(String, i64, i64), &TargetRow> = FxHashMap::default();
            for r in &s.rows {
                row_by_key.insert((r.contig.clone(), r.start, r.end), r);
            }
            let mut tn_vec = Vec::with_capacity(target_keys.len());
            let mut seg_vec = Vec::with_capacity(target_keys.len());

            for k in &target_keys {
                let r = row_by_key[k];
                tn_vec.push(r.tn_log2);
                seg_vec.push(r.seg.as_ref().map(|seg| seg.mean));

                if let Some(c) = &r.best_call {
                    let ck = (
                        s.sample_id.clone(),
                        c.call_type.clone(),
                        c.copy_number,
                        c.filter.clone(),
                        c.start,
                        c.end,
                    );
                    if seen_calls.insert(ck) {
                        let mut exons: Vec<i64> = target_keys
                            .iter()
                            .filter(|tk| tk.1 <= c.end && tk.2 >= c.start)
                            .map(|tk| panel[tk].exon_number)
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect();
                        exons.sort();
                        calls_out.push(CallPayload {
                            sample_id: s.sample_id.clone(),
                            call_type: c.call_type.clone(),
                            copy_number: c.copy_number,
                            filter: c.filter.clone(),
                            qual: c.qual,
                            start: c.start,
                            end: c.end,
                            exons,
                        });
                    }
                }
            }

            sample_traces.push(SampleTracePayload {
                sample_id: s.sample_id.clone(),
                ntc: ntc_re.is_match(&s.sample_id),
                tn_log2: tn_vec,
                segment_mean: seg_vec,
            });
        }

        let _ = &s_unused_suppress(&parsed_samples);
        payloads.push(GeneGroupPayload {
            gene: gene.clone(),
            strand,
            contig,
            targets: t_cols,
            samples: sample_traces,
            calls: calls_out,
        });
    }

    let mut group_scripts = Vec::with_capacity(payloads.len());
    let mut group_ids = BTreeMap::new();
    let mut largest = 0usize;

    for (idx, p) in payloads.iter().enumerate() {
        let json = serde_json::to_string(p)?.replace("</script", "<\\/script");
        largest = largest.max(json.len());
        let tag_id = format!("cnvgrp-{idx}");
        group_ids.insert(p.gene.clone(), tag_id.clone());
        group_scripts.push(format!(
            "<script type=\"application/json\" id=\"{tag_id}\">{json}</script>"
        ));
    }

    let eff_default = if group_ids.contains_key(default_group) {
        default_group.to_string()
    } else {
        payloads[0].gene.clone()
    };

    let config_json = serde_json::to_string(&PageConfig {
        default_group: eff_default.clone(),
        ntc_pattern: ntc_pattern.to_string(),
        group_count: payloads.len(),
        group_ids,
    })?;

    let plotly_tag = if let Some(inline_path) = inline_plotly_path {
        let js = fs::read_to_string(inline_path)?.replace("</script", "<\\/script");
        format!("<script>{js}</script>")
    } else {
        format!("<script src=\"{}\"></script>", html_escape(plotly_src))
    };

    let options_html: String = payloads
        .iter()
        .map(|p| format!("<option value=\"{}\"></option>", html_escape(&p.gene)))
        .collect();

    let title_esc = html_escape(title);
    let def_esc = html_escape(&eff_default);
    let scripts_joined = group_scripts.join("");

    let doc = format!(
        "<!doctype html>
<html lang=\"en\">
<head>
<meta charset=\"utf-8\">
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
<title>{title_esc}</title>
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
body {{ margin: 0; background: var(--wash); color: var(--ink); font-family: -apple-system, BlinkMacSystemFont, \"Segoe UI\", sans-serif; font-size: 14px; }}
main {{ max-width: 1540px; margin: 0 auto; padding: 26px 34px 46px; }}
header {{ display: flex; align-items: flex-end; justify-content: space-between; gap: 24px; margin-bottom: 18px; }}
h1 {{ margin: 0; font-size: clamp(25px, 2.2vw, 34px); letter-spacing: -0.035em; line-height: 1.08; font-weight: 700; }}
.kicker {{ margin: 0 0 8px; color: var(--accent); font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px; letter-spacing: .13em; text-transform: uppercase; }}
.notice {{ margin: 9px 0 0; color: var(--muted); }}
.summary {{ color: var(--muted); font-size: 12px; text-align: right; white-space: nowrap; }}
.toolbar, .plot-panel, .calls-panel {{ background: var(--panel); border: 1px solid var(--line); border-radius: var(--radius); }}
.toolbar {{ display: flex; align-items: center; flex-wrap: wrap; gap: 14px 18px; padding: 13px 16px; margin-bottom: 16px; }}
.control {{ display: inline-flex; align-items: center; gap: 7px; color: var(--muted); }}
.control label {{ font-size: 12px; font-weight: 600; }}
input[type=\"search\"], select {{ border: 1px solid #bcc9d1; border-radius: 7px; color: var(--ink); background: #fff; padding: 8px 10px; font: inherit; min-height: 35px; }}
input[type=\"search\"] {{ width: 180px; }}
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
.disclaimer {{ margin-top: 18px; color: var(--muted); font-size: 12px; text-align: center; }}
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
      <p class=\"kicker\">Unofficial review aid &middot; pure-Rust synthetic demo</p>
      <h1>{title_esc}</h1>
      <p class=\"notice\">Synthetic data only. Coordinates and calls are illustrative, not clinical results. Independent open-source project &mdash; not an official Illumina product.</p>
    </div>
    <div class=\"summary\">{g_len} genes / {s_len} samples / lazy per-gene payloads</div>
  </header>
  <section class=\"toolbar\" aria-label=\"Plot controls\">
    <div class=\"control\"><label for=\"gene-search\">Gene</label><input id=\"gene-search\" type=\"search\" list=\"gene-list\" placeholder=\"Search genes\" autocomplete=\"off\"><datalist id=\"gene-list\">{options_html}</datalist></div>
    <div class=\"control\"><label for=\"metric\">Y metric</label><select id=\"metric\"><option value=\"normalized_ratio\">normalized_ratio</option><option value=\"tn_log2\">tn_log2</option><option value=\"segment_mean\">segment_mean</option></select></div>
    <div class=\"checks\">
      <label class=\"check\"><input id=\"toggle-reference\" type=\"checkbox\" checked> reference lines</label>
      <label class=\"check\"><input id=\"toggle-bands\" type=\"checkbox\" checked> PASS bands</label>
      <label class=\"check\"><input id=\"toggle-strand\" type=\"checkbox\" checked> strand reversal</label>
      <label class=\"check\"><input id=\"toggle-ntc\" type=\"checkbox\" checked> NTC styling</label>
    </div>
  </section>
  <section class=\"plot-panel\" aria-label=\"Exon-level copy-number plot\">
    <div id=\"plot\" aria-live=\"polite\"></div>
    <div id=\"status\" class=\"status\"></div>
  </section>
  <section class=\"calls-panel\" aria-labelledby=\"calls-heading\">
    <div class=\"calls-heading\"><h2 id=\"calls-heading\">Calls in <span id=\"current-gene\">{def_esc}</span></h2><span id=\"strand-chip\" class=\"strand-chip\"></span></div>
    <div id=\"calls-table\" class=\"table-wrap\"></div>
    <div id=\"definitions\" class=\"definitions\"></div>
  </section>
  <p class=\"disclaimer\">Unofficial independent review aid for DRAGEN CNV output formats &mdash; not affiliated with, endorsed by, or sponsored by Illumina, Inc., and not for diagnostic (IVD) use.</p>
</main>
{scripts_joined}
<script id=\"cnv-config\" type=\"application/json\">{config_json}</script>
{plotly_tag}
<script>{PLOT_RUNTIME_JS}</script>
</body>
</html>
",
        g_len = payloads.len(),
        s_len = parsed_samples.len(),
    );

    if let Some(parent) = out_html.parent() {
        fs::create_dir_all(parent)?;
    }
    let bytes = doc.as_bytes();
    fs::write(out_html, bytes)?;
    Ok((bytes.len(), largest))
}

#[inline]
fn s_unused_suppress(samples: &[SampleData]) -> usize {
    samples.iter().map(|s| s.result_dir.as_os_str().len()).sum()
}
