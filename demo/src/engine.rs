//! Single-pass multi-core Rust showcase engine (`rayon` + `MultiGzDecoder` + 100kb bins).
//!
//! Emits the v4 interactive HTML report (`plot-meta` + per-gene `<script id="gN">` lazy JSON
//! payloads + Light/Dark theme toggle) directly from synthetic DRAGEN folders.

use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use regex::Regex;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

pub const V8_STRING_LIMIT: usize = 536_870_888;
const BIN_SIZE: i64 = 100_000;
static PAGE_CSS: &str = include_str!("page.css");
static PAGE_JS: &str = include_str!("page.js");

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
    rows: Vec<TargetRow>,
}

#[derive(Serialize)]
struct PlotMetaJson<'a> {
    has_gene_data: bool,
    group_label: &'a str,
    samples: &'a [String],
    ntc_samples: &'a [String],
    groups: &'a [String],
    default_group: &'a str,
    dropped_unannotated_rows: usize,
}

#[derive(Serialize)]
struct SerializedCall {
    sample_id: String,
    status: &'static str,
    call_type: String,
    call_copy_number: String,
    call_filter: String,
    contig: String,
    call_start: String,
    call_end: String,
    exon_labels: Vec<String>,
    x_start: Option<usize>,
    x_end: Option<usize>,
    target_count: usize,
}

#[derive(Serialize)]
struct SerializedGroup {
    targets: Vec<[String; 5]>,
    s: Vec<usize>,
    x: Vec<usize>,
    v: Vec<Option<f64>>,
    m: Vec<Option<f64>>,
    c: Vec<i32>,
    calls: Vec<SerializedCall>,
    tv: Vec<f64>,
    tt: Vec<String>,
    bp: Vec<f64>,
    strand: String,
}

fn round4(x: f64) -> Option<f64> {
    if !x.is_finite() {
        return None;
    }
    format!("{x:.4}").parse::<f64>().ok()
}

fn build_axis_metadata(labels_by_rank: &[(usize, String)]) -> (Vec<f64>, Vec<String>, Vec<f64>) {
    if labels_by_rank.is_empty() {
        return (Vec::new(), Vec::new(), Vec::new());
    }
    let mut tick_vals = Vec::new();
    let mut tick_text = Vec::new();
    let mut boundaries = Vec::new();
    let mut run_label = labels_by_rank[0].1.clone();
    let mut run_start = labels_by_rank[0].0;
    let mut last = run_start;
    for (rank, label) in &labels_by_rank[1..] {
        if *label != run_label {
            tick_vals.push((run_start + last) as f64 / 2.0);
            tick_text.push(run_label);
            boundaries.push(last as f64 + 0.5);
            run_label = label.clone();
            run_start = *rank;
        }
        last = *rank;
    }
    tick_vals.push((run_start + last) as f64 / 2.0);
    tick_text.push(run_label);
    (tick_vals, tick_text, boundaries)
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

    Ok(SampleData { sample_id, rows })
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

    let sorted_samples: Vec<String> = parsed_samples
        .iter()
        .map(|s| s.sample_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let sample_index: HashMap<String, usize> = sorted_samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();
    let ntc_samples: Vec<String> = sorted_samples
        .iter()
        .filter(|s| ntc_re.is_match(s))
        .cloned()
        .collect();

    let mut group_names_set: BTreeSet<String> = BTreeSet::new();
    for r in &parsed_samples[0].rows {
        if let Some(ann) = panel.get(&(r.contig.clone(), r.start, r.end)) {
            group_names_set.insert(ann.gene.clone());
        }
    }
    let group_names: Vec<String> = group_names_set.into_iter().collect();

    let mut serialized_groups = Vec::with_capacity(group_names.len());
    let mut largest = 0usize;

    for gene in &group_names {
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
        target_keys.sort_by_key(|k| (panel[k].exon_number, k.1, k.2));
        target_keys.dedup();

        let strand = panel[&target_keys[0]].strand.clone();
        let mut targets = Vec::with_capacity(target_keys.len());
        let mut rank_by_key: FxHashMap<(String, i64, i64), usize> = FxHashMap::default();
        let mut labels_by_rank: Vec<(usize, String)> = Vec::with_capacity(target_keys.len());

        for (idx, k) in target_keys.iter().enumerate() {
            let rank = idx + 1;
            let ann = &panel[k];
            let exon_lbl = ann.exon_number.to_string();
            rank_by_key.insert(k.clone(), rank);
            labels_by_rank.push((rank, exon_lbl.clone()));
            targets.push([
                exon_lbl,
                k.0.clone(),
                k.1.to_string(),
                k.2.to_string(),
                ann.target_name.clone(),
            ]);
        }

        let (tv, tt, bp) = build_axis_metadata(&labels_by_rank);

        let mut s_arr = Vec::new();
        let mut x_arr = Vec::new();
        let mut v_arr = Vec::new();
        let mut m_arr = Vec::new();
        let mut c_arr = Vec::new();

        #[derive(Default)]
        struct CallWork {
            sample_id: String,
            call_type: String,
            cn: String,
            filter: String,
            contig: String,
            start: String,
            end: String,
            labels: BTreeSet<i64>,
            ranks: BTreeSet<usize>,
        }
        let mut call_map: FxHashMap<(String, String, i64, String, i64, i64), usize> =
            FxHashMap::default();
        let mut calls_work: Vec<CallWork> = Vec::new();

        for sp in &parsed_samples {
            let si = sample_index[&sp.sample_id];
            for r in &sp.rows {
                let k = (r.contig.clone(), r.start, r.end);
                let Some(&rank) = rank_by_key.get(&k) else {
                    continue;
                };
                let ann = &panel[&k];
                s_arr.push(si);
                x_arr.push(rank);
                v_arr.push(round4(r.tn_log2));
                m_arr.push(r.seg.as_ref().and_then(|seg| round4(seg.mean)));

                if let Some(bc) = &r.best_call {
                    let ck = (
                        sp.sample_id.clone(),
                        bc.call_type.clone(),
                        bc.copy_number,
                        bc.filter.clone(),
                        bc.start,
                        bc.end,
                    );
                    let ci = if let Some(&idx) = call_map.get(&ck) {
                        idx
                    } else {
                        let idx = calls_work.len();
                        call_map.insert(ck, idx);
                        calls_work.push(CallWork {
                            sample_id: sp.sample_id.clone(),
                            call_type: bc.call_type.clone(),
                            cn: bc.copy_number.to_string(),
                            filter: bc.filter.clone(),
                            contig: r.contig.clone(),
                            start: bc.start.to_string(),
                            end: bc.end.to_string(),
                            labels: BTreeSet::new(),
                            ranks: BTreeSet::new(),
                        });
                        idx
                    };
                    calls_work[ci].labels.insert(ann.exon_number);
                    calls_work[ci].ranks.insert(rank);
                    c_arr.push(ci as i32);
                } else {
                    c_arr.push(-1);
                }
            }
        }

        let calls: Vec<SerializedCall> = calls_work
            .into_iter()
            .map(|cw| {
                let x_start = cw.ranks.iter().min().copied();
                let x_end = cw.ranks.iter().max().copied();
                let target_count = cw.ranks.len();
                SerializedCall {
                    sample_id: cw.sample_id,
                    status: if cw.filter == "PASS" { "PASS" } else { "FILTERED" },
                    call_type: cw.call_type,
                    call_copy_number: cw.cn,
                    call_filter: cw.filter,
                    contig: cw.contig,
                    call_start: cw.start,
                    call_end: cw.end,
                    exon_labels: cw.labels.into_iter().map(|n| n.to_string()).collect(),
                    x_start,
                    x_end,
                    target_count,
                }
            })
            .collect();

        let sg = SerializedGroup {
            targets,
            s: s_arr,
            x: x_arr,
            v: v_arr,
            m: m_arr,
            c: c_arr,
            calls,
            tv,
            tt,
            bp,
            strand,
        };
        let payload = serde_json::to_string(&sg)?.replace("</script", "<\\/script");
        largest = largest.max(payload.len());
        serialized_groups.push(payload);
    }

    let eff_default = if group_names.iter().any(|g| g == default_group) {
        default_group.to_string()
    } else {
        group_names[0].clone()
    };

    let meta_json = serde_json::to_string(&PlotMetaJson {
        has_gene_data: true,
        group_label: "gene",
        samples: &sorted_samples,
        ntc_samples: &ntc_samples,
        groups: &group_names,
        default_group: &eff_default,
        dropped_unannotated_rows: 0,
    })?;

    let plotly_tag = if let Some(inline_path) = inline_plotly_path {
        let js = fs::read_to_string(inline_path)?.replace("</script", "<\\/script");
        format!("<script>\n{js}\n</script>")
    } else {
        format!("<script src=\"{}\"></script>", html_escape(plotly_src))
    };

    let escaped_title = html_escape(title);
    if let Some(parent) = out_html.parent() {
        fs::create_dir_all(parent)?;
    }

    let mut out = String::with_capacity(256 * 1024);
    out.push_str(&format!(
        "<!doctype html>
<html lang=\"en\" data-theme=\"light\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>{escaped_title}</title>
  {plotly_tag}
  <style>{PAGE_CSS}</style>
</head>
<body>
  <div class=\"brand-spectrum\" aria-hidden=\"true\"></div>
  <div class=\"wrap\">
    <div class=\"header\">
      <div class=\"header-brand\">
        <div class=\"header-emblem\" aria-hidden=\"true\">
          <svg width=\"26\" height=\"26\" viewBox=\"0 0 32 32\" fill=\"none\">
            <path d=\"M3 10C9 10 11 22 16 22C21 22 23 10 29 10\" stroke=\"#66ddff\" stroke-width=\"2.2\" stroke-linecap=\"round\"/>
            <path d=\"M3 22C9 22 11 10 16 10C21 10 23 22 29 22\" stroke=\"#00a3e0\" stroke-width=\"2.2\" stroke-linecap=\"round\" opacity=\"0.7\"/>
            <line x1=\"6\" y1=\"16\" x2=\"26\" y2=\"16\" stroke=\"#c7ecff\" stroke-width=\"1.3\" stroke-linecap=\"round\"/>
            <rect x=\"11\" y=\"10\" width=\"2.8\" height=\"6\" rx=\"1.4\" fill=\"#ff6a00\"/>
            <rect x=\"14.8\" y=\"8\" width=\"2.8\" height=\"8\" rx=\"1.4\" fill=\"#ff6a00\"/>
            <rect x=\"18.6\" y=\"16\" width=\"2.8\" height=\"6.5\" rx=\"1.4\" fill=\"#2e8bff\"/>
          </svg>
        </div>
        <div>
          <div class=\"eyebrow\">
            <span class=\"eyebrow-word\">UNOFFICIAL AID</span>
            <span class=\"eyebrow-dot\" aria-hidden=\"true\"></span>
            <span class=\"eyebrow-suite\">Exon-Level CNV Review &middot; Synthetic Demo</span>
          </div>
          <h1 class=\"title\">{escaped_title}</h1>
          <p class=\"subtitle\">100% synthetic batch (BRCA1, CFTR, EGFR, PTEN, DMD, TP53). Legend click toggles samples; double-click isolates one sample. Negative-strand genes read 5&#39;-&gt;3&#39;.</p>
        </div>
      </div>
      <div class=\"header-actions\">
        <div class=\"theme-toggle\" role=\"group\" aria-label=\"Color theme\">
          <button type=\"button\" class=\"theme-btn is-active\" id=\"theme-light-btn\">Light</button>
          <button type=\"button\" class=\"theme-btn\" id=\"theme-dark-btn\">Dark</button>
        </div>
      </div>
    </div>

    <div class=\"controls\">
      <div class=\"control\">
        <label id=\"group-label\" for=\"group-input\">Group</label>
        <input id=\"group-input\" list=\"group-list\" placeholder=\"Type a gene or contig\" autocomplete=\"off\">
        <datalist id=\"group-list\"></datalist>
      </div>
      <div class=\"control\">
        <label for=\"metric-select\">Metric</label>
        <select id=\"metric-select\">
          <option value=\"normalized_ratio\">Normalized Ratio</option>
          <option value=\"tn_log2\">TN Log2</option>
          <option value=\"segment_mean\">Segment Mean</option>
        </select>
      </div>
      <div class=\"control\">
        <label for=\"sort-select\">Sort</label>
        <select id=\"sort-select\">
          <option value=\"x\">Target Order</option>
          <option value=\"y_desc\">Y Desc</option>
          <option value=\"y_asc\">Y Asc</option>
        </select>
      </div>
      <button id=\"render-btn\" type=\"button\">Render</button>
      <button id=\"reset-btn\" type=\"button\">Reset</button>
      <div class=\"controls-options\">
        <label class=\"toggle\"><input type=\"checkbox\" id=\"toggle-refs\" checked> Reference lines (0.5 / 1.0 / 1.5)</label>
        <label class=\"toggle\"><input type=\"checkbox\" id=\"toggle-bands\" checked> Highlight PASS CNV regions</label>
        <label class=\"toggle\"><input type=\"checkbox\" id=\"toggle-strand\" checked> Reverse axis for - strand genes</label>
        <label class=\"toggle\"><input type=\"checkbox\" id=\"toggle-ntc\" checked> Style NTC samples</label>
      </div>
    </div>

    <div class=\"meta\">
      <div class=\"chip\" id=\"group-count-chip\"></div>
      <div class=\"chip\" id=\"sample-count-chip\"></div>
      <div class=\"chip\" id=\"ntc-chip\" hidden></div>
      <div class=\"chip\" id=\"mode-chip\"></div>
      <div class=\"chip\" id=\"strand-chip\" hidden></div>
      <div class=\"chip\" id=\"dropped-chip\"></div>
    </div>

    <div class=\"panel\">
      <div id=\"plot\"></div>
      <div class=\"hint\">PASS CNV calls in the selected group are drawn as colored bands behind the points. Toggle them above if they hide markers.</div>
    </div>

    <div class=\"calls\">
      <h2>CNV Calls In Current View</h2>
      <p class=\"calls-desc\" id=\"calls-desc\"></p>
      <div id=\"calls-content\"></div>
      <div id=\"calls-filter-notes\" class=\"filter-notes\"></div>
    </div>
    <p class=\"hint\" style=\"margin-top:18px;text-align:center;\">Unofficial independent review aid for DRAGEN CNV output files &mdash; not affiliated with or endorsed by Illumina, Inc., and not for diagnostic (IVD) use. All coordinates and values on this page are synthetic.</p>
  </div>
<script id=\"plot-meta\" type=\"application/json\">{meta_json}</script>
"
    ));

    let _unused = FxHashSet::<u8>::default();
    for (idx, payload) in serialized_groups.iter().enumerate() {
        out.push_str(&format!(
            "<script id=\"g{idx}\" type=\"application/json\">{payload}</script>\n"
        ));
    }
    out.push_str(&format!("  <script>{PAGE_JS}</script>\n</body>\n</html>\n"));

    let bytes = out.as_bytes();
    fs::write(out_html, bytes)?;
    Ok((bytes.len(), largest))
}
