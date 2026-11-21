mod engine;
mod synth;

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Result};
use engine::{build_showcase_html, V8_STRING_LIMIT};
use synth::{generate_synthetic_run, DEFAULT_SEED};

fn main() -> Result<()> {
    let mut keep_work = false;
    let mut emit_csv = false;
    let mut inline_plotly = false;
    let mut seed = DEFAULT_SEED;
    let mut out_html: Option<PathBuf> = None;

    let args: Vec<String> = env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--keep-work" => keep_work = true,
            "--emit-csv" => emit_csv = true,
            "--inline-plotly" => inline_plotly = true,
            "--seed" => {
                i += 1;
                seed = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(DEFAULT_SEED);
            }
            "-o" | "--output" => {
                i += 1;
                out_html = args.get(i).map(PathBuf::from);
            }
            "-h" | "--help" => {
                println!(
                    "Usage: cargo run --release --manifest-path demo/Cargo.toml -- [--keep-work] [--emit-csv] [--inline-plotly] [--seed N] [-o PATH]"
                );
                return Ok(());
            }
            other => bail!("Unknown flag: {other}"),
        }
        i += 1;
    }

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo_root = manifest_dir
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));

    let work_dir = manifest_dir.join("_work");
    let synth_dir = work_dir.join("synthetic_run");
    let csv_path = work_dir.join("cnv_plot_table.csv");
    let html_path = out_html.unwrap_or_else(|| repo_root.join("docs/demo/index.html"));
    let vendor_plotly = repo_root.join("docs/vendor/plotly-2.35.2.min.js");

    let t0 = Instant::now();
    if synth_dir.exists() {
        fs::remove_dir_all(&synth_dir)?;
    }
    let (n_samples, n_targets) = generate_synthetic_run(&synth_dir, seed)?;
    println!(
        "[synth] generated {n_samples} synthetic samples x {n_targets} capture targets (seed={seed})"
    );

    let mut sample_dirs: Vec<PathBuf> = fs::read_dir(&synth_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().map(|ft| ft.is_dir()).unwrap_or(false))
        .map(|e| e.path())
        .collect();
    sample_dirs.sort();

    let (bytes_written, largest_json) = build_showcase_html(
        &sample_dirs,
        &synth_dir.join("synthetic_panel_targets.tsv"),
        &html_path,
        if emit_csv { Some(csv_path.as_path()) } else { None },
        "CNV Plot Builder \u{2013} synthetic demo",
        "BRCA1",
        r"(?i)\bntc\b",
        "../vendor/plotly-2.35.2.min.js",
        if inline_plotly {
            Some(vendor_plotly.as_path())
        } else {
            None
        },
    )?;

    if !keep_work {
        let _ = fs::remove_dir_all(&synth_dir);
    }

    let pct = (largest_json as f64 / V8_STRING_LIMIT as f64) * 100.0;
    println!(
        "[build] largest embedded per-gene JSON string: {largest_json} chars / {V8_STRING_LIMIT} V8 limit ({pct:.6}% used)"
    );
    println!(
        "[done] wrote {} ({bytes_written} bytes) in {:.1} ms",
        html_path.display(),
        t0.elapsed().as_secs_f64() * 1000.0
    );
    Ok(())
}
