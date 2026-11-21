//! Deterministic synthetic DRAGEN-shaped CNV batch generator.
//!
//! Produces illustrative synthetic files only — zero patient data, zero clinical
//! sample identifiers, and zero proprietary panel coordinates.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use flate2::write::GzEncoder;
use flate2::Compression;

pub const DEFAULT_SEED: u64 = 20260729;

#[derive(Clone, Copy)]
pub struct GeneSpec {
    pub name: &'static str,
    pub contig: &'static str,
    pub strand: &'static str,
    pub exon_count: u32,
}

#[derive(Clone)]
pub struct TargetSpec {
    pub contig: &'static str,
    pub start: i64,
    pub end: i64,
    pub name: String,
    pub gene: &'static str,
    pub exon_number: u32,
    pub strand: &'static str,
}

#[derive(Clone, Copy)]
struct PlantedEvent {
    gene: &'static str,
    sample: &'static str,
    start_exon: u32,
    end_exon: u32,
    ratio: f64,
    call_type: &'static str,
    copy_number: u32,
    filter_tag: &'static str,
    qual: f64,
}

const GENES: &[GeneSpec] = &[
    GeneSpec { name: "BRCA1", contig: "chr17", strand: "-", exon_count: 23 },
    GeneSpec { name: "BRCA2", contig: "chr13", strand: "+", exon_count: 27 },
    GeneSpec { name: "TP53", contig: "chr17", strand: "-", exon_count: 11 },
    GeneSpec { name: "MLH1", contig: "chr3", strand: "+", exon_count: 19 },
    GeneSpec { name: "MSH2", contig: "chr2", strand: "+", exon_count: 16 },
    GeneSpec { name: "PTEN", contig: "chr10", strand: "+", exon_count: 9 },
    GeneSpec { name: "EGFR", contig: "chr7", strand: "+", exon_count: 28 },
    GeneSpec { name: "DMD", contig: "chrX", strand: "-", exon_count: 20 },
    GeneSpec { name: "CFTR", contig: "chr7", strand: "+", exon_count: 27 },
    GeneSpec { name: "STK11", contig: "chr19", strand: "+", exon_count: 9 },
    GeneSpec { name: "PMS2", contig: "chr7", strand: "-", exon_count: 15 },
    GeneSpec { name: "CDKN2A", contig: "chr9", strand: "-", exon_count: 3 },
];

const EVENTS: &[PlantedEvent] = &[
    PlantedEvent { gene: "BRCA1", sample: "DEMO-RUN01-S03", start_exon: 1, end_exon: 23, ratio: 0.50, call_type: "DEL", copy_number: 1, filter_tag: "PASS", qual: 48.0 },
    PlantedEvent { gene: "EGFR", sample: "DEMO-RUN01-S07", start_exon: 18, end_exon: 21, ratio: 1.50, call_type: "DUP", copy_number: 3, filter_tag: "PASS", qual: 43.0 },
    PlantedEvent { gene: "PTEN", sample: "DEMO-RUN01-S05", start_exon: 6, end_exon: 6, ratio: 0.45, call_type: "DEL", copy_number: 1, filter_tag: "cnvLength", qual: 7.5 },
    PlantedEvent { gene: "DMD", sample: "DEMO-RUN01-S09", start_exon: 8, end_exon: 13, ratio: 0.50, call_type: "DEL", copy_number: 1, filter_tag: "PASS", qual: 39.0 },
    PlantedEvent { gene: "TP53", sample: "DEMO-RUN01-S02", start_exon: 4, end_exon: 6, ratio: 0.55, call_type: "DEL", copy_number: 1, filter_tag: "cnvQual", qual: 8.2 },
];

/// Portable SplitMix64 PRNG so synthetic runs are bit-reproducible across OS/architectures.
struct DetRng {
    state: u64,
    spare: Option<f64>,
}

impl DetRng {
    fn new(seed: u64) -> Self {
        Self { state: seed, spare: None }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn next_f64_open(&mut self) -> f64 {
        let u = (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64);
        u.clamp(1e-12, 1.0 - 1e-12)
    }

    fn randint(&mut self, lo: i64, hi: i64) -> i64 {
        let span = (hi - lo + 1) as u64;
        lo + (self.next_u64() % span) as i64
    }

    fn gauss(&mut self, mean: f64, std_dev: f64) -> f64 {
        if let Some(z1) = self.spare.take() {
            return mean + std_dev * z1;
        }
        let u1 = self.next_f64_open();
        let u2 = self.next_f64_open();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        self.spare = Some(r * theta.sin());
        mean + std_dev * (r * theta.cos())
    }
}

fn event_for(gene: &str, sample: &str) -> Option<&'static PlantedEvent> {
    EVENTS.iter().find(|e| e.gene == gene && e.sample == sample)
}

fn target_ratio(gene: &GeneSpec, exon_number: u32, sample: &str, rng: &mut DetRng) -> f64 {
    let is_ntc = sample.contains("NTC");
    let (mut ratio, ev) = if is_ntc {
        (rng.gauss(0.10, 0.04).max(0.02), None)
    } else {
        (rng.gauss(1.00, 0.06).max(0.02), event_for(gene.name, sample))
    };

    if let Some(e) = ev {
        if e.start_exon <= exon_number && exon_number <= e.end_exon {
            ratio = rng.gauss(e.ratio, e.ratio * 0.025).max(0.02);
        }
    }

    // Shared batch capture artifact at CFTR exons 10-11 across all samples.
    if gene.name == "CFTR" && (10..=11).contains(&exon_number) {
        ratio = (ratio * 0.65).max(0.02);
    }
    ratio
}

fn make_targets(rng: &mut DetRng) -> Vec<TargetSpec> {
    let mut cursor_by_contig: BTreeMap<&'static str, i64> = BTreeMap::new();
    for g in GENES {
        cursor_by_contig.entry(g.contig).or_insert(1_000_000);
    }

    let mut targets = Vec::new();
    for gene in GENES {
        let mut cursor = cursor_by_contig[gene.contig] + rng.randint(2_000, 6_000);
        for genomic_idx in 1..=gene.exon_count {
            let exon_number = if gene.strand == "+" {
                genomic_idx
            } else {
                gene.exon_count - genomic_idx + 1
            };
            cursor += rng.randint(450, 1_500);
            let t_count = rng.randint(1, 4);
            for within in 1..=t_count {
                let len = rng.randint(120, 200);
                let start = cursor;
                let end = start + len - 1;
                targets.push(TargetSpec {
                    contig: gene.contig,
                    start,
                    end,
                    name: format!("{}_exon_{:02}_target_{}", gene.name, exon_number, within),
                    gene: gene.name,
                    exon_number,
                    strand: gene.strand,
                });
                cursor = end + rng.randint(90, 260);
            }
            cursor += rng.randint(500, 1_400);
        }
        cursor_by_contig.insert(gene.contig, cursor + rng.randint(2_000, 6_000));
    }
    targets
}

pub fn sample_ids() -> Vec<String> {
    let mut ids: Vec<String> = (1..=12).map(|i| format!("DEMO-RUN01-S{i:02}")).collect();
    ids.push("DEMO-RUN01-NTC-A".to_string());
    ids.push("DEMO-RUN01-NTC-B".to_string());
    ids
}

pub fn generate_synthetic_run(out_dir: &Path, seed: u64) -> Result<(usize, usize)> {
    let mut rng = DetRng::new(seed);
    fs::create_dir_all(out_dir)?;

    let targets = make_targets(&mut rng);
    let panel_path = out_dir.join("synthetic_panel_targets.tsv");
    {
        let mut w = BufWriter::new(File::create(&panel_path)?);
        writeln!(w, "chrom\tstart\tend\tname\tGeneSymbol\texonNumber\tstrand")?;
        for t in &targets {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                t.contig, t.start, t.end, t.name, t.gene, t.exon_number, t.strand
            )?;
        }
        w.flush()?;
    }

    let mut by_gene: BTreeMap<&'static str, Vec<&TargetSpec>> = BTreeMap::new();
    for g in GENES {
        by_gene.insert(g.name, Vec::new());
    }
    for t in &targets {
        by_gene.get_mut(t.gene).unwrap().push(t);
    }

    let samples = sample_ids();
    for sample in &samples {
        let caller_dir = out_dir
            .join(sample)
            .join("vcf")
            .join("dragen")
            .join(sample)
            .join("cnv_exome_varcaller");
        fs::create_dir_all(&caller_dir)?;

        let mut ratios: BTreeMap<&str, f64> = BTreeMap::new();
        let mut seg_rows: Vec<(String, i64, i64, usize, f64)> = Vec::new();

        for gene in GENES {
            let g_targets = &by_gene[gene.name];
            for &t in g_targets {
                let r = target_ratio(gene, t.exon_number, sample, &mut rng);
                ratios.insert(&t.name, r);
            }

            let mut boundaries: BTreeSet<u32> = BTreeSet::from([1, gene.exon_count + 1]);
            if !sample.contains("NTC") {
                if let Some(ev) = event_for(gene.name, sample) {
                    boundaries.insert(ev.start_exon);
                    boundaries.insert(ev.end_exon + 1);
                }
            }
            if gene.name == "CFTR" {
                boundaries.insert(10);
                boundaries.insert(12);
            }
            let b_vec: Vec<u32> = boundaries.into_iter().collect();
            for pair in b_vec.windows(2) {
                let (left, right_excl) = (pair[0], pair[1]);
                let members: Vec<&&TargetSpec> = g_targets
                    .iter()
                    .filter(|t| left <= t.exon_number && t.exon_number < right_excl)
                    .collect();
                if members.is_empty() {
                    continue;
                }
                let mean: f64 =
                    members.iter().map(|t| ratios[t.name.as_str()]).sum::<f64>() / members.len() as f64;
                seg_rows.push((
                    gene.contig.to_string(),
                    members.first().unwrap().start,
                    members.last().unwrap().end,
                    members.len(),
                    mean,
                ));
            }
        }

        // 1. *.dragen.tn.tsv.gz
        {
            let tn_path = caller_dir.join(format!("{sample}.dragen.tn.tsv.gz"));
            let mut gz = GzEncoder::new(BufWriter::new(File::create(tn_path)?), Compression::default());
            writeln!(gz, "#contig\tstart\tend\ttarget_name\ttn_log2\timproper_pairs")?;
            for t in &targets {
                let r = ratios[t.name.as_str()];
                let imp = rng.randint(0, 3);
                writeln!(
                    gz,
                    "{}\t{}\t{}\t{}\t{:.6}\t{}",
                    t.contig,
                    t.start,
                    t.end,
                    t.name,
                    r.log2(),
                    imp
                )?;
            }
            gz.finish()?.flush()?;
        }

        // 2. *.dragen.seg
        {
            let seg_path = caller_dir.join(format!("{sample}.dragen.seg"));
            let mut w = BufWriter::new(File::create(seg_path)?);
            writeln!(w, "contig\tstart\tend\tnum_probes\tsegment_mean")?;
            for (c, s, e, n, m) in &seg_rows {
                writeln!(w, "{c}\t{s}\t{e}\t{n}\t{m:.4}")?;
            }
            w.flush()?;
        }

        // 3. *.dragen.cnv.gff3
        {
            let gff_path = caller_dir.join(format!("{sample}.dragen.cnv.gff3"));
            let mut w = BufWriter::new(File::create(gff_path)?);
            writeln!(w, "##gff-version 3")?;
            writeln!(w, "# Synthetic calls are illustrative and not clinical results.")?;
            for ev in EVENTS.iter().filter(|e| e.sample == sample) {
                let gene = GENES.iter().find(|g| g.name == ev.gene).unwrap();
                let members: Vec<&&TargetSpec> = by_gene[gene.name]
                    .iter()
                    .filter(|t| ev.start_exon <= t.exon_number && t.exon_number <= ev.end_exon)
                    .collect();
                let start = members.first().unwrap().start;
                let end = members.last().unwrap().end;
                writeln!(
                    w,
                    "{}\tsynthetic\t{}\t{}\t{}\t.\t{}\t.\tID={}_{}_{};CN={};Filter={};QUAL={:.1}",
                    gene.contig,
                    ev.call_type,
                    start,
                    end,
                    gene.strand,
                    sample,
                    gene.name,
                    ev.call_type,
                    ev.copy_number,
                    ev.filter_tag,
                    ev.qual
                )?;
            }
            w.flush()?;
        }
    }

    Ok((samples.len(), targets.len()))
}
