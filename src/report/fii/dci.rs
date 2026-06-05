use std::fmt::Write as _;
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::cli::ComputeDciArgs;
use crate::io;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DciConfig {
    pub low_threshold: f64,
    pub high_threshold: f64,
}

impl Default for DciConfig {
    fn default() -> Self {
        Self {
            low_threshold: 0.33,
            high_threshold: 0.66,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Sample {
    pub sample_label: String,
    pub order_rank: usize,
    pub nucleus: Option<f64>,
    pub splice: Option<f64>,
    pub proteostasis: Option<f64>,
    pub mito: Option<f64>,
    pub tme: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct DciRow {
    pub sample_label: String,
    pub order_rank: usize,
    pub dci: f64,
    pub dci_core: f64,
    pub dci_extended: f64,
}

pub fn run_compute_dci(args: &ComputeDciArgs) -> Result<(), String> {
    let cfg = read_config(args.config.as_deref())?;
    let raw = std::fs::read_to_string(&args.input)
        .map_err(|e| format!("READ failed for {}: {e}", args.input.display()))?;
    let root: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", args.input.display()))?;
    let mut samples = read_samples(&root)?;
    samples.sort_by(|a, b| {
        a.order_rank
            .cmp(&b.order_rank)
            .then(a.sample_label.cmp(&b.sample_label))
    });

    let rows = samples
        .iter()
        .map(|s| compute_row(s, &cfg))
        .collect::<Vec<_>>();

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;
    io::write_bytes_atomic(
        &args.out.join("decision_concordance.tsv"),
        render_tsv(&rows).as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for decision_concordance.tsv: {e}"))?;
    let json_value = json!({
        "samples": rows.iter().map(|r| {
            json!({
                "sample_label": r.sample_label,
                "order_rank": r.order_rank,
                "DCI": r.dci,
                "DCI_core": r.dci_core,
                "DCI_extended": r.dci_extended
            })
        }).collect::<Vec<_>>()
    });
    io::write_json_atomic(&args.out.join("decision_concordance.json"), &json_value)
        .map_err(|e| format!("WRITE failed for decision_concordance.json: {e}"))?;
    Ok(())
}

fn read_config(path: Option<&Path>) -> Result<DciConfig, String> {
    let mut cfg = DciConfig::default();
    let Some(path) = path else { return Ok(cfg) };
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for config {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for config {}: {e}", path.display()))?;
    if let Some(x) = v.get("low_threshold").and_then(Value::as_f64) {
        cfg.low_threshold = x.clamp(0.0, 1.0);
    }
    if let Some(x) = v.get("high_threshold").and_then(Value::as_f64) {
        cfg.high_threshold = x.clamp(0.0, 1.0);
    }
    if cfg.high_threshold < cfg.low_threshold {
        std::mem::swap(&mut cfg.high_threshold, &mut cfg.low_threshold);
    }
    Ok(cfg)
}

fn read_samples(root: &Value) -> Result<Vec<Sample>, String> {
    let arr = root
        .get("samples")
        .and_then(Value::as_array)
        .ok_or_else(|| "MISSING samples array".to_string())?;
    let mut out = Vec::new();
    for item in arr {
        let label = item
            .get("sample_label")
            .or_else(|| item.get("label"))
            .and_then(Value::as_str)
            .ok_or_else(|| "MISSING sample label".to_string())?
            .to_string();
        let order_rank =
            item.get("order_rank")
                .or_else(|| item.get("rank"))
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("MISSING order_rank for {label}"))? as usize;
        out.push(Sample {
            sample_label: label,
            order_rank,
            nucleus: pick(item, &["FII_nucleus", "fii_nucleus", "nucleus_fii"]),
            splice: pick(item, &["FII_splice", "fii_splice", "splice_fii"]),
            proteostasis: pick(
                item,
                &["FII_proteostasis", "fii_proteostasis", "proteostasis_fii"],
            ),
            mito: pick(item, &["FII_mito", "fii_mito", "mito_fii"]),
            tme: pick(item, &["FII_tme", "fii_tme", "tme_fii"]),
        });
    }
    Ok(out)
}

fn pick(v: &Value, keys: &[&str]) -> Option<f64> {
    for key in keys {
        let Some(x) = v.get(*key).and_then(Value::as_f64) else {
            continue;
        };
        if x.is_finite() {
            return Some(x.clamp(0.0, 1.0));
        }
    }
    None
}

fn vote(v: Option<f64>, cfg: &DciConfig) -> Option<i32> {
    let v = v?;
    if v < cfg.low_threshold {
        Some(0)
    } else if v < cfg.high_threshold {
        Some(1)
    } else {
        Some(2)
    }
}

pub fn compute_row(sample: &Sample, cfg: &DciConfig) -> DciRow {
    let mut all_votes = Vec::new();
    let mut core_votes = Vec::new();
    let nucleus = vote(sample.nucleus, cfg);
    let splice = vote(sample.splice, cfg);
    let proteostasis = vote(sample.proteostasis, cfg);
    let mito = vote(sample.mito, cfg);
    let tme = vote(sample.tme, cfg);
    for x in [nucleus, splice, proteostasis, mito, tme]
        .into_iter()
        .flatten()
    {
        all_votes.push(x);
    }
    for x in [nucleus, splice, proteostasis].into_iter().flatten() {
        core_votes.push(x);
    }
    let dci = pairwise_mean_agreement(&all_votes);
    let dci_core = pairwise_mean_agreement(&core_votes);
    DciRow {
        sample_label: sample.sample_label.clone(),
        order_rank: sample.order_rank,
        dci,
        dci_core,
        dci_extended: dci,
    }
}

fn pairwise_mean_agreement(votes: &[i32]) -> f64 {
    if votes.len() < 2 {
        return 0.0;
    }
    let mut sum = 0.0;
    let mut n = 0usize;
    for i in 0..votes.len() {
        for j in (i + 1)..votes.len() {
            let a = votes[i];
            let b = votes[j];
            let agr = 1.0 - (a - b).abs() as f64 / 2.0;
            sum += agr;
            n += 1;
        }
    }
    if n == 0 {
        0.0
    } else {
        (sum / n as f64).clamp(0.0, 1.0)
    }
}

fn render_tsv(rows: &[DciRow]) -> String {
    let mut out = String::with_capacity(rows.len() * 64 + 64);
    out.push_str("sample_label\torder_rank\tDCI\tDCI_core\tDCI_extended\n");
    for r in rows {
        let _ = writeln!(
            &mut out,
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}",
            r.sample_label, r.order_rank, r.dci, r.dci_core, r.dci_extended
        );
    }
    out
}
