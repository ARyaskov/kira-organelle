use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::cli::ComputeDciArgs;
use crate::io;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DciConfig {
    low_threshold: f64,
    high_threshold: f64,
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
struct Sample {
    sample_label: String,
    order_rank: usize,
    nucleus: Option<f64>,
    splice: Option<f64>,
    proteostasis: Option<f64>,
    mito: Option<f64>,
    tme: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct DciRow {
    sample_label: String,
    order_rank: usize,
    dci: f64,
    dci_core: f64,
    dci_extended: f64,
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

fn compute_row(sample: &Sample, cfg: &DciConfig) -> DciRow {
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
    let mut out = String::from("sample_label\torder_rank\tDCI\tDCI_core\tDCI_extended\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\n",
            r.sample_label, r.order_rank, r.dci, r.dci_core, r.dci_extended
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fully_concordant_is_high() {
        let cfg = DciConfig::default();
        let s = Sample {
            sample_label: "x".to_string(),
            order_rank: 0,
            nucleus: Some(0.9),
            splice: Some(0.91),
            proteostasis: Some(0.95),
            mito: Some(0.88),
            tme: Some(0.9),
        };
        let row = compute_row(&s, &cfg);
        assert!((row.dci - 1.0).abs() < 1e-12);
    }

    #[test]
    fn fully_conflicting_is_low() {
        let cfg = DciConfig::default();
        let s = Sample {
            sample_label: "x".to_string(),
            order_rank: 0,
            nucleus: Some(0.1),
            splice: Some(0.9),
            proteostasis: Some(0.1),
            mito: Some(0.9),
            tme: Some(0.1),
        };
        let row = compute_row(&s, &cfg);
        assert!(row.dci < 0.5);
    }

    #[test]
    fn mixed_partial_is_mid() {
        let cfg = DciConfig::default();
        let s = Sample {
            sample_label: "x".to_string(),
            order_rank: 0,
            nucleus: Some(0.2),
            splice: Some(0.5),
            proteostasis: Some(0.8),
            mito: Some(0.5),
            tme: Some(0.2),
        };
        let row = compute_row(&s, &cfg);
        assert!(row.dci > 0.2 && row.dci < 0.9);
    }

    #[test]
    fn missing_is_robust() {
        let cfg = DciConfig::default();
        let s = Sample {
            sample_label: "x".to_string(),
            order_rank: 0,
            nucleus: Some(0.2),
            splice: None,
            proteostasis: Some(0.8),
            mito: None,
            tme: Some(0.2),
        };
        let row = compute_row(&s, &cfg);
        assert!((0.0..=1.0).contains(&row.dci));
    }

    #[test]
    fn deterministic() {
        let cfg = DciConfig::default();
        let s = Sample {
            sample_label: "x".to_string(),
            order_rank: 0,
            nucleus: Some(0.2),
            splice: Some(0.5),
            proteostasis: Some(0.8),
            mito: Some(0.4),
            tme: Some(0.9),
        };
        let a = compute_row(&s, &cfg);
        let b = compute_row(&s, &cfg);
        assert!((a.dci - b.dci).abs() < 1e-12);
        assert!((a.dci_core - b.dci_core).abs() < 1e-12);
    }
}
