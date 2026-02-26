use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::cli::ComputeIliArgs;
use crate::io;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct IliConfig {
    min_signal_threshold: f64,
    tie_epsilon: f64,
}

impl Default for IliConfig {
    fn default() -> Self {
        Self {
            min_signal_threshold: 0.02,
            tie_epsilon: 1e-6,
        }
    }
}

#[derive(Debug, Clone)]
struct SampleRow {
    sample_label: String,
    order_rank: usize,
    nucleus: Option<f64>,
    splice: Option<f64>,
    proteostasis: Option<f64>,
    mito: Option<f64>,
    tme: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct IliRow {
    order_rank: usize,
    sample_label: String,
    ili_category: String,
    ili_confidence: f64,
    leading_deltafii_value: f64,
}

pub fn run_compute_ili(args: &ComputeIliArgs) -> Result<(), String> {
    let cfg = read_config(args.config.as_deref())?;
    let mut samples = read_samples(&args.input)?;
    samples.sort_by(|a, b| {
        a.order_rank
            .cmp(&b.order_rank)
            .then(a.sample_label.cmp(&b.sample_label))
    });
    validate_order(&samples)?;

    let mut rows = Vec::new();
    for i in 1..samples.len() {
        rows.push(compute_transition(&samples[i - 1], &samples[i], &cfg));
    }

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;
    io::write_bytes_atomic(
        &args.out.join("irreversibility_localization.tsv"),
        render_tsv(&rows).as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for irreversibility_localization.tsv: {e}"))?;
    let json_value = json!({
        "samples": rows.iter().map(|r| {
            json!({
                "sample_label": r.sample_label,
                "order_rank": r.order_rank,
                "ili": r.ili_category,
                "confidence": r.ili_confidence,
                "leading_ΔFII_value": r.leading_deltafii_value
            })
        }).collect::<Vec<_>>()
    });
    io::write_json_atomic(
        &args.out.join("irreversibility_localization.json"),
        &json_value,
    )
    .map_err(|e| format!("WRITE failed for irreversibility_localization.json: {e}"))?;
    Ok(())
}

fn read_config(path: Option<&Path>) -> Result<IliConfig, String> {
    let mut cfg = IliConfig::default();
    let Some(path) = path else { return Ok(cfg) };
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for config {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for config {}: {e}", path.display()))?;
    if let Some(x) = v.get("min_signal_threshold").and_then(Value::as_f64) {
        cfg.min_signal_threshold = x;
    }
    if let Some(x) = v.get("tie_epsilon").and_then(Value::as_f64) {
        cfg.tie_epsilon = x;
    }
    Ok(cfg)
}

fn read_samples(path: &Path) -> Result<Vec<SampleRow>, String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
    let arr = v
        .get("samples")
        .and_then(Value::as_array)
        .ok_or_else(|| format!("MISSING samples array in {}", path.display()))?;
    let mut out = Vec::new();
    for item in arr {
        let sample_label = item
            .get("label")
            .or_else(|| item.get("sample_label"))
            .and_then(Value::as_str)
            .ok_or_else(|| "MISSING sample label".to_string())?
            .to_string();
        let order_rank = item
            .get("rank")
            .or_else(|| item.get("order_rank"))
            .and_then(Value::as_u64)
            .ok_or_else(|| format!("MISSING order_rank for {}", sample_label))?
            as usize;
        out.push(SampleRow {
            sample_label,
            order_rank,
            nucleus: pick_component(item, &["FII_nucleus", "fii_nucleus", "nucleus_fii"]),
            splice: pick_component(item, &["FII_splice", "fii_splice", "splice_fii"]),
            proteostasis: pick_component(
                item,
                &["FII_proteostasis", "fii_proteostasis", "proteostasis_fii"],
            ),
            mito: pick_component(item, &["FII_mito", "fii_mito", "mito_fii"]),
            tme: pick_component(item, &["FII_tme", "fii_tme", "tme_fii"]),
        });
    }
    Ok(out)
}

fn pick_component(v: &Value, keys: &[&str]) -> Option<f64> {
    for key in keys {
        let Some(x) = v.get(*key).and_then(Value::as_f64) else {
            continue;
        };
        if x.is_finite() {
            return Some(x);
        }
    }
    None
}

fn validate_order(samples: &[SampleRow]) -> Result<(), String> {
    if samples.len() < 2 {
        return Err("MISSING samples: at least 2 ordered samples required for ILI".to_string());
    }
    for i in 1..samples.len() {
        if samples[i].order_rank <= samples[i - 1].order_rank {
            return Err("INVALID order_rank sequence for ILI".to_string());
        }
    }
    Ok(())
}

fn compute_transition(prev: &SampleRow, curr: &SampleRow, cfg: &IliConfig) -> IliRow {
    let delta_order = (curr.order_rank as f64 - prev.order_rank as f64).max(1.0);
    let mut deltas = Vec::<(&str, f64)>::new();
    for (name, a, b) in [
        ("NUCLEUS_DRIVEN", prev.nucleus, curr.nucleus),
        ("SPLICE_DRIVEN", prev.splice, curr.splice),
        ("PROTEOSTASIS_DRIVEN", prev.proteostasis, curr.proteostasis),
        ("METABOLIC_DRIVEN", prev.mito, curr.mito),
        ("TME_DRIVEN", prev.tme, curr.tme),
    ] {
        let (Some(p), Some(c)) = (a, b) else { continue };
        let d = (c - p) / delta_order;
        if d > 0.0 && d.is_finite() {
            deltas.push((name, d));
        }
    }
    if deltas.is_empty() {
        return IliRow {
            order_rank: curr.order_rank,
            sample_label: curr.sample_label.clone(),
            ili_category: "UNRESOLVED".to_string(),
            ili_confidence: 0.0,
            leading_deltafii_value: 0.0,
        };
    }
    deltas.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(b.0)));
    let max_val = deltas[0].1;
    if max_val < cfg.min_signal_threshold {
        return IliRow {
            order_rank: curr.order_rank,
            sample_label: curr.sample_label.clone(),
            ili_category: "UNRESOLVED".to_string(),
            ili_confidence: 0.0,
            leading_deltafii_value: max_val,
        };
    }
    let n_tie = deltas
        .iter()
        .take_while(|(_, v)| (max_val - *v).abs() <= cfg.tie_epsilon)
        .count();
    if n_tie > 1 {
        return IliRow {
            order_rank: curr.order_rank,
            sample_label: curr.sample_label.clone(),
            ili_category: "UNRESOLVED".to_string(),
            ili_confidence: 0.0,
            leading_deltafii_value: max_val,
        };
    }
    let sum_pos = deltas.iter().map(|(_, v)| *v).sum::<f64>();
    let confidence = if sum_pos > 0.0 {
        max_val / sum_pos
    } else {
        0.0
    };
    IliRow {
        order_rank: curr.order_rank,
        sample_label: curr.sample_label.clone(),
        ili_category: deltas[0].0.to_string(),
        ili_confidence: confidence.clamp(0.0, 1.0),
        leading_deltafii_value: max_val,
    }
}

fn render_tsv(rows: &[IliRow]) -> String {
    let mut out = String::from(
        "order_rank\tsample_label\tILI_category\tILI_confidence\tleading_ΔFII_value\n",
    );
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{:.6}\n",
            r.order_rank,
            r.sample_label,
            r.ili_category,
            r.ili_confidence,
            r.leading_deltafii_value
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(
        label: &str,
        rank: usize,
        nucleus: Option<f64>,
        splice: Option<f64>,
        proteostasis: Option<f64>,
        mito: Option<f64>,
        tme: Option<f64>,
    ) -> SampleRow {
        SampleRow {
            sample_label: label.to_string(),
            order_rank: rank,
            nucleus,
            splice,
            proteostasis,
            mito,
            tme,
        }
    }

    #[test]
    fn ili_picks_known_leading_organelle() {
        let cfg = IliConfig::default();
        let prev = sample(
            "baseline",
            0,
            Some(0.20),
            Some(0.20),
            Some(0.20),
            Some(0.20),
            Some(0.20),
        );
        let curr = sample(
            "step1",
            1,
            Some(0.21),
            Some(0.48),
            Some(0.23),
            Some(0.22),
            Some(0.20),
        );
        let row = compute_transition(&prev, &curr, &cfg);
        assert_eq!(row.ili_category, "SPLICE_DRIVEN");
        assert!(row.ili_confidence > 0.5);
    }

    #[test]
    fn ili_tie_and_low_signal_resolve_to_unresolved() {
        let mut cfg = IliConfig::default();
        cfg.tie_epsilon = 1e-5;
        cfg.min_signal_threshold = 0.02;
        let prev = sample(
            "baseline",
            0,
            Some(0.10),
            Some(0.10),
            Some(0.10),
            Some(0.10),
            Some(0.10),
        );
        let tie_curr = sample(
            "tie",
            1,
            Some(0.30),
            Some(0.30),
            Some(0.10),
            Some(0.10),
            Some(0.10),
        );
        let tie_row = compute_transition(&prev, &tie_curr, &cfg);
        assert_eq!(tie_row.ili_category, "UNRESOLVED");

        let low_curr = sample(
            "low",
            1,
            Some(0.105),
            Some(0.103),
            Some(0.102),
            Some(0.101),
            Some(0.100),
        );
        let low_row = compute_transition(&prev, &low_curr, &cfg);
        assert_eq!(low_row.ili_category, "UNRESOLVED");
    }

    #[test]
    fn ili_skips_missing_components_and_still_computes() {
        let cfg = IliConfig::default();
        let prev = sample("baseline", 0, None, Some(0.15), None, Some(0.10), None);
        let curr = sample("step1", 1, None, Some(0.22), None, Some(0.12), None);
        let row = compute_transition(&prev, &curr, &cfg);
        assert_eq!(row.ili_category, "SPLICE_DRIVEN");
        assert!(row.leading_deltafii_value > 0.0);
    }

    #[test]
    fn ili_transition_deterministic() {
        let cfg = IliConfig::default();
        let prev = sample(
            "baseline",
            0,
            Some(0.20),
            Some(0.20),
            Some(0.20),
            Some(0.20),
            Some(0.20),
        );
        let curr = sample(
            "step1",
            2,
            Some(0.24),
            Some(0.42),
            Some(0.23),
            Some(0.28),
            Some(0.21),
        );
        let a = compute_transition(&prev, &curr, &cfg);
        let b = compute_transition(&prev, &curr, &cfg);
        assert_eq!(a.ili_category, b.ili_category);
        assert_eq!(a.order_rank, b.order_rank);
        assert!((a.ili_confidence - b.ili_confidence).abs() < 1e-12);
        assert!((a.leading_deltafii_value - b.leading_deltafii_value).abs() < 1e-12);
    }
}
