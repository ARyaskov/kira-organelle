use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::cli::ComputeCaiArgs;
use crate::io;
use crate::util::tsv::TsvReader;

use super::read::resolve_inputs;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CaiConfig {
    skewness_weight: f64,
    tail_heaviness_weight: f64,
    tail_mass_weight: f64,
    hard_threshold: f64,
}

impl Default for CaiConfig {
    fn default() -> Self {
        Self {
            skewness_weight: 0.33,
            tail_heaviness_weight: 0.33,
            tail_mass_weight: 0.34,
            hard_threshold: 0.66,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
struct CaiRow {
    sample_label: String,
    order_rank: usize,
    cai: f64,
    skewness_component: f64,
    tail_heaviness_component: f64,
    tail_mass: f64,
}

pub fn run_compute_cai(args: &ComputeCaiArgs) -> Result<(), String> {
    if args.inputs.is_empty() && args.manifest.is_none() {
        return Err("INVALID compute-cai args: --inputs or --manifest required".to_string());
    }
    let mut cfg = read_config(args.config.as_deref())?;
    if let Some(spec) = args.cai_weights.as_deref() {
        apply_weight_spec(&mut cfg, spec)?;
    }
    normalize_weights(&mut cfg)?;

    let mut specs = resolve_inputs(&args.inputs, args.manifest.as_deref())?;
    specs.sort_by(|a, b| a.order.cmp(&b.order).then(a.label.cmp(&b.label)));

    let mut rows = Vec::with_capacity(specs.len());
    for spec in &specs {
        let values = read_fii_values(&spec.root.join("functional_irreversibility_index.tsv"))?;
        rows.push(compute_cai(
            &spec.label,
            spec.order,
            &values,
            cfg.hard_threshold,
            cfg.skewness_weight,
            cfg.tail_heaviness_weight,
            cfg.tail_mass_weight,
        ));
    }

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;
    io::write_bytes_atomic(
        &args.out.join("commitment_asymmetry.tsv"),
        render_tsv(&rows).as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for commitment_asymmetry.tsv: {e}"))?;

    let json_value = json!({
        "samples": rows.iter().map(|r| {
            json!({
                "sample_label": r.sample_label,
                "order_rank": r.order_rank,
                "CAI": r.cai,
                "components": {
                    "skewness": r.skewness_component,
                    "tail_heaviness": r.tail_heaviness_component,
                    "tail_mass": r.tail_mass
                }
            })
        }).collect::<Vec<_>>()
    });
    io::write_json_atomic(&args.out.join("commitment_asymmetry.json"), &json_value)
        .map_err(|e| format!("WRITE failed for commitment_asymmetry.json: {e}"))?;
    Ok(())
}

fn read_config(path: Option<&Path>) -> Result<CaiConfig, String> {
    let mut cfg = CaiConfig::default();
    let Some(path) = path else { return Ok(cfg) };
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for config {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for config {}: {e}", path.display()))?;
    if let Some(x) = v.get("skewness_weight").and_then(Value::as_f64) {
        cfg.skewness_weight = x;
    }
    if let Some(x) = v.get("tail_heaviness_weight").and_then(Value::as_f64) {
        cfg.tail_heaviness_weight = x;
    }
    if let Some(x) = v.get("tail_mass_weight").and_then(Value::as_f64) {
        cfg.tail_mass_weight = x;
    }
    if let Some(x) = v.get("hard_threshold").and_then(Value::as_f64) {
        cfg.hard_threshold = x.clamp(0.0, 1.0);
    }
    Ok(cfg)
}

fn apply_weight_spec(cfg: &mut CaiConfig, spec: &str) -> Result<(), String> {
    let mut seen = BTreeMap::<String, f64>::new();
    for chunk in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let (k, v) = chunk
            .split_once(':')
            .ok_or_else(|| format!("INVALID --cai-weights entry '{chunk}', expected key:value"))?;
        let key = k.trim().to_ascii_lowercase();
        let val = v
            .trim()
            .parse::<f64>()
            .map_err(|e| format!("INVALID --cai-weights value '{v}': {e}"))?;
        seen.insert(key, val);
    }
    if let Some(v) = seen.get("skewness").or_else(|| seen.get("skew")) {
        cfg.skewness_weight = *v;
    }
    if let Some(v) = seen.get("tail_heaviness").or_else(|| seen.get("tail")) {
        cfg.tail_heaviness_weight = *v;
    }
    if let Some(v) = seen.get("tail_mass").or_else(|| seen.get("mass")) {
        cfg.tail_mass_weight = *v;
    }
    Ok(())
}

fn normalize_weights(cfg: &mut CaiConfig) -> Result<(), String> {
    let sum = cfg.skewness_weight + cfg.tail_heaviness_weight + cfg.tail_mass_weight;
    if !sum.is_finite() || sum <= 0.0 {
        return Err("INVALID CAI weights: sum must be > 0".to_string());
    }
    cfg.skewness_weight /= sum;
    cfg.tail_heaviness_weight /= sum;
    cfg.tail_mass_weight /= sum;
    Ok(())
}

fn read_fii_values(path: &Path) -> Result<Vec<f64>, String> {
    if !path.is_file() {
        return Err(format!(
            "MISSING functional_irreversibility_index.tsv in {}",
            path.display()
        ));
    }
    let mut reader =
        TsvReader::open(path).map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let mut header = Vec::new();
    if !reader
        .read_record(&mut header)
        .map_err(|e| format!("READ header failed for {}: {e}", path.display()))?
    {
        return Err(format!("MISSING header in {}", path.display()));
    }
    let cols = (0..TsvReader::fields_len(&header))
        .map(|i| {
            reader
                .field(&header, i)
                .unwrap_or_default()
                .trim()
                .to_string()
        })
        .collect::<Vec<_>>();
    let fii_idx = cols
        .iter()
        .position(|c| c == "functional_irreversibility_index")
        .ok_or_else(|| {
            format!(
                "MISSING required column 'functional_irreversibility_index' in {}",
                path.display()
            )
        })?;
    let mut values = Vec::new();
    let mut fields = Vec::new();
    while reader
        .read_record(&mut fields)
        .map_err(|e| format!("READ row failed for {}: {e}", path.display()))?
    {
        let parsed = reader
            .field(&fields, fii_idx)
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .ok();
        if let Some(v) = parsed
            && v.is_finite()
        {
            values.push(v.clamp(0.0, 1.0));
        }
    }
    Ok(values)
}

fn compute_cai(
    sample_label: &str,
    order_rank: usize,
    values: &[f64],
    hard_threshold: f64,
    w_skew: f64,
    w_tail_heavy: f64,
    w_tail_mass: f64,
) -> CaiRow {
    if values.is_empty() {
        return CaiRow {
            sample_label: sample_label.to_string(),
            order_rank,
            cai: 0.0,
            skewness_component: 0.0,
            tail_heaviness_component: 0.0,
            tail_mass: 0.0,
        };
    }

    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let q25 = percentile(&sorted, 0.25);
    let q50 = percentile(&sorted, 0.50);
    let q75 = percentile(&sorted, 0.75);
    let q95 = percentile(&sorted, 0.95);
    let iqr = (q75 - q25).max(0.0);

    // Bowley skewness robustly captures right-tail asymmetry.
    let skew = if iqr > 1e-12 {
        ((q75 + q25 - 2.0 * q50) / iqr).max(0.0)
    } else {
        0.0
    };
    let skew_norm = skew.clamp(0.0, 1.0);

    // Tail-heaviness proxy from upper-spread / IQR.
    let upper_spread = (q95 - q75).max(0.0);
    let ratio = if iqr > 1e-12 {
        upper_spread / iqr
    } else if upper_spread > 0.0 {
        1.0
    } else {
        0.0
    };
    let tail_heavy_norm = (ratio / (ratio + 1.0)).clamp(0.0, 1.0);

    let tail_mass =
        sorted.iter().filter(|v| **v >= hard_threshold).count() as f64 / sorted.len() as f64;
    // Small non-zero resistant tail indicates subpopulation commitment.
    let tail_mass_norm = if tail_mass > 0.0 {
        (1.0 - tail_mass).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let cai = (w_skew * skew_norm + w_tail_heavy * tail_heavy_norm + w_tail_mass * tail_mass_norm)
        .clamp(0.0, 1.0);

    CaiRow {
        sample_label: sample_label.to_string(),
        order_rank,
        cai,
        skewness_component: skew_norm,
        tail_heaviness_component: tail_heavy_norm,
        tail_mass,
    }
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let q = q.clamp(0.0, 1.0);
    let n = sorted.len();
    let rank = if q <= 0.0 {
        1usize
    } else if q >= 1.0 {
        n
    } else {
        (q * n as f64).ceil() as usize
    };
    sorted[rank.saturating_sub(1).min(n - 1)]
}

fn render_tsv(rows: &[CaiRow]) -> String {
    let mut out = String::from(
        "sample_label\torder_rank\tCAI\tskewness_component\ttail_heaviness_component\ttail_mass\n",
    );
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\n",
            r.sample_label,
            r.order_rank,
            r.cai,
            r.skewness_component,
            r.tail_heaviness_component,
            r.tail_mass
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn linspace(n: usize, start: f64, end: f64) -> Vec<f64> {
        if n <= 1 {
            return vec![start];
        }
        (0..n)
            .map(|i| start + (end - start) * (i as f64 / (n as f64 - 1.0)))
            .collect()
    }

    #[test]
    fn synthetic_distributions_show_expected_ordering() {
        let uniform = linspace(1000, 0.2, 0.8);
        let mut skewed = vec![0.05; 900];
        skewed.extend(vec![0.95; 100]);
        let mut bimodal = vec![0.15; 500];
        bimodal.extend(vec![0.85; 500]);

        let a = compute_cai("u", 0, &uniform, 0.66, 0.33, 0.33, 0.34);
        let b = compute_cai("s", 1, &skewed, 0.66, 0.33, 0.33, 0.34);
        let c = compute_cai("b", 2, &bimodal, 0.66, 0.33, 0.33, 0.34);

        assert!(b.cai > a.cai);
        assert!(b.cai > c.cai || c.cai > a.cai);
    }

    #[test]
    fn sensitive_to_small_hard_tail() {
        let diffuse = vec![0.55; 1000];
        let mut small_tail = vec![0.55; 980];
        small_tail.extend(vec![0.95; 20]);

        let a = compute_cai("d", 0, &diffuse, 0.66, 0.33, 0.33, 0.34);
        let b = compute_cai("t", 1, &small_tail, 0.66, 0.33, 0.33, 0.34);
        assert!(b.cai > a.cai);
    }

    #[test]
    fn stable_under_count_scaling() {
        let base = vec![0.2, 0.2, 0.3, 0.4, 0.9];
        let mut scaled = Vec::new();
        for _ in 0..200 {
            scaled.extend(base.iter().copied());
        }
        let a = compute_cai("a", 0, &base, 0.66, 0.33, 0.33, 0.34);
        let b = compute_cai("b", 1, &scaled, 0.66, 0.33, 0.33, 0.34);
        assert!((a.cai - b.cai).abs() < 1e-9);
    }

    #[test]
    fn deterministic_across_runs() {
        let mut values = vec![0.1; 700];
        values.extend(vec![0.95; 40]);
        values.extend(vec![0.5; 260]);
        let a = compute_cai("x", 0, &values, 0.66, 0.33, 0.33, 0.34);
        let b = compute_cai("x", 0, &values, 0.66, 0.33, 0.33, 0.34);
        assert!((a.cai - b.cai).abs() < 1e-12);
        assert!((a.skewness_component - b.skewness_component).abs() < 1e-12);
        assert!((a.tail_heaviness_component - b.tail_heaviness_component).abs() < 1e-12);
        assert!((a.tail_mass - b.tail_mass).abs() < 1e-12);
    }
}
