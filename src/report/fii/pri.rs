use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::cli::ComputePriArgs;
use crate::io;

use super::read::{SampleInput, resolve_inputs};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PriConfig {
    pub w_nuclear: f64,
    pub w_splice: f64,
    pub w_translation: f64,
}

impl Default for PriConfig {
    fn default() -> Self {
        Self {
            w_nuclear: 0.4,
            w_splice: 0.3,
            w_translation: 0.3,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PriRow {
    pub sample_label: String,
    pub order_rank: usize,
    pub pri: f64,
    pub nuclear_plasticity_component: f64,
    pub splicing_integrity_component: f64,
    pub translational_selectivity_component: f64,
    pub low_confidence: bool,
}

pub fn run_compute_pri(args: &ComputePriArgs) -> Result<(), String> {
    if args.inputs.is_empty() && args.manifest.is_none() {
        return Err("INVALID compute-pri args: --inputs or --manifest required".to_string());
    }
    let mut cfg = read_config(args.config.as_deref())?;
    if let Some(spec) = args.pri_weights.as_deref() {
        apply_weight_spec(&mut cfg, spec)?;
    }
    normalize_weights(&mut cfg)?;

    let mut specs = resolve_inputs(&args.inputs, args.manifest.as_deref())?;
    specs.sort_by(|a, b| a.order.cmp(&b.order).then(a.label.cmp(&b.label)));

    let mut rows = Vec::with_capacity(specs.len());
    for spec in &specs {
        rows.push(compute_for_sample(spec, &cfg)?);
    }

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;
    io::write_bytes_atomic(
        &args.out.join("plasticity_reserve.tsv"),
        render_tsv(&rows).as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for plasticity_reserve.tsv: {e}"))?;
    let json_value = json!({
        "samples": rows.iter().map(|r| {
            json!({
                "sample_label": r.sample_label,
                "order_rank": r.order_rank,
                "PRI": r.pri,
                "components": {
                    "nuclear": r.nuclear_plasticity_component,
                    "splice": r.splicing_integrity_component,
                    "translation": r.translational_selectivity_component
                },
                "low_confidence": r.low_confidence
            })
        }).collect::<Vec<_>>()
    });
    io::write_json_atomic(&args.out.join("plasticity_reserve.json"), &json_value)
        .map_err(|e| format!("WRITE failed for plasticity_reserve.json: {e}"))?;
    Ok(())
}

fn compute_for_sample(spec: &SampleInput, cfg: &PriConfig) -> Result<PriRow, String> {
    let jsons = collect_candidate_jsons(&spec.root)?;

    let plastic = jsons.iter().find_map(|v| {
        find_numeric(
            v,
            &["fraction_PlasticAdaptive", "fraction_plastic_adaptive"],
        )
    });
    let committed = jsons
        .iter()
        .find_map(|v| find_numeric(v, &["fraction_CommittedState", "fraction_committed_state"]));
    let npc = plastic
        .map(|v| v.clamp(0.0, 1.0))
        .or_else(|| committed.map(|v| (1.0 - v).clamp(0.0, 1.0)));

    let sic = jsons
        .iter()
        .find_map(|v| find_numeric(v, &["low_splice_noise_fraction"]))
        .map(|v| v.clamp(0.0, 1.0));
    let tsc = jsons
        .iter()
        .find_map(|v| find_numeric(v, &["translation_selectivity_index"]))
        .map(|v| v.clamp(0.0, 1.0));

    let mut weighted_sum = 0.0;
    let mut weight_sum = 0.0;
    if let Some(v) = npc {
        weighted_sum += cfg.w_nuclear * v;
        weight_sum += cfg.w_nuclear;
    }
    if let Some(v) = sic {
        weighted_sum += cfg.w_splice * v;
        weight_sum += cfg.w_splice;
    }
    if let Some(v) = tsc {
        weighted_sum += cfg.w_translation * v;
        weight_sum += cfg.w_translation;
    }

    let pri = if weight_sum > 0.0 {
        (weighted_sum / weight_sum).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let low_confidence = weight_sum < 0.999_999;

    Ok(PriRow {
        sample_label: spec.label.clone(),
        order_rank: spec.order,
        pri,
        nuclear_plasticity_component: npc.unwrap_or(0.0),
        splicing_integrity_component: sic.unwrap_or(0.0),
        translational_selectivity_component: tsc.unwrap_or(0.0),
        low_confidence,
    })
}

fn collect_candidate_jsons(root: &Path) -> Result<Vec<Value>, String> {
    let mut paths = vec![root.join("summary.json"), root.join("state.json")];
    for tool in [
        "kira-nuclearqc",
        "kira-spliceqc",
        "kira-riboqc",
        "kira-translatomeqc",
    ] {
        paths.push(root.join(tool).join("summary.json"));
    }
    let mut out = Vec::new();
    for path in paths {
        if !path.is_file() {
            continue;
        }
        let raw = std::fs::read_to_string(&path)
            .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
        let v: Value = serde_json::from_str(&raw)
            .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
        out.push(v);
    }
    Ok(out)
}

fn find_numeric(root: &Value, keys: &[&str]) -> Option<f64> {
    let target = keys
        .iter()
        .map(|k| k.to_ascii_lowercase())
        .collect::<Vec<_>>();
    fn walk(v: &Value, target: &[String]) -> Option<f64> {
        match v {
            Value::Object(map) => {
                for (k, vv) in map {
                    let kl = k.to_ascii_lowercase();
                    if target.iter().any(|t| t == &kl)
                        && let Some(x) = vv.as_f64()
                        && x.is_finite()
                    {
                        return Some(x);
                    }
                }
                for vv in map.values() {
                    if let Some(x) = walk(vv, target) {
                        return Some(x);
                    }
                }
                None
            }
            Value::Array(arr) => {
                for vv in arr {
                    if let Some(x) = walk(vv, target) {
                        return Some(x);
                    }
                }
                None
            }
            _ => None,
        }
    }
    walk(root, &target)
}

fn read_config(path: Option<&Path>) -> Result<PriConfig, String> {
    let mut cfg = PriConfig::default();
    let Some(path) = path else { return Ok(cfg) };
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for config {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for config {}: {e}", path.display()))?;
    if let Some(x) = v.get("w_nuclear").and_then(Value::as_f64) {
        cfg.w_nuclear = x;
    }
    if let Some(x) = v.get("w_splice").and_then(Value::as_f64) {
        cfg.w_splice = x;
    }
    if let Some(x) = v.get("w_translation").and_then(Value::as_f64) {
        cfg.w_translation = x;
    }
    Ok(cfg)
}

fn apply_weight_spec(cfg: &mut PriConfig, spec: &str) -> Result<(), String> {
    let mut seen = BTreeMap::<String, f64>::new();
    for chunk in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let (k, v) = chunk
            .split_once(':')
            .ok_or_else(|| format!("INVALID --pri-weights entry '{chunk}', expected key:value"))?;
        let key = k.trim().to_ascii_lowercase();
        let val = v
            .trim()
            .parse::<f64>()
            .map_err(|e| format!("INVALID --pri-weights value '{v}': {e}"))?;
        seen.insert(key, val);
    }
    if let Some(v) = seen.get("nuclear").or_else(|| seen.get("n")) {
        cfg.w_nuclear = *v;
    }
    if let Some(v) = seen.get("splice").or_else(|| seen.get("s")) {
        cfg.w_splice = *v;
    }
    if let Some(v) = seen.get("translation").or_else(|| seen.get("t")) {
        cfg.w_translation = *v;
    }
    Ok(())
}

pub fn normalize_weights(cfg: &mut PriConfig) -> Result<(), String> {
    let sum = cfg.w_nuclear + cfg.w_splice + cfg.w_translation;
    if !sum.is_finite() || sum <= 0.0 {
        return Err("INVALID PRI weights: sum must be > 0".to_string());
    }
    cfg.w_nuclear /= sum;
    cfg.w_splice /= sum;
    cfg.w_translation /= sum;
    Ok(())
}

fn render_tsv(rows: &[PriRow]) -> String {
    let mut out = String::with_capacity(rows.len() * 96 + 160);
    out.push_str(
        "sample_label\torder_rank\tPRI\tnuclear_plasticity_component\tsplicing_integrity_component\ttranslational_selectivity_component\n",
    );
    for r in rows {
        let _ = writeln!(
            &mut out,
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            r.sample_label,
            r.order_rank,
            r.pri,
            r.nuclear_plasticity_component,
            r.splicing_integrity_component,
            r.translational_selectivity_component
        );
    }
    out
}
