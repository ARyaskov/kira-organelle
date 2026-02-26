use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::cli::DecisionArgs;
use crate::io;
use crate::util::tsv::TsvReader;

const TOTAL_DEFINED_FLAGS: f64 = 4.0;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DecisionThresholds {
    stable_threshold: f64,
    resistant_threshold: f64,
    velocity_low_threshold: f64,
    velocity_mid_threshold: f64,
    velocity_high_threshold: f64,
    acceleration_mid_threshold: f64,
    acceleration_high_threshold: f64,
}

impl Default for DecisionThresholds {
    fn default() -> Self {
        Self {
            stable_threshold: 0.30,
            resistant_threshold: 0.66,
            velocity_low_threshold: 0.05,
            velocity_mid_threshold: 0.10,
            velocity_high_threshold: 0.15,
            acceleration_mid_threshold: 0.05,
            acceleration_high_threshold: 0.10,
        }
    }
}

#[derive(Debug, Clone)]
struct SampleInput {
    sample_label: String,
    order_rank: usize,
    fii_mean: Option<f64>,
    fii_median: Option<f64>,
    fii_velocity: Option<f64>,
    fii_acceleration: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct SampleDecision {
    sample_label: String,
    order_rank: usize,
    decision_tier: String,
    confidence: f64,
    cai_context: Option<f64>,
    pri_context: Option<f64>,
    cocs_context: Option<f64>,
    dci_context: Option<f64>,
    drivers: Vec<String>,
}

pub fn run_decision(args: &DecisionArgs) -> Result<(), String> {
    let thresholds = read_thresholds(args.config.as_deref())?;
    let samples = read_phase_samples(&args.input)?;
    let flags_path = args
        .input
        .parent()
        .map(|p| p.join("early_warning_flags.tsv"))
        .unwrap_or_else(|| Path::new(".").join("early_warning_flags.tsv"));
    let flags_map = if flags_path.is_file() {
        read_flags(&flags_path)?
    } else {
        BTreeMap::new()
    };
    let cai_path = args
        .input
        .parent()
        .map(|p| p.join("commitment_asymmetry.json"))
        .unwrap_or_else(|| Path::new(".").join("commitment_asymmetry.json"));
    let cai_map = if cai_path.is_file() {
        read_cai(&cai_path)?
    } else {
        BTreeMap::new()
    };
    let pri_path = args
        .input
        .parent()
        .map(|p| p.join("plasticity_reserve.json"))
        .unwrap_or_else(|| Path::new(".").join("plasticity_reserve.json"));
    let pri_map = if pri_path.is_file() {
        read_pri(&pri_path)?
    } else {
        BTreeMap::new()
    };
    let cocs_path = args
        .input
        .parent()
        .map(|p| p.join("cross_organelle_coupling.json"))
        .unwrap_or_else(|| Path::new(".").join("cross_organelle_coupling.json"));
    let cocs_map = if cocs_path.is_file() {
        read_cocs(&cocs_path)?
    } else {
        BTreeMap::new()
    };
    let dci_path = args
        .input
        .parent()
        .map(|p| p.join("decision_concordance.json"))
        .unwrap_or_else(|| Path::new(".").join("decision_concordance.json"));
    let dci_map = if dci_path.is_file() {
        read_dci(&dci_path)?
    } else {
        BTreeMap::new()
    };

    let mut decisions = Vec::with_capacity(samples.len());
    for s in &samples {
        let flags = flags_map.get(&s.sample_label).cloned().unwrap_or_default();
        let cai = cai_map.get(&s.sample_label).copied();
        let pri = pri_map.get(&s.sample_label).copied();
        let cocs = cocs_map.get(&s.sample_label).copied();
        let dci = dci_map.get(&s.sample_label).copied();
        decisions.push(assign_decision(s, &flags, &thresholds, cai, pri, cocs, dci));
    }

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;

    let tsv = render_tsv(&decisions);
    io::write_bytes_atomic(&args.out.join("sample_decision.tsv"), tsv.as_bytes())
        .map_err(|e| format!("WRITE failed for sample_decision.tsv: {e}"))?;

    let json_value = json!({
        "samples": decisions,
        "thresholds": thresholds
    });
    io::write_json_atomic(&args.out.join("sample_decision.json"), &json_value)
        .map_err(|e| format!("WRITE failed for sample_decision.json: {e}"))?;

    Ok(())
}

fn read_thresholds(path: Option<&Path>) -> Result<DecisionThresholds, String> {
    let mut t = DecisionThresholds::default();
    let Some(path) = path else { return Ok(t) };
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for config {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for config {}: {e}", path.display()))?;
    if let Some(x) = v.get("stable_threshold").and_then(Value::as_f64) {
        t.stable_threshold = x;
    }
    if let Some(x) = v.get("resistant_threshold").and_then(Value::as_f64) {
        t.resistant_threshold = x;
    }
    if let Some(x) = v.get("velocity_low_threshold").and_then(Value::as_f64) {
        t.velocity_low_threshold = x;
    }
    if let Some(x) = v.get("velocity_mid_threshold").and_then(Value::as_f64) {
        t.velocity_mid_threshold = x;
    }
    if let Some(x) = v.get("velocity_high_threshold").and_then(Value::as_f64) {
        t.velocity_high_threshold = x;
    }
    if let Some(x) = v.get("acceleration_mid_threshold").and_then(Value::as_f64) {
        t.acceleration_mid_threshold = x;
    }
    if let Some(x) = v.get("acceleration_high_threshold").and_then(Value::as_f64) {
        t.acceleration_high_threshold = x;
    }
    Ok(t)
}

fn read_phase_samples(path: &Path) -> Result<Vec<SampleInput>, String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let root: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
    let arr = root
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
        let order_rank =
            item.get("rank")
                .or_else(|| item.get("order_rank"))
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("MISSING rank for {}", sample_label))? as usize;
        let fii_mean = item.get("fii_mean").and_then(Value::as_f64);
        let fii_median = item.get("fii_median").and_then(Value::as_f64);
        let fii_velocity = item.get("velocity").and_then(Value::as_f64);
        let fii_acceleration = item.get("acceleration").and_then(Value::as_f64);
        out.push(SampleInput {
            sample_label,
            order_rank,
            fii_mean,
            fii_median,
            fii_velocity,
            fii_acceleration,
        });
    }
    out.sort_by(|a, b| {
        a.order_rank
            .cmp(&b.order_rank)
            .then(a.sample_label.cmp(&b.sample_label))
    });
    Ok(out)
}

fn read_flags(path: &Path) -> Result<BTreeMap<String, Vec<String>>, String> {
    let mut reader =
        TsvReader::open(path).map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let mut header = Vec::new();
    if !reader
        .read_record(&mut header)
        .map_err(|e| format!("READ header failed for {}: {e}", path.display()))?
    {
        return Ok(BTreeMap::new());
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
    let idx = cols
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i))
        .collect::<BTreeMap<_, _>>();
    let sample_idx = *idx
        .get("sample_label")
        .ok_or_else(|| format!("MISSING sample_label column in {}", path.display()))?;
    let flag_idx = *idx
        .get("flag")
        .ok_or_else(|| format!("MISSING flag column in {}", path.display()))?;
    let mut out = BTreeMap::<String, Vec<String>>::new();
    let mut fields = Vec::new();
    while reader
        .read_record(&mut fields)
        .map_err(|e| format!("READ row failed for {}: {e}", path.display()))?
    {
        let sample = reader
            .field(&fields, sample_idx)
            .unwrap_or_default()
            .trim()
            .to_string();
        let flag = reader
            .field(&fields, flag_idx)
            .unwrap_or_default()
            .trim()
            .to_string();
        if sample.is_empty() || flag.is_empty() {
            continue;
        }
        out.entry(sample).or_default().push(flag);
    }
    for flags in out.values_mut() {
        flags.sort();
        flags.dedup();
    }
    Ok(out)
}

fn assign_decision(
    sample: &SampleInput,
    flags: &[String],
    t: &DecisionThresholds,
    cai_context: Option<f64>,
    pri_context: Option<f64>,
    cocs_context: Option<f64>,
    dci_context: Option<f64>,
) -> SampleDecision {
    let fii_mean = sample.fii_mean.unwrap_or(0.0);
    let _fii_median = sample.fii_median.unwrap_or(0.0);
    let vel = sample.fii_velocity.unwrap_or(0.0);
    let acc = sample.fii_acceleration.unwrap_or(0.0);
    let mut drivers = BTreeSet::new();
    for f in flags {
        drivers.insert(f.clone());
    }
    if fii_mean >= t.resistant_threshold {
        drivers.insert("HIGH_FII".to_string());
    }
    if vel >= t.velocity_high_threshold {
        drivers.insert("HIGH_VELOCITY".to_string());
    }
    if acc >= t.acceleration_high_threshold {
        drivers.insert("HIGH_ACCELERATION".to_string());
    }

    let decision_tier = if fii_mean >= t.resistant_threshold {
        "FIXED_RESISTANT".to_string()
    } else if vel >= t.velocity_high_threshold
        && acc >= t.acceleration_high_threshold
        && !flags.is_empty()
    {
        "PRE_RESISTANT".to_string()
    } else if fii_mean < t.resistant_threshold
        && (vel >= t.velocity_mid_threshold || acc >= t.acceleration_mid_threshold)
        && !flags.is_empty()
    {
        "TRANSITION_RISK".to_string()
    } else if fii_mean < t.stable_threshold
        && vel.abs() < t.velocity_low_threshold
        && flags.is_empty()
    {
        "STABLE_ADAPTIVE".to_string()
    } else {
        // Missing/ambiguous input degradation path.
        "TRANSITION_RISK".to_string()
    };

    let confidence = compute_confidence(fii_mean, vel, acc, flags.len());

    SampleDecision {
        sample_label: sample.sample_label.clone(),
        order_rank: sample.order_rank,
        decision_tier,
        confidence,
        cai_context,
        pri_context,
        cocs_context,
        dci_context,
        drivers: drivers.into_iter().collect(),
    }
}

fn compute_confidence(fii_mean: f64, vel: f64, acc: f64, n_flags: usize) -> f64 {
    let norm_fii = clamp01(fii_mean);
    let norm_vel = clamp01(vel.abs() / 0.2);
    let norm_acc = clamp01(acc.max(0.0) / 0.2);
    let flag_fraction = clamp01(n_flags as f64 / TOTAL_DEFINED_FLAGS);
    clamp01((norm_fii + norm_vel + norm_acc + flag_fraction) / 4.0)
}

fn clamp01(v: f64) -> f64 {
    if !v.is_finite() {
        0.0
    } else {
        v.clamp(0.0, 1.0)
    }
}

fn render_tsv(rows: &[SampleDecision]) -> String {
    let mut out = String::from(
        "sample_label\torder_rank\tdecision_tier\tconfidence_score\tcai_context\tpri_context\tcocs_context\tdci_context\tdrivers\n",
    );
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\n",
            row.sample_label,
            row.order_rank,
            row.decision_tier,
            row.confidence,
            row.cai_context
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            row.pri_context
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            row.cocs_context
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            row.dci_context
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            row.drivers.join(",")
        ));
    }
    out
}

fn read_cai(path: &Path) -> Result<BTreeMap<String, f64>, String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
    let arr = v
        .get("samples")
        .and_then(Value::as_array)
        .ok_or_else(|| format!("MISSING samples array in {}", path.display()))?;
    let mut out = BTreeMap::new();
    for item in arr {
        let label = item
            .get("sample_label")
            .or_else(|| item.get("label"))
            .and_then(Value::as_str)
            .unwrap_or_default()
            .trim()
            .to_string();
        if label.is_empty() {
            continue;
        }
        let cai = item
            .get("CAI")
            .or_else(|| item.get("cai"))
            .and_then(Value::as_f64)
            .map(|v| v.clamp(0.0, 1.0));
        if let Some(v) = cai {
            out.insert(label, v);
        }
    }
    Ok(out)
}

fn read_pri(path: &Path) -> Result<BTreeMap<String, f64>, String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
    let arr = v
        .get("samples")
        .and_then(Value::as_array)
        .ok_or_else(|| format!("MISSING samples array in {}", path.display()))?;
    let mut out = BTreeMap::new();
    for item in arr {
        let label = item
            .get("sample_label")
            .or_else(|| item.get("label"))
            .and_then(Value::as_str)
            .unwrap_or_default()
            .trim()
            .to_string();
        if label.is_empty() {
            continue;
        }
        let pri = item
            .get("PRI")
            .or_else(|| item.get("pri"))
            .and_then(Value::as_f64)
            .map(|v| v.clamp(0.0, 1.0));
        if let Some(v) = pri {
            out.insert(label, v);
        }
    }
    Ok(out)
}

fn read_cocs(path: &Path) -> Result<BTreeMap<String, f64>, String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
    let arr = v
        .get("samples")
        .and_then(Value::as_array)
        .ok_or_else(|| format!("MISSING samples array in {}", path.display()))?;
    let mut out = BTreeMap::new();
    for item in arr {
        let label = item
            .get("sample_label")
            .or_else(|| item.get("label"))
            .and_then(Value::as_str)
            .unwrap_or_default()
            .trim()
            .to_string();
        if label.is_empty() {
            continue;
        }
        let cocs = item
            .get("COCS_global")
            .or_else(|| item.get("cocs_global"))
            .and_then(Value::as_f64)
            .map(|v| v.clamp(0.0, 1.0));
        if let Some(v) = cocs {
            out.insert(label, v);
        }
    }
    Ok(out)
}

fn read_dci(path: &Path) -> Result<BTreeMap<String, f64>, String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
    let arr = v
        .get("samples")
        .and_then(Value::as_array)
        .ok_or_else(|| format!("MISSING samples array in {}", path.display()))?;
    let mut out = BTreeMap::new();
    for item in arr {
        let label = item
            .get("sample_label")
            .or_else(|| item.get("label"))
            .and_then(Value::as_str)
            .unwrap_or_default()
            .trim()
            .to_string();
        if label.is_empty() {
            continue;
        }
        let dci = item
            .get("DCI")
            .or_else(|| item.get("dci"))
            .and_then(Value::as_f64)
            .map(|v| v.clamp(0.0, 1.0));
        if let Some(v) = dci {
            out.insert(label, v);
        }
    }
    Ok(out)
}
