use std::collections::BTreeMap;
use std::path::Path;

use serde::Serialize;
use serde_json::{Value, json};

use crate::cli::ReportDecisionTimelineArgs;
use crate::io;
use crate::util::tsv::TsvReader;

#[derive(Debug, Clone, Serialize)]
struct TimelineRow {
    order_rank: usize,
    sample_label: String,
    decision_tier: String,
    confidence_score: f64,
    cai_context: Option<f64>,
    pri_context: Option<f64>,
    cocs_context: Option<f64>,
    dci_context: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct Metrics {
    flip_count: u64,
    volatility: f64,
    early_flip_index: f64,
    stability_score: f64,
    persistence: BTreeMap<String, u64>,
}

pub fn run_report_decision_timeline(args: &ReportDecisionTimelineArgs) -> Result<(), String> {
    let mut timeline = read_decision_rows(&args.input)?;
    timeline.sort_by(|a, b| {
        a.order_rank
            .cmp(&b.order_rank)
            .then(a.sample_label.cmp(&b.sample_label))
    });
    let metrics = compute_metrics(&timeline);

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;
    io::write_bytes_atomic(
        &args.out.join("decision_timeline.tsv"),
        render_timeline_tsv(&timeline).as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for decision_timeline.tsv: {e}"))?;
    io::write_bytes_atomic(
        &args.out.join("decision_stability.tsv"),
        render_stability_tsv(&metrics).as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for decision_stability.tsv: {e}"))?;
    let json_value = build_json(&timeline, &metrics);
    io::write_json_atomic(&args.out.join("decision_stability.json"), &json_value)
        .map_err(|e| format!("WRITE failed for decision_stability.json: {e}"))?;
    Ok(())
}

fn read_decision_rows(input_dir: &Path) -> Result<Vec<TimelineRow>, String> {
    let json_path = input_dir.join("sample_decision.json");
    if json_path.is_file() {
        let raw = std::fs::read_to_string(&json_path)
            .map_err(|e| format!("READ failed for {}: {e}", json_path.display()))?;
        let v: Value = serde_json::from_str(&raw)
            .map_err(|e| format!("PARSE failed for {}: {e}", json_path.display()))?;
        let arr = v
            .get("samples")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("MISSING samples array in {}", json_path.display()))?;
        let mut out = Vec::new();
        for item in arr {
            let sample_label = item
                .get("sample_label")
                .and_then(Value::as_str)
                .ok_or_else(|| "MISSING sample_label".to_string())?
                .to_string();
            let order_rank = item
                .get("order_rank")
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("MISSING order_rank for {}", sample_label))?
                as usize;
            let decision_tier = item
                .get("decision_tier")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("MISSING decision_tier for {}", sample_label))?
                .to_string();
            let confidence_score = item
                .get("confidence")
                .or_else(|| item.get("confidence_score"))
                .and_then(Value::as_f64)
                .unwrap_or(0.0)
                .clamp(0.0, 1.0);
            let cai_context = item
                .get("cai_context")
                .or_else(|| item.get("CAI"))
                .or_else(|| item.get("cai"))
                .and_then(Value::as_f64)
                .map(|v| v.clamp(0.0, 1.0));
            let pri_context = item
                .get("pri_context")
                .or_else(|| item.get("PRI"))
                .or_else(|| item.get("pri"))
                .and_then(Value::as_f64)
                .map(|v| v.clamp(0.0, 1.0));
            let cocs_context = item
                .get("cocs_context")
                .or_else(|| item.get("COCS_global"))
                .or_else(|| item.get("cocs_global"))
                .and_then(Value::as_f64)
                .map(|v| v.clamp(0.0, 1.0));
            let dci_context = item
                .get("dci_context")
                .or_else(|| item.get("DCI"))
                .or_else(|| item.get("dci"))
                .and_then(Value::as_f64)
                .map(|v| v.clamp(0.0, 1.0));
            out.push(TimelineRow {
                order_rank,
                sample_label,
                decision_tier,
                confidence_score,
                cai_context,
                pri_context,
                cocs_context,
                dci_context,
            });
        }
        return Ok(out);
    }

    let tsv_path = input_dir.join("sample_decision.tsv");
    if !tsv_path.is_file() {
        return Err(format!(
            "MISSING sample_decision.json/sample_decision.tsv in {}",
            input_dir.display()
        ));
    }
    let mut reader = TsvReader::open(&tsv_path)
        .map_err(|e| format!("READ failed for {}: {e}", tsv_path.display()))?;
    let mut header = Vec::new();
    if !reader
        .read_record(&mut header)
        .map_err(|e| format!("READ header failed for {}: {e}", tsv_path.display()))?
    {
        return Ok(Vec::new());
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
        .ok_or_else(|| format!("MISSING sample_label in {}", tsv_path.display()))?;
    let rank_idx = *idx
        .get("order_rank")
        .ok_or_else(|| format!("MISSING order_rank in {}", tsv_path.display()))?;
    let tier_idx = *idx
        .get("decision_tier")
        .ok_or_else(|| format!("MISSING decision_tier in {}", tsv_path.display()))?;
    let conf_idx = *idx
        .get("confidence_score")
        .ok_or_else(|| format!("MISSING confidence_score in {}", tsv_path.display()))?;
    let cai_idx = idx.get("cai_context").copied();
    let pri_idx = idx.get("pri_context").copied();
    let cocs_idx = idx.get("cocs_context").copied();
    let dci_idx = idx.get("dci_context").copied();

    let mut out = Vec::new();
    let mut fields = Vec::new();
    while reader
        .read_record(&mut fields)
        .map_err(|e| format!("READ row failed for {}: {e}", tsv_path.display()))?
    {
        let sample_label = reader
            .field(&fields, sample_idx)
            .unwrap_or_default()
            .trim()
            .to_string();
        let order_rank = reader
            .field(&fields, rank_idx)
            .unwrap_or_default()
            .trim()
            .parse::<usize>()
            .unwrap_or(0);
        let decision_tier = reader
            .field(&fields, tier_idx)
            .unwrap_or_default()
            .trim()
            .to_string();
        let confidence_score = reader
            .field(&fields, conf_idx)
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .unwrap_or(0.0)
            .clamp(0.0, 1.0);
        let cai_context = cai_idx
            .and_then(|i| reader.field(&fields, i))
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .ok()
            .map(|v| v.clamp(0.0, 1.0));
        let pri_context = pri_idx
            .and_then(|i| reader.field(&fields, i))
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .ok()
            .map(|v| v.clamp(0.0, 1.0));
        let cocs_context = cocs_idx
            .and_then(|i| reader.field(&fields, i))
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .ok()
            .map(|v| v.clamp(0.0, 1.0));
        let dci_context = dci_idx
            .and_then(|i| reader.field(&fields, i))
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .ok()
            .map(|v| v.clamp(0.0, 1.0));
        out.push(TimelineRow {
            order_rank,
            sample_label,
            decision_tier,
            confidence_score,
            cai_context,
            pri_context,
            cocs_context,
            dci_context,
        });
    }
    Ok(out)
}

fn compute_metrics(timeline: &[TimelineRow]) -> Metrics {
    let mut persistence = BTreeMap::new();
    for tier in [
        "STABLE_ADAPTIVE",
        "TRANSITION_RISK",
        "PRE_RESISTANT",
        "FIXED_RESISTANT",
    ] {
        persistence.insert(tier.to_string(), 0u64);
    }

    let n = timeline.len();
    if n == 0 {
        return Metrics {
            flip_count: 0,
            volatility: 0.0,
            early_flip_index: 1.0,
            stability_score: 1.0,
            persistence,
        };
    }

    // longest consecutive run per tier
    let mut i = 0usize;
    while i < n {
        let tier = timeline[i].decision_tier.clone();
        let mut j = i + 1;
        while j < n && timeline[j].decision_tier == tier {
            j += 1;
        }
        let len = (j - i) as u64;
        if let Some(entry) = persistence.get_mut(&tier)
            && len > *entry
        {
            *entry = len;
        }
        i = j;
    }

    let mut flip_count = 0u64;
    let mut weighted_jump_sum = 0.0;
    for i in 1..n {
        let a = tier_rank(&timeline[i - 1].decision_tier);
        let b = tier_rank(&timeline[i].decision_tier);
        if a != b {
            flip_count += 1;
        }
        let avg_conf = (timeline[i - 1].confidence_score + timeline[i].confidence_score) / 2.0;
        weighted_jump_sum += (a - b).abs() as f64 * avg_conf;
    }
    let denom = if n > 1 { (n - 1) as f64 * 3.0 } else { 1.0 };
    let stability_score = (1.0 - (weighted_jump_sum / denom)).clamp(0.0, 1.0);
    let volatility = if n > 1 {
        flip_count as f64 / (n - 1) as f64
    } else {
        0.0
    };

    let early_flip = timeline
        .iter()
        .position(|t| {
            matches!(
                t.decision_tier.as_str(),
                "PRE_RESISTANT" | "FIXED_RESISTANT"
            )
        })
        .map(|idx| idx as f64)
        .unwrap_or((n - 1) as f64);
    let early_flip_index = if n > 1 {
        early_flip / (n - 1) as f64
    } else {
        1.0
    };

    Metrics {
        flip_count,
        volatility,
        early_flip_index,
        stability_score,
        persistence,
    }
}

fn tier_rank(tier: &str) -> i32 {
    match tier {
        "STABLE_ADAPTIVE" => 0,
        "TRANSITION_RISK" => 1,
        "PRE_RESISTANT" => 2,
        "FIXED_RESISTANT" => 3,
        _ => 1,
    }
}

fn render_timeline_tsv(rows: &[TimelineRow]) -> String {
    let mut out = String::from(
        "order_rank\tsample_label\tdecision_tier\tconfidence_score\tcai_context\tpri_context\tcocs_context\tdci_context\n",
    );
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\n",
            r.order_rank,
            r.sample_label,
            r.decision_tier,
            r.confidence_score,
            r.cai_context.map(|v| format!("{v:.6}")).unwrap_or_default(),
            r.pri_context.map(|v| format!("{v:.6}")).unwrap_or_default(),
            r.cocs_context
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            r.dci_context.map(|v| format!("{v:.6}")).unwrap_or_default()
        ));
    }
    out
}

fn render_stability_tsv(m: &Metrics) -> String {
    let mut out = String::from("metric\tvalue\n");
    out.push_str(&format!("flip_count\t{}\n", m.flip_count));
    out.push_str(&format!("volatility\t{:.6}\n", m.volatility));
    out.push_str(&format!("early_flip_index\t{:.6}\n", m.early_flip_index));
    out.push_str(&format!("stability_score\t{:.6}\n", m.stability_score));
    for tier in [
        "STABLE_ADAPTIVE",
        "TRANSITION_RISK",
        "PRE_RESISTANT",
        "FIXED_RESISTANT",
    ] {
        out.push_str(&format!(
            "persistence_length_{}\t{}\n",
            tier,
            m.persistence.get(tier).copied().unwrap_or(0)
        ));
    }
    out
}

fn build_json(timeline: &[TimelineRow], metrics: &Metrics) -> Value {
    json!({
        "timeline": timeline.iter().map(|r| {
            json!({
                "sample": r.sample_label,
                "tier": r.decision_tier,
                "confidence": r.confidence_score,
                "cai_context": r.cai_context,
                "pri_context": r.pri_context,
                "cocs_context": r.cocs_context,
                "dci_context": r.dci_context
            })
        }).collect::<Vec<_>>(),
        "metrics": {
            "flip_count": metrics.flip_count,
            "volatility": metrics.volatility,
            "early_flip_index": metrics.early_flip_index,
            "stability_score": metrics.stability_score,
            "persistence": metrics.persistence
        }
    })
}
