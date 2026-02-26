mod assets;
mod binning;
mod cai;
mod cocs;
mod dci;
mod decision;
mod decision_timeline;
mod dynamics;
mod ili;
mod phase;
mod pri;
mod read;
mod render;

use std::path::Path;

use serde_json::{Value, json};

use crate::cli::ReportFiiArgs;
use crate::io;
use crate::version;

pub use read::SampleReportData;
use read::{read_sample, resolve_inputs};
use render::{ReportPayload, render_html};

pub fn run_report_fii(args: &ReportFiiArgs) -> Result<(), String> {
    if args.inputs.is_empty() && args.manifest.is_none() {
        return Err("INVALID report-fii args: --inputs or --manifest required".to_string());
    }

    let specs = resolve_inputs(&args.inputs, args.manifest.as_deref())?;
    if specs.is_empty() {
        return Err("MISSING samples for report-fii".to_string());
    }

    let mut samples = Vec::with_capacity(specs.len());
    for spec in &specs {
        let data = read_sample(spec)?;
        samples.push(data);
    }
    samples.sort_by(|a, b| a.order.cmp(&b.order).then(a.label.cmp(&b.label)));

    let title = args
        .title
        .clone()
        .unwrap_or_else(|| "Functional Irreversibility Landscape".to_string());
    let payload = ReportPayload {
        title,
        generated_at: io::deterministic_rfc3339_utc(),
        version: version::tool_version().to_string(),
        git_commit: option_env!("KIRA_GIT_COMMIT").map(|s| s.to_string()),
        samples: samples.clone(),
        decision_data: discover_optional_json(
            &specs,
            &args.out,
            &["sample_decision.json", "sample_decision.tsv"],
        )?,
        decision_stability_data: discover_optional_json(
            &specs,
            &args.out,
            &["decision_stability.json"],
        )?,
        phase_portrait_data: discover_optional_json(
            &specs,
            &args.out,
            &["phase_portrait_summary.json"],
        )?,
        ili_data: discover_optional_json(
            &specs,
            &args.out,
            &["irreversibility_localization.json"],
        )?,
        cai_data: discover_optional_json(&specs, &args.out, &["commitment_asymmetry.json"])?,
        pri_data: discover_optional_json(&specs, &args.out, &["plasticity_reserve.json"])?,
        cocs_data: discover_optional_json(&specs, &args.out, &["cross_organelle_coupling.json"])?,
        dci_data: discover_optional_json(&specs, &args.out, &["decision_concordance.json"])?,
    };
    let html = render_html(&payload)?;

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;
    let html_path = args.out.join("fii_landscape.html");
    io::write_bytes_atomic(&html_path, html.as_bytes())
        .map_err(|e| format!("WRITE failed for {}: {e}", html_path.display()))?;

    let index_path = args.out.join("fii_landscape_index.json");
    let index = build_index(&samples, &html_path);
    io::write_json_atomic(&index_path, &index)
        .map_err(|e| format!("WRITE failed for {}: {e}", index_path.display()))?;

    Ok(())
}

fn discover_optional_json(
    specs: &[read::SampleInput],
    out_dir: &std::path::Path,
    names: &[&str],
) -> Result<Value, String> {
    let mut candidates = Vec::new();
    candidates.push(out_dir.to_path_buf());
    if let Some(parent) = out_dir.parent() {
        candidates.push(parent.to_path_buf());
    }
    for spec in specs {
        candidates.push(spec.root.clone());
        if let Some(p) = spec.root.parent() {
            candidates.push(p.to_path_buf());
        }
    }

    for base in candidates {
        for name in names {
            let path = base.join(name);
            if !path.is_file() {
                continue;
            }
            if path.extension().and_then(|e| e.to_str()) == Some("json") {
                let raw = std::fs::read_to_string(&path)
                    .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
                let value: Value = serde_json::from_str(&raw)
                    .map_err(|e| format!("PARSE failed for {}: {e}", path.display()))?;
                return Ok(value);
            }
            // TSV decision fallback: convert to JSON shape with samples list.
            if *name == "sample_decision.tsv" {
                return read_decision_tsv_as_json(&path);
            }
        }
    }
    Ok(json!({}))
}

fn read_decision_tsv_as_json(path: &std::path::Path) -> Result<Value, String> {
    let mut reader = crate::util::tsv::TsvReader::open(path)
        .map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let mut header = Vec::new();
    if !reader
        .read_record(&mut header)
        .map_err(|e| format!("READ header failed for {}: {e}", path.display()))?
    {
        return Ok(json!({"samples":[]}));
    }
    let cols = (0..crate::util::tsv::TsvReader::fields_len(&header))
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
        .collect::<std::collections::BTreeMap<_, _>>();
    let sample_idx = *idx
        .get("sample_label")
        .ok_or_else(|| format!("MISSING sample_label in {}", path.display()))?;
    let rank_idx = *idx
        .get("order_rank")
        .ok_or_else(|| format!("MISSING order_rank in {}", path.display()))?;
    let tier_idx = *idx
        .get("decision_tier")
        .ok_or_else(|| format!("MISSING decision_tier in {}", path.display()))?;
    let conf_idx = *idx
        .get("confidence_score")
        .ok_or_else(|| format!("MISSING confidence_score in {}", path.display()))?;
    let drivers_idx = idx.get("drivers").copied();
    let cai_idx = idx.get("cai_context").copied();
    let pri_idx = idx.get("pri_context").copied();
    let cocs_idx = idx.get("cocs_context").copied();
    let dci_idx = idx.get("dci_context").copied();

    let mut samples = Vec::new();
    let mut fields = Vec::new();
    while reader
        .read_record(&mut fields)
        .map_err(|e| format!("READ row failed for {}: {e}", path.display()))?
    {
        let sample_label = reader
            .field(&fields, sample_idx)
            .unwrap_or_default()
            .trim()
            .to_string();
        if sample_label.is_empty() {
            continue;
        }
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
        let confidence = reader
            .field(&fields, conf_idx)
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .unwrap_or(0.0);
        let drivers = drivers_idx
            .and_then(|i| reader.field(&fields, i))
            .unwrap_or_default()
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(ToOwned::to_owned)
            .collect::<Vec<_>>();
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
        samples.push(json!({
            "sample_label": sample_label,
            "order_rank": order_rank,
            "decision_tier": decision_tier,
            "confidence": confidence,
            "cai_context": cai_context,
            "pri_context": pri_context,
            "cocs_context": cocs_context,
            "dci_context": dci_context,
            "drivers": drivers
        }));
    }
    Ok(json!({ "samples": samples }))
}

pub fn run_compute_state_dynamics(
    args: &crate::cli::ComputeStateDynamicsArgs,
) -> Result<(), String> {
    dynamics::run_compute_state_dynamics(args)
}

pub fn run_report_phase(args: &crate::cli::ReportPhaseArgs) -> Result<(), String> {
    phase::run_report_phase(args)
}

pub fn run_decision(args: &crate::cli::DecisionArgs) -> Result<(), String> {
    decision::run_decision(args)
}

pub fn run_report_decision_timeline(
    args: &crate::cli::ReportDecisionTimelineArgs,
) -> Result<(), String> {
    decision_timeline::run_report_decision_timeline(args)
}

pub fn run_compute_ili(args: &crate::cli::ComputeIliArgs) -> Result<(), String> {
    ili::run_compute_ili(args)
}

pub fn run_compute_cai(args: &crate::cli::ComputeCaiArgs) -> Result<(), String> {
    cai::run_compute_cai(args)
}

pub fn run_compute_pri(args: &crate::cli::ComputePriArgs) -> Result<(), String> {
    pri::run_compute_pri(args)
}

pub fn run_compute_cocs(args: &crate::cli::ComputeCocsArgs) -> Result<(), String> {
    cocs::run_compute_cocs(args)
}

pub fn run_compute_dci(args: &crate::cli::ComputeDciArgs) -> Result<(), String> {
    dci::run_compute_dci(args)
}

fn build_index(samples: &[SampleReportData], html_path: &Path) -> Value {
    json!({
        "artifact": "fii_landscape.html",
        "path": html_path.display().to_string(),
        "samples": samples.iter().map(|s| {
            json!({
                "id": s.id,
                "label": s.label,
                "path": s.path,
                "order": s.order,
                "n_cells": s.n_cells
            })
        }).collect::<Vec<_>>()
    })
}
