pub mod exec;
pub mod experiments;
pub mod plan;
pub mod timeseries;
pub mod tool;

use crate::cli::RunArgs;
use crate::fii::FiiWeights;
use crate::integration_export;
use crate::{AggregateOptions, run_aggregate, set_write_cells_json_enabled};
use serde_json::{Map, Value};
use std::ffi::OsString;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

pub fn run_pipeline(args: &RunArgs) -> Result<(), String> {
    let base_out = args
        .out
        .clone()
        .unwrap_or_else(|| args.input.join("kira-pipeline"));
    let experiments = experiments::discover_experiments(&args.input)?;
    let multiple = experiments.len() > 1;

    if args.dry_run {
        if multiple {
            println!("Detected {} experiments", experiments.len());
        }
        for spec in &experiments {
            let exp_out = if multiple {
                base_out.join(&spec.name)
            } else {
                base_out.clone()
            };
            let exp_input = if spec.uses_original_directory() {
                args.input.clone()
            } else {
                base_out.join(".kira-input").join(&spec.name)
            };
            let exp_args = RunArgs {
                input: exp_input,
                integration_manifest: args.integration_manifest.clone(),
                out: Some(exp_out),
                threads: args.threads,
                no_cache: args.no_cache,
                strict: args.strict,
                dry_run: true,
                fii_weights: args.fii_weights.clone(),
                mitoqc_write_profile_json: args.mitoqc_write_profile_json,
                write_cells_json: args.write_cells_json,
            };
            println!("=== experiment: {} ===", spec.name);
            let plan = plan::build_execution_plan(&exp_args);
            println!("{}", plan::render_dry_run_plan(&plan));
        }
        return Ok(());
    }

    let staging_root = base_out.join(".kira-input");
    let mut completed: Vec<(String, std::path::PathBuf)> = Vec::new();
    for spec in experiments {
        let prepared_input = experiments::materialize_experiment_input(&spec, &staging_root)?;
        let exp_out = if multiple {
            base_out.join(&spec.name)
        } else {
            base_out.clone()
        };
        let exp_args = RunArgs {
            input: prepared_input.clone(),
            integration_manifest: args.integration_manifest.clone(),
            out: Some(exp_out.clone()),
            threads: args.threads,
            no_cache: args.no_cache,
            strict: args.strict,
            dry_run: false,
            fii_weights: args.fii_weights.clone(),
            mitoqc_write_profile_json: args.mitoqc_write_profile_json,
            write_cells_json: args.write_cells_json,
        };
        run_single_pipeline(&exp_args, &spec.name)?;
        completed.push((spec.name.clone(), exp_out));
    }

    if multiple {
        let report_parent = base_out
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| base_out.clone());
        let report_dir = report_parent.join(derive_multi_report_dir_name(&args.input));
        timeseries::write_timeseries_report(&completed, &report_dir)?;
        integration_export::write_multi_integration_artifacts(
            &completed,
            &report_dir,
            args.integration_manifest.as_deref(),
        )?;
    }
    Ok(())
}

fn derive_multi_report_dir_name(input: &Path) -> String {
    let raw_name = input
        .file_name()
        .and_then(|s| s.to_str())
        .filter(|s| !s.is_empty())
        .unwrap_or("timeseries");
    let trimmed = raw_name
        .strip_suffix("_RAW")
        .or_else(|| raw_name.strip_suffix("-RAW"))
        .or_else(|| raw_name.strip_suffix("_raw"))
        .or_else(|| raw_name.strip_suffix("-raw"))
        .unwrap_or(raw_name);
    format!("{trimmed}-report")
}

fn run_single_pipeline(args: &RunArgs, experiment_name: &str) -> Result<(), String> {
    let total_started = Instant::now();
    let plan_started = Instant::now();
    let plan = plan::build_execution_plan(args);
    let plan_ms = plan_started.elapsed().as_millis() as u64;

    let exec_started = Instant::now();
    let issues = exec::execute_plan(&plan, args.strict, args.mitoqc_write_profile_json)?;
    let exec_ms = exec_started.elapsed().as_millis() as u64;
    if !issues.is_empty() {
        for issue in &issues {
            tracing::warn!(tool = ?issue.tool, code = %issue.code, message = %issue.message);
        }
    }

    let aggregate_started = Instant::now();
    let aggregate_opts = AggregateOptions {
        input: plan.out_root.clone(),
        input_b: None,
        out: Some(plan.organelle_out.clone()),
        strict: args.strict,
        json: true,
        validate_only: false,
        fii_weights: parse_fii_weights(args.fii_weights.as_deref())?,
        export_systems_model: None,
    };
    set_write_cells_json_enabled(args.write_cells_json);
    run_aggregate(&aggregate_opts)?;
    run_perturb_simulator(&plan, args)?;
    run_irreversibility(&plan)?;
    mirror_primary_aggregate_outputs(&plan.organelle_out, &plan.out_root)?;
    write_organelle_bundle_zip(&plan.organelle_out, &plan.out_root, experiment_name)?;

    let aggregate_ms = aggregate_started.elapsed().as_millis() as u64;
    let total_ms = total_started.elapsed().as_millis() as u64;
    tracing::info!(
        plan_ms,
        exec_ms,
        aggregate_ms,
        total_ms,
        peak_rss_kb = ?peak_rss_kb(),
        "run stage timing summary"
    );
    Ok(())
}

fn parse_fii_weights(raw: Option<&str>) -> Result<Option<FiiWeights>, String> {
    match raw {
        Some(v) => Ok(Some(FiiWeights::parse(v)?)),
        None => Ok(None),
    }
}

fn write_organelle_bundle_zip(
    organelle_out: &Path,
    out_root: &Path,
    experiment_name: &str,
) -> Result<(), String> {
    let zip_name = format!("{}.zip", experiments::sanitize_name(experiment_name));
    let zip_path = out_root.join(zip_name);
    if zip_path.exists() {
        std::fs::remove_file(&zip_path)
            .map_err(|e| format!("failed replacing {}: {e}", zip_path.display()))?;
    }

    let mut files = Vec::new();
    for name in [
        "interpretation.json",
        "llm_report.md",
        "pipeline_step.json",
        "report.html",
        "state.json",
    ] {
        let src = organelle_out.join(name);
        if src.is_file() {
            files.push(src);
        }
    }
    if files.is_empty() {
        return Ok(());
    }

    let mut cmd = Command::new("zip");
    cmd.arg("-q").arg("-j").arg(&zip_path);
    for file in files {
        cmd.arg(file);
    }
    let status = cmd
        .status()
        .map_err(|e| format!("failed running zip for {}: {e}", zip_path.display()))?;
    if !status.success() {
        return Err(format!(
            "zip command failed for {} with status {}",
            zip_path.display(),
            status
        ));
    }
    Ok(())
}

fn mirror_primary_aggregate_outputs(organelle_out: &Path, out_root: &Path) -> Result<(), String> {
    for name in ["state.json", "cells.json", "interpretation.json"] {
        let src = organelle_out.join(name);
        if !src.is_file() {
            continue;
        }
        let dst = out_root.join(name);
        let bytes =
            std::fs::read(&src).map_err(|e| format!("failed reading {}: {e}", src.display()))?;
        crate::io::write_bytes_atomic(&dst, &bytes)
            .map_err(|e| format!("failed writing {}: {e}", dst.display()))?;
    }
    Ok(())
}

fn peak_rss_kb() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in status.lines() {
            if let Some(rest) = line.strip_prefix("VmHWM:") {
                let value = rest.split_whitespace().next()?;
                return value.parse::<u64>().ok();
            }
        }
        None
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

fn run_irreversibility(plan: &plan::ExecutionPlan) -> Result<(), String> {
    let timeseries_path = plan
        .organelle_out
        .join("integration")
        .join("timeseries.tsv");
    if !timeseries_path.is_file() {
        return Err(format!(
            "cannot run kira-irreversibility: timeseries missing at {}",
            timeseries_path.display()
        ));
    }
    if !timeseries_has_dynamics(&timeseries_path)? {
        tracing::warn!(
            tool = "kira-irreversibility",
            timeseries = %timeseries_path.display(),
            "skipping irreversibility: no entity has >=2 timepoints"
        );
        return Ok(());
    }

    tracing::info!(
        tool = "kira-irreversibility",
        "resolved tool binary by name"
    );
    let mut cmd = Command::new("kira-irreversibility");
    enrich_command_path_with_exe_dir(&mut cmd);
    let weights_path = resolve_or_generate_irreversibility_weights_path(plan, &timeseries_path)?;
    let irreversibility_out = plan.out_root.join("kira-irreversibility");

    cmd.arg("run")
        .arg("--organelle-timeseries")
        .arg(&timeseries_path)
        .arg("--mode")
        .arg("cell")
        .arg("--weights")
        .arg(weights_path)
        .arg("--out")
        .arg(&irreversibility_out);

    let status = cmd.status().map_err(|e| {
        format!(
            "failed running kira-irreversibility for {}: {e}",
            plan.out_root.display()
        )
    })?;
    if !status.success() {
        return Err(format!(
            "kira-irreversibility failed for {} with status {}",
            plan.out_root.display(),
            status
        ));
    }
    Ok(())
}

fn timeseries_has_dynamics(path: &Path) -> Result<bool, String> {
    let file =
        std::fs::File::open(path).map_err(|e| format!("failed opening {}: {e}", path.display()))?;
    let mut reader = std::io::BufReader::new(file);
    let mut line = String::new();

    let header_bytes = reader
        .read_line(&mut line)
        .map_err(|e| format!("failed reading header from {}: {e}", path.display()))?;
    if header_bytes == 0 {
        return Ok(false);
    }

    let mut counts: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
    line.clear();
    loop {
        let bytes = reader
            .read_line(&mut line)
            .map_err(|e| format!("failed reading {}: {e}", path.display()))?;
        if bytes == 0 {
            break;
        }
        let row = line.trim_end_matches(&['\r', '\n'][..]);
        if row.is_empty() {
            line.clear();
            continue;
        }
        let mut parts = row.split('\t');
        let Some(entity_id) = parts.next() else {
            line.clear();
            continue;
        };
        let entry = counts.entry(entity_id.to_string()).or_insert(0);
        *entry += 1;
        if *entry >= 2 {
            return Ok(true);
        }
        line.clear();
    }
    Ok(false)
}

fn resolve_or_generate_irreversibility_weights_path(
    plan: &plan::ExecutionPlan,
    timeseries_path: &Path,
) -> Result<PathBuf, String> {
    if let Ok(path) = std::env::var("KIRA_IRREVERSIBILITY_WEIGHTS") {
        let pb = PathBuf::from(path);
        if pb.is_file() {
            return Ok(pb);
        }
        return Err(format!(
            "KIRA_IRREVERSIBILITY_WEIGHTS points to missing file: {}",
            pb.display()
        ));
    }

    let axis_names = read_timeseries_axis_names(timeseries_path)?;
    if axis_names.is_empty() {
        return Err(format!(
            "cannot generate irreversibility weights: no axis columns in {}",
            timeseries_path.display()
        ));
    }

    let weights_dir = plan.out_root.join("kira-irreversibility-input");
    std::fs::create_dir_all(&weights_dir).map_err(|e| {
        format!(
            "failed creating irreversibility input directory {}: {e}",
            weights_dir.display()
        )
    })?;
    let weights_path = weights_dir.join("weights.generated.toml");
    write_generated_irreversibility_weights(&weights_path, &axis_names)?;
    Ok(weights_path)
}

fn read_timeseries_axis_names(path: &Path) -> Result<Vec<String>, String> {
    let file =
        std::fs::File::open(path).map_err(|e| format!("failed opening {}: {e}", path.display()))?;
    let mut reader = std::io::BufReader::new(file);
    let mut header = String::new();
    let bytes = reader
        .read_line(&mut header)
        .map_err(|e| format!("failed reading header from {}: {e}", path.display()))?;
    if bytes == 0 {
        return Err(format!("empty timeseries file: {}", path.display()));
    }

    if header.ends_with('\n') {
        header.pop();
        if header.ends_with('\r') {
            header.pop();
        }
    }

    let cols: Vec<&str> = header.split('\t').collect();
    if cols.len() < 3 {
        return Ok(Vec::new());
    }
    Ok(cols.iter().skip(2).map(|s| (*s).to_string()).collect())
}

fn write_generated_irreversibility_weights(
    path: &Path,
    axis_names: &[String],
) -> Result<(), String> {
    let mut body = String::new();
    body.push_str("lambda_entropy = 0.20\n\n[weights]\n");
    for axis in axis_names {
        let escaped = axis.replace('\\', "\\\\").replace('"', "\\\"");
        body.push('"');
        body.push_str(&escaped);
        body.push_str("\" = 1.0\n");
    }
    std::fs::write(path, body).map_err(|e| format!("failed writing {}: {e}", path.display()))
}

fn run_perturb_simulator(plan: &plan::ExecutionPlan, args: &RunArgs) -> Result<(), String> {
    let state_path = if plan.out_root.join("state.json").is_file() {
        plan.out_root.join("state.json")
    } else {
        plan.organelle_out.join("state.json")
    };
    let baseline_axes = normalize_axes_for_perturb(extract_baseline_axes(&state_path)?);
    if baseline_axes.is_empty() {
        return Err(format!(
            "cannot run kira-perturb-simulator: no baseline axes in {}",
            state_path.display()
        ));
    }

    let perturb_input_dir = plan.out_root.join("kira-perturb-simulator-input");
    std::fs::create_dir_all(&perturb_input_dir).map_err(|e| {
        format!(
            "failed creating perturb input directory {}: {e}",
            perturb_input_dir.display()
        )
    })?;

    let summary_path = perturb_input_dir.join("organelle_summary.json");
    let expression_path = perturb_input_dir.join("expression.json");
    let perturbation_path =
        resolve_or_create_perturbation_config(&args.input, &perturb_input_dir, &baseline_axes)?;
    write_perturb_summary(&summary_path, &baseline_axes)?;
    write_perturb_expression(&expression_path, &baseline_axes)?;

    let perturb_out = plan.out_root.join("kira-perturb-simulator");
    ensure_clean_perturb_output_dir(&perturb_out)?;
    tracing::info!(
        tool = "kira-perturb-simulator",
        "resolved tool binary by name"
    );
    let mut cmd = Command::new("kira-perturb-simulator");
    enrich_command_path_with_exe_dir(&mut cmd);
    cmd.arg("run");
    cmd.arg("--organelle-summary")
        .arg(&summary_path)
        .arg("--expression")
        .arg(&expression_path)
        .arg("--perturbation")
        .arg(&perturbation_path)
        .arg("--mode")
        .arg("axis")
        .arg("--out")
        .arg(&perturb_out);

    let status = cmd.status().map_err(|e| {
        format!(
            "failed running kira-perturb-simulator for {}: {e}",
            plan.out_root.display()
        )
    })?;
    if !status.success() {
        return Err(format!(
            "kira-perturb-simulator failed for {} with status {}",
            plan.out_root.display(),
            status
        ));
    }

    Ok(())
}

fn ensure_clean_perturb_output_dir(path: &Path) -> Result<(), String> {
    if !path.exists() {
        return Ok(());
    }
    if !path.is_dir() {
        return Err(format!(
            "kira-perturb-simulator output path exists and is not a directory: {}",
            path.display()
        ));
    }
    std::fs::remove_dir_all(path).map_err(|e| {
        format!(
            "failed cleaning previous perturb output {}: {e}",
            path.display()
        )
    })
}

fn enrich_command_path_with_exe_dir(cmd: &mut Command) {
    let Ok(exe) = std::env::current_exe() else {
        return;
    };
    let Some(parent) = exe.parent() else {
        return;
    };

    let mut new_path = OsString::from(parent.as_os_str());
    if let Some(existing) = std::env::var_os("PATH") {
        new_path.push(":");
        new_path.push(existing);
    }
    cmd.env("PATH", new_path);
}

fn extract_baseline_axes(state_path: &Path) -> Result<Vec<(String, f64)>, String> {
    let bytes = std::fs::read(state_path)
        .map_err(|e| format!("failed reading {}: {e}", state_path.display()))?;
    let value: Value = serde_json::from_slice(&bytes)
        .map_err(|e| format!("failed parsing {}: {e}", state_path.display()))?;
    let organelles = value
        .get("organelle_states")
        .and_then(Value::as_array)
        .ok_or_else(|| {
            format!(
                "invalid state.json: missing organelle_states array in {}",
                state_path.display()
            )
        })?;

    let mut axes = Vec::new();
    for org in organelles {
        let Some(organelle) = org.get("organelle").and_then(Value::as_str) else {
            continue;
        };
        let Some(org_axes) = org.get("axes").and_then(Value::as_array) else {
            continue;
        };
        for axis in org_axes {
            let Some(name) = axis.get("name").and_then(Value::as_str) else {
                continue;
            };
            let Some(median) = axis.get("median").and_then(Value::as_f64) else {
                continue;
            };
            if !median.is_finite() {
                return Err(format!(
                    "invalid non-finite median for axis {}::{} in {}",
                    organelle,
                    name,
                    state_path.display()
                ));
            }
            axes.push((format!("{organelle}::{name}"), median));
        }
    }

    axes.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(axes)
}

fn write_perturb_summary(path: &Path, axes: &[(String, f64)]) -> Result<(), String> {
    let mut axes_map = Map::new();
    for (name, value) in axes {
        axes_map.insert(name.clone(), Value::from(*value));
    }

    let mut root = Map::new();
    root.insert("axes".to_string(), Value::Object(axes_map));
    crate::io::write_json_atomic(path, &Value::Object(root))
        .map_err(|e| format!("failed writing {}: {e}", path.display()))
}

fn write_perturb_expression(path: &Path, axes: &[(String, f64)]) -> Result<(), String> {
    let mut sample_values = Map::new();
    for (name, value) in axes {
        sample_values.insert(name.clone(), Value::from(*value));
    }

    let mut samples = Map::new();
    samples.insert("baseline".to_string(), Value::Object(sample_values));

    let mut root = Map::new();
    root.insert("samples".to_string(), Value::Object(samples));
    crate::io::write_json_atomic(path, &Value::Object(root))
        .map_err(|e| format!("failed writing {}: {e}", path.display()))
}

fn resolve_or_create_perturbation_config(
    input_root: &Path,
    perturb_input_dir: &Path,
    axes: &[(String, f64)],
) -> Result<std::path::PathBuf, String> {
    for candidate in [
        input_root.join("kira-perturbation.toml"),
        input_root.join("perturbation.toml"),
    ] {
        if candidate.is_file() {
            return Ok(candidate);
        }
    }

    let first_axis = axes
        .first()
        .map(|(name, _)| name.clone())
        .ok_or_else(|| "cannot create default perturbation: no axes available".to_string())?;
    let path = perturb_input_dir.join("perturbation.toml");
    let body = format!(
        "[[targets]]\n\
type = \"axis\"\n\
axis = \"{first_axis}\"\n\
delta = 0.1\n\
mode = \"additive\"\n\
min = 0.0\n\
max = 1.0\n\
\n\
[[coupling.rows]]\n\
from = \"{first_axis}\"\n\
[[coupling.rows.targets]]\n\
to = \"{first_axis}\"\n\
alpha = 0.25\n"
    );
    std::fs::write(&path, body).map_err(|e| {
        format!(
            "failed writing default perturbation config {}: {e}",
            path.display()
        )
    })?;
    Ok(path)
}

fn normalize_axes_for_perturb(axes: Vec<(String, f64)>) -> Vec<(String, f64)> {
    let mut normalized = Vec::with_capacity(axes.len());
    for (name, value) in axes {
        let clamped = value.clamp(0.0, 1.0);
        if (clamped - value).abs() > f64::EPSILON {
            tracing::warn!(
                axis = %name,
                original = value,
                clamped,
                "axis value out of [0,1] for perturb simulator input; clamping"
            );
        }
        normalized.push((name, clamped));
    }
    normalized
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::derive_multi_report_dir_name;

    #[test]
    fn report_dir_name_from_input_dataset() {
        assert_eq!(
            derive_multi_report_dir_name(Path::new("/tmp/GSE264586_RAW")),
            "GSE264586-report"
        );
        assert_eq!(
            derive_multi_report_dir_name(Path::new("/tmp/GSE127465-RAW")),
            "GSE127465-report"
        );
        assert_eq!(
            derive_multi_report_dir_name(Path::new("/tmp/experiment_set")),
            "experiment_set-report"
        );
    }
}
