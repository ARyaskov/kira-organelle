use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::process::Stdio;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use flate2::read::GzDecoder;
use kira_mitoqc::io::mtx::load_mtx_dir;
use tracing::{error, info, warn};

use crate::contracts::types::{Issue, Severity};

use super::plan::{ExecutionPlan, ToolStep};
use super::tool::{ToolInvocationMode, resolve_tool_invocation};

#[derive(Debug, Clone)]
pub struct ToolResult {
    pub tool: String,
    pub success: bool,
    pub message: String,
}

pub fn prepare_layout(plan: &ExecutionPlan) -> Result<(), String> {
    std::fs::create_dir_all(&plan.out_root)
        .map_err(|e| format!("failed creating out root {}: {e}", plan.out_root.display()))?;
    std::fs::create_dir_all(&plan.organelle_out).map_err(|e| {
        format!(
            "failed creating organelle out {}: {e}",
            plan.organelle_out.display()
        )
    })?;

    for step in &plan.steps {
        std::fs::create_dir_all(&step.out_dir)
            .map_err(|e| format!("failed creating step out {}: {e}", step.out_dir.display()))?;
    }

    Ok(())
}

pub fn execute_plan(plan: &ExecutionPlan, strict: bool) -> Result<Vec<Issue>, String> {
    prepare_layout(plan)?;
    let mut issues = Vec::new();
    let mut shared_cache_bin_for_downstream = prepare_pipeline_shared_cache(plan)?;
    let mut expr_cache_for_microenvironment: Option<PathBuf> = None;

    for step in &plan.steps {
        let mut exec_step = step.clone();
        if tool_uses_pipeline_root_out(&exec_step.tool) {
            exec_step.out_dir = plan.out_root.clone();
        }
        if step.tool == "kira-microenvironment" {
            exec_step.cache_path = expr_cache_for_microenvironment.clone();
        } else {
            exec_step.cache_path = if supports_cache_arg(&step.tool) {
                shared_cache_bin_for_downstream.clone()
            } else {
                None
            };
        }

        let started = Instant::now();
        info!(tool = %exec_step.tool, mode = %mode_string(&exec_step.mode), "tool step started");
        let result = execute_step(&exec_step);
        let elapsed = started.elapsed();

        if result.success {
            info!(tool = %exec_step.tool, duration_ms = elapsed.as_millis() as u64, "tool step succeeded");
            if exec_step.tool == "kira-mitoqc" {
                let expr_cache = resolve_expression_cache_path(
                    &exec_step.input,
                    &plan.out_root,
                    &exec_step.out_dir,
                )?;
                if let Some(resolved) = shared_cache_bin_for_downstream.as_ref() {
                    info!(
                        tool = %exec_step.tool,
                        shared_cache = %resolved.display(),
                        "using shared cache for downstream pipeline tools"
                    );
                }
                info!(
                    tool = %exec_step.tool,
                    expression_cache = %expr_cache.display(),
                    "resolved expression cache for microenvironment stage"
                );
                expr_cache_for_microenvironment = Some(expr_cache);
            }
            continue;
        }

        warn!(
            tool = %step.tool,
            duration_ms = elapsed.as_millis() as u64,
            error = %result.message,
            "tool step failed"
        );
        issues.push(Issue {
            severity: if strict {
                Severity::Error
            } else {
                Severity::Warn
            },
            tool: Some(step.tool.clone()),
            code: "TOOL_EXECUTION_FAILED".to_string(),
            message: result.message.clone(),
            path: Some(step.out_dir.display().to_string()),
        });

        if strict {
            return Err(format!("strict mode: tool '{}' failed", step.tool));
        }
    }

    Ok(issues)
}

fn prepare_pipeline_shared_cache(plan: &ExecutionPlan) -> Result<Option<PathBuf>, String> {
    let first_cache = plan.steps.first().and_then(|s| s.cache_path.clone());
    let Some(target_cache) = first_cache else {
        return Ok(None);
    };

    if target_cache.is_file() {
        info!(
            shared_cache = %target_cache.display(),
            "pipeline shared cache already exists"
        );
        return Ok(Some(target_cache));
    }

    let input_root = if plan.input.is_file() {
        plan.input.parent().unwrap_or(&plan.input)
    } else {
        &plan.input
    };
    let detected_prefix = detect_dataset_prefix(input_root)?;
    let source_cache = input_root.join(kira_shared_sc_cache::resolve_shared_cache_filename(
        detected_prefix.as_deref(),
    ));
    if source_cache.is_file() {
        std::fs::copy(&source_cache, &target_cache).map_err(|e| {
            format!(
                "failed to copy shared cache {} to {}: {e}",
                source_cache.display(),
                target_cache.display()
            )
        })?;
        info!(
            source = %source_cache.display(),
            shared_cache = %target_cache.display(),
            "copied shared cache into pipeline out root"
        );
        return Ok(Some(target_cache));
    }

    let mtx = load_mtx_dir(&plan.input, None)
        .map_err(|e| format!("failed loading MTX to build shared cache: {e}"))?;
    let col_ptr: Vec<u64> = mtx
        .matrix
        .indptr()
        .raw_storage()
        .iter()
        .map(|v| *v as u64)
        .collect();
    let mut row_idx = Vec::with_capacity(mtx.matrix.nnz());
    let mut values_u32 = Vec::with_capacity(mtx.matrix.nnz());
    for col in mtx.matrix.outer_iterator() {
        for (row, value) in col.iter() {
            if !value.is_finite()
                || *value < 0.0
                || value.fract() != 0.0
                || *value > u32::MAX as f32
            {
                return Err(
                    "failed building shared cache: matrix contains non-integer count values"
                        .to_string(),
                );
            }
            row_idx.push(row as u32);
            values_u32.push(*value as u32);
        }
    }
    let input = kira_shared_sc_cache::SharedCacheWriteInput {
        genes: &mtx.features,
        barcodes: &mtx.barcodes,
        col_ptr: &col_ptr,
        row_idx: &row_idx,
        values_u32: &values_u32,
    };
    kira_shared_sc_cache::write_shared_cache(&target_cache, &input).map_err(|e| {
        format!(
            "failed writing shared cache {}: {e}",
            target_cache.display()
        )
    })?;
    info!(
        shared_cache = %target_cache.display(),
        "created pipeline shared cache before kira-mitoqc step"
    );
    Ok(Some(target_cache))
}

fn detect_dataset_prefix(dir: &Path) -> Result<Option<String>, String> {
    if !dir.is_dir() {
        return Ok(None);
    }
    let mut prefixes = Vec::new();
    let rd = std::fs::read_dir(dir)
        .map_err(|e| format!("failed reading directory {}: {e}", dir.display()))?;
    for entry in rd {
        let entry =
            entry.map_err(|e| format!("failed reading directory {}: {e}", dir.display()))?;
        if !entry
            .file_type()
            .map_err(|e| format!("failed reading file type in {}: {e}", dir.display()))?
            .is_file()
        {
            continue;
        }
        let name = entry.file_name();
        let Some(name) = name.to_str() else {
            continue;
        };
        for suffix in [
            "_matrix.mtx",
            "_matrix.mtx.gz",
            "_features.tsv",
            "_features.tsv.gz",
            "_genes.tsv",
            "_genes.tsv.gz",
            "_barcodes.tsv",
            "_barcodes.tsv.gz",
        ] {
            if let Some(prefix) = name.strip_suffix(suffix)
                && !prefix.is_empty()
            {
                prefixes.push(prefix.to_string());
            }
        }
    }
    prefixes.sort();
    prefixes.dedup();
    if prefixes.len() > 1 {
        return Err("multiple prefixed datasets found in input directory".to_string());
    }
    Ok(prefixes.into_iter().next())
}

fn tool_uses_pipeline_root_out(tool: &str) -> bool {
    matches!(
        tool,
        "kira-nuclearqc"
            | "kira-spliceqc"
            | "kira-proteoqc"
            | "kira-autolys"
            | "kira-secretion"
            | "kira-energetics"
            | "kira-microenvironment"
    )
}

fn resolve_shared_cache_path(
    input: &Path,
    out_root: &Path,
    mito_out_dir: &Path,
) -> Result<PathBuf, String> {
    let input_root = if input.is_file() {
        input.parent().unwrap_or(input)
    } else {
        input
    };
    let exact_input = input_root.join("kira-organelle.bin");
    if let Some(p) = resolve_cache_file_candidate(&exact_input) {
        return Ok(p);
    }

    let exact_out = out_root.join("kira-organelle.bin");
    if let Some(p) = resolve_cache_file_candidate(&exact_out) {
        return Ok(p);
    }

    let exact_mito_out = mito_out_dir.join("kira-organelle.bin");
    if let Some(p) = resolve_cache_file_candidate(&exact_mito_out) {
        return Ok(p);
    }

    // Some tools produce a direct expr.bin at the pipeline root.
    let root_expr = out_root.join("expr.bin");
    if root_expr.is_file() {
        return Ok(root_expr);
    }

    let mut candidates = collect_shared_cache_candidates(input_root)?;
    candidates.extend(collect_shared_cache_candidates(out_root)?);
    candidates.extend(collect_shared_cache_candidates(mito_out_dir)?);
    candidates.sort();

    match candidates.first() {
        Some(path) => Ok(path.clone()),
        None => Err(format!(
            "shared cache file not found; expected {}, {} or {} (or *.kira-organelle.bin file in those dirs)",
            exact_input.display(),
            exact_out.display(),
            exact_mito_out.display()
        )),
    }
}

fn resolve_shared_cache_bin_path(
    input: &Path,
    out_root: &Path,
    mito_out_dir: &Path,
) -> Result<PathBuf, String> {
    let input_root = if input.is_file() {
        input.parent().unwrap_or(input)
    } else {
        input
    };
    let exact_input = input_root.join("kira-organelle.bin");
    if exact_input.is_file() {
        return Ok(exact_input);
    }

    let exact_out = out_root.join("kira-organelle.bin");
    if exact_out.is_file() {
        return Ok(exact_out);
    }

    let exact_mito_out = mito_out_dir.join("kira-organelle.bin");
    if exact_mito_out.is_file() {
        return Ok(exact_mito_out);
    }

    let mut candidates = collect_shared_cache_bin_candidates(input_root)?;
    candidates.extend(collect_shared_cache_bin_candidates(out_root)?);
    candidates.extend(collect_shared_cache_bin_candidates(mito_out_dir)?);
    candidates.sort();

    match candidates.first() {
        Some(path) => Ok(path.clone()),
        None => Err(format!(
            "shared cache file not found; expected {}, {} or {} (or *.kira-organelle.bin file in those dirs)",
            exact_input.display(),
            exact_out.display(),
            exact_mito_out.display()
        )),
    }
}

fn collect_shared_cache_bin_candidates(dir: &Path) -> Result<Vec<PathBuf>, String> {
    if dir.is_file() {
        return Ok(Vec::new());
    }
    let rd = std::fs::read_dir(dir)
        .map_err(|e| format!("failed reading directory {}: {e}", dir.display()))?;
    let mut candidates = Vec::new();
    for entry in rd {
        let entry =
            entry.map_err(|e| format!("failed reading directory {}: {e}", dir.display()))?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name == "kira-organelle.bin" || name.ends_with(".kira-organelle.bin") {
            candidates.push(path);
        }
    }
    Ok(candidates)
}

fn resolve_expression_cache_path(
    input: &Path,
    out_root: &Path,
    mito_out_dir: &Path,
) -> Result<PathBuf, String> {
    resolve_shared_cache_path(input, out_root, mito_out_dir)
}

fn refresh_mito_shared_cache(step: &ToolStep, out_root: &Path) -> Result<(), String> {
    let mut retry_step = step.clone();
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_millis();
    let isolated_cache = out_root.join(format!(".kira-mitoqc-refresh-cache.{ts}"));
    retry_step.cache_path = Some(isolated_cache);
    let result = execute_step(&retry_step);
    if result.success {
        Ok(())
    } else {
        Err(format!(
            "failed to refresh kira-mitoqc shared cache: {}",
            result.message
        ))
    }
}

fn collect_shared_cache_candidates(dir: &Path) -> Result<Vec<PathBuf>, String> {
    if dir.is_file() {
        return Ok(Vec::new());
    }
    let rd = std::fs::read_dir(dir)
        .map_err(|e| format!("failed reading directory {}: {e}", dir.display()))?;
    let mut candidates = Vec::new();
    for entry in rd {
        let entry =
            entry.map_err(|e| format!("failed reading directory {}: {e}", dir.display()))?;
        let path = entry.path();
        if !path.is_file() && !path.is_dir() {
            continue;
        }
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name.ends_with(".kira-organelle.bin") {
            if let Some(cache_file) = resolve_cache_file_candidate(&path) {
                candidates.push(cache_file);
            }
            continue;
        }
        if name == "expr.bin" && path.is_file() {
            candidates.push(path);
        }
    }
    Ok(candidates)
}

fn resolve_cache_file_candidate(path: &Path) -> Option<PathBuf> {
    if path.is_file() {
        return Some(path.to_path_buf());
    }
    if path.is_dir() {
        let expr = path.join("expr.bin");
        if expr.is_file() {
            return Some(expr);
        }
    }
    None
}

pub fn execute_step(step: &ToolStep) -> ToolResult {
    match &step.mode {
        ToolInvocationMode::Library => ToolResult {
            tool: step.tool.clone(),
            success: false,
            message: "library integration not implemented for this build".to_string(),
        },
        ToolInvocationMode::Unavailable => ToolResult {
            tool: step.tool.clone(),
            success: false,
            message: "tool invocation unavailable (no library, no binary)".to_string(),
        },
        ToolInvocationMode::Binary(path) => execute_binary_step(step, path),
    }
}

fn execute_binary_step(step: &ToolStep, path: &std::path::Path) -> ToolResult {
    let mut cmd = Command::new(path);
    if let Some(dir) = path.parent() {
        // Some tools load assets via relative paths, so run from the binary directory.
        cmd.current_dir(dir);
    }
    let args = match build_tool_args(step, path) {
        Ok(v) => v,
        Err(e) => {
            return ToolResult {
                tool: step.tool.clone(),
                success: false,
                message: e,
            };
        }
    };
    for arg in args {
        cmd.arg(arg);
    }
    cmd.stdout(Stdio::inherit());
    cmd.stderr(Stdio::inherit());

    match cmd.status() {
        Ok(status) => {
            if status.success() {
                ToolResult {
                    tool: step.tool.clone(),
                    success: true,
                    message: "ok".to_string(),
                }
            } else {
                ToolResult {
                    tool: step.tool.clone(),
                    success: false,
                    message: format!("exit {status}"),
                }
            }
        }
        Err(e) => {
            error!(tool = %step.tool, error = %e, "failed spawning binary");
            ToolResult {
                tool: step.tool.clone(),
                success: false,
                message: format!("spawn error: {e}"),
            }
        }
    }
}

fn build_tool_args(step: &ToolStep, bin_path: &Path) -> Result<Vec<String>, String> {
    let cwd = std::env::current_dir().ok();
    let input = absolutize_from_cwd(&step.input, cwd.as_deref());
    let out_dir = absolutize_from_cwd(&step.out_dir, cwd.as_deref());
    let tool_input = if step.tool == "kira-mitoqc" {
        input.clone()
    } else if input.is_file() {
        input
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| input.clone())
    } else {
        input.clone()
    };

    if step.tool == "kira-energetics" {
        return Ok(vec![
            "run".to_string(),
            "--deps".to_string(),
            out_dir.display().to_string(),
            "--out".to_string(),
            out_dir.display().to_string(),
        ]);
    }

    if step.tool == "kira-microenvironment" {
        let cache = step
            .cache_path
            .as_ref()
            .ok_or_else(|| "kira-microenvironment requires resolved shared cache".to_string())?;
        let cache = absolutize_from_cwd(cache, cwd.as_deref());
        let expr_path = resolve_microenvironment_expr_path(&input, &out_dir, &cache);
        let barcodes = resolve_barcodes_path(&input, &out_dir)?;
        ensure_microenvironment_features_sidecar(&expr_path, &input, &barcodes)?;
        let resources = resolve_microenvironment_resource_paths()?;

        let mut args = vec![
            "run".to_string(),
            "--expr".to_string(),
            expr_path.display().to_string(),
            "--barcodes".to_string(),
            barcodes.display().to_string(),
            "--resources".to_string(),
            resources.resources_dir.display().to_string(),
            "--lr-profile".to_string(),
            resources.lr_profile.clone(),
            "--out".to_string(),
            out_dir.display().to_string(),
        ];
        if let Some(auto) = resources.auto_groups {
            args.push("--auto-groups-mode".to_string());
            args.push("hierarchical".to_string());
            args.push("--auto-groups-coarse".to_string());
            args.push(auto.coarse.display().to_string());
            args.push("--auto-groups-fine".to_string());
            args.push(auto.fine.display().to_string());
            if let Some(anti) = auto.anti {
                args.push("--auto-groups-anti".to_string());
                args.push(anti.display().to_string());
            }
            args.push("--auto-groups-unknown".to_string());
            args.push("Unclassified".to_string());
            args.push("--auto-groups-emit-scores".to_string());
        } else {
            let groups = ensure_microenvironment_groups(&out_dir)?;
            args.push("--groups".to_string());
            args.push(groups.display().to_string());
        }
        let secretion_dir = out_dir.join("kira-secretion");
        if secretion_dir.is_dir() {
            args.push("--secretion".to_string());
            args.push(secretion_dir.display().to_string());
        }
        return Ok(args);
    }

    let mut args = vec!["run".to_string()];
    args.push("--input".to_string());
    args.push(tool_input.display().to_string());
    args.push("--out".to_string());
    args.push(out_dir.display().to_string());

    // Keep mode explicit for compatibility across tools where it is required.
    args.push("--mode".to_string());
    args.push("cell".to_string());
    args.push("--run-mode".to_string());
    args.push("pipeline".to_string());

    match step.tool.as_str() {
        "kira-mitoqc" => {
            if let Some(n) = step.threads {
                args.push("--threads".to_string());
                args.push(n.to_string());
            }
            if let Some(assets_dir) = resolve_mitoqc_assets_dir(bin_path) {
                args.push("--assets".to_string());
                args.push(assets_dir.display().to_string());
            }
            if step.cache_path.is_some() {
                // kira-mitoqc interprets --cache as a directory for expr.bin, not a shared-cache file path.
                args.push("--cache".to_string());
                args.push(out_dir.display().to_string());
            }
        }
        "kira-spliceqc" | "kira-proteoqc" => {
            if let Some(n) = step.threads {
                args.push("--threads".to_string());
                args.push(n.to_string());
            }
        }
        "kira-nuclearqc" | "kira-autolys" | "kira-secretion" => {}
        _ => {
            if let Some(n) = step.threads {
                args.push("--threads".to_string());
                args.push(n.to_string());
            }
        }
    }

    if supports_cache_arg(&step.tool) && step.cache_path.is_some() && step.tool != "kira-mitoqc" {
        let cache = step.cache_path.as_ref().expect("checked is_some");
        args.push("--cache".to_string());
        args.push(
            absolutize_from_cwd(cache, cwd.as_deref())
                .display()
                .to_string(),
        );
    }

    Ok(args)
}

fn resolve_mitoqc_assets_dir(bin_path: &Path) -> Option<PathBuf> {
    if let Some(p) = std::env::var_os("KIRA_ORGANELLE_MITOQC_ASSETS").map(PathBuf::from)
        && p.is_dir()
    {
        return Some(p);
    }

    if let Some(bin_dir) = bin_path.parent() {
        let direct = bin_dir.join("assets");
        if direct.is_dir() {
            return Some(direct);
        }
        for ancestor in bin_dir.ancestors() {
            let candidate = ancestor.join("assets");
            if candidate.is_dir() {
                return Some(candidate);
            }
        }
    }

    if let Ok(cwd) = std::env::current_dir() {
        let cwd_assets = cwd.join("assets");
        if cwd_assets.is_dir() {
            return Some(cwd_assets);
        }
        for ancestor in cwd.ancestors() {
            let candidate = ancestor.join("kira-mitoqc").join("assets");
            if candidate.is_dir() {
                return Some(candidate);
            }
        }
    }

    None
}

fn resolve_microenvironment_expr_path(input: &Path, out_dir: &Path, cache: &Path) -> PathBuf {
    let input_root = if input.is_file() {
        input.parent().unwrap_or(input)
    } else {
        input
    };

    let preferred = input_root.join("expr.bin");
    if preferred.is_file() {
        return preferred;
    }

    let root_expr = out_dir.join("expr.bin");
    if root_expr.is_file() {
        return root_expr;
    }

    let cache_expr = out_dir.join("kira-organelle.bin").join("expr.bin");
    if cache_expr.is_file() {
        return cache_expr;
    }

    cache.to_path_buf()
}

fn ensure_microenvironment_features_sidecar(
    expr_path: &Path,
    input: &Path,
    barcodes: &Path,
) -> Result<(), String> {
    let Some(expr_dir) = expr_path.parent() else {
        return Ok(());
    };
    let existing_plain = expr_dir.join("features.tsv");
    let existing_gz = expr_dir.join("features.tsv.gz");
    if existing_plain.is_file() || existing_gz.is_file() {
        return Ok(());
    }

    let input_root = if input.is_file() {
        input.parent().unwrap_or(input)
    } else {
        input
    };

    let mut candidates = Vec::new();
    if let Some(parent) = barcodes.parent() {
        candidates.push(parent.join("features.tsv.gz"));
        candidates.push(parent.join("features.tsv"));
    }

    let barcode_name = barcodes
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_string();
    let stem = barcode_name
        .strip_suffix(".barcodes.tsv.gz")
        .or_else(|| barcode_name.strip_suffix(".barcodes.tsv"))
        .or_else(|| barcode_name.strip_suffix("_barcodes.tsv.gz"))
        .or_else(|| barcode_name.strip_suffix("_barcodes.tsv"))
        .map(|s| s.to_string());
    if let Some(stem) = stem {
        candidates.push(input_root.join(format!("{stem}.features.tsv.gz")));
        candidates.push(input_root.join(format!("{stem}.features.tsv")));
        candidates.push(input_root.join(format!("{stem}_features.tsv.gz")));
        candidates.push(input_root.join(format!("{stem}_features.tsv")));
    }

    candidates.push(input_root.join("features.tsv.gz"));
    candidates.push(input_root.join("features.tsv"));

    let Some(source) = candidates.into_iter().find(|p| p.is_file()) else {
        tracing::warn!(
            tool = "kira-microenvironment",
            input = %input_root.display(),
            barcodes = %barcodes.display(),
            "features sidecar was not resolved; gene lookup may degrade to fallback symbols"
        );
        return Ok(());
    };

    let target = if source
        .extension()
        .and_then(|x| x.to_str())
        .map(|x| x.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        existing_gz
    } else {
        existing_plain
    };

    std::fs::copy(&source, &target).map_err(|e| {
        format!(
            "failed copying features sidecar {} -> {}: {e}",
            source.display(),
            target.display()
        )
    })?;
    Ok(())
}

fn supports_cache_arg(tool: &str) -> bool {
    matches!(
        tool,
        "kira-mitoqc"
            | "kira-nuclearqc"
            | "kira-spliceqc"
            | "kira-proteoqc"
            | "kira-autolys"
            | "kira-secretion"
    )
}

pub fn resolve_barcodes_path(input: &Path, out_dir: &Path) -> Result<PathBuf, String> {
    let input_root = if input.is_file() {
        input.parent().unwrap_or(input)
    } else {
        input
    };

    let gz = input_root.join("barcodes.tsv.gz");
    if gz.is_file() {
        return Ok(gz);
    }
    let plain = input_root.join("barcodes.tsv");
    if plain.is_file() {
        return Ok(plain);
    }

    let mut candidates = Vec::new();
    let rd = std::fs::read_dir(input_root).map_err(|e| {
        format!(
            "failed reading input directory {}: {e}",
            input_root.display()
        )
    })?;
    for entry in rd {
        let entry = entry.map_err(|e| {
            format!(
                "failed reading input directory {}: {e}",
                input_root.display()
            )
        })?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name.ends_with("_barcodes.tsv.gz")
            || name.ends_with("_barcodes.tsv")
            || name.ends_with(".barcodes.tsv.gz")
            || name.ends_with(".barcodes.tsv")
        {
            candidates.push(path);
        }
    }
    candidates.sort();
    if let Some(found) = candidates.first() {
        return Ok(found.clone());
    }

    if let Some(raw_counts) = find_raw_counts_path(input)? {
        let target_dir = out_dir.join("kira-microenvironment");
        std::fs::create_dir_all(&target_dir)
            .map_err(|e| format!("failed creating {}: {e}", target_dir.display()))?;
        let target = target_dir.join("barcodes.from_raw_counts.tsv");
        write_barcodes_from_raw_counts(&raw_counts, &target)?;
        return Ok(target);
    }

    Err(format!(
        "barcodes file not found in {} (expected barcodes.tsv(.gz), *_barcodes.tsv(.gz), *.barcodes.tsv(.gz), or raw_counts header)",
        input_root.display()
    ))
}

fn find_raw_counts_path(input: &Path) -> Result<Option<PathBuf>, String> {
    if input.is_file() {
        let name = input
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default();
        if name.eq_ignore_ascii_case("raw_counts.tsv")
            || name.eq_ignore_ascii_case("raw_counts.tsv.gz")
            || name.ends_with("_raw_counts.tsv")
            || name.ends_with("_raw_counts.tsv.gz")
        {
            return Ok(Some(input.to_path_buf()));
        }
    }
    let input_root = if input.is_file() {
        input.parent().unwrap_or(input)
    } else {
        input
    };
    let mut candidates = Vec::new();
    let rd = std::fs::read_dir(input_root).map_err(|e| {
        format!(
            "failed reading input directory {}: {e}",
            input_root.display()
        )
    })?;
    for entry in rd {
        let entry = entry.map_err(|e| {
            format!(
                "failed reading input directory {}: {e}",
                input_root.display()
            )
        })?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name.eq_ignore_ascii_case("raw_counts.tsv")
            || name.eq_ignore_ascii_case("raw_counts.tsv.gz")
            || name.ends_with("_raw_counts.tsv")
            || name.ends_with("_raw_counts.tsv.gz")
        {
            candidates.push(path);
        }
    }
    candidates.sort();
    Ok(candidates.into_iter().next())
}

fn write_barcodes_from_raw_counts(raw_counts: &Path, target: &Path) -> Result<(), String> {
    let reader = open_maybe_gz(raw_counts)?;
    let mut header = None;
    for line in reader.lines() {
        let line = line.map_err(|e| format!("failed reading {}: {e}", raw_counts.display()))?;
        let cleaned = line.trim_end_matches(&['\r', '\n'][..]);
        if cleaned.trim().is_empty() || cleaned.trim_start().starts_with('#') {
            continue;
        }
        header = Some(cleaned.to_string());
        break;
    }
    let header =
        header.ok_or_else(|| format!("raw_counts file is empty: {}", raw_counts.display()))?;
    let cols = header
        .split('\t')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>();
    if cols.len() <= 1 {
        return Err(format!(
            "raw_counts header has insufficient columns in {}",
            raw_counts.display()
        ));
    }
    let barcodes = if first_col_looks_like_gene_header(cols.first().map(String::as_str)) {
        cols.into_iter().skip(1).collect::<Vec<_>>()
    } else {
        cols
    };
    if barcodes.is_empty() {
        return Err(format!(
            "no barcodes resolved from raw_counts header in {}",
            raw_counts.display()
        ));
    }
    let mut out = String::new();
    for bc in barcodes {
        out.push_str(&bc);
        out.push('\n');
    }
    crate::io::write_bytes_atomic(target, out.as_bytes())
        .map_err(|e| format!("failed writing {}: {e}", target.display()))
}

fn open_maybe_gz(path: &Path) -> Result<BufReader<Box<dyn Read>>, String> {
    if path
        .extension()
        .and_then(|v| v.to_str())
        .map(|v| v.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        let file =
            File::open(path).map_err(|e| format!("failed opening {}: {e}", path.display()))?;
        let gz = GzDecoder::new(file);
        Ok(BufReader::new(Box::new(gz)))
    } else {
        let file =
            File::open(path).map_err(|e| format!("failed opening {}: {e}", path.display()))?;
        Ok(BufReader::new(Box::new(file)))
    }
}

fn first_col_looks_like_gene_header(v: Option<&str>) -> bool {
    let Some(v) = v else { return false };
    if v.trim().is_empty() {
        return true;
    }
    matches!(
        v.trim().to_ascii_lowercase().as_str(),
        "gene" | "genes" | "gene_symbol" | "genesymbol" | "symbol" | "feature" | "features"
    )
}

fn ensure_microenvironment_groups(out_root: &Path) -> Result<PathBuf, String> {
    let classify = out_root.join("kira-secretion").join("classify.tsv");
    if !classify.is_file() {
        return Err(format!(
            "kira-microenvironment requires secretion classify.tsv at {}",
            classify.display()
        ));
    }
    let raw = std::fs::read_to_string(&classify)
        .map_err(|e| format!("failed reading {}: {e}", classify.display()))?;
    let mut lines = raw.lines();
    let header = lines.next().unwrap_or_default();
    let cols = header.split('\t').collect::<Vec<_>>();
    let cell_idx = cols
        .iter()
        .position(|v| v.eq_ignore_ascii_case("cell_id"))
        .ok_or_else(|| format!("missing cell_id column in {}", classify.display()))?;
    let regime_idx = cols
        .iter()
        .position(|v| v.eq_ignore_ascii_case("regime"))
        .ok_or_else(|| format!("missing regime column in {}", classify.display()))?;

    let target_dir = out_root.join("kira-microenvironment");
    std::fs::create_dir_all(&target_dir)
        .map_err(|e| format!("failed creating {}: {e}", target_dir.display()))?;
    let target = target_dir.join("groups.from_secretion.tsv");

    let mut out = String::from("cell_id\tgroup\n");
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() <= regime_idx || fields.len() <= cell_idx {
            continue;
        }
        out.push_str(fields[cell_idx]);
        out.push('\t');
        out.push_str(fields[regime_idx]);
        out.push('\n');
    }
    crate::io::write_bytes_atomic(&target, out.as_bytes())
        .map_err(|e| format!("failed writing {}: {e}", target.display()))?;
    Ok(target)
}

struct MicroenvironmentAutoGroupsPaths {
    coarse: PathBuf,
    fine: PathBuf,
    anti: Option<PathBuf>,
}

struct MicroenvironmentResources {
    resources_dir: PathBuf,
    lr_profile: String,
    auto_groups: Option<MicroenvironmentAutoGroupsPaths>,
}

fn resolve_microenvironment_resource_paths() -> Result<MicroenvironmentResources, String> {
    let root = std::env::var_os("KIRA_MICROENV_ROOT")
        .map(PathBuf::from)
        .or_else(resolve_microenvironment_root_from_binary)
        .ok_or_else(|| {
            "kira-microenvironment root is not configured; set KIRA_MICROENV_ROOT or ensure kira-microenvironment binary is discoverable"
                .to_string()
        })?;

    let resources_candidates = [root.join("resources_onco"), root.join("resources")];
    let resources_dir = resources_candidates
        .iter()
        .find(|p| p.is_dir())
        .cloned()
        .unwrap_or_else(|| root.join("resources"));

    let lr_profile = {
        let onco = root.join("resources").join("lr_pairs_mvp_onco.tsv");
        let full = root.join("resources").join("lr_pairs.tsv");
        if onco.is_file() {
            onco.display().to_string()
        } else if full.is_file() {
            full.display().to_string()
        } else {
            "mvp".to_string()
        }
    };

    let force_secretion_groups = std::env::var_os("KIRA_MICROENV_FORCE_SECRETION_GROUPS")
        .map(|v| {
            let s = v.to_string_lossy();
            let t = s.trim();
            !t.is_empty() && t != "0" && !t.eq_ignore_ascii_case("false")
        })
        .unwrap_or(false);
    let auto_groups = if force_secretion_groups {
        None
    } else {
        let coarse = root.join("marker_panels_coarse.tsv");
        let fine = root.join("marker_panels.tsv");
        let anti = root.join("anti_markers.tsv");
        if coarse.is_file() && fine.is_file() {
            Some(MicroenvironmentAutoGroupsPaths {
                coarse,
                fine,
                anti: anti.is_file().then_some(anti),
            })
        } else {
            None
        }
    };

    Ok(MicroenvironmentResources {
        resources_dir,
        lr_profile,
        auto_groups,
    })
}

fn resolve_microenvironment_root_from_binary() -> Option<PathBuf> {
    match resolve_tool_invocation("kira-microenvironment") {
        ToolInvocationMode::Binary(path) => {
            let parent = path.parent()?;
            if parent.file_name().and_then(|v| v.to_str()) == Some("release")
                || parent.file_name().and_then(|v| v.to_str()) == Some("debug")
            {
                let target_dir = parent.parent()?;
                return target_dir.parent().map(Path::to_path_buf);
            }
            parent.parent().map(Path::to_path_buf)
        }
        _ => None,
    }
}

fn absolutize_from_cwd(path: &Path, cwd: Option<&Path>) -> PathBuf {
    if path.is_absolute() {
        return path.to_path_buf();
    }
    match cwd {
        Some(base) => base.join(path),
        None => path.to_path_buf(),
    }
}

fn mode_string(mode: &ToolInvocationMode) -> &'static str {
    match mode {
        ToolInvocationMode::Library => "library",
        ToolInvocationMode::Binary(_) => "binary",
        ToolInvocationMode::Unavailable => "unavailable",
    }
}
