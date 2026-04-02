pub mod cells;
pub mod cli;
pub mod compare;
pub mod contracts;
pub mod error;
pub mod fii;
pub mod integration_export;
pub mod interpret;
pub mod io;
pub mod logging;
pub mod model;
pub mod normalize;
pub mod registry;
pub mod render;
pub mod report;
pub mod run;
pub mod simd;
pub mod state;
pub mod stress_localization;
pub mod systems;
pub mod util;
pub mod version;

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};

use cells::types::CellsState;
use contracts::types::{Issue, Severity, ToolContracts};
use model::organelle::OrganelleId;
use serde_json::{Map, Value};
use state::{AggregatedState, OrganelleState, STATE_SCHEMA_V1};
use tracing::info;

#[derive(Debug, Clone)]
pub struct AggregateOptions {
    pub input: PathBuf,
    pub input_b: Option<PathBuf>,
    pub out: Option<PathBuf>,
    pub strict: bool,
    pub json: bool,
    pub validate_only: bool,
    pub fii_weights: Option<fii::FiiWeights>,
    pub export_systems_model: Option<PathBuf>,
}

const MIN_LLM_REPORT_CONFIDENCE: f64 = 0.6;
static PROFILE_STAGES: AtomicBool = AtomicBool::new(false);
static WRITE_CELLS_JSON: AtomicBool = AtomicBool::new(true);

#[derive(Debug, Clone, Copy, Default)]
struct BuildBundleProfile {
    normalize_ms: f64,
}

pub fn set_profile_stages_enabled(enabled: bool) {
    PROFILE_STAGES.store(enabled, Ordering::Relaxed);
}

fn profile_stages_enabled() -> bool {
    PROFILE_STAGES.load(Ordering::Relaxed)
}

pub fn set_write_cells_json_enabled(enabled: bool) {
    WRITE_CELLS_JSON.store(enabled, Ordering::Relaxed);
}

fn write_cells_json_enabled() -> bool {
    WRITE_CELLS_JSON.load(Ordering::Relaxed)
}

pub fn run_aggregate(opts: &AggregateOptions) -> Result<AggregatedState, String> {
    validate_input_dir(&opts.input)?;
    if let Some(input_b) = &opts.input_b {
        validate_input_dir(input_b)?;
    }

    let profile_stages = profile_stages_enabled();
    let plan_t0 = std::time::Instant::now();
    let (mut state, mut cells_state, bundle_profile) =
        build_bundle_from_raw(&opts.input, opts.input_b.as_deref(), opts.strict)?;
    let plan_ms = plan_t0.elapsed().as_secs_f64() * 1000.0;

    let fii = if let Some(weights) = opts.fii_weights {
        Some(fii::compute_fii(&cells_state, weights)?)
    } else {
        None
    };
    if let Some(fii) = &fii {
        state.functional_irreversibility = Some(fii.summary.clone());
    }

    let invalid_reason = if !opts.validate_only && is_pipeline_contract_layout(&opts.input) {
        validate_pipeline_contract(&opts.input, &mut state, &mut cells_state)
    } else {
        None
    };

    let out_dir = opts
        .out
        .clone()
        .unwrap_or_else(|| opts.input.join("kira-organelle"));

    if opts.validate_only {
        info!("validate-only mode, skipping file writes");
        return Ok(state);
    }

    std::fs::create_dir_all(&out_dir).map_err(|e| {
        format!(
            "failed to create output directory {}: {e}",
            out_dir.display()
        )
    })?;

    if !opts.json {
        tracing::debug!("--json not set; Stage 3 still emits state.json");
    }
    let state_value = state.to_json_value();
    io::write_json_atomic(&out_dir.join("state.json"), &state_value)
        .map_err(|e| format!("failed writing state.json: {e}"))?;

    let write_cells = write_cells_json_enabled();
    let cells_path = out_dir.join("cells.json");
    if write_cells {
        let cells_value = cells_state.to_json_value();
        io::write_json_atomic(&cells_path, &cells_value)
            .map_err(|e| format!("failed writing cells.json: {e}"))?;
    } else if let Err(e) = std::fs::remove_file(&cells_path) {
        if e.kind() != std::io::ErrorKind::NotFound {
            return Err(format!("failed removing cells.json: {e}"));
        }
    }
    if let Some(fii) = &fii {
        let tsv = fii::render_fii_tsv(&fii.rows);
        io::write_bytes_atomic(
            &out_dir.join("functional_irreversibility_index.tsv"),
            tsv.as_bytes(),
        )
        .map_err(|e| format!("failed writing functional_irreversibility_index.tsv: {e}"))?;
    }
    let mut systems_profile = systems::composites::SystemsStageProfile::default();
    integration_export::write_single_integration_artifacts(
        &state,
        &cells_state,
        &out_dir,
        opts.export_systems_model.as_deref(),
        Some(&mut systems_profile),
    )?;
    if profile_stages {
        tracing::info!(
            stage = "plan",
            ms = format_args!("{:.3}", plan_ms),
            "stage_profile"
        );
        tracing::info!(
            stage = "normalize",
            ms = format_args!("{:.3}", bundle_profile.normalize_ms),
            "stage_profile"
        );
        tracing::info!(
            stage = "composite",
            ms = format_args!("{:.3}", systems_profile.composite_ms),
            "stage_profile"
        );
        tracing::info!(
            stage = "aggregate",
            ms = format_args!("{:.3}", systems_profile.aggregate_ms),
            "stage_profile"
        );
        tracing::info!(
            stage = "fragility",
            ms = format_args!("{:.3}", systems_profile.fragility_ms),
            "stage_profile"
        );
    }

    if let Some(reason) = invalid_reason {
        return Err(reason);
    }

    let comparison = if let Some(input_b) = &opts.input_b {
        let mut comparison_issues = Vec::new();

        let (a_state_cmp, a_cells_cmp) =
            compare::load_side_from_output(&opts.input, &mut comparison_issues)
                .unwrap_or_else(|| (state.clone(), cells_state.clone()));

        let (b_state_cmp, b_cells_cmp) =
            match compare::load_side_from_output(input_b, &mut comparison_issues) {
                Some(v) => v,
                None => {
                    let (s, c, _) = build_bundle_from_raw(input_b, None, opts.strict)?;
                    (s, c)
                }
            };

        Some(compare::build_comparison(
            &a_state_cmp,
            &a_cells_cmp,
            &b_state_cmp,
            &b_cells_cmp,
            &opts.input,
            input_b,
            comparison_issues,
        ))
    } else {
        None
    };
    let interpretation_issues = comparison
        .as_ref()
        .map(|c| c.issues.clone())
        .unwrap_or_else(|| state.issues.clone());
    let interpretation =
        interpret::engine::build_interpretation(&state, comparison.as_ref(), interpretation_issues);

    if let Some(comparison) = comparison.as_ref() {
        io::write_json_atomic(
            &out_dir.join("comparison.json"),
            &comparison.to_json_value(),
        )
        .map_err(|e| format!("failed writing comparison.json: {e}"))?;
    }
    io::write_json_atomic(
        &out_dir.join("interpretation.json"),
        &interpretation.to_json_value(),
    )
    .map_err(|e| format!("failed writing interpretation.json: {e}"))?;

    let llm_report_enabled = interpretation.confidence >= MIN_LLM_REPORT_CONFIDENCE;
    if llm_report_enabled {
        let llm_report = render::llm::build_llm_report_markdown(&state, &interpretation);
        io::write_bytes_atomic(&out_dir.join("llm_report.md"), llm_report.as_bytes())
            .map_err(|e| format!("failed writing llm_report.md: {e}"))?;
    } else if let Err(e) = std::fs::remove_file(out_dir.join("llm_report.md")) {
        if e.kind() != std::io::ErrorKind::NotFound {
            return Err(format!("failed removing llm_report.md: {e}"));
        }
    }

    let report_html = render::html::build_report_html(&state, comparison.as_ref());
    io::write_bytes_atomic(&out_dir.join("report.html"), report_html.as_bytes())
        .map_err(|e| format!("failed writing report.html: {e}"))?;

    let pipeline = build_pipeline_step_json(
        &opts.input,
        opts.input_b.as_deref(),
        opts.input_b.is_some(),
        write_cells,
        llm_report_enabled,
        fii.is_some(),
    );
    io::write_json_atomic(&out_dir.join("pipeline_step.json"), &pipeline)
        .map_err(|e| format!("failed writing pipeline_step.json: {e}"))?;

    Ok(state)
}

pub fn run_run_command(args: &cli::RunArgs) -> Result<(), String> {
    run::run_pipeline(args)
}

pub fn run_report_fii_command(args: &cli::ReportFiiArgs) -> Result<(), String> {
    report::fii::run_report_fii(args)
}

pub fn run_compute_state_dynamics_command(
    args: &cli::ComputeStateDynamicsArgs,
) -> Result<(), String> {
    report::fii::run_compute_state_dynamics(args)
}

pub fn run_report_phase_command(args: &cli::ReportPhaseArgs) -> Result<(), String> {
    report::fii::run_report_phase(args)
}

pub fn run_decision_command(args: &cli::DecisionArgs) -> Result<(), String> {
    report::fii::run_decision(args)
}

pub fn run_report_decision_timeline_command(
    args: &cli::ReportDecisionTimelineArgs,
) -> Result<(), String> {
    report::fii::run_report_decision_timeline(args)
}

pub fn run_compute_ili_command(args: &cli::ComputeIliArgs) -> Result<(), String> {
    report::fii::run_compute_ili(args)
}

pub fn run_compute_cai_command(args: &cli::ComputeCaiArgs) -> Result<(), String> {
    report::fii::run_compute_cai(args)
}

pub fn run_compute_pri_command(args: &cli::ComputePriArgs) -> Result<(), String> {
    report::fii::run_compute_pri(args)
}

pub fn run_compute_cocs_command(args: &cli::ComputeCocsArgs) -> Result<(), String> {
    report::fii::run_compute_cocs(args)
}

pub fn run_compute_dci_command(args: &cli::ComputeDciArgs) -> Result<(), String> {
    report::fii::run_compute_dci(args)
}

fn validate_input_dir(path: &Path) -> Result<(), String> {
    if !path.exists() {
        return Err(format!(
            "input directory does not exist: {}",
            path.display()
        ));
    }
    if !path.is_dir() {
        return Err(format!("input path is not a directory: {}", path.display()));
    }
    Ok(())
}

fn build_bundle_from_raw(
    input: &Path,
    input_b: Option<&Path>,
    strict: bool,
) -> Result<(AggregatedState, CellsState, BuildBundleProfile), String> {
    let mut issues = Vec::new();
    let tools = contracts::discover_and_read(input, strict, &mut issues)?;

    let mut state = build_state(input, input_b, tools.clone(), issues);
    let mut cells = cells::join::build_cells_state(input, &tools, &state.issues);
    let normalize_t0 = std::time::Instant::now();
    let normalization = normalize::global_robust::compute_global_robust_normalization(&cells);
    let normalize_ms = normalize_t0.elapsed().as_secs_f64() * 1000.0;
    state.global_normalization = Some(normalization.summary.clone());
    cells.normalization_context = Some(normalization);
    stress_localization::compute_and_attach(&mut state, &cells);

    Ok((state, cells, BuildBundleProfile { normalize_ms }))
}

fn build_state(
    input: &Path,
    input_b: Option<&Path>,
    tools: Vec<ToolContracts>,
    mut issues: Vec<Issue>,
) -> AggregatedState {
    let mut organelle_states = tools
        .iter()
        .filter_map(|tool| {
            model::extract::extract_organelle_state(&tool.name, &tool.summary_json, &mut issues)
        })
        .collect::<Vec<OrganelleState>>();

    organelle_states.sort_by_key(|s| s.organelle);

    AggregatedState {
        schema: STATE_SCHEMA_V1.to_string(),
        tool: state::ToolMeta {
            name: "kira-organelle".to_string(),
            version: version::tool_version().to_string(),
        },
        inputs: state::Inputs {
            a: input.to_string_lossy().to_string(),
            b: input_b.map(|p| p.to_string_lossy().to_string()),
        },
        organelle_states,
        stress_localization: None,
        functional_irreversibility: None,
        systems_state_model: None,
        global_normalization: None,
        issues,
        timestamp: io::deterministic_rfc3339_utc(),
    }
}

pub fn build_pipeline_step_json(
    input: &Path,
    input_b: Option<&Path>,
    has_comparison: bool,
    has_cells: bool,
    has_llm_report: bool,
    has_fii: bool,
) -> Value {
    let mut root = Map::new();

    let mut tool = Map::new();
    tool.insert(
        "name".to_string(),
        Value::String("kira-organelle".to_string()),
    );
    tool.insert(
        "version".to_string(),
        Value::String(version::tool_version().to_string()),
    );
    root.insert("tool".to_string(), Value::Object(tool));
    root.insert("mode".to_string(), Value::String("aggregate".to_string()));

    let mut inputs = Map::new();
    inputs.insert(
        "a".to_string(),
        Value::String(input.to_string_lossy().to_string()),
    );
    match input_b {
        Some(p) => {
            inputs.insert(
                "b".to_string(),
                Value::String(p.to_string_lossy().to_string()),
            );
        }
        None => {
            inputs.insert("b".to_string(), Value::Null);
        }
    }
    root.insert("inputs".to_string(), Value::Object(inputs));

    let mut artifacts = Map::new();
    artifacts.insert("state".to_string(), Value::String("state.json".to_string()));
    if has_cells {
        artifacts.insert("cells".to_string(), Value::String("cells.json".to_string()));
    }
    artifacts.insert(
        "report".to_string(),
        Value::String("report.html".to_string()),
    );
    artifacts.insert(
        "interpretation".to_string(),
        Value::String("interpretation.json".to_string()),
    );
    artifacts.insert(
        "integration_timeseries".to_string(),
        Value::String("integration/timeseries.tsv".to_string()),
    );
    artifacts.insert(
        "integration_summary".to_string(),
        Value::String("integration/summary.json".to_string()),
    );
    artifacts.insert(
        "integration_expression_aggregated".to_string(),
        Value::String("integration/expression_aggregated.tsv".to_string()),
    );
    artifacts.insert(
        "integration_metrics".to_string(),
        Value::String("integration/metrics.tsv".to_string()),
    );
    artifacts.insert(
        "integration_metrics_normalized".to_string(),
        Value::String("integration/metrics_normalized.tsv".to_string()),
    );
    artifacts.insert(
        "integration_systems_model".to_string(),
        Value::String("integration/systems_model.json".to_string()),
    );
    artifacts.insert(
        "integration_landscape_report".to_string(),
        Value::String("integration/landscape.html".to_string()),
    );
    if has_llm_report {
        artifacts.insert(
            "llm_report".to_string(),
            Value::String("llm_report.md".to_string()),
        );
    }
    if has_comparison {
        artifacts.insert(
            "comparison".to_string(),
            Value::String("comparison.json".to_string()),
        );
    }
    if has_fii {
        artifacts.insert(
            "functional_irreversibility_index".to_string(),
            Value::String("functional_irreversibility_index.tsv".to_string()),
        );
    }
    root.insert("artifacts".to_string(), Value::Object(artifacts));

    root.insert(
        "timestamp".to_string(),
        Value::String(io::now_rfc3339_utc()),
    );

    Value::Object(root)
}

fn is_pipeline_contract_layout(input: &Path) -> bool {
    const MIN_EXPECTED_TOOL_DIRS_FOR_FAIL_FAST: usize = 3;

    contracts::EXPECTED_TOOLS
        .iter()
        .filter(|tool| input.join(*tool).is_dir())
        .count()
        >= MIN_EXPECTED_TOOL_DIRS_FOR_FAIL_FAST
}

fn validate_pipeline_contract(
    input: &Path,
    state: &mut AggregatedState,
    cells: &mut CellsState,
) -> Option<String> {
    promote_malformed_issues_to_error(&mut state.issues);
    promote_malformed_issues_to_error(&mut cells.issues);

    let mut violations = Vec::new();

    if cells.n_cells == 0 {
        violations.push("n_cells == 0".to_string());
    }

    let valid_organelle_count = state
        .organelle_states
        .iter()
        .filter(|s| !s.axes.is_empty() && s.regimes.is_some())
        .count();
    if valid_organelle_count < 2 {
        violations.push(format!(
            "valid organelles with axes+regimes < 2 (got {valid_organelle_count})"
        ));
    }

    let has_nuclear = has_valid_organelle(state, OrganelleId::Nucleus);
    let has_mito = has_valid_organelle(state, OrganelleId::Mitochondria);
    let has_proteo = has_valid_organelle(state, OrganelleId::Proteostasis);
    if !(has_nuclear && (has_mito || has_proteo)) {
        violations.push(
            "required organelle set missing: need nuclearqc + (mitoqc or proteoqc)".to_string(),
        );
    }

    let missing_summary_count = state
        .issues
        .iter()
        .filter(|i| i.code == "MISSING_SUMMARY")
        .count();
    if missing_summary_count >= 3 {
        violations.push(format!(
            "critical modules missing summary.json >= 3 (got {missing_summary_count})"
        ));
    }

    let malformed_errors = state
        .issues
        .iter()
        .filter(|i| i.severity == Severity::Error && is_malformed_contract_code(&i.code))
        .count();
    if malformed_errors > 0 {
        violations.push(format!(
            "malformed contract schema errors present (count {malformed_errors})"
        ));
    }

    if violations.is_empty() {
        return None;
    }

    let reason = format!("RUN_INVALID: {}", violations.join("; "));
    let issue = Issue {
        severity: Severity::Error,
        tool: None,
        code: "RUN_INVALID".to_string(),
        message: reason.clone(),
        path: Some(input.display().to_string()),
    };
    state.issues.push(issue.clone());
    cells.issues.push(issue);
    Some(reason)
}

fn has_valid_organelle(state: &AggregatedState, organelle: OrganelleId) -> bool {
    state
        .organelle_states
        .iter()
        .any(|s| s.organelle == organelle && !s.axes.is_empty() && s.regimes.is_some())
}

fn promote_malformed_issues_to_error(issues: &mut [Issue]) {
    for issue in issues {
        if is_malformed_contract_code(&issue.code) {
            issue.severity = Severity::Error;
        }
    }
}

fn is_malformed_contract_code(code: &str) -> bool {
    code.starts_with("MALFORMED_")
        || matches!(code, "MISSING_REGIME_COUNTS" | "MISSING_REGIME_FRACTIONS")
}

pub fn strict_missing_issue(
    tool: Option<&str>,
    code: &str,
    msg: String,
    path: Option<String>,
) -> Issue {
    Issue {
        severity: Severity::Error,
        tool: tool.map(ToOwned::to_owned),
        code: code.to_string(),
        message: msg,
        path,
    }
}

pub fn warn_missing_issue(
    tool: Option<&str>,
    code: &str,
    msg: String,
    path: Option<String>,
) -> Issue {
    Issue {
        severity: Severity::Warn,
        tool: tool.map(ToOwned::to_owned),
        code: code.to_string(),
        message: msg,
        path,
    }
}
