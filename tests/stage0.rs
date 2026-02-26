use std::fs;

use kira_organelle::contracts::discover_and_read;
use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn discover_tools_partial() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    fs::create_dir_all(input.join("kira-mitoqc")).expect("mkdir");
    fs::create_dir_all(input.join("kira-secretion")).expect("mkdir");
    fs::write(
        input.join("kira-mitoqc/summary.json"),
        r#"{"status":"ok","pass":true}"#,
    )
    .expect("write");
    fs::write(
        input.join("kira-secretion/summary.json"),
        r#"{"status":"ok","n_reads":42}"#,
    )
    .expect("write");

    let mut issues = Vec::new();
    let tools = discover_and_read(&input, false, &mut issues).expect("discover ok");

    assert_eq!(tools.len(), 2);
    let missing_dir_warnings = issues
        .iter()
        .filter(|i| {
            i.code == "MISSING_TOOL_DIR"
                && i.severity == kira_organelle::contracts::types::Severity::Warn
        })
        .count();
    assert_eq!(missing_dir_warnings, 6);
}

#[test]
fn parse_pipeline_step_optional() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    fs::create_dir_all(input.join("kira-mitoqc")).expect("mkdir");
    fs::create_dir_all(input.join("kira-secretion")).expect("mkdir");

    fs::write(input.join("kira-mitoqc/summary.json"), r#"{"status":"ok"}"#).expect("write");
    fs::write(
        input.join("kira-secretion/summary.json"),
        r#"{"status":"ok"}"#,
    )
    .expect("write");
    fs::write(
        input.join("kira-secretion/pipeline_step.json"),
        r#"{"artifacts":{"primary_metrics":"metrics.tsv","summary":"summary.json"}}"#,
    )
    .expect("write");

    let mut issues = Vec::new();
    let tools = discover_and_read(&input, false, &mut issues).expect("discover ok");

    let mito = tools
        .iter()
        .find(|t| t.name == "kira-mitoqc")
        .expect("mito");
    assert!(mito.pipeline_step_path.is_none());

    let sec = tools
        .iter()
        .find(|t| t.name == "kira-secretion")
        .expect("sec");
    assert!(sec.pipeline_step_path.is_some());
    assert_eq!(sec.primary_metrics_path.as_deref(), Some("metrics.tsv"));
}

#[test]
fn state_json_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out = dir.path().join("out");
    fs::create_dir_all(input.join("kira-mitoqc")).expect("mkdir");
    fs::write(
        input.join("kira-mitoqc/summary.json"),
        r#"{"status":"ok","pass":true,"n_reads":100}"#,
    )
    .expect("write");

    let opts = AggregateOptions {
        input: input.clone(),
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    };

    run_aggregate(&opts).expect("first run");
    let first = fs::read(out.join("state.json")).expect("read first");

    run_aggregate(&opts).expect("second run");
    let second = fs::read(out.join("state.json")).expect("read second");

    assert_eq!(first, second);
}
