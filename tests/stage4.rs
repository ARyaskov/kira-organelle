use std::fs;

use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn html_single_input_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0,"p90":2.0}},"qc":{"low_confidence_fraction":0.05}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );

    let opts = AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    };

    run_aggregate(&opts).expect("first");
    let first = fs::read(out.join("report.html")).expect("read first");

    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("report.html")).expect("read second");

    assert_eq!(first, second);
}

#[test]
fn html_comparison_deterministic() {
    let dir = tempdir().expect("tempdir");
    let a = dir.path().join("a");
    let b = dir.path().join("b");
    let out = dir.path().join("out");

    write_tool(
        &a,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}},"regimes":{"counts":{"r1":1},"fractions":{"r1":0.1}}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );
    write_tool(
        &b,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":2.0}},"regimes":{"counts":{"r2":1},"fractions":{"r2":0.3}}}"#,
        "barcode\taxis_x\nC1\t2.0\n",
    );

    let opts = AggregateOptions {
        input: a,
        input_b: Some(b),
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    };

    run_aggregate(&opts).expect("first");
    let first = fs::read(out.join("report.html")).expect("read first");

    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("report.html")).expect("read second");

    assert_eq!(first, second);
}

#[test]
fn html_contains_expected_sections() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );

    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect("aggregate");

    let html = fs::read_to_string(out.join("report.html")).expect("read");

    for section in [
        "<h2>Header</h2>",
        "<h2>Global Summary</h2>",
        "<h2>Organelle Sections</h2>",
        "<h2>Issues &amp; Warnings</h2>",
        "<h2>Footer / Metadata</h2>",
    ] {
        assert!(html.contains(section), "missing section: {section}");
    }
}

#[test]
fn llm_report_deterministic_and_contains_tools() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"decay_score":{"median":0.9},"bioenergetics":{"median":0.1}},"regimes":{"counts":{"fragile":1},"fractions":{"fragile":1.0}},"qc":{"low_confidence_fraction":0.1}}"#,
        "barcode\tdecay_score\tbioenergetics\nC1\t0.9\t0.1\n",
    );
    write_tool(
        &input,
        "kira-nuclearqc",
        r#"{"distributions":{"nucleus_stress":{"median":0.4}},"regimes":{"counts":{"stable":1},"fractions":{"stable":1.0}},"qc":{"low_confidence_fraction":0.1}}"#,
        "barcode\tnucleus_stress\nC1\t0.4\n",
    );
    write_tool(
        &input,
        "kira-proteoqc",
        r#"{"distributions":{"proteostasis_load":{"median":0.3}},"qc":{"low_confidence_fraction":0.1}}"#,
        "barcode\tproteostasis_load\nC1\t0.3\n",
    );
    write_tool(
        &input,
        "kira-secretion",
        r#"{"distributions":{"secretory_load":{"median":0.2}},"qc":{"low_confidence_fraction":0.1}}"#,
        "barcode\tsecretory_load\nC1\t0.2\n",
    );

    let opts = AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    };

    run_aggregate(&opts).expect("first");
    let first = fs::read(out.join("llm_report.md")).expect("read first");
    let pipeline: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("pipeline_step.json")).expect("read pipeline"))
            .expect("parse pipeline");

    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("llm_report.md")).expect("read second");

    assert_eq!(first, second);

    let text = String::from_utf8(first).expect("utf8");
    assert!(text.contains("## Utility Goals"));
    assert!(text.contains("### kira-mitoqc"));
    assert!(text.contains("### kira-secretion"));
    assert_eq!(
        pipeline
            .get("artifacts")
            .and_then(|v| v.get("llm_report"))
            .and_then(|v| v.as_str()),
        Some("llm_report.md")
    );
}

#[test]
fn llm_report_suppressed_when_confidence_low() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}},"qc":{"low_confidence_fraction":1.0}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );

    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect("aggregate");

    assert!(!out.join("llm_report.md").is_file());
    let pipeline: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("pipeline_step.json")).expect("read pipeline"))
            .expect("parse pipeline");
    assert!(
        pipeline
            .get("artifacts")
            .and_then(|v| v.get("llm_report"))
            .is_none()
    );
}

#[test]
fn full_pipeline_layout_fail_fast_on_invalid_bundle() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    for tool in [
        "kira-mitoqc",
        "kira-nuclearqc",
        "kira-spliceqc",
        "kira-proteoqc",
        "kira-autolys",
        "kira-secretion",
    ] {
        fs::create_dir_all(input.join(tool)).expect("mkdir tool dir");
    }

    fs::write(
        input.join("kira-nuclearqc/summary.json"),
        r#"{"distributions":{"axes":1},"regimes":{}}"#,
    )
    .expect("write nuclear summary");

    let err = run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect_err("must fail-fast");

    assert!(err.contains("RUN_INVALID"));
    assert!(out.join("state.json").is_file());
    assert!(out.join("cells.json").is_file());
    assert!(!out.join("report.html").is_file());
    assert!(!out.join("interpretation.json").is_file());
    assert!(!out.join("llm_report.md").is_file());

    let state: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("state.json")).expect("read state"))
            .expect("parse state");

    assert!(
        state
            .get("issues")
            .and_then(|v| v.as_array())
            .expect("issues")
            .iter()
            .any(|i| i.get("code").and_then(|v| v.as_str()) == Some("RUN_INVALID"))
    );
}

#[test]
fn partial_pipeline_layout_still_fail_fast_when_contract_invalid() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    for tool in ["kira-mitoqc", "kira-nuclearqc", "kira-secretion"] {
        fs::create_dir_all(input.join(tool)).expect("mkdir tool dir");
    }

    fs::write(
        input.join("kira-nuclearqc/summary.json"),
        r#"{"distributions":{"axes":1},"regimes":{}}"#,
    )
    .expect("write nuclear summary");

    let err = run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect_err("must fail-fast");

    assert!(err.contains("RUN_INVALID"));
    assert!(out.join("state.json").is_file());
    assert!(out.join("cells.json").is_file());
    assert!(!out.join("report.html").is_file());
    assert!(!out.join("interpretation.json").is_file());
    assert!(!out.join("llm_report.md").is_file());
}

fn write_tool(input_root: &std::path::Path, tool: &str, summary_json: &str, tsv: &str) {
    let tool_dir = input_root.join(tool);
    fs::create_dir_all(&tool_dir).expect("mkdir tool");

    fs::write(tool_dir.join("summary.json"), summary_json).expect("summary");
    fs::write(
        tool_dir.join("pipeline_step.json"),
        r#"{"artifacts":{"primary_metrics":"primary_metrics.tsv","summary":"summary.json"}}"#,
    )
    .expect("pipeline");
    fs::write(tool_dir.join("primary_metrics.tsv"), tsv).expect("tsv");
}
