use std::fs;

use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn extract_axes_partial_quantiles() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");

    fs::create_dir_all(input.join("kira-mitoqc")).expect("mkdir");
    fs::write(
        input.join("kira-mitoqc/summary.json"),
        r#"{"distributions":{"mt_ratio":{"median":0.42}},"qc":{"low_confidence_fraction":0.12}}"#,
    )
    .expect("write");

    let state = run_aggregate(&AggregateOptions {
        input: input.clone(),
        input_b: None,
        out: Some(dir.path().join("out")),
        strict: false,
        json: true,
        validate_only: true,
        fii_weights: None,
    })
    .expect("aggregate");

    assert_eq!(state.organelle_states.len(), 1);
    let axis = &state.organelle_states[0].axes[0];
    assert_eq!(axis.name, "mt_ratio");
    assert_eq!(axis.median, 0.42);
    assert_eq!(axis.p90, None);
    assert_eq!(axis.p99, None);
}

#[test]
fn organelle_order_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");

    let tools = ["kira-secretion", "kira-mitoqc", "kira-autolys"];
    for tool in tools {
        fs::create_dir_all(input.join(tool)).expect("mkdir");
        fs::write(
            input.join(tool).join("summary.json"),
            r#"{"distributions":{"x":{"median":1.0}}}"#,
        )
        .expect("write");
    }

    let state = run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(dir.path().join("out")),
        strict: false,
        json: true,
        validate_only: true,
        fii_weights: None,
    })
    .expect("aggregate");

    let names = state
        .organelle_states
        .iter()
        .map(|s| serde_json::to_string(&s.organelle).expect("serialize"))
        .collect::<Vec<_>>();

    assert_eq!(
        names,
        vec!["\"mitochondria\"", "\"autophagy\"", "\"secretion\""]
    );
}

#[test]
fn missing_distributions_warns() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");

    fs::create_dir_all(input.join("kira-mitoqc")).expect("mkdir");
    fs::write(
        input.join("kira-mitoqc/summary.json"),
        r#"{"qc":{"low_confidence_fraction":0.1}}"#,
    )
    .expect("write");

    let state = run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(dir.path().join("out")),
        strict: false,
        json: true,
        validate_only: true,
        fii_weights: None,
    })
    .expect("aggregate");

    assert_eq!(state.organelle_states.len(), 1);
    assert!(
        state
            .issues
            .iter()
            .any(|i| i.code == "MISSING_DISTRIBUTIONS"
                && i.severity == kira_organelle::contracts::types::Severity::Warn)
    );
}

#[test]
fn state_json_schema_v1_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    fs::create_dir_all(input.join("kira-mitoqc")).expect("mkdir");
    fs::write(
        input.join("kira-mitoqc/summary.json"),
        r#"{"distributions":{"mt_ratio":{"median":0.5,"p90":0.8}},"qc":{"low_confidence_fraction":0.01}}"#,
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

    run_aggregate(&opts).expect("first");
    let first = fs::read(out.join("state.json")).expect("read first");

    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("state.json")).expect("read second");

    assert_eq!(first, second);

    let parsed: serde_json::Value = serde_json::from_slice(&first).expect("parse");
    assert_eq!(
        parsed.get("schema").and_then(|v| v.as_str()),
        Some("kira-organelle-state-v1")
    );
}
