use std::collections::BTreeSet;
use std::fs;

use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn single_mode_interpretation_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"decay_score":{"median":0.8},"bioenergetics":{"median":0.2}},"qc":{"low_confidence_fraction":0.1}}"#,
    );

    let opts = AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
        export_systems_model: None,
    };

    run_aggregate(&opts).expect("first");
    let first = fs::read(out.join("interpretation.json")).expect("read first");

    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("interpretation.json")).expect("read second");

    assert_eq!(first, second);
}

#[test]
fn comparison_mode_signal_trigger() {
    let dir = tempdir().expect("tempdir");
    let a = dir.path().join("a");
    let b = dir.path().join("b");
    let out = dir.path().join("out");

    write_tool(
        &a,
        "kira-secretion",
        r#"{"distributions":{"er_golgi_pressure":{"median":0.2},"stress_secretion_index":{"median":0.1}}}"#,
    );
    write_tool(
        &b,
        "kira-secretion",
        r#"{"distributions":{"er_golgi_pressure":{"median":0.7},"stress_secretion_index":{"median":0.4}}}"#,
    );

    run_aggregate(&AggregateOptions {
        input: a,
        input_b: Some(b),
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
        export_systems_model: None,
    })
    .expect("aggregate");

    let v: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("interpretation.json")).expect("read"))
            .expect("parse");

    let ids = v
        .get("signals")
        .and_then(|x| x.as_array())
        .expect("signals")
        .iter()
        .filter_map(|s| s.get("id").and_then(|x| x.as_str()))
        .collect::<Vec<_>>();

    assert!(ids.contains(&"secretory_pressure_increase"));
}

#[test]
fn confidence_computation_stable() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"x":{"median":0.2}},"qc":{"low_confidence_fraction":0.2}}"#,
    );

    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
        export_systems_model: None,
    })
    .expect("aggregate");

    let v: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("interpretation.json")).expect("read"))
            .expect("parse");

    let confidence = v
        .get("confidence")
        .and_then(|x| x.as_f64())
        .expect("confidence");

    let expected = (1.0 - 0.2) * 0.5 + (1.0 / 7.0) * 0.3;
    assert!((confidence - expected).abs() < 1e-12);
}

#[test]
fn no_duplicate_signals() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"decay_score":{"median":0.9},"bioenergetics":{"median":0.1},"stress_index":{"median":0.9}},"regimes":{"counts":{"fragile":1},"fractions":{"fragile":1.0}},"qc":{"low_confidence_fraction":0.1}}"#,
    );
    write_tool(
        &input,
        "kira-nuclearqc",
        r#"{"distributions":{"nucleus_stress":{"median":0.8}},"regimes":{"counts":{"stable":1},"fractions":{"stable":1.0}}}"#,
    );
    write_tool(
        &input,
        "kira-proteoqc",
        r#"{"distributions":{"proteostasis_load":{"median":0.9},"stress_index":{"median":0.9}}}"#,
    );
    write_tool(
        &input,
        "kira-autolys",
        r#"{"distributions":{"stress_index":{"median":0.9}}}"#,
    );
    write_tool(
        &input,
        "kira-secretion",
        r#"{"distributions":{"er_golgi_pressure":{"median":0.9},"stress_secretion_index":{"median":0.9}}}"#,
    );

    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
        export_systems_model: None,
    })
    .expect("aggregate");

    let v: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("interpretation.json")).expect("read"))
            .expect("parse");

    let ids = v
        .get("signals")
        .and_then(|x| x.as_array())
        .expect("signals")
        .iter()
        .filter_map(|s| s.get("id").and_then(|x| x.as_str()).map(ToOwned::to_owned))
        .collect::<Vec<_>>();

    let unique = ids.iter().cloned().collect::<BTreeSet<_>>();
    assert_eq!(ids.len(), unique.len());
}

fn write_tool(input_root: &std::path::Path, tool: &str, summary_json: &str) {
    let tool_dir = input_root.join(tool);
    fs::create_dir_all(&tool_dir).expect("mkdir tool");

    fs::write(tool_dir.join("summary.json"), summary_json).expect("summary");
    fs::write(
        tool_dir.join("pipeline_step.json"),
        r#"{"artifacts":{"primary_metrics":"primary_metrics.tsv","summary":"summary.json"}}"#,
    )
    .expect("pipeline");
    fs::write(tool_dir.join("primary_metrics.tsv"), "barcode\tx\nC1\t1\n").expect("tsv");
}
