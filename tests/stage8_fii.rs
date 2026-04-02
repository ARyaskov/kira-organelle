use std::fs;

use kira_organelle::fii::FiiWeights;
use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn fii_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");
    write_fii_tools(&input);

    let opts = AggregateOptions {
        input: input.clone(),
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: Some(FiiWeights::default()),
        export_systems_model: None,
    };
    run_aggregate(&opts).expect("first");
    let first = fs::read(out.join("functional_irreversibility_index.tsv")).expect("read first");
    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("functional_irreversibility_index.tsv")).expect("read second");
    assert_eq!(first, second);
}

#[test]
fn fii_weight_variation_changes_score() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out_a = dir.path().join("out_a");
    let out_b = dir.path().join("out_b");
    write_fii_tools(&input);

    run_aggregate(&AggregateOptions {
        input: input.clone(),
        input_b: None,
        out: Some(out_a.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: Some(FiiWeights {
            mitochondrial: 1.0,
            translation: 0.0,
            splice: 0.0,
        }),
        export_systems_model: None,
    })
    .expect("aggregate a");
    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out_b.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: Some(FiiWeights {
            mitochondrial: 0.0,
            translation: 1.0,
            splice: 0.0,
        }),
        export_systems_model: None,
    })
    .expect("aggregate b");

    let a = read_fii_for_cell(&out_a, "BC1");
    let b = read_fii_for_cell(&out_b, "BC1");
    assert!(a > b);
}

#[test]
fn fii_regime_boundaries() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        "barcode\tmitochondrial_stress_adaptation_score\nBC_A\t0.0\nBC_B\t0.33\nBC_C\t0.66\n",
    );
    write_tool(
        &input,
        "kira-riboqc",
        "barcode\ttranslation_commitment_score\nBC_A\t0.0\nBC_B\t0.33\nBC_C\t0.66\n",
    );
    write_tool(
        &input,
        "kira-spliceqc",
        "barcode\tsplice_irreversibility_index\nBC_A\t0.0\nBC_B\t0.33\nBC_C\t0.66\n",
    );

    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: Some(FiiWeights::default()),
        export_systems_model: None,
    })
    .expect("aggregate");

    let tsv = fs::read_to_string(out.join("functional_irreversibility_index.tsv")).expect("tsv");
    assert!(tsv.contains("\nBC_A\t0.000000\t0.000000\t0.000000\t0.000000\tAdaptive\tfalse\n"));
    assert!(tsv.contains("\nBC_B\t0.330000\t0.330000\t0.330000\t0.330000\tTransition\tfalse\n"));
    assert!(tsv.contains("\nBC_C\t0.660000\t0.660000\t0.660000\t0.660000\tResistant\tfalse\n"));
}

#[test]
fn fii_low_confidence_when_inputs_missing() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        "barcode\tmitochondrial_stress_adaptation_score\nBC1\t0.8\n",
    );
    write_tool(
        &input,
        "kira-spliceqc",
        "barcode\tsplice_irreversibility_index\nBC1\t0.6\n",
    );
    // No riboqc row -> low confidence expected.

    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: Some(FiiWeights::default()),
        export_systems_model: None,
    })
    .expect("aggregate");

    let tsv = fs::read_to_string(out.join("functional_irreversibility_index.tsv")).expect("tsv");
    assert!(tsv.contains("\ttrue\n"));

    let state: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("state.json")).expect("state")).expect("parse");
    let frac = state
        .get("functional_irreversibility")
        .and_then(|v| v.get("low_confidence_fraction"))
        .and_then(|v| v.as_f64())
        .expect("low_confidence_fraction");
    assert!((frac - 1.0).abs() < 1e-9);
}

#[test]
fn fii_disabled_keeps_legacy_artifacts() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");
    write_fii_tools(&input);

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

    assert!(!out.join("functional_irreversibility_index.tsv").exists());
    let state: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("state.json")).expect("state")).expect("parse");
    assert!(
        state
            .get("functional_irreversibility")
            .is_some_and(serde_json::Value::is_null)
    );
}

fn read_fii_for_cell(out: &std::path::Path, cell_id: &str) -> f64 {
    let tsv = fs::read_to_string(out.join("functional_irreversibility_index.tsv")).expect("tsv");
    let line = tsv
        .lines()
        .find(|l| l.starts_with(&format!("{cell_id}\t")))
        .expect("cell row");
    let cols = line.split('\t').collect::<Vec<_>>();
    cols[4].parse::<f64>().expect("fii parse")
}

fn write_fii_tools(input_root: &std::path::Path) {
    write_tool(
        input_root,
        "kira-mitoqc",
        "barcode\tmitochondrial_stress_adaptation_score\nBC1\t0.90\nBC2\t0.20\n",
    );
    write_tool(
        input_root,
        "kira-riboqc",
        "barcode\ttranslation_commitment_score\nBC1\t0.10\nBC2\t0.70\n",
    );
    write_tool(
        input_root,
        "kira-spliceqc",
        "barcode\tsplice_irreversibility_index\nBC1\t0.50\nBC2\t0.40\n",
    );
}

fn write_tool(input_root: &std::path::Path, tool: &str, metrics_body: &str) {
    let tool_dir = input_root.join(tool);
    fs::create_dir_all(&tool_dir).expect("mkdir");
    fs::write(
        tool_dir.join("summary.json"),
        r#"{"distributions":{"x":{"median":1.0}}}"#,
    )
    .expect("summary");
    fs::write(
        tool_dir.join("pipeline_step.json"),
        r#"{"artifacts":{"primary_metrics":"primary_metrics.tsv","summary":"summary.json"}}"#,
    )
    .expect("pipeline");
    fs::write(tool_dir.join("primary_metrics.tsv"), metrics_body).expect("metrics");
}
