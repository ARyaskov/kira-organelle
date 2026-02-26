use std::fs;

use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn axis_delta_correctness() {
    let dir = tempdir().expect("tempdir");
    let a = dir.path().join("a");
    let b = dir.path().join("b");
    let out = dir.path().join("out");

    write_tool(
        &a,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0,"p90":2.0,"p99":3.0}}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );
    write_tool(
        &b,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":2.5,"p90":4.0,"p99":7.0}}}"#,
        "barcode\taxis_x\nC1\t2.0\n",
    );

    run_aggregate(&AggregateOptions {
        input: a,
        input_b: Some(b),
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect("aggregate");

    let cmp: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("comparison.json")).expect("read comparison"))
            .expect("parse comparison");

    let axis = find_axis_delta(&cmp, "mitochondria", "axis_x").expect("axis delta");
    assert_eq!(
        axis.pointer("/delta/median").and_then(|v| v.as_f64()),
        Some(1.5)
    );
    assert_eq!(
        axis.pointer("/delta/p90").and_then(|v| v.as_f64()),
        Some(2.0)
    );
    assert_eq!(
        axis.pointer("/delta/p99").and_then(|v| v.as_f64()),
        Some(4.0)
    );
}

#[test]
fn regime_union_delta() {
    let dir = tempdir().expect("tempdir");
    let a = dir.path().join("a");
    let b = dir.path().join("b");
    let out = dir.path().join("out");

    write_tool(
        &a,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}},"regimes":{"counts":{"r1":10},"fractions":{"r1":0.4}}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );
    write_tool(
        &b,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}},"regimes":{"counts":{"r2":10},"fractions":{"r2":0.7}}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );

    run_aggregate(&AggregateOptions {
        input: a,
        input_b: Some(b),
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect("aggregate");

    let cmp: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("comparison.json")).expect("read comparison"))
            .expect("parse comparison");

    let reg = find_organelle_delta(&cmp, "mitochondria")
        .and_then(|x| x.get("regimes"))
        .expect("regimes");

    assert_eq!(
        reg.pointer("/delta/r1").and_then(|v| v.as_f64()),
        Some(-0.4)
    );
    assert_eq!(reg.pointer("/delta/r2").and_then(|v| v.as_f64()), Some(0.7));
}

#[test]
fn cell_delta_common_only() {
    let dir = tempdir().expect("tempdir");
    let a = dir.path().join("a");
    let b = dir.path().join("b");
    let out = dir.path().join("out");

    write_tool(
        &a,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}}}"#,
        "barcode\taxis_x\nC1\t1.0\nC2\t5.0\n",
    );
    write_tool(
        &b,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}}}"#,
        "barcode\taxis_x\nC1\t4.0\nC3\t9.0\n",
    );

    run_aggregate(&AggregateOptions {
        input: a,
        input_b: Some(b),
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect("aggregate");

    let cmp: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("comparison.json")).expect("read comparison"))
            .expect("parse comparison");

    assert_eq!(
        cmp.pointer("/cell_deltas/n_common_cells")
            .and_then(|v| v.as_u64()),
        Some(1)
    );

    let axis = cmp
        .pointer("/cell_deltas/axes")
        .and_then(|v| v.as_array())
        .and_then(|v| v.first())
        .expect("cell axis summary");

    assert_eq!(axis.get("median_delta").and_then(|v| v.as_f64()), Some(3.0));
    assert_eq!(
        axis.get("p90_abs_delta").and_then(|v| v.as_f64()),
        Some(3.0)
    );
}

#[test]
fn comparison_deterministic() {
    let dir = tempdir().expect("tempdir");
    let a = dir.path().join("a");
    let b = dir.path().join("b");
    let out = dir.path().join("out");

    write_tool(
        &a,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":1.0}}}"#,
        "barcode\taxis_x\nC1\t1.0\n",
    );
    write_tool(
        &b,
        "kira-mitoqc",
        r#"{"distributions":{"axis_x":{"median":2.0}}}"#,
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
    let first = fs::read(out.join("comparison.json")).expect("read first");

    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("comparison.json")).expect("read second");

    assert_eq!(first, second);

    let parsed: serde_json::Value = serde_json::from_slice(&first).expect("parse");
    assert_eq!(
        parsed.get("schema").and_then(|v| v.as_str()),
        Some("kira-organelle-comparison-v1")
    );
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

fn find_organelle_delta<'a>(
    cmp: &'a serde_json::Value,
    organelle: &str,
) -> Option<&'a serde_json::Value> {
    cmp.get("organelle_deltas")
        .and_then(|v| v.as_array())
        .and_then(|arr| {
            arr.iter()
                .find(|o| o.get("organelle").and_then(|v| v.as_str()) == Some(organelle))
        })
}

fn find_axis_delta<'a>(
    cmp: &'a serde_json::Value,
    organelle: &str,
    axis_name: &str,
) -> Option<&'a serde_json::Value> {
    find_organelle_delta(cmp, organelle)
        .and_then(|o| o.get("axes"))
        .and_then(|v| v.as_array())
        .and_then(|arr| {
            arr.iter()
                .find(|a| a.get("name").and_then(|v| v.as_str()) == Some(axis_name))
        })
}
