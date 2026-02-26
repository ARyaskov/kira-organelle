use std::fs;

use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn join_partial_cells() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        "primary_metrics.tsv",
        "barcode\taxis_a\nBC01\t1.0\nBC02\t2.0\n",
    );
    write_tool(
        &input,
        "kira-riboqc",
        "primary_metrics.tsv",
        "barcode\taxis_b\nBC02\t3.0\nBC03\t4.0\n",
    );

    run_aggregate(&AggregateOptions {
        input: input.clone(),
        input_b: None,
        out: Some(out.clone()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
    })
    .expect("aggregate");

    let cells: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("cells.json")).expect("read cells"))
            .expect("parse");

    assert_eq!(cells.get("n_cells").and_then(|v| v.as_u64()), Some(3));
    let ids = cells
        .get("cells")
        .and_then(|v| v.as_array())
        .expect("cells array")
        .iter()
        .filter_map(|c| c.get("id").and_then(|v| v.as_str()).map(ToOwned::to_owned))
        .collect::<Vec<_>>();
    assert_eq!(ids, vec!["BC01", "BC02", "BC03"]);
}

#[test]
fn cell_key_detection() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        "primary_metrics.tsv",
        "sample\taxis_a\nS1\t1.0\n",
    );
    write_tool(
        &input,
        "kira-riboqc",
        "primary_metrics.tsv",
        "barcode\taxis_b\nBC1\t2.0\n",
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

    let cells: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("cells.json")).expect("read cells"))
            .expect("parse");
    assert_eq!(
        cells.get("cell_key").and_then(|v| v.as_str()),
        Some("barcode")
    );
}

#[test]
fn axis_sorting_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        "primary_metrics.tsv",
        "barcode\tz_axis\ta_axis\tm_axis\nBC1\t9\t1\t5\n",
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

    let cells: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("cells.json")).expect("read cells"))
            .expect("parse");

    let axes = cells
        .get("organelle_axes")
        .and_then(|v| v.as_array())
        .and_then(|v| v.first())
        .and_then(|v| v.get("axes"))
        .and_then(|v| v.as_array())
        .expect("axes");

    let names = axes
        .iter()
        .filter_map(|v| v.as_str().map(ToOwned::to_owned))
        .collect::<Vec<_>>();
    assert_eq!(names, vec!["a_axis", "m_axis", "z_axis"]);
}

#[test]
fn cells_json_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        "primary_metrics.tsv",
        "barcode\taxis_a\tregime\tconfidence\tflags\nBC1\t1.0\talpha\t0.9\tf1;f2\n",
    );

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
    let first = fs::read(out.join("cells.json")).expect("read first");

    run_aggregate(&opts).expect("second");
    let second = fs::read(out.join("cells.json")).expect("read second");

    assert_eq!(first, second);

    let parsed: serde_json::Value = serde_json::from_slice(&first).expect("parse");
    assert_eq!(
        parsed.get("schema").and_then(|v| v.as_str()),
        Some("kira-organelle-cells-v1")
    );
}

fn write_tool(input_root: &std::path::Path, tool: &str, metrics_name: &str, metrics_body: &str) {
    let tool_dir = input_root.join(tool);
    fs::create_dir_all(&tool_dir).expect("mkdir tool");

    fs::write(
        tool_dir.join("summary.json"),
        r#"{"distributions":{"x":{"median":1.0}}}"#,
    )
    .expect("write summary");

    fs::write(
        tool_dir.join("pipeline_step.json"),
        format!(
            "{{\"artifacts\":{{\"primary_metrics\":\"{}\",\"summary\":\"summary.json\"}}}}",
            metrics_name
        ),
    )
    .expect("write pipeline");

    fs::write(tool_dir.join(metrics_name), metrics_body).expect("write tsv");
}
