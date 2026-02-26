use std::fs;

use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn aggregate_writes_integration_artifacts() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");

    write_tool(
        &input,
        "kira-mitoqc",
        r#"{"distributions":{"mito_axis":{"median":0.70}}}"#,
        "barcode\tmito_axis\nC1\t0.7\nC2\t0.5\n",
    );
    write_tool(
        &input,
        "kira-proteoqc",
        r#"{"distributions":{"proteo_axis":{"median":0.40}}}"#,
        "barcode\tproteo_axis\nC1\t0.4\nC2\t0.3\n",
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

    let timeseries =
        fs::read_to_string(out.join("integration/timeseries.tsv")).expect("timeseries");
    let mut lines = timeseries.lines();
    assert_eq!(
        lines.next().expect("header"),
        "entity_id\ttimepoint\tmito\tproteostasis\tsplice\tsecretion\tenergetics\tautophagy"
    );
    assert!(lines.next().is_some());

    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("integration/summary.json")).expect("summary"))
            .expect("summary parse");
    assert_eq!(
        summary
            .get("schema")
            .and_then(|v| v.as_str())
            .expect("schema"),
        "kira-organelle-integration-v1"
    );
    assert_eq!(
        summary
            .get("expression_available")
            .and_then(|v| v.as_bool())
            .expect("expression_available"),
        true
    );

    let expression =
        fs::read_to_string(out.join("integration/expression_aggregated.tsv")).expect("expression");
    assert!(expression.contains("mitochondria::mito_axis"));
    assert!(expression.contains("proteostasis::proteo_axis"));
}

fn write_tool(input: &std::path::Path, tool: &str, summary_json: &str, primary_tsv: &str) {
    let tool_dir = input.join(tool);
    fs::create_dir_all(&tool_dir).expect("mkdir tool");
    fs::write(tool_dir.join("summary.json"), summary_json).expect("summary");
    fs::write(tool_dir.join("primary_metrics.tsv"), primary_tsv).expect("metrics");
    fs::write(
        tool_dir.join("pipeline_step.json"),
        r#"{"artifacts":{"summary":"summary.json","primary_metrics":"primary_metrics.tsv"}}"#,
    )
    .expect("pipeline");
}
