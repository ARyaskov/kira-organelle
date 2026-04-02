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
        "barcode\tmito_axis\tRSS\tSII\tOSL\tTSM\tPCP\tASM\tMSM\tTPI\tCCI\tPCI\tLCI\tAPB\tIMSI\tMCB\tHSI\nC1\t0.7\t1.0\t0.2\t0.9\t0.8\t1.2\t0.7\t0.9\t1.1\t0.6\t0.5\t0.7\t0.4\t0.8\t0.6\t0.9\nC2\t0.5\t0.4\t0.1\t0.3\t0.2\t0.4\t0.2\t0.3\t0.4\t0.5\t0.6\t0.5\t0.6\t0.2\t0.5\t0.3\n",
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
        export_systems_model: Some(out.join("systems_model.json")),
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
    let model = summary
        .get("model")
        .and_then(|v| v.as_object())
        .expect("model");
    assert_eq!(
        model
            .get("model_version")
            .and_then(|v| v.as_str())
            .expect("model_version"),
        "2.0.0"
    );

    let expression =
        fs::read_to_string(out.join("integration/expression_aggregated.tsv")).expect("expression");
    assert!(expression.contains("mitochondria::mito_axis"));
    assert!(expression.contains("proteostasis::proteo_axis"));

    let metrics = fs::read_to_string(out.join("integration/metrics.tsv")).expect("metrics");
    let mut metric_lines = metrics.lines();
    assert_eq!(
        metric_lines.next().expect("metrics header"),
        "cell_id\tcluster\tStressVector\tCompensationDeficit\tCPI\tAFS\tIMSC\tRegimeClass\tRareState\tMahalanobisDistance\tPotential\tStabilityGradient\tTPI_landscape\tBasinId\tTransitionCandidate\tDominantVulnerableAxis\tLSI\tMaxTrajectory\tBEE"
    );
    assert!(metric_lines.next().is_some());

    let normalized = fs::read_to_string(out.join("integration/metrics_normalized.tsv"))
        .expect("metrics_normalized");
    assert!(normalized.contains("\tCCI\t"));

    let landscape =
        fs::read_to_string(out.join("integration/landscape.html")).expect("landscape html");
    assert!(landscape.contains("Systems Landscape And Transition Geometry"));
    assert!(landscape.contains("Regime Transition Matrix"));

    let integration_systems_model: serde_json::Value = serde_json::from_slice(
        &fs::read(out.join("integration/systems_model.json")).expect("integration systems model"),
    )
    .expect("integration systems model parse");
    assert!(integration_systems_model.get("state_graph").is_some());
    assert!(
        integration_systems_model
            .get("state_transition_matrix")
            .is_some()
    );

    let systems = summary
        .get("systems_state_model")
        .and_then(|v| v.as_object())
        .expect("systems_state_model");
    assert!(systems.get("global_stats").is_some());
    assert!(systems.get("cluster_stats").is_some());
    assert!(systems.get("regime_distribution").is_some());
    assert!(systems.get("global_cpi_p50").is_some());
    assert!(systems.get("global_cpi_p90").is_some());
    assert!(systems.get("fragility").is_some());
    assert!(systems.get("therapeutic_projection").is_some());
    assert!(systems.get("dynamic_stability").is_some());
    assert!(systems.get("compatibility").is_some());
    assert!(systems.get("landscape").is_some());

    let normalization = summary
        .get("global_normalization")
        .and_then(|v| v.as_object())
        .expect("global_normalization");
    assert_eq!(
        normalization
            .get("method")
            .and_then(|v| v.as_str())
            .expect("method"),
        "median_mad"
    );
    assert!(
        normalization
            .get("per_metric")
            .and_then(|v| v.as_array())
            .map(|arr| !arr.is_empty())
            .unwrap_or(false)
    );

    let ingestion = summary
        .get("ingestion_diagnostics")
        .and_then(|v| v.as_array())
        .expect("ingestion_diagnostics");
    assert_eq!(ingestion.len(), 2);
    assert!(
        ingestion
            .iter()
            .any(|x| x.get("tool").and_then(|v| v.as_str()) == Some("kira-mitoqc"))
    );

    let formal_model: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("systems_model.json")).expect("systems model"))
            .expect("systems model parse");
    assert_eq!(
        formal_model
            .get("model_version")
            .and_then(|v| v.as_str())
            .expect("model_version"),
        "2.0.0"
    );
    assert!(formal_model.get("state_graph").is_some());
    assert!(formal_model.get("ode_proxy").is_some());
    assert!(formal_model.get("state_transition_matrix").is_some());
    assert!(formal_model.get("basin_connectivity").is_some());
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
