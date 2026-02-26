use std::fs;

use kira_organelle::cli::DecisionArgs;
use kira_organelle::{AggregateOptions, run_aggregate, run_decision_command};
use tempfile::tempdir;

#[test]
fn decision_tier_assignment_correctness() {
    let dir = tempdir().expect("tempdir");
    let out = dir.path().join("out");
    write_phase_summary(
        &out,
        r#"{
          "samples":[
            {"label":"baseline","rank":0,"fii_mean":0.22,"fii_median":0.21,"velocity":0.01,"acceleration":0.00},
            {"label":"c2","rank":1,"fii_mean":0.55,"fii_median":0.54,"velocity":0.16,"acceleration":0.11},
            {"label":"eot","rank":2,"fii_mean":0.72,"fii_median":0.71,"velocity":0.05,"acceleration":-0.02}
          ]
        }"#,
    );
    write_flags(
        &out,
        "sample_label\torder_rank\tflag\tmetric\tvalue\tthreshold\tnote\nc2\t1\tVELOCITY_SPIKE\tfii_velocity\t0.16\t0.10\tn\nc2\t1\tACCELERATION_PEAK\tfii_acceleration\t0.11\t0.08\tn\n",
    );

    run_decision_command(&DecisionArgs {
        input: out.join("phase_portrait_summary.json"),
        config: None,
        out: out.clone(),
    })
    .expect("decision");

    let tsv = fs::read_to_string(out.join("sample_decision.tsv")).expect("tsv");
    assert!(tsv.contains("baseline\t0\tSTABLE_ADAPTIVE"));
    assert!(tsv.contains("c2\t1\tPRE_RESISTANT"));
    assert!(tsv.contains("eot\t2\tFIXED_RESISTANT"));
}

#[test]
fn decision_priority_pre_resistant_over_transition() {
    let dir = tempdir().expect("tempdir");
    let out = dir.path().join("out");
    write_phase_summary(
        &out,
        r#"{"samples":[{"label":"s1","rank":1,"fii_mean":0.50,"fii_median":0.50,"velocity":0.17,"acceleration":0.12}]}"#,
    );
    write_flags(
        &out,
        "sample_label\torder_rank\tflag\tmetric\tvalue\tthreshold\tnote\ns1\t1\tVELOCITY_SPIKE\tfii_velocity\t0.17\t0.10\tn\n",
    );
    run_decision_command(&DecisionArgs {
        input: out.join("phase_portrait_summary.json"),
        config: None,
        out: out.clone(),
    })
    .expect("decision");
    let tsv = fs::read_to_string(out.join("sample_decision.tsv")).expect("tsv");
    assert!(tsv.contains("s1\t1\tPRE_RESISTANT"));
}

#[test]
fn decision_confidence_deterministic() {
    let dir = tempdir().expect("tempdir");
    let out = dir.path().join("out");
    write_phase_summary(
        &out,
        r#"{"samples":[{"label":"x","rank":0,"fii_mean":0.40,"fii_median":0.39,"velocity":0.10,"acceleration":0.10}]}"#,
    );
    write_flags(
        &out,
        "sample_label\torder_rank\tflag\tmetric\tvalue\tthreshold\tnote\nx\t0\tACCELERATION_PEAK\tfii_acceleration\t0.10\t0.08\tn\n",
    );
    let args = DecisionArgs {
        input: out.join("phase_portrait_summary.json"),
        config: None,
        out: out.clone(),
    };
    run_decision_command(&args).expect("first");
    let a = fs::read(out.join("sample_decision.tsv")).expect("a");
    run_decision_command(&args).expect("second");
    let b = fs::read(out.join("sample_decision.tsv")).expect("b");
    assert_eq!(a, b);
}

#[test]
fn decision_threshold_override_behavior() {
    let dir = tempdir().expect("tempdir");
    let out = dir.path().join("out");
    write_phase_summary(
        &out,
        r#"{"samples":[{"label":"s","rank":0,"fii_mean":0.58,"fii_median":0.58,"velocity":0.09,"acceleration":0.04}]}"#,
    );
    write_flags(
        &out,
        "sample_label\torder_rank\tflag\tmetric\tvalue\tthreshold\tnote\n",
    );
    let cfg = out.join("decision_config.json");
    fs::write(
        &cfg,
        r#"{"stable_threshold":0.3,"resistant_threshold":0.55,"velocity_mid_threshold":0.08}"#,
    )
    .expect("cfg");

    run_decision_command(&DecisionArgs {
        input: out.join("phase_portrait_summary.json"),
        config: Some(cfg),
        out: out.clone(),
    })
    .expect("decision");
    let tsv = fs::read_to_string(out.join("sample_decision.tsv")).expect("tsv");
    assert!(tsv.contains("s\t0\tFIXED_RESISTANT"));
}

#[test]
fn decision_missing_input_degradation() {
    let dir = tempdir().expect("tempdir");
    let out = dir.path().join("out");
    write_phase_summary(&out, r#"{"samples":[{"label":"m","rank":0}]}"#);
    write_flags(
        &out,
        "sample_label\torder_rank\tflag\tmetric\tvalue\tthreshold\tnote\n",
    );
    run_decision_command(&DecisionArgs {
        input: out.join("phase_portrait_summary.json"),
        config: None,
        out: out.clone(),
    })
    .expect("decision");
    let tsv = fs::read_to_string(out.join("sample_decision.tsv")).expect("tsv");
    assert!(tsv.contains("m\t0\tTRANSITION_RISK") || tsv.contains("m\t0\tSTABLE_ADAPTIVE"));
}

#[test]
fn backward_compat_without_decision_invoked() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in");
    let out = dir.path().join("out");
    write_tool(
        &input,
        "kira-mitoqc",
        "barcode\tmitochondrial_stress_adaptation_score\nBC1\t0.3\n",
    );
    write_tool(
        &input,
        "kira-riboqc",
        "barcode\ttranslation_commitment_score\nBC1\t0.4\n",
    );
    write_tool(
        &input,
        "kira-spliceqc",
        "barcode\tsplice_irreversibility_index\nBC1\t0.5\n",
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
    assert!(!out.join("sample_decision.tsv").exists());
    assert!(!out.join("sample_decision.json").exists());
}

fn write_phase_summary(out: &std::path::Path, json_body: &str) {
    fs::create_dir_all(out).expect("mkdir");
    fs::write(out.join("phase_portrait_summary.json"), json_body).expect("summary");
}

fn write_flags(out: &std::path::Path, body: &str) {
    fs::write(out.join("early_warning_flags.tsv"), body).expect("flags");
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
