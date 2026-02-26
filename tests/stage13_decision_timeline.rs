use std::fs;

use kira_organelle::cli::ReportDecisionTimelineArgs;
use kira_organelle::{AggregateOptions, run_aggregate, run_report_decision_timeline_command};
use tempfile::tempdir;

#[test]
fn known_sequence_metrics() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("decision");
    let out = dir.path().join("out");
    write_decision_json(
        &input,
        r#"{
          "samples":[
            {"sample_label":"s0","order_rank":0,"decision_tier":"STABLE_ADAPTIVE","confidence":0.2,"drivers":[]},
            {"sample_label":"s1","order_rank":1,"decision_tier":"TRANSITION_RISK","confidence":0.8,"drivers":[]},
            {"sample_label":"s2","order_rank":2,"decision_tier":"PRE_RESISTANT","confidence":0.9,"drivers":[]},
            {"sample_label":"s3","order_rank":3,"decision_tier":"FIXED_RESISTANT","confidence":0.9,"drivers":[]}
          ]
        }"#,
    );
    run_report_decision_timeline_command(&ReportDecisionTimelineArgs {
        input: input.clone(),
        out: out.clone(),
    })
    .expect("report");
    let tsv = fs::read_to_string(out.join("decision_stability.tsv")).expect("stability");
    assert!(tsv.contains("flip_count\t3"));
    assert!(tsv.contains("volatility\t1.000000"));
    assert!(tsv.contains("early_flip_index\t0.666667"));
}

#[test]
fn edge_cases_single_and_two_samples() {
    let dir = tempdir().expect("tempdir");
    let single = dir.path().join("single");
    let out_single = dir.path().join("out_single");
    write_decision_json(
        &single,
        r#"{"samples":[{"sample_label":"s0","order_rank":0,"decision_tier":"STABLE_ADAPTIVE","confidence":0.2,"drivers":[]}]}"#,
    );
    run_report_decision_timeline_command(&ReportDecisionTimelineArgs {
        input: single,
        out: out_single.clone(),
    })
    .expect("single");
    let s = fs::read_to_string(out_single.join("decision_stability.tsv")).expect("s");
    assert!(s.contains("flip_count\t0"));
    assert!(s.contains("volatility\t0.000000"));

    let two = dir.path().join("two");
    let out_two = dir.path().join("out_two");
    write_decision_json(
        &two,
        r#"{"samples":[
          {"sample_label":"a","order_rank":0,"decision_tier":"STABLE_ADAPTIVE","confidence":0.9,"drivers":[]},
          {"sample_label":"b","order_rank":1,"decision_tier":"FIXED_RESISTANT","confidence":0.9,"drivers":[]}
        ]}"#,
    );
    run_report_decision_timeline_command(&ReportDecisionTimelineArgs {
        input: two,
        out: out_two.clone(),
    })
    .expect("two");
    let t = fs::read_to_string(out_two.join("decision_stability.tsv")).expect("t");
    assert!(t.contains("flip_count\t1"));
    assert!(t.contains("volatility\t1.000000"));
}

#[test]
fn confidence_weighted_stability_is_correct() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("decision");
    let out = dir.path().join("out");
    write_decision_json(
        &input,
        r#"{"samples":[
          {"sample_label":"a","order_rank":0,"decision_tier":"STABLE_ADAPTIVE","confidence":1.0,"drivers":[]},
          {"sample_label":"b","order_rank":1,"decision_tier":"FIXED_RESISTANT","confidence":1.0,"drivers":[]}
        ]}"#,
    );
    run_report_decision_timeline_command(&ReportDecisionTimelineArgs {
        input,
        out: out.clone(),
    })
    .expect("run");
    let tsv = fs::read_to_string(out.join("decision_stability.tsv")).expect("tsv");
    assert!(tsv.contains("stability_score\t0.000000"));
}

#[test]
fn timeline_report_deterministic() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("decision");
    let out = dir.path().join("out");
    write_decision_json(
        &input,
        r#"{"samples":[
          {"sample_label":"a","order_rank":0,"decision_tier":"STABLE_ADAPTIVE","confidence":0.1,"drivers":[]},
          {"sample_label":"b","order_rank":1,"decision_tier":"TRANSITION_RISK","confidence":0.2,"drivers":[]}
        ]}"#,
    );
    let args = ReportDecisionTimelineArgs {
        input: input.clone(),
        out: out.clone(),
    };
    run_report_decision_timeline_command(&args).expect("first");
    let a1 = fs::read(out.join("decision_timeline.tsv")).expect("a1");
    let b1 = fs::read(out.join("decision_stability.tsv")).expect("b1");
    let c1 = fs::read(out.join("decision_stability.json")).expect("c1");
    run_report_decision_timeline_command(&args).expect("second");
    let a2 = fs::read(out.join("decision_timeline.tsv")).expect("a2");
    let b2 = fs::read(out.join("decision_stability.tsv")).expect("b2");
    let c2 = fs::read(out.join("decision_stability.json")).expect("c2");
    assert_eq!(a1, a2);
    assert_eq!(b1, b2);
    assert_eq!(c1, c2);
}

#[test]
fn backward_compatibility_not_invoked() {
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
    assert!(!out.join("decision_timeline.tsv").exists());
    assert!(!out.join("decision_stability.tsv").exists());
    assert!(!out.join("decision_stability.json").exists());
}

fn write_decision_json(dir: &std::path::Path, body: &str) {
    fs::create_dir_all(dir).expect("mkdir");
    fs::write(dir.join("sample_decision.json"), body).expect("write");
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
