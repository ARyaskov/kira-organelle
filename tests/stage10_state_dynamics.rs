use std::fs;
use std::path::Path;

use kira_organelle::cli::ComputeStateDynamicsArgs;
use kira_organelle::{AggregateOptions, run_aggregate, run_compute_state_dynamics_command};
use tempfile::tempdir;

#[test]
fn dynamics_deterministic_for_fixed_fixtures() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("samples");
    let out = dir.path().join("out");
    create_sample(&input.join("baseline"), &[0.2, 0.3, 0.4], 0.2, 0.3, 0.4);
    create_sample(&input.join("c2"), &[0.4, 0.5, 0.6], 0.4, 0.5, 0.6);
    create_sample(&input.join("eot"), &[0.7, 0.8, 0.9], 0.7, 0.8, 0.9);

    let args = ComputeStateDynamicsArgs {
        inputs: vec![input.join("baseline"), input.join("c2"), input.join("eot")],
        manifest: None,
        out: out.clone(),
    };
    run_compute_state_dynamics_command(&args).expect("first");
    let v1 = fs::read(out.join("fii_state_velocity.tsv")).expect("v1");
    let a1 = fs::read(out.join("fii_state_acceleration.tsv")).expect("a1");
    run_compute_state_dynamics_command(&args).expect("second");
    let v2 = fs::read(out.join("fii_state_velocity.tsv")).expect("v2");
    let a2 = fs::read(out.join("fii_state_acceleration.tsv")).expect("a2");
    assert_eq!(v1, v2);
    assert_eq!(a1, a2);
}

#[test]
fn dynamics_two_vs_three_samples_behavior() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("samples");
    let out2 = dir.path().join("out2");
    let out3 = dir.path().join("out3");
    create_sample(&input.join("s1"), &[0.1, 0.1], 0.1, 0.1, 0.1);
    create_sample(&input.join("s2"), &[0.4, 0.4], 0.4, 0.4, 0.4);
    create_sample(&input.join("s3"), &[0.9, 0.9], 0.9, 0.9, 0.9);

    run_compute_state_dynamics_command(&ComputeStateDynamicsArgs {
        inputs: vec![input.join("s1"), input.join("s2")],
        manifest: None,
        out: out2.clone(),
    })
    .expect("run2");
    let acc2 = fs::read_to_string(out2.join("fii_state_acceleration.tsv")).expect("acc2");
    assert!(acc2.contains("s1\t0\t0.000000"));
    assert!(acc2.contains("s2\t1\t0.000000"));

    run_compute_state_dynamics_command(&ComputeStateDynamicsArgs {
        inputs: vec![input.join("s1"), input.join("s2"), input.join("s3")],
        manifest: None,
        out: out3.clone(),
    })
    .expect("run3");
    let acc3 = fs::read_to_string(out3.join("fii_state_acceleration.tsv")).expect("acc3");
    assert!(acc3.contains("s3\t2\t0.200000"));
}

#[test]
fn dynamics_manifest_validation_failures() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("samples");
    create_sample(&input.join("a"), &[0.2], 0.2, 0.2, 0.2);
    create_sample(&input.join("b"), &[0.3], 0.3, 0.3, 0.3);

    let manifest = dir.path().join("manifest.csv");
    fs::write(
        &manifest,
        format!(
            "sample_label,sample_dir,order_rank\nA,{},0\nB,{},0\n",
            input.join("a").display(),
            input.join("b").display()
        ),
    )
    .expect("manifest");

    let err = run_compute_state_dynamics_command(&ComputeStateDynamicsArgs {
        inputs: Vec::new(),
        manifest: Some(manifest),
        out: dir.path().join("out"),
    })
    .expect_err("must fail");
    assert!(err.contains("duplicated order_rank"));
}

#[test]
fn dynamics_missing_data_propagation() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("samples");
    let out = dir.path().join("out");
    create_sample(&input.join("s1"), &[0.2, 0.2], 0.2, 0.2, 0.2);
    create_invalid_sample(&input.join("s2"));

    run_compute_state_dynamics_command(&ComputeStateDynamicsArgs {
        inputs: vec![input.join("s1"), input.join("s2")],
        manifest: None,
        out: out.clone(),
    })
    .expect("run");
    let vel = fs::read_to_string(out.join("fii_state_velocity.tsv")).expect("vel");
    assert!(vel.contains("s2\t1\t0.000000\t0.000000\t-0.200000"));
    assert!(vel.contains("\ttrue\n"));
}

#[test]
fn backward_compatibility_without_dynamics_is_unchanged() {
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
        export_systems_model: None,
    })
    .expect("aggregate");
    assert!(!out.join("fii_state_velocity.tsv").exists());
    assert!(!out.join("fii_state_acceleration.tsv").exists());
}

fn create_sample(root: &Path, fii_values: &[f64], mito: f64, tr: f64, sp: f64) {
    fs::create_dir_all(root).expect("mkdir");
    let mut tsv = String::from(
        "cell_id\tmitochondrial_stress_adaptation_score\ttranslation_commitment_score\tsplice_irreversibility_index\tfunctional_irreversibility_index\tfii_regime\tfii_low_confidence\n",
    );
    for (i, v) in fii_values.iter().enumerate() {
        tsv.push_str(&format!(
            "c{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\tTransition\tfalse\n",
            i + 1,
            mito,
            tr,
            sp,
            v
        ));
    }
    fs::write(root.join("functional_irreversibility_index.tsv"), tsv).expect("tsv");
    fs::write(root.join("summary.json"), "{}").expect("summary");
}

fn create_invalid_sample(root: &Path) {
    fs::create_dir_all(root).expect("mkdir");
    fs::write(
        root.join("functional_irreversibility_index.tsv"),
        "cell_id\tfunctional_irreversibility_index\tfii_regime\tfii_low_confidence\nc1\tNaN\tTransition\ttrue\n",
    )
    .expect("tsv");
    fs::write(root.join("summary.json"), "{}").expect("summary");
}

fn write_tool(input_root: &Path, tool: &str, metrics_body: &str) {
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
