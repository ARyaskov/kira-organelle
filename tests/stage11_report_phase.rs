use std::fs;
use std::path::Path;

use kira_organelle::cli::ReportPhaseArgs;
use kira_organelle::{AggregateOptions, run_aggregate, run_report_phase_command};
use tempfile::tempdir;

#[test]
fn phase_binning_contains_all_three_projections() {
    let dir = tempdir().expect("tempdir");
    let samples = dir.path().join("samples");
    let out = dir.path().join("out");
    sample(
        &samples.join("s0"),
        &[0.2, 0.3],
        &["Adaptive", "Transition"],
        &[false, false],
    );
    sample(
        &samples.join("s1"),
        &[0.4, 0.5],
        &["Transition", "Transition"],
        &[false, true],
    );
    sample(
        &samples.join("s2"),
        &[0.8, 0.9],
        &["Resistant", "Resistant"],
        &[false, false],
    );

    run_report_phase_command(&ReportPhaseArgs {
        inputs: vec![samples.join("s0"), samples.join("s1"), samples.join("s2")],
        manifest: None,
        config: None,
        out: out.clone(),
    })
    .expect("report");

    let bins = fs::read_to_string(out.join("phase_portrait_bins.tsv")).expect("bins");
    assert!(bins.contains("fii_vs_velocity"));
    assert!(bins.contains("fii_vs_acceleration"));
    assert!(bins.contains("velocity_vs_acceleration"));
    assert!(bins.contains("\tALL\t"));
}

#[test]
fn phase_flags_trigger_all_types() {
    let dir = tempdir().expect("tempdir");
    let samples = dir.path().join("samples");
    let out = dir.path().join("out");

    sample(
        &samples.join("b0"),
        &[0.20, 0.21],
        &["Adaptive", "Adaptive"],
        &[false, false],
    );
    sample(
        &samples.join("b1"),
        &[0.32, 0.33],
        &["Transition", "Transition"],
        &[false, false],
    );
    sample(
        &samples.join("b2"),
        &[0.55, 0.56],
        &["Transition", "Resistant"],
        &[false, false],
    );
    sample(
        &samples.join("b3"),
        &[0.42, 0.43],
        &["Transition", "Transition"],
        &[false, false],
    );
    sample(
        &samples.join("b4"),
        &[0.60, 0.62],
        &["Resistant", "Resistant"],
        &[false, false],
    );

    run_report_phase_command(&ReportPhaseArgs {
        inputs: vec![
            samples.join("b0"),
            samples.join("b1"),
            samples.join("b2"),
            samples.join("b3"),
            samples.join("b4"),
        ],
        manifest: None,
        config: None,
        out: out.clone(),
    })
    .expect("report");

    let flags = fs::read_to_string(out.join("early_warning_flags.tsv")).expect("flags");
    assert!(flags.contains("VELOCITY_SPIKE"));
    assert!(flags.contains("ACCELERATION_PEAK"));
    assert!(flags.contains("PHASE_DRIFT_TO_RESISTANT"));
    assert!(flags.contains("TRANSITION_VOLATILITY"));
}

#[test]
fn phase_deterministic_outputs() {
    let dir = tempdir().expect("tempdir");
    let samples = dir.path().join("samples");
    let out = dir.path().join("out");
    sample(
        &samples.join("s0"),
        &[0.2, 0.3],
        &["Adaptive", "Transition"],
        &[false, false],
    );
    sample(
        &samples.join("s1"),
        &[0.4, 0.5],
        &["Transition", "Transition"],
        &[false, true],
    );
    sample(
        &samples.join("s2"),
        &[0.8, 0.9],
        &["Resistant", "Resistant"],
        &[false, false],
    );

    let args = ReportPhaseArgs {
        inputs: vec![samples.join("s0"), samples.join("s1"), samples.join("s2")],
        manifest: None,
        config: None,
        out: out.clone(),
    };
    run_report_phase_command(&args).expect("first");
    let a1 = fs::read(out.join("phase_portrait_points.tsv")).expect("a1");
    let b1 = fs::read(out.join("phase_portrait_bins.tsv")).expect("b1");
    let c1 = fs::read(out.join("phase_portrait_summary.json")).expect("c1");
    let d1 = fs::read(out.join("early_warning_flags.tsv")).expect("d1");
    run_report_phase_command(&args).expect("second");
    let a2 = fs::read(out.join("phase_portrait_points.tsv")).expect("a2");
    let b2 = fs::read(out.join("phase_portrait_bins.tsv")).expect("b2");
    let c2 = fs::read(out.join("phase_portrait_summary.json")).expect("c2");
    let d2 = fs::read(out.join("early_warning_flags.tsv")).expect("d2");
    assert_eq!(a1, a2);
    assert_eq!(b1, b2);
    assert_eq!(c1, c2);
    assert_eq!(d1, d2);
}

#[test]
fn phase_manifest_validation_errors() {
    let dir = tempdir().expect("tempdir");
    let samples = dir.path().join("samples");
    sample(&samples.join("a"), &[0.2], &["Adaptive"], &[false]);
    sample(&samples.join("b"), &[0.3], &["Transition"], &[false]);
    sample(&samples.join("c"), &[0.4], &["Transition"], &[false]);

    let dup = dir.path().join("dup.csv");
    fs::write(
        &dup,
        format!(
            "sample_label,sample_dir,order_rank\nA,{},0\nB,{},0\n",
            samples.join("a").display(),
            samples.join("b").display()
        ),
    )
    .expect("dup");
    let err = run_report_phase_command(&ReportPhaseArgs {
        inputs: Vec::new(),
        manifest: Some(dup),
        config: None,
        out: dir.path().join("out_dup"),
    })
    .expect_err("dup fail");
    assert!(err.contains("non-monotonic") || err.contains("duplicated"));

    let missing = dir.path().join("missing.csv");
    fs::write(
        &missing,
        format!(
            "sample_label,sample_dir,order_rank\nA,{},0\nC,{},2\n",
            samples.join("a").display(),
            samples.join("c").display()
        ),
    )
    .expect("missing");
    let err = run_report_phase_command(&ReportPhaseArgs {
        inputs: Vec::new(),
        manifest: Some(missing),
        config: None,
        out: dir.path().join("out_missing"),
    })
    .expect_err("missing fail");
    assert!(err.contains("missing/non-monotonic") || err.contains("non-monotonic"));

    let nonmono = dir.path().join("nonmono.csv");
    fs::write(
        &nonmono,
        format!(
            "sample_label,sample_dir,order_rank\nA,{},0\nC,{},2\nB,{},1\n",
            samples.join("a").display(),
            samples.join("c").display(),
            samples.join("b").display()
        ),
    )
    .expect("nonmono");
    let err = run_report_phase_command(&ReportPhaseArgs {
        inputs: Vec::new(),
        manifest: Some(nonmono),
        config: None,
        out: dir.path().join("out_nonmono"),
    })
    .expect_err("nonmono fail");
    assert!(err.contains("non-monotonic"));
}

#[test]
fn phase_backward_compat_when_not_invoked() {
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
    assert!(!out.join("phase_portrait_points.tsv").exists());
    assert!(!out.join("phase_portrait_bins.tsv").exists());
    assert!(!out.join("phase_portrait_summary.json").exists());
    assert!(!out.join("early_warning_flags.tsv").exists());
}

fn sample(root: &Path, fii: &[f64], regimes: &[&str], low_conf: &[bool]) {
    fs::create_dir_all(root).expect("mkdir");
    let mut tsv = String::from(
        "cell_id\tmitochondrial_stress_adaptation_score\ttranslation_commitment_score\tsplice_irreversibility_index\tfunctional_irreversibility_index\tfii_regime\tfii_low_confidence\n",
    );
    for i in 0..fii.len() {
        tsv.push_str(&format!(
            "c{}\t0.4\t0.4\t0.4\t{:.3}\t{}\t{}\n",
            i + 1,
            fii[i],
            regimes.get(i).copied().unwrap_or("Transition"),
            if *low_conf.get(i).unwrap_or(&false) {
                "true"
            } else {
                "false"
            }
        ));
    }
    fs::write(root.join("functional_irreversibility_index.tsv"), tsv).expect("tsv");
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
