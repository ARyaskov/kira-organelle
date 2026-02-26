use std::fs;
use std::path::Path;

use kira_organelle::cli::ReportFiiArgs;
use kira_organelle::error::{classify_error_message, stable_error_code};
use kira_organelle::run_report_fii_command;
use tempfile::tempdir;

#[test]
fn report_fii_deterministic() {
    let dir = tempdir().expect("tempdir");
    let samples = dir.path().join("samples");
    let out = dir.path().join("out");
    create_sample(&samples.join("baseline"), 0.21, 0.30, 0.49);
    create_sample(&samples.join("eot"), 0.72, 0.18, 0.10);

    let args = ReportFiiArgs {
        inputs: vec![samples.join("baseline"), samples.join("eot")],
        out: out.clone(),
        manifest: None,
        title: Some("FII Determinism".to_string()),
    };
    run_report_fii_command(&args).expect("first");
    let first = fs::read(out.join("fii_landscape.html")).expect("read first");
    run_report_fii_command(&args).expect("second");
    let second = fs::read(out.join("fii_landscape.html")).expect("read second");
    assert_eq!(first, second);
}

#[test]
fn report_fii_golden_sections_present() {
    let dir = tempdir().expect("tempdir");
    let samples = dir.path().join("samples");
    let out = dir.path().join("out");
    create_sample(&samples.join("baseline"), 0.21, 0.30, 0.49);
    create_sample(&samples.join("c2"), 0.38, 0.31, 0.31);

    run_report_fii_command(&ReportFiiArgs {
        inputs: vec![samples.join("baseline"), samples.join("c2")],
        out: out.clone(),
        manifest: None,
        title: Some("Golden".to_string()),
    })
    .expect("report");

    let html = fs::read_to_string(out.join("fii_landscape.html")).expect("html");
    assert!(html.contains("<h2>Distribution View</h2>"));
    assert!(html.contains("<h2>Regime Composition</h2>"));
    assert!(html.contains("<h2>Component Coupling View</h2>"));
    assert!(html.contains("<h2>Sample Trajectory Summary</h2>"));
    assert!(html.contains("\"samples\":["));
}

#[test]
fn report_fii_manifest_ordering_is_respected() {
    let dir = tempdir().expect("tempdir");
    let samples = dir.path().join("samples");
    let out = dir.path().join("out");
    create_sample(&samples.join("s1"), 0.20, 0.30, 0.50);
    create_sample(&samples.join("s2"), 0.80, 0.10, 0.10);

    let manifest = dir.path().join("manifest.tsv");
    fs::write(
        &manifest,
        format!(
            "path\tlabel\torder\ttimepoint\n{}\tEOT\t2\tEOT\n{}\tBaseline\t1\tBL\n",
            samples.join("s2").display(),
            samples.join("s1").display()
        ),
    )
    .expect("manifest");

    run_report_fii_command(&ReportFiiArgs {
        inputs: Vec::new(),
        out: out.clone(),
        manifest: Some(manifest),
        title: None,
    })
    .expect("report");

    let index: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("fii_landscape_index.json")).expect("index"))
            .expect("parse");
    let labels = index
        .get("samples")
        .and_then(|v| v.as_array())
        .expect("samples")
        .iter()
        .filter_map(|v| v.get("label").and_then(|s| s.as_str()))
        .collect::<Vec<_>>();
    assert_eq!(labels, vec!["Baseline", "EOT"]);
}

#[test]
fn report_fii_missing_file_has_stable_error_code() {
    let dir = tempdir().expect("tempdir");
    let sample = dir.path().join("sample");
    fs::create_dir_all(&sample).expect("mkdir");

    let err = run_report_fii_command(&ReportFiiArgs {
        inputs: vec![sample],
        out: dir.path().join("out"),
        manifest: None,
        title: None,
    })
    .expect_err("must fail");
    let code = stable_error_code(classify_error_message(&err));
    assert_eq!(code, "E_CONTRACT");
}

fn create_sample(root: &Path, fii: f64, trans: f64, splice: f64) {
    fs::create_dir_all(root).expect("mkdir");
    fs::write(
        root.join("functional_irreversibility_index.tsv"),
        format!(
            "cell_id\tmitochondrial_stress_adaptation_score\ttranslation_commitment_score\tsplice_irreversibility_index\tfunctional_irreversibility_index\tfii_regime\tfii_low_confidence\n\
            c1\t{fii:.3}\t{trans:.3}\t{splice:.3}\t{fii:.3}\tAdaptive\tfalse\n\
            c2\t{fii:.3}\t{trans:.3}\t{splice:.3}\t{fii:.3}\tTransition\tfalse\n\
            c3\t{fii:.3}\t{trans:.3}\t{splice:.3}\t{fii:.3}\tResistant\ttrue\n"
        ),
    )
    .expect("tsv");
    fs::write(
        root.join("state.json"),
        format!(
            "{{\"functional_irreversibility\":{{\"weights\":{{\"mitochondrial\":0.333333,\"translation\":0.333333,\"splice\":0.333333}},\"distribution\":{{\"mean\":{fii:.6},\"median\":{fii:.6},\"p25\":{fii:.6},\"p75\":{fii:.6}}},\"regime_fractions\":{{\"adaptive\":0.333333,\"transition\":0.333333,\"resistant\":0.333333}},\"low_confidence_fraction\":0.333333}}}}"
        ),
    )
    .expect("state");
}
