use std::fs;

use kira_organelle::cli::ComputeDciArgs;
use kira_organelle::run_compute_dci_command;
use tempfile::tempdir;

#[test]
fn fully_concordant_sample() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in.json");
    let out = dir.path().join("out");
    fs::write(
        &input,
        r#"{"samples":[{"label":"s1","rank":0,"FII_nucleus":0.9,"FII_splice":0.91,"FII_proteostasis":0.95,"FII_mito":0.88,"FII_tme":0.9}]}"#,
    )
    .expect("write");
    run_compute_dci_command(&ComputeDciArgs {
        input,
        config: None,
        out: out.clone(),
    })
    .expect("run");
    let tsv = fs::read_to_string(out.join("decision_concordance.tsv")).expect("tsv");
    assert!(tsv.contains("s1\t0\t1.000000"));
}

#[test]
fn fully_conflicting_sample() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in.json");
    let out = dir.path().join("out");
    fs::write(
        &input,
        r#"{"samples":[{"label":"s1","rank":0,"FII_nucleus":0.1,"FII_splice":0.9,"FII_proteostasis":0.1,"FII_mito":0.9,"FII_tme":0.1}]}"#,
    )
    .expect("write");
    run_compute_dci_command(&ComputeDciArgs {
        input,
        config: None,
        out: out.clone(),
    })
    .expect("run");
    let d = first_dci(&fs::read_to_string(out.join("decision_concordance.tsv")).expect("tsv"));
    assert!(d < 0.6);
}

#[test]
fn mixed_partial_agreement() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in.json");
    let out = dir.path().join("out");
    fs::write(
        &input,
        r#"{"samples":[{"label":"s1","rank":0,"FII_nucleus":0.2,"FII_splice":0.5,"FII_proteostasis":0.8,"FII_mito":0.5,"FII_tme":0.2}]}"#,
    )
    .expect("write");
    run_compute_dci_command(&ComputeDciArgs {
        input,
        config: None,
        out: out.clone(),
    })
    .expect("run");
    let d = first_dci(&fs::read_to_string(out.join("decision_concordance.tsv")).expect("tsv"));
    assert!(d > 0.2 && d < 0.9);
}

#[test]
fn missing_organelle_robustness() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in.json");
    let out = dir.path().join("out");
    fs::write(
        &input,
        r#"{"samples":[{"label":"s1","rank":0,"FII_nucleus":0.2,"FII_proteostasis":0.8,"FII_tme":0.2}]}"#,
    )
    .expect("write");
    run_compute_dci_command(&ComputeDciArgs {
        input,
        config: None,
        out: out.clone(),
    })
    .expect("run");
    let d = first_dci(&fs::read_to_string(out.join("decision_concordance.tsv")).expect("tsv"));
    assert!((0.0..=1.0).contains(&d));
}

#[test]
fn deterministic_across_runs() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in.json");
    let out = dir.path().join("out");
    fs::write(
        &input,
        r#"{"samples":[{"label":"s1","rank":0,"FII_nucleus":0.2,"FII_splice":0.5,"FII_proteostasis":0.8,"FII_mito":0.4,"FII_tme":0.9}]}"#,
    )
    .expect("write");
    let args = ComputeDciArgs {
        input,
        config: None,
        out: out.clone(),
    };
    run_compute_dci_command(&args).expect("first");
    let a = fs::read(out.join("decision_concordance.tsv")).expect("a");
    let b = fs::read(out.join("decision_concordance.json")).expect("b");
    run_compute_dci_command(&args).expect("second");
    let a2 = fs::read(out.join("decision_concordance.tsv")).expect("a2");
    let b2 = fs::read(out.join("decision_concordance.json")).expect("b2");
    assert_eq!(a, a2);
    assert_eq!(b, b2);
}

fn first_dci(tsv: &str) -> f64 {
    for (i, line) in tsv.lines().enumerate() {
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let p = line.split('\t').collect::<Vec<_>>();
        if p.len() >= 3 {
            return p[2].parse::<f64>().unwrap_or(0.0);
        }
    }
    0.0
}
