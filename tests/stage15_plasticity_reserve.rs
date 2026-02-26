use std::fs;

use kira_organelle::cli::ComputePriArgs;
use kira_organelle::run_compute_pri_command;
use tempfile::tempdir;

#[test]
fn synthetic_component_combinations() {
    let dir = tempdir().expect("tempdir");
    let s_hi = dir.path().join("s_hi");
    let s_lo = dir.path().join("s_lo");
    let out = dir.path().join("out");

    write_sample(&s_hi, 0.90, 0.85, 0.80);
    write_sample(&s_lo, 0.10, 0.20, 0.15);

    run_compute_pri_command(&ComputePriArgs {
        inputs: vec![s_hi.clone(), s_lo.clone()],
        manifest: None,
        config: None,
        pri_weights: None,
        out: out.clone(),
    })
    .expect("compute-pri");

    let tsv = fs::read_to_string(out.join("plasticity_reserve.tsv")).expect("tsv");
    let rows = parse_rows(&tsv);
    assert!(rows["s_hi"] > rows["s_lo"]);
}

#[test]
fn missing_component_degradation_behavior() {
    let dir = tempdir().expect("tempdir");
    let s_full = dir.path().join("full");
    let s_missing = dir.path().join("missing");
    let out = dir.path().join("out");
    write_sample(&s_full, 0.7, 0.7, 0.7);
    write_sample_missing_translation(&s_missing, 0.7, 0.7);

    run_compute_pri_command(&ComputePriArgs {
        inputs: vec![s_full.clone(), s_missing.clone()],
        manifest: None,
        config: None,
        pri_weights: None,
        out: out.clone(),
    })
    .expect("compute-pri");

    let json = fs::read_to_string(out.join("plasticity_reserve.json")).expect("json");
    let v: serde_json::Value = serde_json::from_str(&json).expect("parse");
    let arr = v["samples"].as_array().expect("samples");
    let item = arr
        .iter()
        .find(|x| {
            x["sample_label"]
                .as_str()
                .unwrap_or_default()
                .contains("missing")
        })
        .expect("missing sample row");
    assert_eq!(item["low_confidence"].as_bool(), Some(true));
}

#[test]
fn stable_under_sample_reordering_with_manifest() {
    let dir = tempdir().expect("tempdir");
    let s1 = dir.path().join("s1");
    let s2 = dir.path().join("s2");
    let out = dir.path().join("out");
    write_sample(&s1, 0.2, 0.3, 0.4);
    write_sample(&s2, 0.8, 0.9, 0.7);

    let manifest = dir.path().join("manifest.csv");
    fs::write(
        &manifest,
        format!(
            "path,label,order\n{},b,1\n{},a,0\n",
            s2.display(),
            s1.display(),
        ),
    )
    .expect("manifest");

    run_compute_pri_command(&ComputePriArgs {
        inputs: vec![],
        manifest: Some(manifest),
        config: None,
        pri_weights: None,
        out: out.clone(),
    })
    .expect("compute-pri");

    let tsv = fs::read_to_string(out.join("plasticity_reserve.tsv")).expect("tsv");
    let mut lines = tsv.lines();
    let _h = lines.next();
    let first = lines.next().unwrap_or_default();
    assert!(first.starts_with("a\t0\t"));
}

#[test]
fn deterministic_across_runs() {
    let dir = tempdir().expect("tempdir");
    let s = dir.path().join("s");
    let out = dir.path().join("out");
    write_sample(&s, 0.45, 0.62, 0.57);

    let args = ComputePriArgs {
        inputs: vec![s.clone()],
        manifest: None,
        config: None,
        pri_weights: None,
        out: out.clone(),
    };
    run_compute_pri_command(&args).expect("first");
    let a = fs::read(out.join("plasticity_reserve.tsv")).expect("a");
    let b = fs::read(out.join("plasticity_reserve.json")).expect("b");
    run_compute_pri_command(&args).expect("second");
    let a2 = fs::read(out.join("plasticity_reserve.tsv")).expect("a2");
    let b2 = fs::read(out.join("plasticity_reserve.json")).expect("b2");
    assert_eq!(a, a2);
    assert_eq!(b, b2);
}

fn parse_rows(tsv: &str) -> std::collections::BTreeMap<String, f64> {
    let mut out = std::collections::BTreeMap::new();
    for (i, line) in tsv.lines().enumerate() {
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let parts = line.split('\t').collect::<Vec<_>>();
        if parts.len() < 3 {
            continue;
        }
        out.insert(parts[0].to_string(), parts[2].parse::<f64>().unwrap_or(0.0));
    }
    out
}

fn write_sample(root: &std::path::Path, plastic: f64, splice: f64, trans: f64) {
    fs::create_dir_all(root.join("kira-nuclearqc")).expect("mkdir n");
    fs::create_dir_all(root.join("kira-spliceqc")).expect("mkdir s");
    fs::create_dir_all(root.join("kira-riboqc")).expect("mkdir t");
    fs::write(
        root.join("kira-nuclearqc/summary.json"),
        format!("{{\"fraction_PlasticAdaptive\":{plastic}}}"),
    )
    .expect("n");
    fs::write(
        root.join("kira-spliceqc/summary.json"),
        format!("{{\"low_splice_noise_fraction\":{splice}}}"),
    )
    .expect("s");
    fs::write(
        root.join("kira-riboqc/summary.json"),
        format!("{{\"translation_selectivity_index\":{trans}}}"),
    )
    .expect("t");
}

fn write_sample_missing_translation(root: &std::path::Path, plastic: f64, splice: f64) {
    fs::create_dir_all(root.join("kira-nuclearqc")).expect("mkdir n");
    fs::create_dir_all(root.join("kira-spliceqc")).expect("mkdir s");
    fs::write(
        root.join("kira-nuclearqc/summary.json"),
        format!("{{\"fraction_PlasticAdaptive\":{plastic}}}"),
    )
    .expect("n");
    fs::write(
        root.join("kira-spliceqc/summary.json"),
        format!("{{\"low_splice_noise_fraction\":{splice}}}"),
    )
    .expect("s");
}
