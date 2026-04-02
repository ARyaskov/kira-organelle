use std::fs;

use kira_organelle::cli::ComputeCaiArgs;
use kira_organelle::{AggregateOptions, run_aggregate, run_compute_cai_command};
use tempfile::tempdir;

#[test]
fn synthetic_distributions_uniform_skewed_bimodal() {
    let dir = tempdir().expect("tempdir");
    let s_uniform = dir.path().join("uniform");
    let s_skewed = dir.path().join("skewed");
    let s_bimodal = dir.path().join("bimodal");
    let out = dir.path().join("out");

    write_fii(&s_uniform, &linspace(1000, 0.2, 0.8));
    let mut skewed = vec![0.05; 900];
    skewed.extend(vec![0.95; 100]);
    write_fii(&s_skewed, &skewed);
    let mut bimodal = vec![0.2; 500];
    bimodal.extend(vec![0.85; 500]);
    write_fii(&s_bimodal, &bimodal);

    run_compute_cai_command(&ComputeCaiArgs {
        inputs: vec![s_uniform.clone(), s_skewed.clone(), s_bimodal.clone()],
        manifest: None,
        config: None,
        cai_weights: None,
        out: out.clone(),
    })
    .expect("compute-cai");

    let tsv = fs::read_to_string(out.join("commitment_asymmetry.tsv")).expect("tsv");
    assert!(tsv.contains("sample_label\torder_rank\tCAI"));
    assert!(tsv.contains("uniform\t0"));
    assert!(tsv.contains("skewed\t1"));
    assert!(tsv.contains("bimodal\t2"));

    let rows = parse_tsv_rows(&tsv);
    let u = rows["uniform"].0;
    let s = rows["skewed"].0;
    assert!(s > u);
}

#[test]
fn sensitive_to_small_hard_tail() {
    let dir = tempdir().expect("tempdir");
    let s_diffuse = dir.path().join("diffuse");
    let s_tail = dir.path().join("tail");
    let out = dir.path().join("out");

    write_fii(&s_diffuse, &vec![0.55; 1000]);
    let mut with_tail = vec![0.55; 980];
    with_tail.extend(vec![0.95; 20]);
    write_fii(&s_tail, &with_tail);

    run_compute_cai_command(&ComputeCaiArgs {
        inputs: vec![s_diffuse.clone(), s_tail.clone()],
        manifest: None,
        config: None,
        cai_weights: None,
        out: out.clone(),
    })
    .expect("compute-cai");
    let tsv = fs::read_to_string(out.join("commitment_asymmetry.tsv")).expect("tsv");
    let rows = parse_tsv_rows(&tsv);
    assert!(rows["tail"].0 > rows["diffuse"].0);
}

#[test]
fn stable_under_cell_count_scaling() {
    let dir = tempdir().expect("tempdir");
    let s_small = dir.path().join("small");
    let s_scaled = dir.path().join("scaled");
    let out = dir.path().join("out");

    let base = vec![0.2, 0.2, 0.3, 0.4, 0.9];
    write_fii(&s_small, &base);
    let mut expanded = Vec::new();
    for _ in 0..200 {
        expanded.extend(base.iter().copied());
    }
    write_fii(&s_scaled, &expanded);

    run_compute_cai_command(&ComputeCaiArgs {
        inputs: vec![s_small.clone(), s_scaled.clone()],
        manifest: None,
        config: None,
        cai_weights: None,
        out: out.clone(),
    })
    .expect("compute-cai");
    let tsv = fs::read_to_string(out.join("commitment_asymmetry.tsv")).expect("tsv");
    let rows = parse_tsv_rows(&tsv);
    assert!((rows["small"].0 - rows["scaled"].0).abs() < 1e-9);
}

#[test]
fn deterministic_across_runs() {
    let dir = tempdir().expect("tempdir");
    let s = dir.path().join("sample");
    let out = dir.path().join("out");
    let mut values = vec![0.1; 700];
    values.extend(vec![0.95; 40]);
    values.extend(vec![0.5; 260]);
    write_fii(&s, &values);

    let args = ComputeCaiArgs {
        inputs: vec![s.clone()],
        manifest: None,
        config: None,
        cai_weights: None,
        out: out.clone(),
    };
    run_compute_cai_command(&args).expect("first");
    let a = fs::read(out.join("commitment_asymmetry.tsv")).expect("a");
    let b = fs::read(out.join("commitment_asymmetry.json")).expect("b");
    run_compute_cai_command(&args).expect("second");
    let a2 = fs::read(out.join("commitment_asymmetry.tsv")).expect("a2");
    let b2 = fs::read(out.join("commitment_asymmetry.json")).expect("b2");
    assert_eq!(a, a2);
    assert_eq!(b, b2);
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
        export_systems_model: None,
    })
    .expect("aggregate");
    assert!(!out.join("commitment_asymmetry.tsv").exists());
    assert!(!out.join("commitment_asymmetry.json").exists());
}

fn write_fii(dir: &std::path::Path, values: &[f64]) {
    fs::create_dir_all(dir).expect("mkdir");
    let mut body = String::from(
        "cell_id\tmitochondrial_stress_adaptation_score\ttranslation_commitment_score\tsplice_irreversibility_index\tfunctional_irreversibility_index\tfii_regime\tfii_low_confidence\n",
    );
    for (i, v) in values.iter().enumerate() {
        let regime = if *v < 0.33 {
            "Adaptive"
        } else if *v < 0.66 {
            "Transition"
        } else {
            "Resistant"
        };
        body.push_str(&format!(
            "c{idx}\t0.0\t0.0\t0.0\t{fii:.6}\t{regime}\tfalse\n",
            idx = i,
            fii = v
        ));
    }
    fs::write(dir.join("functional_irreversibility_index.tsv"), body).expect("write fii");
    fs::write(
        dir.join("summary.json"),
        r#"{"functional_irreversibility":{"distribution":{"mean":0.0,"median":0.0,"p25":0.0,"p75":0.0},"regime_fractions":{"adaptive":0.0,"transition":1.0,"resistant":0.0},"low_confidence_fraction":0.0}}"#,
    )
    .expect("write summary");
}

fn parse_tsv_rows(tsv: &str) -> std::collections::BTreeMap<String, (f64, usize)> {
    let mut out = std::collections::BTreeMap::new();
    for (i, line) in tsv.lines().enumerate() {
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let parts = line.split('\t').collect::<Vec<_>>();
        if parts.len() < 3 {
            continue;
        }
        let label = parts[0].to_string();
        let order_rank = parts[1].parse::<usize>().unwrap_or(0);
        let cai = parts[2].parse::<f64>().unwrap_or(0.0);
        out.insert(label, (cai, order_rank));
    }
    out
}

fn linspace(n: usize, start: f64, end: f64) -> Vec<f64> {
    if n <= 1 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (end - start) * (i as f64 / (n as f64 - 1.0)))
        .collect()
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
