use std::fs;

use kira_organelle::cli::ComputeCocsArgs;
use kira_organelle::run_compute_cocs_command;
use tempfile::tempdir;

#[test]
fn uncoupled_vs_coupled_transitions() {
    let dir = tempdir().expect("tempdir");
    let in_unc = dir.path().join("uncoupled.json");
    let in_cou = dir.path().join("coupled.json");
    let out_unc = dir.path().join("out_unc");
    let out_cou = dir.path().join("out_cou");

    fs::write(&in_unc, uncoupled_json()).expect("unc json");
    fs::write(&in_cou, coupled_json()).expect("cou json");

    run_compute_cocs_command(&ComputeCocsArgs {
        input: in_unc,
        out: out_unc.clone(),
    })
    .expect("unc run");
    run_compute_cocs_command(&ComputeCocsArgs {
        input: in_cou,
        out: out_cou.clone(),
    })
    .expect("cou run");

    let u = last_cocs_global(
        &fs::read_to_string(out_unc.join("cross_organelle_coupling.tsv")).expect("u"),
    );
    let c = last_cocs_global(
        &fs::read_to_string(out_cou.join("cross_organelle_coupling.tsv")).expect("c"),
    );
    assert!(c > u);
}

#[test]
fn monotonic_more_coupled_than_noisy() {
    let dir = tempdir().expect("tempdir");
    let in_m = dir.path().join("mono.json");
    let in_n = dir.path().join("noisy.json");
    let out_m = dir.path().join("out_m");
    let out_n = dir.path().join("out_n");
    fs::write(&in_m, monotonic_json()).expect("m");
    fs::write(&in_n, noisy_json()).expect("n");

    run_compute_cocs_command(&ComputeCocsArgs {
        input: in_m,
        out: out_m.clone(),
    })
    .expect("m run");
    run_compute_cocs_command(&ComputeCocsArgs {
        input: in_n,
        out: out_n.clone(),
    })
    .expect("n run");

    let m = last_cocs_global(
        &fs::read_to_string(out_m.join("cross_organelle_coupling.tsv")).expect("m tsv"),
    );
    let n = last_cocs_global(
        &fs::read_to_string(out_n.join("cross_organelle_coupling.tsv")).expect("n tsv"),
    );
    assert!(m >= n);
}

#[test]
fn missing_organelle_is_robust() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("missing.json");
    let out = dir.path().join("out");
    fs::write(&input, missing_json()).expect("write");
    run_compute_cocs_command(&ComputeCocsArgs {
        input,
        out: out.clone(),
    })
    .expect("run");
    let tsv = fs::read_to_string(out.join("cross_organelle_coupling.tsv")).expect("tsv");
    assert!(tsv.contains("COCS_global"));
    let v = last_cocs_global(&tsv);
    assert!((0.0..=1.0).contains(&v));
}

#[test]
fn deterministic_across_runs() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("in.json");
    let out = dir.path().join("out");
    fs::write(&input, coupled_json()).expect("write");
    let args = ComputeCocsArgs {
        input: input.clone(),
        out: out.clone(),
    };
    run_compute_cocs_command(&args).expect("first");
    let a = fs::read(out.join("cross_organelle_coupling.tsv")).expect("a");
    let b = fs::read(out.join("cross_organelle_coupling.json")).expect("b");
    run_compute_cocs_command(&args).expect("second");
    let a2 = fs::read(out.join("cross_organelle_coupling.tsv")).expect("a2");
    let b2 = fs::read(out.join("cross_organelle_coupling.json")).expect("b2");
    assert_eq!(a, a2);
    assert_eq!(b, b2);
}

fn last_cocs_global(tsv: &str) -> f64 {
    let mut last = 0.0;
    for (i, line) in tsv.lines().enumerate() {
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let p = line.split('\t').collect::<Vec<_>>();
        if p.len() >= 3 {
            last = p[2].parse::<f64>().unwrap_or(0.0);
        }
    }
    last
}

fn coupled_json() -> &'static str {
    r#"{"samples":[
      {"label":"s0","rank":0,"FII_nucleus":0.10,"FII_splice":0.10,"FII_proteostasis":0.10,"FII_mito":0.10,"FII_tme":0.10},
      {"label":"s1","rank":1,"FII_nucleus":0.20,"FII_splice":0.20,"FII_proteostasis":0.20,"FII_mito":0.20,"FII_tme":0.20},
      {"label":"s2","rank":2,"FII_nucleus":0.35,"FII_splice":0.35,"FII_proteostasis":0.35,"FII_mito":0.35,"FII_tme":0.35},
      {"label":"s3","rank":3,"FII_nucleus":0.55,"FII_splice":0.55,"FII_proteostasis":0.55,"FII_mito":0.55,"FII_tme":0.55}
    ]}"#
}

fn uncoupled_json() -> &'static str {
    r#"{"samples":[
      {"label":"s0","rank":0,"FII_nucleus":0.10,"FII_splice":0.10,"FII_proteostasis":0.10,"FII_mito":0.10,"FII_tme":0.10},
      {"label":"s1","rank":1,"FII_nucleus":0.20,"FII_splice":0.05,"FII_proteostasis":0.13,"FII_mito":0.11,"FII_tme":0.07},
      {"label":"s2","rank":2,"FII_nucleus":0.25,"FII_splice":0.20,"FII_proteostasis":0.11,"FII_mito":0.18,"FII_tme":0.09},
      {"label":"s3","rank":3,"FII_nucleus":0.30,"FII_splice":0.10,"FII_proteostasis":0.20,"FII_mito":0.15,"FII_tme":0.14}
    ]}"#
}

fn monotonic_json() -> &'static str {
    r#"{"samples":[
      {"label":"s0","rank":0,"FII_nucleus":0.10,"FII_splice":0.10,"FII_proteostasis":0.10,"FII_mito":0.10,"FII_tme":0.10},
      {"label":"s1","rank":1,"FII_nucleus":0.14,"FII_splice":0.13,"FII_proteostasis":0.15,"FII_mito":0.145,"FII_tme":0.135},
      {"label":"s2","rank":2,"FII_nucleus":0.23,"FII_splice":0.21,"FII_proteostasis":0.25,"FII_mito":0.235,"FII_tme":0.225},
      {"label":"s3","rank":3,"FII_nucleus":0.26,"FII_splice":0.24,"FII_proteostasis":0.29,"FII_mito":0.27,"FII_tme":0.255},
      {"label":"s4","rank":4,"FII_nucleus":0.38,"FII_splice":0.35,"FII_proteostasis":0.42,"FII_mito":0.39,"FII_tme":0.37}
    ]}"#
}

fn noisy_json() -> &'static str {
    r#"{"samples":[
      {"label":"s0","rank":0,"FII_nucleus":0.10,"FII_splice":0.10,"FII_proteostasis":0.20,"FII_mito":0.20,"FII_tme":0.20},
      {"label":"s1","rank":1,"FII_nucleus":0.25,"FII_splice":0.13,"FII_proteostasis":0.12,"FII_mito":0.31,"FII_tme":0.16},
      {"label":"s2","rank":2,"FII_nucleus":0.13,"FII_splice":0.33,"FII_proteostasis":0.26,"FII_mito":0.26,"FII_tme":0.23},
      {"label":"s3","rank":3,"FII_nucleus":0.35,"FII_splice":0.16,"FII_proteostasis":0.28,"FII_mito":0.10,"FII_tme":0.41},
      {"label":"s4","rank":4,"FII_nucleus":0.17,"FII_splice":0.21,"FII_proteostasis":0.07,"FII_mito":0.34,"FII_tme":0.32},
      {"label":"s5","rank":5,"FII_nucleus":0.36,"FII_splice":0.10,"FII_proteostasis":0.16,"FII_mito":0.30,"FII_tme":0.19}
    ]}"#
}

fn missing_json() -> &'static str {
    r#"{"samples":[
      {"label":"s0","rank":0,"FII_nucleus":0.10,"FII_splice":0.10,"FII_proteostasis":0.10,"FII_mito":0.10},
      {"label":"s1","rank":1,"FII_nucleus":0.22,"FII_splice":0.19,"FII_proteostasis":0.20,"FII_tme":0.18},
      {"label":"s2","rank":2,"FII_nucleus":0.33,"FII_splice":0.25,"FII_proteostasis":0.26,"FII_mito":0.24,"FII_tme":0.23}
    ]}"#
}
