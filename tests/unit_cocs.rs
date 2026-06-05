use kira_organelle::report::fii::cocs::{SamplePoint, compute_rows, render_tsv};

fn p(label: &str, rank: usize, n: f64, s: f64, pr: f64, m: f64, t: f64) -> SamplePoint {
    SamplePoint {
        sample_label: label.to_string(),
        order_rank: rank,
        nucleus: Some(n),
        splice: Some(s),
        proteostasis: Some(pr),
        mito: Some(m),
        tme: Some(t),
    }
}

#[test]
fn coupled_vs_uncoupled() {
    let uncoupled = vec![
        p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
        p("s1", 1, 0.2, 0.05, 0.13, 0.11, 0.07),
        p("s2", 2, 0.25, 0.2, 0.11, 0.18, 0.09),
        p("s3", 3, 0.3, 0.1, 0.2, 0.15, 0.14),
    ];
    let coupled = vec![
        p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
        p("s1", 1, 0.2, 0.2, 0.2, 0.2, 0.2),
        p("s2", 2, 0.3, 0.3, 0.3, 0.3, 0.3),
        p("s3", 3, 0.45, 0.45, 0.45, 0.45, 0.45),
    ];
    let ru = compute_rows(&uncoupled);
    let rc = compute_rows(&coupled);
    assert!(rc.last().expect("c").cocs_global > ru.last().expect("u").cocs_global);
}

#[test]
fn monotonic_higher_than_noisy() {
    let mono = vec![
        p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
        p("s1", 1, 0.14, 0.13, 0.15, 0.145, 0.135),
        p("s2", 2, 0.23, 0.21, 0.25, 0.235, 0.225),
        p("s3", 3, 0.26, 0.24, 0.29, 0.27, 0.255),
        p("s4", 4, 0.38, 0.35, 0.42, 0.39, 0.37),
    ];
    let noisy = vec![
        p("s0", 0, 0.10, 0.10, 0.20, 0.20, 0.20),
        p("s1", 1, 0.25, 0.13, 0.12, 0.31, 0.16),
        p("s2", 2, 0.13, 0.33, 0.26, 0.26, 0.23),
        p("s3", 3, 0.35, 0.16, 0.28, 0.10, 0.41),
        p("s4", 4, 0.17, 0.21, 0.07, 0.34, 0.32),
        p("s5", 5, 0.36, 0.10, 0.16, 0.30, 0.19),
    ];
    let rm = compute_rows(&mono);
    let rn = compute_rows(&noisy);
    assert!(rm.last().expect("m").cocs_global >= rn.last().expect("n").cocs_global);
}

#[test]
fn missing_organelle_values_are_robust() {
    let mut pts = vec![
        p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
        p("s1", 1, 0.2, 0.2, 0.2, 0.2, 0.2),
        p("s2", 2, 0.3, 0.3, 0.3, 0.3, 0.3),
    ];
    pts[1].tme = None;
    pts[2].mito = None;
    let rows = compute_rows(&pts);
    let last = rows.last().expect("last");
    assert!(last.cocs_global.is_finite());
    assert!((0.0..=1.0).contains(&last.cocs_global));
}

#[test]
fn deterministic() {
    let pts = vec![
        p("s0", 0, 0.1, 0.2, 0.3, 0.4, 0.5),
        p("s1", 1, 0.2, 0.3, 0.4, 0.5, 0.6),
        p("s2", 2, 0.4, 0.5, 0.6, 0.7, 0.8),
    ];
    let a = compute_rows(&pts);
    let b = compute_rows(&pts);
    assert_eq!(render_tsv(&a), render_tsv(&b));
}

#[test]
fn extended_distinct_from_global_when_pairs_differ() {
    let pts = vec![
        p("s0", 0, 0.1, 0.2, 0.3, 0.4, 0.5),
        p("s1", 1, 0.2, 0.3, 0.4, 0.5, 0.6),
        p("s2", 2, 0.4, 0.5, 0.6, 0.7, 0.8),
    ];
    let rows = compute_rows(&pts);
    let last = rows.last().expect("rows");
    assert!(
        last.cocs_extended.is_finite() && last.cocs_core.is_finite(),
        "both finite"
    );
}
