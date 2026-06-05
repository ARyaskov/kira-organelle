use kira_organelle::report::fii::dci::{DciConfig, Sample, compute_row};

fn sample(
    nuc: Option<f64>,
    spl: Option<f64>,
    pr: Option<f64>,
    mt: Option<f64>,
    tme: Option<f64>,
) -> Sample {
    Sample {
        sample_label: "x".to_string(),
        order_rank: 0,
        nucleus: nuc,
        splice: spl,
        proteostasis: pr,
        mito: mt,
        tme,
    }
}

#[test]
fn fully_concordant_is_high() {
    let cfg = DciConfig::default();
    let s = sample(Some(0.9), Some(0.91), Some(0.95), Some(0.88), Some(0.9));
    let row = compute_row(&s, &cfg);
    assert!((row.dci - 1.0).abs() < 1e-12);
}

#[test]
fn fully_conflicting_is_low() {
    let cfg = DciConfig::default();
    let s = sample(Some(0.1), Some(0.9), Some(0.1), Some(0.9), Some(0.1));
    let row = compute_row(&s, &cfg);
    assert!(row.dci < 0.5);
}

#[test]
fn mixed_partial_is_mid() {
    let cfg = DciConfig::default();
    let s = sample(Some(0.2), Some(0.5), Some(0.8), Some(0.5), Some(0.2));
    let row = compute_row(&s, &cfg);
    assert!(row.dci > 0.2 && row.dci < 0.9);
}

#[test]
fn missing_is_robust() {
    let cfg = DciConfig::default();
    let s = sample(Some(0.2), None, Some(0.8), None, Some(0.2));
    let row = compute_row(&s, &cfg);
    assert!((0.0..=1.0).contains(&row.dci));
}

#[test]
fn deterministic() {
    let cfg = DciConfig::default();
    let s = sample(Some(0.2), Some(0.5), Some(0.8), Some(0.4), Some(0.9));
    let a = compute_row(&s, &cfg);
    let b = compute_row(&s, &cfg);
    assert!((a.dci - b.dci).abs() < 1e-12);
    assert!((a.dci_core - b.dci_core).abs() < 1e-12);
}
