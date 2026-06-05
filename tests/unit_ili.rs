use kira_organelle::report::fii::ili::{IliConfig, SampleRow, compute_transition};

fn sample(
    label: &str,
    rank: usize,
    nucleus: Option<f64>,
    splice: Option<f64>,
    proteostasis: Option<f64>,
    mito: Option<f64>,
    tme: Option<f64>,
) -> SampleRow {
    SampleRow {
        sample_label: label.to_string(),
        order_rank: rank,
        nucleus,
        splice,
        proteostasis,
        mito,
        tme,
    }
}

#[test]
fn ili_picks_known_leading_organelle() {
    let cfg = IliConfig::default();
    let prev = sample(
        "baseline",
        0,
        Some(0.20),
        Some(0.20),
        Some(0.20),
        Some(0.20),
        Some(0.20),
    );
    let curr = sample(
        "step1",
        1,
        Some(0.21),
        Some(0.48),
        Some(0.23),
        Some(0.22),
        Some(0.20),
    );
    let row = compute_transition(&prev, &curr, &cfg);
    assert_eq!(row.ili_category, "SPLICE_DRIVEN");
    assert!(row.ili_confidence > 0.5);
}

#[test]
fn ili_tie_and_low_signal_resolve_to_unresolved() {
    let cfg = IliConfig {
        tie_epsilon: 1e-5,
        min_signal_threshold: 0.02,
    };
    let prev = sample(
        "baseline",
        0,
        Some(0.10),
        Some(0.10),
        Some(0.10),
        Some(0.10),
        Some(0.10),
    );
    let tie_curr = sample(
        "tie",
        1,
        Some(0.30),
        Some(0.30),
        Some(0.10),
        Some(0.10),
        Some(0.10),
    );
    assert_eq!(
        compute_transition(&prev, &tie_curr, &cfg).ili_category,
        "UNRESOLVED"
    );

    let low_curr = sample(
        "low",
        1,
        Some(0.105),
        Some(0.103),
        Some(0.102),
        Some(0.101),
        Some(0.100),
    );
    assert_eq!(
        compute_transition(&prev, &low_curr, &cfg).ili_category,
        "UNRESOLVED"
    );
}

#[test]
fn ili_skips_missing_components_and_still_computes() {
    let cfg = IliConfig::default();
    let prev = sample("baseline", 0, None, Some(0.15), None, Some(0.10), None);
    let curr = sample("step1", 1, None, Some(0.22), None, Some(0.12), None);
    let row = compute_transition(&prev, &curr, &cfg);
    assert_eq!(row.ili_category, "SPLICE_DRIVEN");
    assert!(row.leading_deltafii_value > 0.0);
}

#[test]
fn ili_transition_deterministic() {
    let cfg = IliConfig::default();
    let prev = sample(
        "baseline",
        0,
        Some(0.20),
        Some(0.20),
        Some(0.20),
        Some(0.20),
        Some(0.20),
    );
    let curr = sample(
        "step1",
        2,
        Some(0.24),
        Some(0.42),
        Some(0.23),
        Some(0.28),
        Some(0.21),
    );
    let a = compute_transition(&prev, &curr, &cfg);
    let b = compute_transition(&prev, &curr, &cfg);
    assert_eq!(a.ili_category, b.ili_category);
    assert_eq!(a.order_rank, b.order_rank);
    assert!((a.ili_confidence - b.ili_confidence).abs() < 1e-12);
    assert!((a.leading_deltafii_value - b.leading_deltafii_value).abs() < 1e-12);
}
