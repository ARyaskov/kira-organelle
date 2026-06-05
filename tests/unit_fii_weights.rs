use kira_organelle::fii::{FiiRegime, FiiWeights, classify_regime};

#[test]
fn parse_weights() {
    let w = FiiWeights::parse("mito:0.34,translation:0.33,splice:0.33").expect("parse");
    assert!((w.mitochondrial - 0.34).abs() < 1e-9);
    assert!((w.translation - 0.33).abs() < 1e-9);
    assert!((w.splice - 0.33).abs() < 1e-9);
}

#[test]
fn reject_bad_weight_sum() {
    let err = FiiWeights::parse("mito:0.5,translation:0.5,splice:0.5").expect_err("err");
    assert!(err.contains("sum"));
}

#[test]
fn regime_boundaries() {
    assert_eq!(classify_regime(0.0), FiiRegime::Adaptive);
    assert_eq!(classify_regime(0.329999), FiiRegime::Adaptive);
    assert_eq!(classify_regime(0.33), FiiRegime::Transition);
    assert_eq!(classify_regime(0.659999), FiiRegime::Transition);
    assert_eq!(classify_regime(0.66), FiiRegime::Resistant);
    assert_eq!(classify_regime(1.0), FiiRegime::Resistant);
}
