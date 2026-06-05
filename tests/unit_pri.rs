use kira_organelle::report::fii::pri::{PriConfig, PriRow, normalize_weights};

#[test]
fn synthetic_component_combinations() {
    let cfg = PriConfig::default();
    let a = PriRow {
        sample_label: "a".to_string(),
        order_rank: 0,
        pri: 0.0,
        nuclear_plasticity_component: 1.0,
        splicing_integrity_component: 1.0,
        translational_selectivity_component: 1.0,
        low_confidence: false,
    };
    let b = PriRow {
        sample_label: "b".to_string(),
        order_rank: 1,
        pri: 0.0,
        nuclear_plasticity_component: 0.1,
        splicing_integrity_component: 0.2,
        translational_selectivity_component: 0.1,
        low_confidence: false,
    };
    let pa = (cfg.w_nuclear * a.nuclear_plasticity_component
        + cfg.w_splice * a.splicing_integrity_component
        + cfg.w_translation * a.translational_selectivity_component)
        .clamp(0.0, 1.0);
    let pb = (cfg.w_nuclear * b.nuclear_plasticity_component
        + cfg.w_splice * b.splicing_integrity_component
        + cfg.w_translation * b.translational_selectivity_component)
        .clamp(0.0, 1.0);
    assert!(pa > pb);
}

#[test]
fn weight_normalization_sums_to_one() {
    let mut cfg = PriConfig {
        w_nuclear: 4.0,
        w_splice: 3.0,
        w_translation: 3.0,
    };
    normalize_weights(&mut cfg).expect("normalize");
    assert!((cfg.w_nuclear + cfg.w_splice + cfg.w_translation - 1.0).abs() < 1e-12);
}
