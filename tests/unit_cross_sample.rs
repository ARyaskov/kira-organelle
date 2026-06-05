use kira_organelle::systems::cross_sample::{basin_overlap_for_test, js_divergence_for_test};

#[test]
fn js_zero_for_identical_distribution() {
    let p = vec![0.2, 0.3, 0.5];
    let q = vec![0.2, 0.3, 0.5];
    assert_eq!(js_divergence_for_test(&p, &q), 0.0);
}

#[test]
fn basin_overlap_is_one_if_both_empty() {
    assert_eq!(basin_overlap_for_test(), 1.0);
}
