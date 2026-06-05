use kira_organelle::registry::metrics::{MetricId, ToolId, find_metric_spec, metric_spec_by_id};

#[test]
fn canonical_order_and_count_are_stable() {
    assert_eq!(MetricId::COUNT, 42);
    assert_eq!(metric_spec_by_id(MetricId::Rss).canonical_name, "RSS");
    assert_eq!(metric_spec_by_id(MetricId::Msm).canonical_name, "MSM");
}

#[test]
fn alias_lookup_is_supported() {
    let spec = find_metric_spec(ToolId::RiboQc, "translation_stress_index").expect("alias");
    assert_eq!(spec.canonical_name, "TSM");
    let exact = find_metric_spec(ToolId::RiboQc, "mTOR_P").expect("exact");
    assert_eq!(exact.canonical_name, "mTOR_P");
}
