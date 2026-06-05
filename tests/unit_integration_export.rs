use kira_organelle::integration_export::{order_experiments, parse_manifest_delimited};
use tempfile::tempdir;

#[test]
fn manifest_tsv_order_is_applied() {
    let raw = "path\torder\ttimepoint\nB\t1\t10\nA\t0\t5\n";
    let parsed = parse_manifest_delimited(raw, '\t').expect("manifest parse");
    let a = parsed.get("A").expect("A");
    let b = parsed.get("B").expect("B");
    assert_eq!(a.order_rank, 0);
    assert_eq!(a.timepoint, Some(5));
    assert_eq!(b.order_rank, 1);
    assert_eq!(b.timepoint, Some(10));
}

#[test]
fn order_falls_back_to_input_for_missing_manifest_entries() {
    let dir = tempdir().expect("tempdir");
    let manifest = dir.path().join("manifest.tsv");
    std::fs::write(&manifest, "label\torder\nB\t0\n").expect("manifest write");

    let exps = vec![
        ("A".to_string(), dir.path().join("out-a")),
        ("B".to_string(), dir.path().join("out-b")),
    ];
    let ordered = order_experiments(&exps, Some(&manifest)).expect("order");
    assert_eq!(ordered[0].name, "B");
    assert_eq!(ordered[0].timepoint, Some(0));
    assert_eq!(ordered[1].name, "A");
    assert_eq!(ordered[1].timepoint, None);
}
