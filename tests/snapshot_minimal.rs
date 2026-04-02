use std::fs;

use kira_organelle::{AggregateOptions, run_aggregate};
use tempfile::tempdir;

#[test]
fn snapshot_minimal_is_byte_identical() {
    let repo = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let input = repo.join("tests/fixtures/snapshot_input");
    let expected_metrics = fs::read(repo.join("tests/fixtures/expected_metrics.tsv"))
        .expect("read expected_metrics.tsv");
    let expected_summary = fs::read(repo.join("tests/fixtures/expected_summary.json"))
        .expect("read expected_summary.json");

    let first = tempdir().expect("tempdir first");
    let second = tempdir().expect("tempdir second");

    run_aggregate(&AggregateOptions {
        input: input.clone(),
        input_b: None,
        out: Some(first.path().to_path_buf()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
        export_systems_model: None,
    })
    .expect("aggregate first");
    run_aggregate(&AggregateOptions {
        input,
        input_b: None,
        out: Some(second.path().to_path_buf()),
        strict: false,
        json: true,
        validate_only: false,
        fii_weights: None,
        export_systems_model: None,
    })
    .expect("aggregate second");

    let first_metrics =
        fs::read(first.path().join("integration/metrics.tsv")).expect("first metrics");
    let first_summary =
        fs::read(first.path().join("integration/summary.json")).expect("first summary");
    let second_metrics =
        fs::read(second.path().join("integration/metrics.tsv")).expect("second metrics");
    let second_summary =
        fs::read(second.path().join("integration/summary.json")).expect("second summary");

    assert_eq!(
        first_metrics, second_metrics,
        "metrics.tsv should be stable"
    );
    assert_eq!(
        first_summary, second_summary,
        "summary.json should be stable"
    );
    assert_eq!(
        first_metrics, expected_metrics,
        "metrics.tsv should match golden fixture"
    );
    assert_eq!(
        first_summary, expected_summary,
        "summary.json should match golden fixture"
    );
}
