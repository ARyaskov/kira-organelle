use std::fs;

use kira_organelle::run::experiments::{ExperimentSource, discover_experiments};
use tempfile::tempdir;

#[test]
fn discover_multi_raw_counts() {
    let dir = tempdir().expect("tempdir");
    fs::write(dir.path().join("GSM1_raw_counts.tsv.gz"), b"").expect("write");
    fs::write(dir.path().join("GSM2_raw_counts.tsv.gz"), b"").expect("write");

    let specs = discover_experiments(dir.path()).expect("discover");
    assert_eq!(specs.len(), 2);
    assert_eq!(specs[0].name, "GSM1_raw_counts");
    assert_eq!(specs[1].name, "GSM2_raw_counts");
}

#[test]
fn discover_multi_tenx_prefixed_dot() {
    let dir = tempdir().expect("tempdir");
    fs::write(dir.path().join("A.matrix.mtx.gz"), b"").expect("write");
    fs::write(dir.path().join("A.features.tsv.gz"), b"").expect("write");
    fs::write(dir.path().join("A.barcodes.tsv.gz"), b"").expect("write");
    fs::write(dir.path().join("B.matrix.mtx.gz"), b"").expect("write");
    fs::write(dir.path().join("B.features.tsv.gz"), b"").expect("write");
    fs::write(dir.path().join("B.barcodes.tsv.gz"), b"").expect("write");

    let specs = discover_experiments(dir.path()).expect("discover");
    assert_eq!(specs.len(), 2);
    assert_eq!(specs[0].name, "A");
    assert_eq!(specs[1].name, "B");
}

#[test]
fn fall_back_to_single_directory_mode() {
    let dir = tempdir().expect("tempdir");
    fs::write(dir.path().join("matrix.mtx.gz"), b"").expect("write");
    fs::write(dir.path().join("features.tsv.gz"), b"").expect("write");
    fs::write(dir.path().join("barcodes.tsv.gz"), b"").expect("write");

    let specs = discover_experiments(dir.path()).expect("discover");
    assert_eq!(specs.len(), 1);
    assert!(matches!(specs[0].source, ExperimentSource::Directory(_)));
}
