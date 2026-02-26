use std::fs;
use std::io::Write;

use flate2::Compression;
use flate2::write::GzEncoder;
use kira_organelle::run::exec::resolve_barcodes_path;
use tempfile::tempdir;

#[test]
fn resolve_barcodes_accepts_prefixed_names() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path();
    let out = dir.path().join("out");
    fs::create_dir_all(&out).expect("mkdir out");
    let target = input.join("SOME_PREFIX_HERE_barcodes.tsv.gz");
    fs::write(&target, b"dummy").expect("write");

    let resolved = resolve_barcodes_path(input, &out).expect("resolve");
    assert_eq!(resolved, target);
}

#[test]
fn resolve_barcodes_falls_back_to_raw_counts_gz_header() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out = dir.path().join("out");
    fs::create_dir_all(&input).expect("mkdir input");
    fs::create_dir_all(&out).expect("mkdir out");

    let raw_path = input.join("GSM123_raw_counts.tsv.gz");
    let file = fs::File::create(&raw_path).expect("create raw.gz");
    let mut enc = GzEncoder::new(file, Compression::default());
    enc.write_all(b"gene\tcellA\tcellB\nGENE1\t1\t2\n")
        .expect("write gz");
    enc.finish().expect("finish gz");

    let resolved = resolve_barcodes_path(&input, &out).expect("resolve");
    let text = fs::read_to_string(&resolved).expect("read out");
    assert!(resolved.ends_with("barcodes.from_raw_counts.tsv"));
    assert_eq!(text, "cellA\ncellB\n");
}
