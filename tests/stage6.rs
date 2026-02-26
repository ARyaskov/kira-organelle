use std::fs;

use kira_organelle::error::{ErrorKind, classify_issue_code, stable_error_code};
use kira_organelle::util::tsv::TsvReader;
use tempfile::tempdir;

#[test]
fn outputs_unchanged_after_optimization() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("metrics.tsv");
    fs::write(&path, "barcode\ta\tb\nBC1\t1.0\t2\nBC2\t\t4\nBC3\t7\t\n").expect("write");

    let mut optimized_rows = Vec::<Vec<String>>::new();
    let mut reader = TsvReader::open(&path).expect("open");
    let mut fields = Vec::new();
    while reader.read_record(&mut fields).expect("read") {
        let mut row = Vec::new();
        for i in 0..fields.len() {
            row.push(reader.field(&fields, i).unwrap_or_default().to_string());
        }
        optimized_rows.push(row);
    }

    let legacy_rows = fs::read_to_string(&path)
        .expect("read string")
        .lines()
        .map(|line| line.split('\t').map(ToOwned::to_owned).collect::<Vec<_>>())
        .collect::<Vec<_>>();

    let optimized_bytes = serde_json::to_vec_pretty(&optimized_rows).expect("json optimized");
    let legacy_bytes = serde_json::to_vec_pretty(&legacy_rows).expect("json legacy");

    assert_eq!(optimized_bytes, legacy_bytes);
}

#[test]
fn error_kind_mapping() {
    assert_eq!(
        classify_issue_code("PRIMARY_METRICS_READ_ERROR"),
        ErrorKind::Io
    );
    assert_eq!(classify_issue_code("SUMMARY_PARSE_ERROR"), ErrorKind::Parse);
    assert_eq!(
        classify_issue_code("MISSING_SUMMARY"),
        ErrorKind::ContractViolation
    );
    assert_eq!(
        classify_issue_code("TOOL_EXECUTION_FAILED"),
        ErrorKind::ToolFailure
    );

    assert_eq!(stable_error_code(ErrorKind::Io), "E_IO");
    assert_eq!(stable_error_code(ErrorKind::Parse), "E_PARSE");
    assert_eq!(
        stable_error_code(ErrorKind::ContractViolation),
        "E_CONTRACT"
    );
    assert_eq!(stable_error_code(ErrorKind::ToolFailure), "E_TOOL");
    assert_eq!(stable_error_code(ErrorKind::Internal), "E_INTERNAL");
}
