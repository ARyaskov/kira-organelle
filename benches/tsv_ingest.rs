use std::fs;

use criterion::{Criterion, criterion_group, criterion_main};
use kira_organelle::util::tsv::TsvReader;
use tempfile::tempdir;

fn bench_tsv_ingest(c: &mut Criterion) {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("synthetic.tsv");

    let mut buf = String::new();
    buf.push_str("barcode\ta\tb\tc\n");
    for i in 0..50_000u32 {
        buf.push_str(&format!("BC{:05}\t{}\t{}\t{}\n", i, i % 7, i % 13, i % 17));
    }
    fs::write(&path, buf).expect("write");

    c.bench_function("tsv_ingest_50k", |b| {
        b.iter(|| {
            let mut reader = TsvReader::open(&path).expect("open");
            let mut fields = Vec::new();
            let mut rows = 0usize;
            while reader.read_record(&mut fields).expect("read") {
                rows += 1;
            }
            rows
        })
    });
}

criterion_group!(benches, bench_tsv_ingest);
criterion_main!(benches);
