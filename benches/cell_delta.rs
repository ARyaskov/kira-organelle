use criterion::{Criterion, criterion_group, criterion_main};
use kira_organelle::util::select::{median_in_place, percentile_nearest_rank_in_place};

fn bench_cell_delta(c: &mut Criterion) {
    let baseline = (0..100_000u32)
        .map(|i| (i as f64 * 0.001).sin())
        .collect::<Vec<_>>();

    c.bench_function("cell_delta_quantiles_100k", |b| {
        b.iter(|| {
            let mut values = baseline.clone();
            let median = median_in_place(&mut values).expect("median");
            for v in &mut values {
                *v = v.abs();
            }
            let p90 = percentile_nearest_rank_in_place(&mut values, 0.90).expect("p90");
            (median, p90)
        })
    });
}

criterion_group!(benches, bench_cell_delta);
criterion_main!(benches);
