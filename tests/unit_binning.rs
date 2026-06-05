use kira_organelle::report::fii::binning::{BinSpec, Heatmap2D, Histogram1D, clamp01};

#[test]
fn histogram_counts_are_stable() {
    let mut h = Histogram1D::new(BinSpec {
        bins: 4,
        min: 0.0,
        max: 1.0,
    });
    h.add(0.0);
    h.add(0.24);
    h.add(0.25);
    h.add(0.9);
    h.add(1.0);
    assert_eq!(h.counts, vec![2, 1, 0, 2]);
}

#[test]
fn heatmap_indices_are_stable() {
    let mut hm = Heatmap2D::new(
        BinSpec {
            bins: 2,
            min: 0.0,
            max: 1.0,
        },
        BinSpec {
            bins: 2,
            min: 0.0,
            max: 1.0,
        },
    );
    hm.add(0.1, 0.1);
    hm.add(0.9, 0.1);
    hm.add(0.9, 0.9);
    assert_eq!(hm.counts, vec![1, 1, 0, 1]);
}

#[test]
fn clamp01_bounds_values() {
    assert_eq!(clamp01(-2.0), 0.0);
    assert_eq!(clamp01(2.0), 1.0);
    assert_eq!(clamp01(f64::NAN), 0.0);
}
