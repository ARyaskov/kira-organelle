#[derive(Debug, Clone, Copy)]
pub struct BinSpec {
    pub bins: usize,
    pub min: f64,
    pub max: f64,
}

#[derive(Debug, Clone)]
pub struct Histogram1D {
    pub spec: BinSpec,
    pub counts: Vec<u64>,
}

impl Histogram1D {
    pub fn new(spec: BinSpec) -> Self {
        Self {
            spec,
            counts: vec![0; spec.bins],
        }
    }

    pub fn add(&mut self, value: f64) {
        let idx = index_1d(self.spec, value);
        self.counts[idx] += 1;
    }
}

#[derive(Debug, Clone)]
pub struct Heatmap2D {
    pub x: BinSpec,
    pub y: BinSpec,
    pub counts: Vec<u64>,
    pub sum_x: Vec<f64>,
    pub sum_y: Vec<f64>,
}

impl Heatmap2D {
    pub fn new(x: BinSpec, y: BinSpec) -> Self {
        let len = x.bins * y.bins;
        Self {
            x,
            y,
            counts: vec![0; len],
            sum_x: vec![0.0; len],
            sum_y: vec![0.0; len],
        }
    }

    pub fn add(&mut self, x: f64, y: f64) {
        let ix = index_1d(self.x, x);
        let iy = index_1d(self.y, y);
        let idx = iy * self.x.bins + ix;
        self.counts[idx] += 1;
        self.sum_x[idx] += x;
        self.sum_y[idx] += y;
    }
}

pub fn clamp01(v: f64) -> f64 {
    if !v.is_finite() {
        0.0
    } else {
        v.clamp(0.0, 1.0)
    }
}

fn index_1d(spec: BinSpec, value: f64) -> usize {
    let v = if value < spec.min {
        spec.min
    } else if value > spec.max {
        spec.max
    } else {
        value
    };
    let ratio = (v - spec.min) / (spec.max - spec.min);
    let mut idx = (ratio * spec.bins as f64).floor() as usize;
    if idx >= spec.bins {
        idx = spec.bins - 1;
    }
    idx
}

#[cfg(test)]
mod tests {
    use super::{BinSpec, Heatmap2D, Histogram1D, clamp01};

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
}
