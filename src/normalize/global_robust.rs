use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::cells::types::{CellsState, MetricStore};
use crate::io::canonicalize_f32;
use crate::registry::metrics::{METRIC_SPECS, MetricId, metric_spec_by_id};
use crate::util::select::median_in_place;

pub const EPS: f32 = 1e-9;
pub const MIN_VALID_CELLS: usize = 20;
const GAUSSIAN_MAD_SCALE: f32 = 1.4826;

#[derive(Debug, Clone)]
pub struct NormStore {
    pub z: Vec<f32>,
    pub present: Vec<bool>,
    pub n_cells: usize,
    pub n_metrics: usize,
}

impl NormStore {
    pub fn new(n_cells: usize) -> Self {
        let n_metrics = MetricId::COUNT;
        Self {
            z: vec![f32::NAN; n_cells * n_metrics],
            present: vec![false; n_cells * n_metrics],
            n_cells,
            n_metrics,
        }
    }

    pub fn get(&self, metric_id: MetricId, cell_idx: usize) -> f32 {
        let idx = self.offset(metric_id, cell_idx);
        match idx {
            Some(i) if self.present[i] => self.z[i],
            _ => f32::NAN,
        }
    }

    pub fn set(&mut self, metric_id: MetricId, cell_idx: usize, value: f32) {
        let Some(idx) = self.offset(metric_id, cell_idx) else {
            return;
        };
        self.z[idx] = value;
        self.present[idx] = !value.is_nan();
    }

    pub fn metric_present(&self, metric_id: MetricId) -> bool {
        let metric_idx = metric_id.as_index();
        if metric_idx >= self.n_metrics {
            return false;
        }
        let start = metric_idx * self.n_cells;
        let end = start + self.n_cells;
        self.present[start..end].iter().any(|v| *v)
    }

    fn offset(&self, metric_id: MetricId, cell_idx: usize) -> Option<usize> {
        if cell_idx >= self.n_cells {
            return None;
        }
        let metric_idx = metric_id.as_index();
        if metric_idx >= self.n_metrics {
            return None;
        }
        Some(metric_idx * self.n_cells + cell_idx)
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct RobustStats {
    pub median: f32,
    pub mad: f32,
    pub scale: f32,
    pub n_valid: u32,
    pub reliable: bool,
    pub scale_degenerate: bool,
    pub missing_fraction: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlobalNormalizationMetricSummary {
    pub metric: String,
    pub median: f32,
    pub mad: f32,
    pub scale: f32,
    pub n_valid: u32,
    pub reliable: bool,
    pub scale_degenerate: bool,
    pub missing_fraction: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlobalNormalizationSummary {
    pub method: String,
    pub eps: f32,
    pub min_valid_cells: usize,
    pub per_metric: Vec<GlobalNormalizationMetricSummary>,
}

#[derive(Debug, Clone)]
pub struct GlobalNormalizationContext {
    pub norm_store: NormStore,
    pub stats: Vec<RobustStats>,
    pub summary: GlobalNormalizationSummary,
}

pub fn compute_global_robust_normalization(cells: &CellsState) -> GlobalNormalizationContext {
    let n_cells = cells.n_cells;
    let mut raw_store = MetricStore::new(n_cells);
    let mut norm_store = NormStore::new(n_cells);
    let mut stats = vec![
        RobustStats {
            median: 0.0,
            mad: 0.0,
            scale: 1.0,
            n_valid: 0,
            reliable: false,
            scale_degenerate: false,
            missing_fraction: 1.0,
        };
        MetricId::COUNT
    ];

    let canonical_to_metric = METRIC_SPECS
        .iter()
        .map(|spec| (spec.canonical_name, spec.id))
        .collect::<BTreeMap<_, _>>();

    for (cell_idx, cell) in cells.cells.iter().enumerate() {
        for org_state in cell.per_organelle.values() {
            for (axis, value) in &org_state.axes {
                if let Some(metric_id) = canonical_to_metric.get(axis.as_str()) {
                    raw_store.set(*metric_id, cell_idx, *value as f32);
                }
            }
        }
    }

    let mut buf = Vec::<f64>::with_capacity(n_cells);
    let mut dev = Vec::<f64>::with_capacity(n_cells);
    let mut per_metric_summary = Vec::<GlobalNormalizationMetricSummary>::new();

    for spec in METRIC_SPECS {
        let metric_id = spec.id;
        buf.clear();
        for cell_idx in 0..n_cells {
            let v = raw_store.get(metric_id, cell_idx);
            if v.is_finite() {
                buf.push(v as f64);
            }
        }

        let n_valid = buf.len();
        let missing_fraction = if n_cells == 0 {
            1.0
        } else {
            1.0 - (n_valid as f32 / n_cells as f32)
        };

        let idx = metric_id.as_index();
        if n_valid < MIN_VALID_CELLS {
            stats[idx] = RobustStats {
                median: 0.0,
                mad: 0.0,
                scale: 1.0,
                n_valid: n_valid as u32,
                reliable: false,
                scale_degenerate: false,
                missing_fraction: canonicalize_f32(missing_fraction),
            };
            for cell_idx in 0..n_cells {
                norm_store.set(metric_id, cell_idx, f32::NAN);
            }
            if n_valid > 0 {
                per_metric_summary.push(GlobalNormalizationMetricSummary {
                    metric: spec.canonical_name.to_string(),
                    median: 0.0,
                    mad: 0.0,
                    scale: 1.0,
                    n_valid: n_valid as u32,
                    reliable: false,
                    scale_degenerate: false,
                    missing_fraction: canonicalize_f32(missing_fraction),
                });
            }
            continue;
        }

        let median = median_in_place(&mut buf).unwrap_or(0.0) as f32;

        dev.clear();
        dev.extend(buf.iter().map(|v| (*v as f32 - median).abs() as f64));
        let mad = median_in_place(&mut dev).unwrap_or(0.0) as f32;

        if mad == 0.0 {
            stats[idx] = RobustStats {
                median,
                mad,
                scale: 1.0,
                n_valid: n_valid as u32,
                reliable: true,
                scale_degenerate: true,
                missing_fraction: canonicalize_f32(missing_fraction),
            };
            for cell_idx in 0..n_cells {
                let raw = raw_store.get(metric_id, cell_idx);
                if raw.is_finite() {
                    norm_store.set(metric_id, cell_idx, 0.0);
                } else {
                    norm_store.set(metric_id, cell_idx, f32::NAN);
                }
            }
            per_metric_summary.push(GlobalNormalizationMetricSummary {
                metric: spec.canonical_name.to_string(),
                median,
                mad,
                scale: 1.0,
                n_valid: n_valid as u32,
                reliable: true,
                scale_degenerate: true,
                missing_fraction: canonicalize_f32(missing_fraction),
            });
            continue;
        }

        let scale = canonicalize_f32(GAUSSIAN_MAD_SCALE * mad);
        stats[idx] = RobustStats {
            median: canonicalize_f32(median),
            mad: canonicalize_f32(mad),
            scale,
            n_valid: n_valid as u32,
            reliable: true,
            scale_degenerate: false,
            missing_fraction: canonicalize_f32(missing_fraction),
        };

        for cell_idx in 0..n_cells {
            let raw = raw_store.get(metric_id, cell_idx);
            if raw.is_finite() {
                let z = (raw - median) / (scale + EPS);
                norm_store.set(metric_id, cell_idx, canonicalize_f32(z));
            } else {
                norm_store.set(metric_id, cell_idx, f32::NAN);
            }
        }

        per_metric_summary.push(GlobalNormalizationMetricSummary {
            metric: spec.canonical_name.to_string(),
            median: canonicalize_f32(median),
            mad: canonicalize_f32(mad),
            scale,
            n_valid: n_valid as u32,
            reliable: true,
            scale_degenerate: false,
            missing_fraction: canonicalize_f32(missing_fraction),
        });
    }

    per_metric_summary.sort_by(|a, b| {
        let ida = canonical_to_metric
            .get(a.metric.as_str())
            .map(|id| id.as_index())
            .unwrap_or(usize::MAX);
        let idb = canonical_to_metric
            .get(b.metric.as_str())
            .map(|id| id.as_index())
            .unwrap_or(usize::MAX);
        ida.cmp(&idb).then(a.metric.cmp(&b.metric))
    });

    GlobalNormalizationContext {
        norm_store,
        stats,
        summary: GlobalNormalizationSummary {
            method: "median_mad".to_string(),
            eps: EPS,
            min_valid_cells: MIN_VALID_CELLS,
            per_metric: per_metric_summary,
        },
    }
}

pub fn robust_stats_for_metric(
    ctx: &GlobalNormalizationContext,
    metric_id: MetricId,
) -> Option<RobustStats> {
    ctx.stats.get(metric_id.as_index()).copied()
}

pub fn metric_name(metric_id: MetricId) -> &'static str {
    metric_spec_by_id(metric_id).canonical_name
}
