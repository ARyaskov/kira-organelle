use std::collections::BTreeMap;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::cells::types::{CellsState, MetricStore};
use crate::io::canonicalize_f32;
use crate::registry::metrics::{METRIC_SPECS, MetricId, metric_spec_by_id};
use crate::util::select::median_in_place;

pub const EPS: f32 = 1e-9;
pub const MIN_VALID_CELLS: usize = 20;
const GAUSSIAN_MAD_SCALE: f32 = 1.4826;

pub type NormStore = MetricStore;

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

struct MetricSlice {
    z: Vec<f32>,
    present: Vec<bool>,
}

struct PerMetricOutcome {
    metric_idx: usize,
    stats: RobustStats,
    slice: MetricSlice,
    summary: Option<GlobalNormalizationMetricSummary>,
}

pub fn compute_global_robust_normalization(cells: &CellsState) -> GlobalNormalizationContext {
    let n_cells = cells.n_cells;
    let mut raw_store = MetricStore::new(n_cells);

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

    let outcomes: Vec<PerMetricOutcome> = METRIC_SPECS
        .par_iter()
        .map(|spec| compute_one_metric(spec, &raw_store, n_cells))
        .collect();

    let mut norm_store = NormStore::new(n_cells);
    let mut stats = vec![default_stats(); MetricId::COUNT];
    let mut per_metric_summary =
        Vec::<GlobalNormalizationMetricSummary>::with_capacity(MetricId::COUNT);
    for outcome in outcomes {
        stats[outcome.metric_idx] = outcome.stats;
        let metric_id = METRIC_SPECS[outcome.metric_idx].id;
        for (cell_idx, (&z, &present)) in outcome
            .slice
            .z
            .iter()
            .zip(outcome.slice.present.iter())
            .enumerate()
        {
            if present {
                norm_store.set(metric_id, cell_idx, z);
            } else {
                norm_store.set(metric_id, cell_idx, f32::NAN);
            }
        }
        if let Some(s) = outcome.summary {
            per_metric_summary.push(s);
        }
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

fn default_stats() -> RobustStats {
    RobustStats {
        median: 0.0,
        mad: 0.0,
        scale: 1.0,
        n_valid: 0,
        reliable: false,
        scale_degenerate: false,
        missing_fraction: 1.0,
    }
}

fn compute_one_metric(
    spec: &crate::registry::metrics::MetricSpec,
    raw_store: &MetricStore,
    n_cells: usize,
) -> PerMetricOutcome {
    let metric_id = spec.id;
    let metric_idx = metric_id.as_index();

    let mut buf = Vec::<f64>::with_capacity(n_cells);
    let mut raw_cache = Vec::<f32>::with_capacity(n_cells);
    for cell_idx in 0..n_cells {
        let v = raw_store.get(metric_id, cell_idx);
        raw_cache.push(v);
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

    let mut slice = MetricSlice {
        z: vec![f32::NAN; n_cells],
        present: vec![false; n_cells],
    };

    if n_valid < MIN_VALID_CELLS {
        let stats = RobustStats {
            n_valid: n_valid as u32,
            missing_fraction: canonicalize_f32(missing_fraction),
            ..default_stats()
        };
        let summary = (n_valid > 0).then(|| GlobalNormalizationMetricSummary {
            metric: spec.canonical_name.to_string(),
            median: 0.0,
            mad: 0.0,
            scale: 1.0,
            n_valid: n_valid as u32,
            reliable: false,
            scale_degenerate: false,
            missing_fraction: canonicalize_f32(missing_fraction),
        });
        return PerMetricOutcome {
            metric_idx,
            stats,
            slice,
            summary,
        };
    }

    let median_f64 = median_in_place(&mut buf).unwrap_or(0.0);
    let median = median_f64 as f32;

    let mut dev = Vec::<f64>::with_capacity(buf.len());
    dev.extend(buf.iter().map(|v| (*v - median_f64).abs()));
    let mad = median_in_place(&mut dev).unwrap_or(0.0) as f32;

    if mad == 0.0 {
        let stats = RobustStats {
            median,
            mad,
            scale: 1.0,
            n_valid: n_valid as u32,
            reliable: true,
            scale_degenerate: true,
            missing_fraction: canonicalize_f32(missing_fraction),
        };
        for (cell_idx, raw) in raw_cache.iter().enumerate() {
            if raw.is_finite() {
                slice.z[cell_idx] = 0.0;
                slice.present[cell_idx] = true;
            }
        }
        let summary = Some(GlobalNormalizationMetricSummary {
            metric: spec.canonical_name.to_string(),
            median,
            mad,
            scale: 1.0,
            n_valid: n_valid as u32,
            reliable: true,
            scale_degenerate: true,
            missing_fraction: canonicalize_f32(missing_fraction),
        });
        return PerMetricOutcome {
            metric_idx,
            stats,
            slice,
            summary,
        };
    }

    let scale = canonicalize_f32(GAUSSIAN_MAD_SCALE * mad);
    let stats = RobustStats {
        median: canonicalize_f32(median),
        mad: canonicalize_f32(mad),
        scale,
        n_valid: n_valid as u32,
        reliable: true,
        scale_degenerate: false,
        missing_fraction: canonicalize_f32(missing_fraction),
    };

    let denom = scale + EPS;
    for (cell_idx, raw) in raw_cache.iter().enumerate() {
        if raw.is_finite() {
            let z = (raw - median) / denom;
            slice.z[cell_idx] = canonicalize_f32(z);
            slice.present[cell_idx] = true;
        }
    }

    let summary = Some(GlobalNormalizationMetricSummary {
        metric: spec.canonical_name.to_string(),
        median: canonicalize_f32(median),
        mad: canonicalize_f32(mad),
        scale,
        n_valid: n_valid as u32,
        reliable: true,
        scale_degenerate: false,
        missing_fraction: canonicalize_f32(missing_fraction),
    });
    PerMetricOutcome {
        metric_idx,
        stats,
        slice,
        summary,
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
