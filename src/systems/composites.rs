use std::collections::BTreeMap;
use std::time::Instant;

use crate::cells::types::CellsState;
use crate::io::canonicalize_f64;
use crate::registry::metrics::{METRIC_SPECS, MetricId};
use crate::util::select::{median_in_place, percentile_nearest_rank_in_place};
use serde::{Deserialize, Serialize};

const RARE_MIN_CLUSTER_SIZE: usize = 10;
const COV_REG_SCALE: f64 = 1e-3;
const SENS_DELTA: f32 = 0.25;
const LANDSCAPE_K: usize = 10;
const BASIN_MIN_SIZE: usize = 5;
const THERAPEUTIC_TOP_AXES: usize = 10;
const THERAPEUTIC_TOP_COMBINATIONS: usize = 5;
const PTP_DELTA: f32 = 0.5;
const SENSITIVITY_AXIS_ORDER: [MetricId; 15] = [
    MetricId::Rss,
    MetricId::Sii,
    MetricId::Osl,
    MetricId::Tsm,
    MetricId::Pcp,
    MetricId::Asm,
    MetricId::Msm,
    MetricId::Tpi,
    MetricId::Cci,
    MetricId::Pci,
    MetricId::Lci,
    MetricId::Hsi,
    MetricId::Apb,
    MetricId::Mcb,
    MetricId::Imsi,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum RegimeClass {
    BaselineCompensated = 0,
    ImmuneSuppressiveEmbedded = 1,
    AdaptiveHighStress = 2,
    HighStressFailedCompensation = 3,
    SystemicCollapseRisk = 4,
}

impl RegimeClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::BaselineCompensated => "BaselineCompensated",
            Self::ImmuneSuppressiveEmbedded => "ImmuneSuppressiveEmbedded",
            Self::AdaptiveHighStress => "AdaptiveHighStress",
            Self::HighStressFailedCompensation => "HighStress_FailedCompensation",
            Self::SystemicCollapseRisk => "SystemicCollapseRisk",
        }
    }
}

#[derive(Debug, Clone)]
pub struct SystemsMetricRow {
    pub cell_id: String,
    pub cluster: String,
    pub stress_vector: f64,
    pub compensation_deficit: f64,
    pub cpi: f64,
    pub afs: f64,
    pub imsc: f64,
    pub regime_class: String,
    pub regime_code: u8,
    pub rare_state: bool,
    pub mahalanobis_distance: f64,
    pub dominant_axis: String,
    pub dominant_axis_sensitivity: f64,
    pub dominant_vulnerable_axis: String,
    pub potential: f64,
    pub stability_gradient: f64,
    pub tpi_landscape: f64,
    pub basin_id: String,
    pub transition_candidate: bool,
    pub lsi: f64,
    pub max_trajectory: f64,
    pub bee: f64,
}

#[derive(Debug, Clone)]
pub struct SystemsNormalizedRow {
    pub cell_id: String,
    pub cluster: String,
    pub metric: String,
    pub raw_value: f64,
    pub z_global: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RobustSummary {
    pub median: f64,
    pub p10: f64,
    pub p90: f64,
    pub mad: f64,
}

#[derive(Debug, Clone)]
pub struct ClusterStats {
    pub cluster: String,
    pub n_cells: usize,
    pub metrics: BTreeMap<String, RobustSummary>,
    pub tail_fraction_cpi: f64,
    pub tail_fraction_stress: f64,
    pub heterogeneity_index: f64,
    pub rare_fraction: f64,
    pub regime_fraction: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AxisFragilityScore {
    pub axis: String,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragilitySummary {
    pub delta: f32,
    pub global_axis_ranking: Vec<AxisFragilityScore>,
    pub cluster_axis_ranking: BTreeMap<String, Vec<AxisFragilityScore>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TherapeuticAxisRanking {
    pub axis: String,
    pub toi: f64,
    pub cri: f64,
    pub tps: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TherapeuticGlobalAxisRanking {
    pub axis: String,
    pub tps: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TherapeuticCombination {
    pub axis_a: String,
    pub axis_b: String,
    pub csp: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TherapeuticProjection {
    pub cluster_rankings: BTreeMap<String, Vec<TherapeuticAxisRanking>>,
    pub global_axis_ranking: Vec<TherapeuticGlobalAxisRanking>,
    pub top_combinations: BTreeMap<String, Vec<TherapeuticCombination>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BasinSummary {
    pub basin_id: String,
    pub size: usize,
    pub depth: f64,
    pub width: f64,
    pub stability: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LandscapeSummary {
    pub num_basins: usize,
    pub top_basins: Vec<BasinSummary>,
    pub all_basins: Vec<BasinSummary>,
    pub basin_threshold: f64,
    pub transition_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BasinRobustnessSummary {
    pub basin_id: String,
    pub ars: f64,
    pub ars_norm: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DynamicStabilitySummary {
    pub mean_lsi: f64,
    pub global_sci_system: f64,
    pub basin_robustness: Vec<BasinRobustnessSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InvariantReport {
    pub total_checks: u64,
    pub violations: u64,
    pub violation_examples: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiologicalConstraintViolation {
    pub constraint: String,
    pub correlation: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiologicalConstraintsReport {
    pub correlations: BTreeMap<String, f64>,
    pub violations: Vec<BiologicalConstraintViolation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationSummary {
    pub invariants: InvariantReport,
    pub robustness_index: f64,
    pub effective_degrees_of_freedom: usize,
    pub condition_number: f64,
    pub biological_constraints: BiologicalConstraintsReport,
    pub sensitivity_stability_index: f64,
    pub overparameterization_flag: bool,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct SystemsStageProfile {
    pub composite_ms: f64,
    pub aggregate_ms: f64,
    pub fragility_ms: f64,
}

#[derive(Debug, Clone)]
pub struct SystemsStateModel {
    pub metrics_rows: Vec<SystemsMetricRow>,
    pub normalized_rows: Vec<SystemsNormalizedRow>,
    pub cluster_stats: Vec<ClusterStats>,
    pub global_regime_counts: BTreeMap<String, u64>,
    pub global_regime_fraction: BTreeMap<String, f64>,
    pub cpi_p50: f64,
    pub cpi_p90: f64,
    pub stress_p50: f64,
    pub stress_p90: f64,
    pub rare_fraction_global: f64,
    pub system_entropy: f64,
    pub normalized_entropy: f64,
    pub fragility: FragilitySummary,
    pub therapeutic_projection: TherapeuticProjection,
    pub landscape: LandscapeSummary,
    pub dynamic_stability: DynamicStabilitySummary,
    pub validation: ValidationSummary,
    pub compatibility_missing_metrics: Vec<String>,
    pub compatibility_degraded_axes: Vec<String>,
    pub missing_metric_counts: BTreeMap<String, u64>,
}

pub fn compute_stress_vector(
    s_nuclear: Option<f32>,
    s_splice: Option<f32>,
    s_mito: Option<f32>,
    s_translate: Option<f32>,
    s_proteo: Option<f32>,
    s_autolys: Option<f32>,
    s_microenv: Option<f32>,
) -> f32 {
    let x = [
        pos(s_nuclear),
        pos(s_splice),
        pos(s_mito),
        pos(s_translate),
        pos(s_proteo),
        pos(s_autolys),
        pos(s_microenv),
    ];
    let mut sum_sq = 0.0f32;
    for v in x {
        sum_sq += v * v;
    }
    sum_sq.sqrt()
}

pub fn compute_compensation_deficit(
    z_tpi: Option<f32>,
    z_cci: Option<f32>,
    z_pci: Option<f32>,
    z_osl: Option<f32>,
    z_lci: Option<f32>,
    z_rss: Option<f32>,
    z_sii: Option<f32>,
    z_hsi: Option<f32>,
    z_apb: Option<f32>,
) -> f32 {
    let gap_tp = match (z_tpi, z_cci, z_pci) {
        (Some(tpi), Some(cci), Some(pci)) => Some(tpi - ((cci + pci) / 2.0)),
        _ => None,
    };
    let gap_ma = match (z_osl, z_lci) {
        (Some(osl), Some(lci)) => Some(osl - lci),
        _ => None,
    };
    let gap_ns = match (z_rss, z_sii) {
        (Some(rss), Some(sii)) => Some(rss - sii),
        _ => None,
    };
    let gap_he = match (z_hsi, z_apb) {
        (Some(hsi), Some(apb)) => Some(hsi - apb),
        _ => None,
    };
    mean_positive([gap_tp, gap_ma, gap_ns, gap_he].iter().copied())
}

pub fn compute_cpi(stress_vector: f32, compensation_deficit: f32, z_pcp: Option<f32>) -> f32 {
    let mut weighted_sum = 0.0f32;
    let mut weight_sum = 0.0f32;

    weighted_sum += 0.3 * stress_vector.max(0.0);
    weight_sum += 0.3;
    weighted_sum += 0.3 * compensation_deficit.max(0.0);
    weight_sum += 0.3;
    if let Some(pcp) = z_pcp {
        weighted_sum += 0.4 * pcp.max(0.0);
        weight_sum += 0.4;
    }

    if weight_sum == 0.0 {
        0.0
    } else {
        (weighted_sum / weight_sum).max(0.0)
    }
}

pub fn compute_afs(
    compensation_deficit: f32,
    z_cci: Option<f32>,
    z_pci: Option<f32>,
    z_lci: Option<f32>,
    z_mcb: Option<f32>,
) -> f32 {
    let adaptation = mean_positive([z_cci, z_pci, z_lci, z_mcb].iter().copied());
    adaptation - compensation_deficit
}

pub fn compute_imsc(z_imsi: Option<f32>, z_msm: Option<f32>) -> f32 {
    let mut weighted_sum = 0.0f32;
    let mut weight_sum = 0.0f32;
    if let Some(v) = z_imsi {
        weighted_sum += 0.5 * v.max(0.0);
        weight_sum += 0.5;
    }
    if let Some(v) = z_msm {
        weighted_sum += 0.5 * v.max(0.0);
        weight_sum += 0.5;
    }
    if weight_sum == 0.0 {
        0.0
    } else {
        weighted_sum / weight_sum
    }
}

pub fn classify_regime(cpi: f32, stress_vector: f32, afs: f32, imsc: f32) -> RegimeClass {
    if cpi >= 3.0 {
        RegimeClass::SystemicCollapseRisk
    } else if stress_vector >= 3.0 && afs < 0.0 {
        RegimeClass::HighStressFailedCompensation
    } else if stress_vector >= 3.0 && afs >= 0.0 {
        RegimeClass::AdaptiveHighStress
    } else if imsc >= 2.0 {
        RegimeClass::ImmuneSuppressiveEmbedded
    } else {
        RegimeClass::BaselineCompensated
    }
}

pub fn compute_systems_state_model(cells: &CellsState) -> Option<SystemsStateModel> {
    compute_systems_state_model_profiled(cells).0
}

pub fn compute_systems_state_model_profiled(
    cells: &CellsState,
) -> (Option<SystemsStateModel>, SystemsStageProfile) {
    let Some(norm) = cells.normalization_context.as_ref() else {
        return (None, SystemsStageProfile::default());
    };
    if cells.cells.is_empty() {
        return (None, SystemsStageProfile::default());
    }
    let mut profile = SystemsStageProfile::default();
    let mut composite_ns: u128 = 0;
    let mut aggregate_ns: u128 = 0;
    let mut fragility_ns: u128 = 0;

    let canonical_to_metric = METRIC_SPECS
        .iter()
        .map(|spec| (spec.canonical_name, spec.id))
        .collect::<BTreeMap<_, _>>();

    let mut metrics_rows = Vec::<SystemsMetricRow>::with_capacity(cells.cells.len());
    let mut normalized_rows = Vec::<SystemsNormalizedRow>::new();
    let mut cluster_indices = BTreeMap::<String, Vec<usize>>::new();
    let mut regime_counts = BTreeMap::<String, u64>::new();
    let mut missing_metric_counts = BTreeMap::<String, u64>::new();
    let active_sensitivity_axes = SENSITIVITY_AXIS_ORDER
        .iter()
        .copied()
        .filter(|m| metric_reliable_and_present(norm, *m))
        .collect::<Vec<_>>();
    let compatibility_missing_metrics = SENSITIVITY_AXIS_ORDER
        .iter()
        .copied()
        .filter(|m| metric_missing(norm, *m))
        .map(|m| metric_name(m).to_string())
        .collect::<Vec<_>>();
    let compatibility_degraded_axes = SENSITIVITY_AXIS_ORDER
        .iter()
        .copied()
        .filter(|m| metric_degraded(norm, *m))
        .map(|m| metric_name(m).to_string())
        .collect::<Vec<_>>();
    let mut per_axis_abs_sensitivity =
        vec![vec![0.0f32; cells.cells.len()]; active_sensitivity_axes.len()];
    let mut per_axis_avs = vec![vec![0.0f32; cells.cells.len()]; active_sensitivity_axes.len()];

    for (cell_idx, cell) in cells.cells.iter().enumerate() {
        let cell_composite_t0 = Instant::now();
        let cluster = "all_cells".to_string();

        for org_state in cell.per_organelle.values() {
            for (axis, value) in &org_state.axes {
                let Some(metric_id) = canonical_to_metric.get(axis.as_str()) else {
                    continue;
                };
                let z = norm.norm_store.get(*metric_id, cell_idx);
                normalized_rows.push(SystemsNormalizedRow {
                    cell_id: cell.id.clone(),
                    cluster: cluster.clone(),
                    metric: axis.clone(),
                    raw_value: *value,
                    z_global: if z.is_finite() { z as f64 } else { f64::NAN },
                });
            }
        }

        let z_rss = z_at(norm, MetricId::Rss, cell_idx, &mut missing_metric_counts);
        let z_sii = z_at(norm, MetricId::Sii, cell_idx, &mut missing_metric_counts);
        let z_osl = z_at(norm, MetricId::Osl, cell_idx, &mut missing_metric_counts);
        let z_tsm = z_at(norm, MetricId::Tsm, cell_idx, &mut missing_metric_counts);
        let z_pcp = z_at(norm, MetricId::Pcp, cell_idx, &mut missing_metric_counts);
        let z_asm = z_at(norm, MetricId::Asm, cell_idx, &mut missing_metric_counts);
        let z_msm = z_at(norm, MetricId::Msm, cell_idx, &mut missing_metric_counts);
        let stress_vector = compute_stress_vector(z_rss, z_sii, z_osl, z_tsm, z_pcp, z_asm, z_msm);

        let z_tpi = z_at(norm, MetricId::Tpi, cell_idx, &mut missing_metric_counts);
        let z_cci = z_at(norm, MetricId::Cci, cell_idx, &mut missing_metric_counts);
        let z_pci = z_at(norm, MetricId::Pci, cell_idx, &mut missing_metric_counts);
        let z_lci = z_at(norm, MetricId::Lci, cell_idx, &mut missing_metric_counts);
        let z_hsi = z_at(norm, MetricId::Hsi, cell_idx, &mut missing_metric_counts);
        let z_apb = z_at(norm, MetricId::Apb, cell_idx, &mut missing_metric_counts);
        let compensation_deficit = compute_compensation_deficit(
            z_tpi, z_cci, z_pci, z_osl, z_lci, z_rss, z_sii, z_hsi, z_apb,
        );

        let cpi = compute_cpi(stress_vector, compensation_deficit, z_pcp);
        let z_mcb = z_at(norm, MetricId::Mcb, cell_idx, &mut missing_metric_counts);
        let afs = compute_afs(compensation_deficit, z_cci, z_pci, z_lci, z_mcb);
        let z_imsi = z_at(norm, MetricId::Imsi, cell_idx, &mut missing_metric_counts);
        let imsc = compute_imsc(z_imsi, z_msm);
        let regime = classify_regime(cpi, stress_vector, afs, imsc);

        let regime_name = regime.as_str().to_string();
        *regime_counts.entry(regime_name.clone()).or_insert(0) += 1;

        let idx = metrics_rows.len();
        metrics_rows.push(SystemsMetricRow {
            cell_id: cell.id.clone(),
            cluster: cluster.clone(),
            stress_vector: canonicalize_f64(stress_vector as f64),
            compensation_deficit: canonicalize_f64(compensation_deficit as f64),
            cpi: canonicalize_f64(cpi as f64),
            afs: canonicalize_f64(afs as f64),
            imsc: canonicalize_f64(imsc as f64),
            regime_class: regime_name,
            regime_code: regime as u8,
            rare_state: false,
            mahalanobis_distance: f64::NAN,
            dominant_axis: "NA".to_string(),
            dominant_axis_sensitivity: 0.0,
            dominant_vulnerable_axis: "NA".to_string(),
            potential: 0.0,
            stability_gradient: 0.0,
            tpi_landscape: 0.0,
            basin_id: "NA".to_string(),
            transition_candidate: false,
            lsi: 0.0,
            max_trajectory: 0.0,
            bee: 0.0,
        });
        composite_ns += cell_composite_t0.elapsed().as_nanos();
        if !active_sensitivity_axes.is_empty() {
            let frag_t0 = Instant::now();
            let mut best_axis: Option<MetricId> = None;
            let mut best_abs = 0.0f32;
            let mut best_vul_axis: Option<MetricId> = None;
            let mut best_avs = 0.0f32;
            for (axis_idx, axis_id) in active_sensitivity_axes.iter().copied().enumerate() {
                let x = norm.norm_store.get(axis_id, cell_idx);
                if !x.is_finite() {
                    per_axis_abs_sensitivity[axis_idx][cell_idx] = 0.0;
                    per_axis_avs[axis_idx][cell_idx] = 0.0;
                    continue;
                }
                let cpi_plus = cpi_from_z_for_cell(norm, cell_idx, Some((axis_id, x + SENS_DELTA)));
                let cpi_minus =
                    cpi_from_z_for_cell(norm, cell_idx, Some((axis_id, x - SENS_DELTA)));
                let sens = (cpi_plus - cpi_minus) / (2.0 * SENS_DELTA);
                let abs_sens = sens.abs();
                per_axis_abs_sensitivity[axis_idx][cell_idx] = abs_sens;
                let comp_capacity = compensatory_capacity_for_axis(norm, cell_idx, axis_id);
                let avs = compute_avs(x, abs_sens, comp_capacity);
                per_axis_avs[axis_idx][cell_idx] = avs;

                let replace = match best_axis {
                    None => true,
                    Some(prev) => {
                        abs_sens > best_abs
                            || (abs_sens == best_abs && axis_id.as_index() < prev.as_index())
                    }
                };
                if replace {
                    best_axis = Some(axis_id);
                    best_abs = abs_sens;
                }
                let replace_vul = match best_vul_axis {
                    None => true,
                    Some(prev) => {
                        avs > best_avs || (avs == best_avs && axis_id.as_index() < prev.as_index())
                    }
                };
                if replace_vul {
                    best_vul_axis = Some(axis_id);
                    best_avs = avs;
                }
            }
            if let Some(row) = metrics_rows.last_mut() {
                if let Some(axis) = best_axis {
                    row.dominant_axis = metric_name(axis).to_string();
                    row.dominant_axis_sensitivity = canonicalize_f64(best_abs as f64);
                }
                if let Some(axis) = best_vul_axis {
                    row.dominant_vulnerable_axis = metric_name(axis).to_string();
                }
            }
            fragility_ns += frag_t0.elapsed().as_nanos();
        }
        cluster_indices.entry(cluster).or_default().push(idx);
    }

    let landscape_t0 = Instant::now();
    let landscape = compute_landscape(&mut metrics_rows);
    aggregate_ns += landscape_t0.elapsed().as_nanos();

    normalized_rows.sort_by(|a, b| {
        a.cell_id
            .cmp(&b.cell_id)
            .then(a.metric.cmp(&b.metric))
            .then(a.cluster.cmp(&b.cluster))
    });

    let mut cpi_all = metrics_rows.iter().map(|r| r.cpi).collect::<Vec<_>>();
    let mut stress_all = metrics_rows
        .iter()
        .map(|r| r.stress_vector)
        .collect::<Vec<_>>();
    let cpi_global = robust_summary(&mut cpi_all);
    let stress_global = robust_summary(&mut stress_all);
    let t_cpi = canonicalize_f64(cpi_global.median + 2.0 * cpi_global.mad);
    let t_stress = canonicalize_f64(stress_global.median + 2.0 * stress_global.mad);
    let global_frag_t0 = Instant::now();
    let global_axis_ranking =
        rank_axis_fragility_global(&active_sensitivity_axes, &per_axis_abs_sensitivity, 10);
    fragility_ns += global_frag_t0.elapsed().as_nanos();

    let mut cluster_stats = Vec::<ClusterStats>::new();
    let mut cluster_axis_ranking = BTreeMap::<String, Vec<AxisFragilityScore>>::new();
    let mut cluster_therapeutic_rankings = BTreeMap::<String, Vec<TherapeuticAxisRanking>>::new();
    let mut cluster_top_combinations = BTreeMap::<String, Vec<TherapeuticCombination>>::new();
    let mut global_tps_per_axis = vec![Vec::<f64>::new(); active_sensitivity_axes.len()];
    let mut rare_cells = 0usize;

    let aggregate_t0 = Instant::now();
    for (cluster, indices) in &cluster_indices {
        let n = indices.len();
        let mut cpi_vals = Vec::<f64>::with_capacity(n);
        let mut stress_vals = Vec::<f64>::with_capacity(n);
        let mut comp_vals = Vec::<f64>::with_capacity(n);
        let mut afs_vals = Vec::<f64>::with_capacity(n);
        let mut imsc_vals = Vec::<f64>::with_capacity(n);
        let mut regime_local = BTreeMap::<String, u64>::new();
        let mut tail_cpi_n = 0usize;
        let mut tail_stress_n = 0usize;

        let mut welford_stress = Welford1::default();
        let mut welford_cpi = Welford1::default();
        let mut welford_afs = Welford1::default();

        for idx in indices {
            let row = &metrics_rows[*idx];
            cpi_vals.push(row.cpi);
            stress_vals.push(row.stress_vector);
            comp_vals.push(row.compensation_deficit);
            afs_vals.push(row.afs);
            imsc_vals.push(row.imsc);
            *regime_local.entry(row.regime_class.clone()).or_insert(0) += 1;
            if row.cpi >= t_cpi {
                tail_cpi_n += 1;
            }
            if row.stress_vector >= t_stress {
                tail_stress_n += 1;
            }
            welford_stress.push(row.stress_vector);
            welford_cpi.push(row.cpi);
            welford_afs.push(row.afs);
        }

        let mut metrics = BTreeMap::<String, RobustSummary>::new();
        metrics.insert("StressVector".to_string(), robust_summary(&mut stress_vals));
        metrics.insert("CPI".to_string(), robust_summary(&mut cpi_vals));
        metrics.insert(
            "CompensationDeficit".to_string(),
            robust_summary(&mut comp_vals),
        );
        metrics.insert("AFS".to_string(), robust_summary(&mut afs_vals));
        metrics.insert("IMSC".to_string(), robust_summary(&mut imsc_vals));

        let tail_fraction_cpi = if n == 0 {
            0.0
        } else {
            tail_cpi_n as f64 / n as f64
        };
        let tail_fraction_stress = if n == 0 {
            0.0
        } else {
            tail_stress_n as f64 / n as f64
        };

        let var_axes =
            (welford_stress.variance() + welford_cpi.variance() + welford_afs.variance()) / 3.0;
        let heterogeneity_index = canonicalize_f64(var_axes + tail_fraction_cpi);

        let (rare_fraction, rare_n) = compute_rare_state_for_cluster(indices, &mut metrics_rows);
        rare_cells += rare_n;
        let rank_t0 = Instant::now();
        cluster_axis_ranking.insert(
            cluster.clone(),
            rank_axis_fragility_cluster(
                indices,
                &active_sensitivity_axes,
                &per_axis_abs_sensitivity,
                5,
            ),
        );
        let therapeutic_rankings = compute_cluster_therapeutic_rankings(
            indices,
            &active_sensitivity_axes,
            &per_axis_avs,
            &per_axis_abs_sensitivity,
            &mut global_tps_per_axis,
        );
        let top_combinations = compute_cluster_csp(
            indices,
            &active_sensitivity_axes,
            &per_axis_avs,
            &per_axis_abs_sensitivity,
            THERAPEUTIC_TOP_COMBINATIONS,
        );
        cluster_therapeutic_rankings.insert(cluster.clone(), therapeutic_rankings);
        cluster_top_combinations.insert(cluster.clone(), top_combinations);
        fragility_ns += rank_t0.elapsed().as_nanos();

        let mut regime_fraction = BTreeMap::<String, f64>::new();
        for (k, v) in regime_local {
            regime_fraction.insert(k, canonicalize_f64(v as f64 / n as f64));
        }

        cluster_stats.push(ClusterStats {
            cluster: cluster.clone(),
            n_cells: n,
            metrics,
            tail_fraction_cpi,
            tail_fraction_stress,
            heterogeneity_index,
            rare_fraction,
            regime_fraction,
        });
    }
    aggregate_ns += aggregate_t0.elapsed().as_nanos();

    cluster_stats.sort_by(|a, b| cluster_key_cmp(&a.cluster, &b.cluster));

    let total = metrics_rows.len() as f64;
    let mut global_regime_fraction = BTreeMap::<String, f64>::new();
    for (k, v) in &regime_counts {
        global_regime_fraction.insert(k.clone(), canonicalize_f64(*v as f64 / total));
    }
    let system_entropy = canonicalize_f64(entropy(&global_regime_fraction));
    let normalized_entropy = if RegimeClass::SystemicCollapseRisk as usize + 1 > 1 {
        let denom = ((RegimeClass::SystemicCollapseRisk as usize + 1) as f64).ln();
        if denom > 0.0 {
            canonicalize_f64(system_entropy / denom)
        } else {
            0.0
        }
    } else {
        0.0
    };
    metrics_rows.sort_by(|a, b| a.cell_id.cmp(&b.cell_id).then(a.cluster.cmp(&b.cluster)));

    profile.composite_ms = composite_ns as f64 / 1_000_000.0;
    profile.aggregate_ms = aggregate_ns as f64 / 1_000_000.0;
    profile.fragility_ms = fragility_ns as f64 / 1_000_000.0;
    let global_therapeutic_ranking = compute_global_therapeutic_ranking(
        &active_sensitivity_axes,
        &global_tps_per_axis,
        THERAPEUTIC_TOP_AXES,
    );
    let dynamic_stability = compute_dynamic_stability(
        &mut metrics_rows,
        norm,
        &landscape,
        &cluster_stats,
        if total == 0.0 {
            0.0
        } else {
            canonicalize_f64(rare_cells as f64 / total)
        },
        system_entropy,
    );
    let validation = compute_validation_summary(
        &metrics_rows,
        norm,
        &landscape,
        &active_sensitivity_axes,
        &per_axis_abs_sensitivity,
    );
    (
        Some(SystemsStateModel {
            metrics_rows,
            normalized_rows,
            cluster_stats,
            global_regime_counts: regime_counts,
            global_regime_fraction,
            cpi_p50: canonicalize_f64(cpi_global.median),
            cpi_p90: canonicalize_f64(cpi_global.p90),
            stress_p50: canonicalize_f64(stress_global.median),
            stress_p90: canonicalize_f64(stress_global.p90),
            rare_fraction_global: if total == 0.0 {
                0.0
            } else {
                canonicalize_f64(rare_cells as f64 / total)
            },
            system_entropy,
            normalized_entropy,
            fragility: FragilitySummary {
                delta: SENS_DELTA,
                global_axis_ranking,
                cluster_axis_ranking,
            },
            therapeutic_projection: TherapeuticProjection {
                cluster_rankings: cluster_therapeutic_rankings,
                global_axis_ranking: global_therapeutic_ranking,
                top_combinations: cluster_top_combinations,
            },
            landscape,
            dynamic_stability,
            validation,
            compatibility_missing_metrics,
            compatibility_degraded_axes,
            missing_metric_counts,
        }),
        profile,
    )
}

fn compute_rare_state_for_cluster(
    indices: &[usize],
    metrics_rows: &mut [SystemsMetricRow],
) -> (f64, usize) {
    if indices.len() < RARE_MIN_CLUSTER_SIZE {
        for idx in indices {
            metrics_rows[*idx].rare_state = false;
            metrics_rows[*idx].mahalanobis_distance = f64::NAN;
        }
        return (0.0, 0);
    }

    let mut w = Welford4::default();
    for idx in indices {
        let r = &metrics_rows[*idx];
        w.push([r.stress_vector, r.cpi, r.afs, r.imsc]);
    }
    let mean = w.mean;
    let mut cov = w.covariance();
    let inv = match invert_4x4(cov) {
        Some(m) => m,
        None => {
            let tr = cov[0][0] + cov[1][1] + cov[2][2] + cov[3][3];
            let lambda = COV_REG_SCALE * (tr / 4.0);
            for (i, row) in cov.iter_mut().enumerate() {
                row[i] += lambda;
            }
            invert_4x4(cov).unwrap_or(identity_4())
        }
    };

    let mut distances = Vec::<f64>::with_capacity(indices.len());
    for idx in indices {
        let r = &metrics_rows[*idx];
        let d = mahalanobis([r.stress_vector, r.cpi, r.afs, r.imsc], mean, inv);
        distances.push(d);
    }

    let mut dist_copy = distances.clone();
    let dist_summary = robust_summary(&mut dist_copy);
    let threshold = canonicalize_f64(dist_summary.median + 2.0 * dist_summary.mad);

    let mut rare_n = 0usize;
    for (local_idx, idx) in indices.iter().enumerate() {
        let d = distances[local_idx];
        metrics_rows[*idx].mahalanobis_distance = canonicalize_f64(d);
        let rare = d.is_finite() && d >= threshold;
        metrics_rows[*idx].rare_state = rare;
        if rare {
            rare_n += 1;
        }
    }

    (
        canonicalize_f64(rare_n as f64 / indices.len() as f64),
        rare_n,
    )
}

fn compute_landscape(rows: &mut [SystemsMetricRow]) -> LandscapeSummary {
    let n = rows.len();
    if n == 0 {
        return LandscapeSummary {
            num_basins: 0,
            top_basins: Vec::new(),
            all_basins: Vec::new(),
            basin_threshold: 0.0,
            transition_fraction: 0.0,
        };
    }

    let mut features = Vec::<[f64; 4]>::with_capacity(n);
    let mut potentials = Vec::<f64>::with_capacity(n);
    for row in rows.iter_mut() {
        let f = [
            finite_or_zero(row.stress_vector),
            finite_or_zero(row.cpi),
            finite_or_zero(row.afs),
            finite_or_zero(row.imsc),
        ];
        let potential = canonicalize_f64((f[1] - f[2]).max(0.0));
        row.potential = potential;
        features.push(f);
        potentials.push(potential);
    }

    let k = LANDSCAPE_K.min(n.saturating_sub(1));
    if k == 0 {
        return LandscapeSummary {
            num_basins: 0,
            top_basins: Vec::new(),
            all_basins: Vec::new(),
            basin_threshold: 0.0,
            transition_fraction: 0.0,
        };
    }

    let mut order = (0..n).collect::<Vec<_>>();
    order.sort_by(|a, b| {
        feature_cmp(features[*a], features[*b])
            .then(rows[*a].cell_id.cmp(&rows[*b].cell_id))
            .then(a.cmp(b))
    });
    let mut pos = vec![0usize; n];
    for (rank, idx) in order.iter().copied().enumerate() {
        pos[idx] = rank;
    }

    let mut knn = vec![Vec::<(usize, f64)>::new(); n];
    let mut knn_edges = Vec::<f64>::with_capacity(n * k);
    for idx in 0..n {
        let neighbors = knn_for_idx(idx, &features, &order, &pos, k);
        for (_, d) in &neighbors {
            knn_edges.push(*d);
        }
        knn[idx] = neighbors;
    }

    for idx in 0..n {
        let neigh = &knn[idx];
        if neigh.is_empty() {
            rows[idx].stability_gradient = 0.0;
            continue;
        }
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for (j, d) in neigh {
            num += (rows[idx].potential - rows[*j].potential).abs();
            den += *d;
        }
        let grad = if den <= 0.0 { 0.0 } else { num / den };
        rows[idx].stability_gradient = canonicalize_f64(grad);
    }

    let mut pot_copy = potentials.clone();
    let t_basin = median_in_place(&mut pot_copy).unwrap_or(0.0);
    let mut edge_copy = knn_edges.clone();
    let eps = median_in_place(&mut edge_copy).unwrap_or(0.0);
    let low_energy = potentials.iter().map(|p| *p <= t_basin).collect::<Vec<_>>();

    let mut adjacency = vec![Vec::<usize>::new(); n];
    if eps > 0.0 {
        for i in 0..n {
            if !low_energy[i] {
                continue;
            }
            for (j, d) in &knn[i] {
                if !low_energy[*j] || *d >= eps {
                    continue;
                }
                adjacency[i].push(*j);
                adjacency[*j].push(i);
            }
        }
    }
    for edges in &mut adjacency {
        edges.sort_unstable();
        edges.dedup();
    }

    let mut visited = vec![false; n];
    let mut components = Vec::<Vec<usize>>::new();
    for i in 0..n {
        if !low_energy[i] || visited[i] {
            continue;
        }
        let mut stack = vec![i];
        visited[i] = true;
        let mut comp = Vec::<usize>::new();
        while let Some(cur) = stack.pop() {
            comp.push(cur);
            let nbrs = &adjacency[cur];
            for ni in (0..nbrs.len()).rev() {
                let nxt = nbrs[ni];
                if visited[nxt] {
                    continue;
                }
                visited[nxt] = true;
                stack.push(nxt);
            }
        }
        comp.sort_unstable();
        if comp.len() >= BASIN_MIN_SIZE {
            components.push(comp);
        }
    }

    let mut basin_infos = Vec::<(Vec<usize>, f64, f64, f64, usize)>::new();
    for comp in components {
        let mut comp_p = comp
            .iter()
            .map(|idx| rows[*idx].potential)
            .collect::<Vec<_>>();
        let med_p = median_in_place(&mut comp_p).unwrap_or(0.0);
        let depth = canonicalize_f64((t_basin - med_p).max(0.0));
        let width = basin_width_from_knn(&comp, &knn);
        let stability = canonicalize_f64(depth / (1.0 + width));
        let min_idx = comp[0];
        basin_infos.push((comp, depth, width, stability, min_idx));
    }

    basin_infos.sort_by(|a, b| {
        b.3.partial_cmp(&a.3)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.0.len().cmp(&a.0.len()))
            .then(a.4.cmp(&b.4))
    });

    let mut top_basins = Vec::<BasinSummary>::new();
    let mut all_basins = Vec::<BasinSummary>::new();
    for (rank, (cells, depth, width, stability, _)) in basin_infos.iter().enumerate() {
        let basin_id = format!("Basin_{}", rank + 1);
        for idx in cells {
            rows[*idx].basin_id = basin_id.clone();
        }
        let basin_summary = BasinSummary {
            basin_id: basin_id.clone(),
            size: cells.len(),
            depth: *depth,
            width: *width,
            stability: *stability,
        };
        all_basins.push(basin_summary.clone());
        if rank < 10 {
            top_basins.push(basin_summary);
        }
    }

    let mut tpi_vals = Vec::<f64>::with_capacity(n);
    for row in rows.iter_mut() {
        let modulator = 1.0 / (1.0 + (row.potential - t_basin).abs());
        let tpi = canonicalize_f64(row.stability_gradient * modulator);
        row.tpi_landscape = tpi;
        tpi_vals.push(tpi);
    }
    let mut tpi_copy = tpi_vals.clone();
    let tpi_median = median_in_place(&mut tpi_copy).unwrap_or(0.0);
    let mut dev = tpi_vals
        .iter()
        .map(|v| (*v - tpi_median).abs())
        .collect::<Vec<_>>();
    let tpi_mad = median_in_place(&mut dev).unwrap_or(0.0);
    let tpi_thr = canonicalize_f64(tpi_median + 2.0 * tpi_mad);

    let mut transition_n = 0usize;
    for row in rows.iter_mut() {
        let is_transition = row.tpi_landscape >= tpi_thr && row.tpi_landscape.is_finite();
        row.transition_candidate = is_transition;
        if is_transition {
            transition_n += 1;
        }
    }

    LandscapeSummary {
        num_basins: basin_infos.len(),
        top_basins,
        all_basins,
        basin_threshold: canonicalize_f64(t_basin),
        transition_fraction: canonicalize_f64(transition_n as f64 / n as f64),
    }
}

fn compute_dynamic_stability(
    rows: &mut [SystemsMetricRow],
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    landscape: &LandscapeSummary,
    cluster_stats: &[ClusterStats],
    rare_fraction_global: f64,
    system_entropy: f64,
) -> DynamicStabilitySummary {
    let n = rows.len();
    if n == 0 {
        return DynamicStabilitySummary {
            mean_lsi: 0.0,
            global_sci_system: 0.0,
            basin_robustness: Vec::new(),
        };
    }

    let mut features = Vec::<[f64; 4]>::with_capacity(n);
    for row in rows.iter() {
        features.push([
            finite_or_zero(row.stress_vector),
            finite_or_zero(row.cpi),
            finite_or_zero(row.afs),
            finite_or_zero(row.imsc),
        ]);
    }

    let k = LANDSCAPE_K.min(n.saturating_sub(1));
    let mut knn = vec![Vec::<(usize, f64)>::new(); n];
    if k > 0 {
        let mut order = (0..n).collect::<Vec<_>>();
        order.sort_by(|a, b| {
            feature_cmp(features[*a], features[*b])
                .then(rows[*a].cell_id.cmp(&rows[*b].cell_id))
                .then(a.cmp(b))
        });
        let mut pos = vec![0usize; n];
        for (rank, idx) in order.iter().copied().enumerate() {
            pos[idx] = rank;
        }
        for idx in 0..n {
            knn[idx] = knn_for_idx(idx, &features, &order, &pos, k);
        }
    }

    let stress_axes = [
        MetricId::Rss,
        MetricId::Sii,
        MetricId::Osl,
        MetricId::Tsm,
        MetricId::Pcp,
        MetricId::Asm,
        MetricId::Msm,
    ];
    let mut lsi_sum = 0.0f64;
    let t_basin = landscape.basin_threshold;
    for idx in 0..n {
        let neigh = &knn[idx];
        let mut disp_sum = 0.0f64;
        let mut potential_mean = 0.0f64;
        let mut potential_m2 = 0.0f64;
        let mut potential_n = 0usize;

        for (j, d) in neigh {
            disp_sum += *d;
            let p = rows[*j].potential;
            potential_n += 1;
            let delta = p - potential_mean;
            potential_mean += delta / potential_n as f64;
            potential_m2 += delta * (p - potential_mean);
        }
        let local_dispersion = if neigh.is_empty() {
            0.0
        } else {
            canonicalize_f64(disp_sum / neigh.len() as f64)
        };
        let local_potential_variance = if potential_n < 2 {
            0.0
        } else {
            canonicalize_f64(potential_m2 / (potential_n - 1) as f64)
        };
        let lsi = if neigh.is_empty() {
            1.0
        } else {
            canonicalize_f64(1.0 / (1.0 + local_dispersion * (1.0 + local_potential_variance)))
                .clamp(0.0, 1.0)
        };
        rows[idx].lsi = lsi;
        lsi_sum += lsi;

        let baseline = [
            rows[idx].stress_vector,
            rows[idx].cpi,
            rows[idx].afs,
            rows[idx].imsc,
        ];
        let mut max_traj = 0.0f64;
        for axis in stress_axes {
            let x = norm.norm_store.get(axis, idx);
            if !x.is_finite() {
                continue;
            }
            let perturbed = systems_from_z_for_cell(norm, idx, Some((axis, x + PTP_DELTA)));
            let traj = euclidean4(
                baseline,
                [
                    perturbed.0 as f64,
                    perturbed.1 as f64,
                    perturbed.2 as f64,
                    perturbed.3 as f64,
                ],
            );
            if traj > max_traj {
                max_traj = traj;
            }
        }
        rows[idx].max_trajectory = canonicalize_f64(max_traj);
        rows[idx].bee = if rows[idx].basin_id != "NA" && rows[idx].potential < t_basin {
            canonicalize_f64((t_basin - rows[idx].potential).max(0.0))
        } else {
            0.0
        };
    }

    let mean_lsi = canonicalize_f64(lsi_sum / n as f64);

    let mut gradients = rows
        .iter()
        .map(|r| finite_or_zero(r.stability_gradient))
        .collect::<Vec<_>>();
    let gradient_median = median_in_place(&mut gradients).unwrap_or(0.0);
    let gradient_mean = rows
        .iter()
        .map(|r| finite_or_zero(r.stability_gradient))
        .sum::<f64>()
        / n as f64;
    let norm_gradient = normalize_by_median(gradient_mean, gradient_median);

    let mut rare_vals = rows
        .iter()
        .map(|r| if r.rare_state { 1.0 } else { 0.0 })
        .collect::<Vec<_>>();
    let rare_median = median_in_place(&mut rare_vals).unwrap_or(0.0);
    let norm_rare = normalize_by_median(rare_fraction_global, rare_median);

    let mut cluster_entropy = Vec::<f64>::new();
    for cluster in cluster_stats {
        cluster_entropy.push(entropy(&cluster.regime_fraction));
    }
    let cluster_entropy_median = median_in_place(&mut cluster_entropy).unwrap_or(0.0);
    let norm_entropy = normalize_by_median(system_entropy, cluster_entropy_median);

    let global_sci_system =
        canonicalize_f64(0.4 * norm_gradient + 0.3 * norm_rare + 0.3 * norm_entropy)
            .clamp(0.0, 1.0);

    let max_ars = landscape
        .all_basins
        .iter()
        .map(|b| b.stability)
        .fold(0.0f64, f64::max);
    let basin_robustness = landscape
        .all_basins
        .iter()
        .map(|b| BasinRobustnessSummary {
            basin_id: b.basin_id.clone(),
            ars: canonicalize_f64(b.depth / (1.0 + b.width)),
            ars_norm: if max_ars > 0.0 {
                canonicalize_f64((b.stability / max_ars).clamp(0.0, 1.0))
            } else {
                0.0
            },
        })
        .collect::<Vec<_>>();

    DynamicStabilitySummary {
        mean_lsi,
        global_sci_system,
        basin_robustness,
    }
}

fn compute_validation_summary(
    rows: &[SystemsMetricRow],
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    landscape: &LandscapeSummary,
    active_axes: &[MetricId],
    per_axis_abs_sensitivity: &[Vec<f32>],
) -> ValidationSummary {
    let invariants = compute_invariants(rows, landscape.basin_threshold);
    let landscape_delta = compute_landscape_delta_proxy(rows);
    let fragility_delta = compute_fragility_delta(rows, norm, active_axes);
    let robustness_index =
        canonicalize_f64(1.0 / (1.0 + landscape_delta + fragility_delta)).clamp(0.0, 1.0);

    let (rank, condition_number) = covariance_rank_condition(norm, active_axes);
    let biological_constraints = biological_constraints_report(norm);

    let sensitivity_stability_index =
        compute_sensitivity_stability(rows, norm, active_axes, per_axis_abs_sensitivity);

    let composite_metrics = 5usize;
    let model_complexity = composite_metrics
        .saturating_mul(active_axes.len())
        .saturating_mul(landscape.num_basins.max(1));
    let overparameterization_flag = model_complexity > rank.max(1) * 10;

    ValidationSummary {
        invariants,
        robustness_index,
        effective_degrees_of_freedom: rank,
        condition_number,
        biological_constraints,
        sensitivity_stability_index,
        overparameterization_flag,
    }
}

fn compute_invariants(rows: &[SystemsMetricRow], t_basin: f64) -> InvariantReport {
    let mut total_checks = 0u64;
    let mut violations = 0u64;
    let mut examples = Vec::<String>::new();
    for row in rows {
        total_checks += 1;
        if row.cpi < 0.0 {
            violations += 1;
            if examples.len() < 10 {
                examples.push(row.cell_id.clone());
            }
        }
        total_checks += 1;
        if row.potential < t_basin && row.bee <= 0.0 {
            violations += 1;
            if examples.len() < 10 {
                examples.push(row.cell_id.clone());
            }
        }
        total_checks += 1;
        if row.regime_code == RegimeClass::SystemicCollapseRisk as u8 && row.cpi < 3.0 {
            violations += 1;
            if examples.len() < 10 {
                examples.push(row.cell_id.clone());
            }
        }
    }
    InvariantReport {
        total_checks,
        violations,
        violation_examples: examples,
    }
}

fn compute_landscape_delta_proxy(rows: &[SystemsMetricRow]) -> f64 {
    if rows.len() < 2 {
        return 0.0;
    }
    let features = rows
        .iter()
        .map(|r| [r.stress_vector, r.cpi, r.afs, r.imsc])
        .collect::<Vec<_>>();
    let n = rows.len();
    let mut order = (0..n).collect::<Vec<_>>();
    order.sort_by(|a, b| {
        feature_cmp(features[*a], features[*b])
            .then(rows[*a].cell_id.cmp(&rows[*b].cell_id))
            .then(a.cmp(b))
    });
    let mut pos = vec![0usize; n];
    for (rank, idx) in order.iter().copied().enumerate() {
        pos[idx] = rank;
    }

    let mut sum = 0.0f64;
    let mut m = 0usize;
    for idx in 0..n {
        let cpi8 = neighbor_mean_cpi(rows, &features, &order, &pos, idx, 8);
        let cpi10 = neighbor_mean_cpi(rows, &features, &order, &pos, idx, 10);
        sum += (cpi10 - cpi8).abs();
        m += 1;
    }
    if m == 0 {
        0.0
    } else {
        canonicalize_f64(sum / m as f64)
    }
}

fn neighbor_mean_cpi(
    rows: &[SystemsMetricRow],
    features: &[[f64; 4]],
    order: &[usize],
    pos: &[usize],
    idx: usize,
    k: usize,
) -> f64 {
    let kk = k.min(rows.len().saturating_sub(1));
    if kk == 0 {
        return rows[idx].cpi;
    }
    let nn = knn_for_idx(idx, features, order, pos, kk);
    let mut sum = rows[idx].cpi;
    let mut n = 1usize;
    for (j, _) in nn {
        sum += rows[j].cpi;
        n += 1;
    }
    canonicalize_f64(sum / n as f64)
}

fn compute_fragility_delta(
    rows: &[SystemsMetricRow],
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    active_axes: &[MetricId],
) -> f64 {
    if rows.is_empty() || active_axes.is_empty() {
        return 0.0;
    }
    let baseline = rank_axis_by_delta(rows, norm, active_axes, 0.25);
    let alt = rank_axis_by_delta(rows, norm, active_axes, 0.2);
    if baseline.is_empty() {
        return 0.0;
    }
    let mut shift_sum = 0.0f64;
    for (rank, axis) in baseline.iter().enumerate() {
        let alt_rank = alt.iter().position(|a| a == axis).unwrap_or(rank);
        shift_sum += (rank as i64 - alt_rank as i64).unsigned_abs() as f64;
    }
    canonicalize_f64(shift_sum / baseline.len() as f64)
}

fn rank_axis_by_delta(
    rows: &[SystemsMetricRow],
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    active_axes: &[MetricId],
    delta: f32,
) -> Vec<MetricId> {
    let mut scores = Vec::<(MetricId, f64)>::new();
    for axis in active_axes {
        let mut vals = Vec::<f64>::with_capacity(rows.len());
        for idx in 0..rows.len() {
            let x = norm.norm_store.get(*axis, idx);
            if !x.is_finite() {
                vals.push(0.0);
                continue;
            }
            let cpi_plus = cpi_from_z_for_cell(norm, idx, Some((*axis, x + delta)));
            let cpi_minus = cpi_from_z_for_cell(norm, idx, Some((*axis, x - delta)));
            vals.push(((cpi_plus - cpi_minus) / (2.0 * delta)).abs() as f64);
        }
        scores.push((*axis, median_in_place(&mut vals).unwrap_or(0.0)));
    }
    scores.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.as_index().cmp(&b.0.as_index()))
    });
    scores.into_iter().map(|(axis, _)| axis).collect()
}

fn covariance_rank_condition(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    active_axes: &[MetricId],
) -> (usize, f64) {
    let p = active_axes.len();
    if p == 0 {
        return (0, 0.0);
    }
    let n_cells = norm.norm_store.n_cells;
    let mut means = vec![0.0f64; p];
    let mut counts = vec![0usize; p];
    for (j, axis) in active_axes.iter().enumerate() {
        for c in 0..n_cells {
            let z = norm.norm_store.get(*axis, c);
            if z.is_finite() {
                means[j] += z as f64;
                counts[j] += 1;
            }
        }
        means[j] = if counts[j] == 0 {
            0.0
        } else {
            means[j] / counts[j] as f64
        };
    }
    let mut cov = vec![vec![0.0f64; p]; p];
    for c in 0..n_cells {
        let mut row = vec![0.0f64; p];
        let mut valid = true;
        for (j, axis) in active_axes.iter().enumerate() {
            let z = norm.norm_store.get(*axis, c);
            if !z.is_finite() {
                valid = false;
                break;
            }
            row[j] = z as f64 - means[j];
        }
        if !valid {
            continue;
        }
        for i in 0..p {
            for j in 0..p {
                cov[i][j] += row[i] * row[j];
            }
        }
    }
    let denom = (n_cells.saturating_sub(1)).max(1) as f64;
    for row in cov.iter_mut().take(p) {
        for v in row.iter_mut().take(p) {
            *v = canonicalize_f64(*v / denom);
        }
    }
    let eigen = jacobi_eigenvalues(cov);
    let mut max_e = 0.0f64;
    let mut min_pos = f64::INFINITY;
    let mut rank = 0usize;
    for e in eigen {
        if e > max_e {
            max_e = e;
        }
        if e > 1e-6 {
            rank += 1;
            if e < min_pos {
                min_pos = e;
            }
        }
    }
    let condition = if !min_pos.is_finite() || min_pos <= 0.0 {
        0.0
    } else {
        canonicalize_f64(max_e / min_pos)
    };
    (rank, condition)
}

fn jacobi_eigenvalues(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    if n == 0 {
        return Vec::new();
    }
    for _ in 0..64 {
        let mut p = 0usize;
        let mut q = 1usize.min(n - 1);
        let mut max = 0.0f64;
        for (i, row) in a.iter().enumerate().take(n) {
            for (j, val_ref) in row.iter().enumerate().skip(i + 1) {
                let val = val_ref.abs();
                if val > max {
                    max = val;
                    p = i;
                    q = j;
                }
            }
        }
        if max < 1e-12 {
            break;
        }
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let phi = 0.5 * (2.0 * apq).atan2(aqq - app);
        let c = phi.cos();
        let s = phi.sin();
        for row in a.iter_mut().take(n) {
            let aip = row[p];
            let aiq = row[q];
            row[p] = c * aip - s * aiq;
            row[q] = s * aip + c * aiq;
        }
        for j in 0..n {
            let apj = a[p][j];
            let aqj = a[q][j];
            a[p][j] = c * apj - s * aqj;
            a[q][j] = s * apj + c * aqj;
        }
        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
    }
    (0..n).map(|i| canonicalize_f64(a[i][i].max(0.0))).collect()
}

fn biological_constraints_report(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
) -> BiologicalConstraintsReport {
    let constraints = [
        ("C1_OSL_PCP_nonnegative", MetricId::Osl, MetricId::Pcp),
        ("C2_HSI_ISS_nonnegative", MetricId::Hsi, MetricId::Iss),
        ("C3_TSM_CCI_nonnegative", MetricId::Tsm, MetricId::Cci),
    ];
    let mut correlations = BTreeMap::<String, f64>::new();
    let mut violations = Vec::<BiologicalConstraintViolation>::new();
    for (name, a, b) in constraints {
        let corr = pearson_corr(norm, a, b);
        correlations.insert(name.to_string(), corr);
        if corr < -0.2 {
            violations.push(BiologicalConstraintViolation {
                constraint: name.to_string(),
                correlation: corr,
            });
        }
    }
    BiologicalConstraintsReport {
        correlations,
        violations,
    }
}

fn pearson_corr(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    a: MetricId,
    b: MetricId,
) -> f64 {
    let n = norm.norm_store.n_cells;
    let mut cnt = 0usize;
    let mut mean_a = 0.0f64;
    let mut mean_b = 0.0f64;
    let mut cov = 0.0f64;
    let mut m2a = 0.0f64;
    let mut m2b = 0.0f64;
    for i in 0..n {
        let za = norm.norm_store.get(a, i);
        let zb = norm.norm_store.get(b, i);
        if !za.is_finite() || !zb.is_finite() {
            continue;
        }
        let x = za as f64;
        let y = zb as f64;
        cnt += 1;
        let dx = x - mean_a;
        mean_a += dx / cnt as f64;
        let dy = y - mean_b;
        mean_b += dy / cnt as f64;
        cov += dx * (y - mean_b);
        m2a += dx * (x - mean_a);
        m2b += dy * (y - mean_b);
    }
    if cnt < 2 || m2a <= 0.0 || m2b <= 0.0 {
        0.0
    } else {
        canonicalize_f64((cov / (m2a.sqrt() * m2b.sqrt())).clamp(-1.0, 1.0))
    }
}

fn compute_sensitivity_stability(
    rows: &[SystemsMetricRow],
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    active_axes: &[MetricId],
    baseline_abs_sens: &[Vec<f32>],
) -> f64 {
    if rows.is_empty() || active_axes.is_empty() {
        return 1.0;
    }
    let mut order = (0..rows.len()).collect::<Vec<_>>();
    order.sort_by(|a, b| {
        rows[*b]
            .cpi
            .partial_cmp(&rows[*a].cpi)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(b))
    });
    let remove_n = ((rows.len() as f64) * 0.05).floor() as usize;
    let mut keep = vec![true; rows.len()];
    for idx in order.into_iter().take(remove_n) {
        keep[idx] = false;
    }

    let baseline = rank_axis_from_abs(active_axes, baseline_abs_sens, None);
    let reduced = rank_axis_from_abs(active_axes, baseline_abs_sens, Some(&keep));
    if baseline.is_empty() {
        return 1.0;
    }
    let mut shift_sum = 0.0f64;
    for (rank, axis) in baseline.iter().enumerate() {
        let r2 = reduced.iter().position(|a| a == axis).unwrap_or(rank);
        shift_sum += (rank as i64 - r2 as i64).unsigned_abs() as f64;
    }
    let mean_shift = shift_sum / baseline.len() as f64;
    let _ = norm;
    canonicalize_f64(1.0 / (1.0 + mean_shift)).clamp(0.0, 1.0)
}

fn rank_axis_from_abs(
    active_axes: &[MetricId],
    values: &[Vec<f32>],
    keep_mask: Option<&[bool]>,
) -> Vec<MetricId> {
    let mut scored = Vec::<(MetricId, f64)>::new();
    for (axis_idx, axis) in active_axes.iter().copied().enumerate() {
        let mut vals = Vec::<f64>::new();
        for cell_idx in 0..values[axis_idx].len() {
            if let Some(mask) = keep_mask
                && !mask[cell_idx]
            {
                continue;
            }
            vals.push(values[axis_idx][cell_idx] as f64);
        }
        scored.push((axis, median_in_place(&mut vals).unwrap_or(0.0)));
    }
    scored.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.as_index().cmp(&b.0.as_index()))
    });
    scored.into_iter().map(|(axis, _)| axis).collect()
}

fn knn_for_idx(
    idx: usize,
    features: &[[f64; 4]],
    order: &[usize],
    pos: &[usize],
    k: usize,
) -> Vec<(usize, f64)> {
    let n = order.len();
    let p = pos[idx];
    let target = features[idx];
    let mut left = p as isize - 1;
    let mut right = p + 1;
    let mut cand = Vec::<(f64, usize)>::with_capacity(k + 4);

    while left >= 0 || right < n {
        let lb_left = if left >= 0 {
            (target[0] - features[order[left as usize]][0]).abs()
        } else {
            f64::INFINITY
        };
        let lb_right = if right < n {
            (target[0] - features[order[right]][0]).abs()
        } else {
            f64::INFINITY
        };

        let worst = if cand.len() < k {
            f64::INFINITY
        } else {
            cand.last().map(|x| x.0).unwrap_or(f64::INFINITY)
        };
        if lb_left > worst && lb_right > worst && cand.len() >= k {
            break;
        }

        let take_left = lb_left <= lb_right;
        let next_rank = if take_left {
            let r = left as usize;
            left -= 1;
            r
        } else {
            let r = right;
            right += 1;
            r
        };
        let j = order[next_rank];
        if j == idx {
            continue;
        }
        let dist = euclidean4(target, features[j]);
        insert_top_k(&mut cand, (dist, j), k);
    }

    cand.into_iter().map(|(d, j)| (j, d)).collect()
}

fn insert_top_k(cand: &mut Vec<(f64, usize)>, item: (f64, usize), k: usize) {
    cand.push(item);
    cand.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.cmp(&b.1))
    });
    if cand.len() > k {
        cand.truncate(k);
    }
}

fn basin_width_from_knn(comp: &[usize], knn: &[Vec<(usize, f64)>]) -> f64 {
    let mut in_comp = vec![false; knn.len()];
    for idx in comp {
        in_comp[*idx] = true;
    }
    let mut sum = 0.0f64;
    let mut n = 0usize;
    for i in comp {
        for (j, d) in &knn[*i] {
            if in_comp[*j] {
                sum += *d;
                n += 1;
            }
        }
    }
    if n == 0 {
        0.0
    } else {
        canonicalize_f64(sum / n as f64)
    }
}

fn feature_cmp(a: [f64; 4], b: [f64; 4]) -> std::cmp::Ordering {
    a[0].partial_cmp(&b[0])
        .unwrap_or(std::cmp::Ordering::Equal)
        .then(a[1].partial_cmp(&b[1]).unwrap_or(std::cmp::Ordering::Equal))
        .then(a[2].partial_cmp(&b[2]).unwrap_or(std::cmp::Ordering::Equal))
        .then(a[3].partial_cmp(&b[3]).unwrap_or(std::cmp::Ordering::Equal))
}

fn euclidean4(a: [f64; 4], b: [f64; 4]) -> f64 {
    let d0 = a[0] - b[0];
    let d1 = a[1] - b[1];
    let d2 = a[2] - b[2];
    let d3 = a[3] - b[3];
    canonicalize_f64((d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3).sqrt())
}

fn finite_or_zero(v: f64) -> f64 {
    if v.is_finite() { v } else { 0.0 }
}

fn cpi_from_z_for_cell(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    cell_idx: usize,
    override_metric: Option<(MetricId, f32)>,
) -> f32 {
    systems_from_z_for_cell(norm, cell_idx, override_metric).1
}

fn systems_from_z_for_cell(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    cell_idx: usize,
    override_metric: Option<(MetricId, f32)>,
) -> (f32, f32, f32, f32) {
    let z = |metric: MetricId| -> Option<f32> {
        if let Some((override_id, override_value)) = override_metric {
            if override_id == metric {
                return if override_value.is_finite() {
                    Some(override_value)
                } else {
                    None
                };
            }
        }
        let v = norm.norm_store.get(metric, cell_idx);
        if v.is_finite() { Some(v) } else { None }
    };

    let z_rss = z(MetricId::Rss);
    let z_sii = z(MetricId::Sii);
    let z_osl = z(MetricId::Osl);
    let z_tsm = z(MetricId::Tsm);
    let z_pcp = z(MetricId::Pcp);
    let z_asm = z(MetricId::Asm);
    let z_msm = z(MetricId::Msm);
    let stress_vector = compute_stress_vector(z_rss, z_sii, z_osl, z_tsm, z_pcp, z_asm, z_msm);

    let z_tpi = z(MetricId::Tpi);
    let z_cci = z(MetricId::Cci);
    let z_pci = z(MetricId::Pci);
    let z_lci = z(MetricId::Lci);
    let z_hsi = z(MetricId::Hsi);
    let z_apb = z(MetricId::Apb);
    let compensation_deficit = compute_compensation_deficit(
        z_tpi, z_cci, z_pci, z_osl, z_lci, z_rss, z_sii, z_hsi, z_apb,
    );
    let cpi = compute_cpi(stress_vector, compensation_deficit, z_pcp);
    let z_mcb = z(MetricId::Mcb);
    let afs = compute_afs(compensation_deficit, z_cci, z_pci, z_lci, z_mcb);
    let z_imsi = z(MetricId::Imsi);
    let imsc = compute_imsc(z_imsi, z_msm);
    (stress_vector, cpi, afs, imsc)
}

fn metric_reliable_and_present(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    metric: MetricId,
) -> bool {
    let reliable = norm
        .stats
        .get(metric.as_index())
        .map(|s| s.reliable)
        .unwrap_or(false);
    reliable && norm.norm_store.metric_present(metric)
}

fn metric_missing(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    metric: MetricId,
) -> bool {
    norm.stats
        .get(metric.as_index())
        .map(|s| s.n_valid == 0)
        .unwrap_or(true)
}

fn metric_degraded(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    metric: MetricId,
) -> bool {
    norm.stats
        .get(metric.as_index())
        .map(|s| s.n_valid > 0 && !s.reliable)
        .unwrap_or(false)
}

fn rank_axis_fragility_global(
    axes: &[MetricId],
    per_axis_abs_sensitivity: &[Vec<f32>],
    top_n: usize,
) -> Vec<AxisFragilityScore> {
    let mut scored = Vec::<(MetricId, f64)>::with_capacity(axes.len());
    for (axis_idx, axis) in axes.iter().copied().enumerate() {
        let mut vals = per_axis_abs_sensitivity[axis_idx]
            .iter()
            .map(|v| *v as f64)
            .collect::<Vec<_>>();
        let score = median_in_place(&mut vals).unwrap_or(0.0);
        scored.push((axis, score));
    }
    sort_axis_scores(&mut scored);
    scored
        .into_iter()
        .take(top_n)
        .map(|(axis, score)| AxisFragilityScore {
            axis: metric_name(axis).to_string(),
            score: canonicalize_f64(score),
        })
        .collect()
}

fn rank_axis_fragility_cluster(
    indices: &[usize],
    axes: &[MetricId],
    per_axis_abs_sensitivity: &[Vec<f32>],
    top_n: usize,
) -> Vec<AxisFragilityScore> {
    let mut scored = Vec::<(MetricId, f64)>::with_capacity(axes.len());
    for (axis_idx, axis) in axes.iter().copied().enumerate() {
        let mut vals = Vec::<f64>::with_capacity(indices.len());
        for idx in indices {
            vals.push(per_axis_abs_sensitivity[axis_idx][*idx] as f64);
        }
        let score = median_in_place(&mut vals).unwrap_or(0.0);
        scored.push((axis, score));
    }
    sort_axis_scores(&mut scored);
    scored
        .into_iter()
        .take(top_n)
        .map(|(axis, score)| AxisFragilityScore {
            axis: metric_name(axis).to_string(),
            score: canonicalize_f64(score),
        })
        .collect()
}

fn compute_cluster_therapeutic_rankings(
    indices: &[usize],
    axes: &[MetricId],
    per_axis_avs: &[Vec<f32>],
    per_axis_abs_sensitivity: &[Vec<f32>],
    global_tps_per_axis: &mut [Vec<f64>],
) -> Vec<TherapeuticAxisRanking> {
    if axes.is_empty() || indices.is_empty() {
        return Vec::new();
    }
    let mut toi = vec![0.0f64; axes.len()];
    for axis_idx in 0..axes.len() {
        let mut vals = Vec::<f64>::with_capacity(indices.len());
        for idx in indices {
            vals.push(per_axis_avs[axis_idx][*idx] as f64);
        }
        toi[axis_idx] = canonicalize_f64(median_in_place(&mut vals).unwrap_or(0.0));
    }
    let mut toi_buf = toi.clone();
    let global_median_toi = median_in_place(&mut toi_buf).unwrap_or(0.0);

    let mut out = Vec::<(MetricId, f64, f64, f64)>::with_capacity(axes.len());
    for axis_idx in 0..axes.len() {
        let toi_norm = canonicalize_f64(toi[axis_idx] / (1.0 + global_median_toi));
        let cri = compute_cluster_cri(axis_idx, axes, indices, per_axis_abs_sensitivity);
        let tps = canonicalize_f64((toi_norm * (1.0 - cri)).max(0.0));
        out.push((axes[axis_idx], toi_norm, cri, tps));
        global_tps_per_axis[axis_idx].push(tps);
    }

    out.sort_by(|a, b| {
        b.3.partial_cmp(&a.3)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.as_index().cmp(&b.0.as_index()))
    });
    out.into_iter()
        .map(|(axis, toi_v, cri_v, tps_v)| TherapeuticAxisRanking {
            axis: metric_name(axis).to_string(),
            toi: toi_v,
            cri: cri_v,
            tps: tps_v,
        })
        .collect()
}

fn compute_global_therapeutic_ranking(
    axes: &[MetricId],
    global_tps_per_axis: &[Vec<f64>],
    top_n: usize,
) -> Vec<TherapeuticGlobalAxisRanking> {
    let mut scored = Vec::<(MetricId, f64)>::with_capacity(axes.len());
    for (axis_idx, axis) in axes.iter().copied().enumerate() {
        let mut vals = global_tps_per_axis
            .get(axis_idx)
            .cloned()
            .unwrap_or_default();
        let score = median_in_place(&mut vals).unwrap_or(0.0);
        scored.push((axis, canonicalize_f64(score)));
    }
    scored.sort_by(|(a_id, a), (b_id, b)| {
        b.partial_cmp(a)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a_id.as_index().cmp(&b_id.as_index()))
    });
    scored
        .into_iter()
        .take(top_n)
        .map(|(axis, tps)| TherapeuticGlobalAxisRanking {
            axis: metric_name(axis).to_string(),
            tps,
        })
        .collect()
}

fn compute_cluster_csp(
    indices: &[usize],
    axes: &[MetricId],
    per_axis_avs: &[Vec<f32>],
    per_axis_abs_sensitivity: &[Vec<f32>],
    top_n: usize,
) -> Vec<TherapeuticCombination> {
    if axes.len() < 2 || indices.is_empty() {
        return Vec::new();
    }
    let mut toi = vec![0.0f64; axes.len()];
    for axis_idx in 0..axes.len() {
        let mut vals = Vec::<f64>::with_capacity(indices.len());
        for idx in indices {
            vals.push(per_axis_avs[axis_idx][*idx] as f64);
        }
        toi[axis_idx] = canonicalize_f64(median_in_place(&mut vals).unwrap_or(0.0));
    }

    let mut combos = Vec::<(MetricId, MetricId, f64)>::new();
    for i in 0..axes.len() {
        for j in (i + 1)..axes.len() {
            let corr = sensitivity_corr(
                indices,
                &per_axis_abs_sensitivity[i],
                &per_axis_abs_sensitivity[j],
            )
            .abs();
            let csp = canonicalize_f64((toi[i] + toi[j]) * corr);
            combos.push((axes[i], axes[j], csp));
        }
    }
    combos.sort_by(|a, b| {
        b.2.partial_cmp(&a.2)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.as_index().cmp(&b.0.as_index()))
            .then(a.1.as_index().cmp(&b.1.as_index()))
    });
    combos
        .into_iter()
        .take(top_n)
        .map(|(a, b, csp)| TherapeuticCombination {
            axis_a: metric_name(a).to_string(),
            axis_b: metric_name(b).to_string(),
            csp,
        })
        .collect()
}

fn compute_cluster_cri(
    axis_idx: usize,
    axes: &[MetricId],
    indices: &[usize],
    per_axis_abs_sensitivity: &[Vec<f32>],
) -> f64 {
    let mut sum = 0.0f64;
    let mut n = 0usize;
    for comp in compensatory_axes_for_metric(axes[axis_idx]) {
        let Some(comp_idx) = axes.iter().position(|a| *a == *comp) else {
            continue;
        };
        let corr = sensitivity_corr(
            indices,
            &per_axis_abs_sensitivity[axis_idx],
            &per_axis_abs_sensitivity[comp_idx],
        );
        if corr.is_finite() {
            sum += corr.max(0.0);
            n += 1;
        }
    }
    if n == 0 {
        0.0
    } else {
        canonicalize_f64((sum / n as f64).clamp(0.0, 1.0))
    }
}

fn sensitivity_corr(indices: &[usize], a: &[f32], b: &[f32]) -> f64 {
    if indices.len() < 2 {
        return 0.0;
    }
    let mut n = 0usize;
    let mut mean_a = 0.0f64;
    let mut mean_b = 0.0f64;
    let mut c = 0.0f64;
    let mut m2_a = 0.0f64;
    let mut m2_b = 0.0f64;
    for idx in indices {
        let x = a[*idx] as f64;
        let y = b[*idx] as f64;
        n += 1;
        let dx = x - mean_a;
        mean_a += dx / n as f64;
        let dy = y - mean_b;
        mean_b += dy / n as f64;
        c += dx * (y - mean_b);
        m2_a += dx * (x - mean_a);
        m2_b += dy * (y - mean_b);
    }
    if n < 2 || m2_a <= 0.0 || m2_b <= 0.0 {
        return 0.0;
    }
    canonicalize_f64((c / (m2_a.sqrt() * m2_b.sqrt())).clamp(-1.0, 1.0))
}

fn compute_avs(stress_z: f32, sensitivity: f32, comp_capacity: Option<f32>) -> f32 {
    if !stress_z.is_finite() || !sensitivity.is_finite() {
        return 0.0;
    }
    let stress = stress_z.max(0.0);
    let sens = sensitivity.max(0.0);
    if stress == 0.0 || sens == 0.0 {
        return 0.0;
    }
    let base = stress * sens;
    match comp_capacity {
        Some(c) if c.is_finite() => base * (1.0 / (1.0 + c.max(0.0))),
        _ => base,
    }
}

fn compensatory_capacity_for_axis(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    cell_idx: usize,
    axis: MetricId,
) -> Option<f32> {
    match axis {
        MetricId::Osl => z_pos(norm, MetricId::Lci, cell_idx),
        MetricId::Tsm | MetricId::Pcp | MetricId::Tpi => {
            mean_z_pos2(norm, MetricId::Cci, MetricId::Pci, cell_idx)
        }
        MetricId::Rss => z_pos(norm, MetricId::Sii, cell_idx),
        MetricId::Sii => z_pos(norm, MetricId::Rss, cell_idx),
        MetricId::Asm | MetricId::Mcb => z_pos(norm, MetricId::Lci, cell_idx),
        MetricId::Msm | MetricId::Imsi => z_pos(norm, MetricId::Apb, cell_idx),
        MetricId::Hsi => z_pos(norm, MetricId::Apb, cell_idx),
        MetricId::Apb => z_pos(norm, MetricId::Hsi, cell_idx),
        MetricId::Cci => z_pos(norm, MetricId::Pci, cell_idx),
        MetricId::Pci => z_pos(norm, MetricId::Cci, cell_idx),
        MetricId::Lci => z_pos(norm, MetricId::Mcb, cell_idx),
        _ => None,
    }
}

fn compensatory_axes_for_metric(axis: MetricId) -> &'static [MetricId] {
    match axis {
        MetricId::Osl => &[MetricId::Lci, MetricId::Cci, MetricId::Pci],
        MetricId::Tsm | MetricId::Pcp | MetricId::Tpi => &[MetricId::Cci, MetricId::Pci],
        MetricId::Rss => &[MetricId::Sii],
        MetricId::Sii => &[MetricId::Rss],
        MetricId::Asm | MetricId::Mcb => &[MetricId::Lci],
        MetricId::Msm | MetricId::Imsi => &[MetricId::Apb],
        MetricId::Hsi => &[MetricId::Apb],
        MetricId::Apb => &[MetricId::Hsi],
        MetricId::Cci => &[MetricId::Pci],
        MetricId::Pci => &[MetricId::Cci],
        MetricId::Lci => &[MetricId::Mcb],
        _ => &[],
    }
}

fn z_pos(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    metric: MetricId,
    cell_idx: usize,
) -> Option<f32> {
    let v = norm.norm_store.get(metric, cell_idx);
    if v.is_finite() {
        Some(v.max(0.0))
    } else {
        None
    }
}

fn mean_z_pos2(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    a: MetricId,
    b: MetricId,
    cell_idx: usize,
) -> Option<f32> {
    match (z_pos(norm, a, cell_idx), z_pos(norm, b, cell_idx)) {
        (Some(x), Some(y)) => Some((x + y) * 0.5),
        (Some(x), None) => Some(x),
        (None, Some(y)) => Some(y),
        (None, None) => None,
    }
}

fn normalize_by_median(value: f64, median: f64) -> f64 {
    let v = finite_or_zero(value);
    let m = finite_or_zero(median);
    canonicalize_f64((v / (v + m + 1e-9)).clamp(0.0, 1.0))
}

fn sort_axis_scores(scored: &mut [(MetricId, f64)]) {
    scored.sort_by(|(a_id, a), (b_id, b)| {
        b.partial_cmp(a)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a_id.as_index().cmp(&b_id.as_index()))
    });
}

fn cluster_key_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    match (a.parse::<i64>(), b.parse::<i64>()) {
        (Ok(ia), Ok(ib)) => ia.cmp(&ib).then(a.cmp(b)),
        _ => a.cmp(b),
    }
}

fn z_at(
    norm: &crate::normalize::global_robust::GlobalNormalizationContext,
    metric: MetricId,
    cell_idx: usize,
    missing_counts: &mut BTreeMap<String, u64>,
) -> Option<f32> {
    let z = norm.norm_store.get(metric, cell_idx);
    if z.is_finite() {
        Some(z)
    } else {
        *missing_counts
            .entry(metric_name(metric).to_string())
            .or_insert(0) += 1;
        None
    }
}

fn metric_name(metric: MetricId) -> &'static str {
    match metric {
        MetricId::Rss => "RSS",
        MetricId::Sii => "SII",
        MetricId::Osl => "OSL",
        MetricId::Tsm => "TSM",
        MetricId::Pcp => "PCP",
        MetricId::Asm => "ASM",
        MetricId::Msm => "MSM",
        MetricId::Tpi => "TPI",
        MetricId::Cci => "CCI",
        MetricId::Pci => "PCI",
        MetricId::Lci => "LCI",
        MetricId::Hsi => "HSI",
        MetricId::Apb => "APB",
        MetricId::Mcb => "MCB",
        MetricId::Imsi => "IMSI",
        _ => "UNKNOWN",
    }
}

fn pos(v: Option<f32>) -> f32 {
    match v {
        Some(x) if x.is_finite() => x.max(0.0),
        _ => 0.0,
    }
}

fn mean_positive<'a>(values: impl Iterator<Item = Option<f32>> + 'a) -> f32 {
    let mut sum = 0.0f32;
    let mut n = 0usize;
    for v in values {
        if let Some(x) = v {
            sum += x.max(0.0);
            n += 1;
        }
    }
    if n == 0 { 0.0 } else { sum / n as f32 }
}

fn robust_summary(values: &mut [f64]) -> RobustSummary {
    if values.is_empty() {
        return RobustSummary {
            median: 0.0,
            p10: 0.0,
            p90: 0.0,
            mad: 0.0,
        };
    }
    let mut med_buf = values.to_vec();
    let median = median_in_place(&mut med_buf).unwrap_or(0.0);

    let mut p10_buf = values.to_vec();
    let p10 = percentile_nearest_rank_in_place(&mut p10_buf, 0.10).unwrap_or(0.0);
    let mut p90_buf = values.to_vec();
    let p90 = percentile_nearest_rank_in_place(&mut p90_buf, 0.90).unwrap_or(0.0);

    let mut dev = values
        .iter()
        .map(|v| (v - median).abs())
        .collect::<Vec<_>>();
    let mad = median_in_place(&mut dev).unwrap_or(0.0);

    RobustSummary {
        median: canonicalize_f64(median),
        p10: canonicalize_f64(p10),
        p90: canonicalize_f64(p90),
        mad: canonicalize_f64(mad),
    }
}

#[derive(Default)]
struct Welford1 {
    n: usize,
    mean: f64,
    m2: f64,
}

impl Welford1 {
    fn push(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        let delta2 = x - self.mean;
        self.m2 += delta * delta2;
    }

    fn variance(&self) -> f64 {
        if self.n < 2 {
            0.0
        } else {
            self.m2 / (self.n as f64 - 1.0)
        }
    }
}

#[derive(Default)]
struct Welford4 {
    n: usize,
    mean: [f64; 4],
    m2: [[f64; 4]; 4],
}

impl Welford4 {
    fn push(&mut self, x: [f64; 4]) {
        self.n += 1;
        let mut delta = [0.0f64; 4];
        for i in 0..4 {
            delta[i] = x[i] - self.mean[i];
            self.mean[i] += delta[i] / self.n as f64;
        }
        let mut delta2 = [0.0f64; 4];
        for i in 0..4 {
            delta2[i] = x[i] - self.mean[i];
        }
        for i in 0..4 {
            for (j, val) in delta2.iter().enumerate() {
                self.m2[i][j] += delta[i] * *val;
            }
        }
    }

    fn covariance(&self) -> [[f64; 4]; 4] {
        if self.n < 2 {
            return identity_4();
        }
        let mut out = [[0.0f64; 4]; 4];
        let denom = self.n as f64 - 1.0;
        for (i, row) in out.iter_mut().enumerate() {
            for (j, item) in row.iter_mut().enumerate() {
                *item = self.m2[i][j] / denom;
            }
        }
        out
    }
}

fn mahalanobis(x: [f64; 4], mean: [f64; 4], inv_cov: [[f64; 4]; 4]) -> f64 {
    let mut d = [0.0f64; 4];
    for i in 0..4 {
        d[i] = x[i] - mean[i];
    }
    let mut tmp = [0.0f64; 4];
    for i in 0..4 {
        for (j, dv) in d.iter().enumerate() {
            tmp[i] += inv_cov[i][j] * *dv;
        }
    }
    let mut quad = 0.0f64;
    for i in 0..4 {
        quad += d[i] * tmp[i];
    }
    quad.max(0.0).sqrt()
}

fn invert_4x4(a: [[f64; 4]; 4]) -> Option<[[f64; 4]; 4]> {
    let mut aug = [[0.0f64; 8]; 4];
    for i in 0..4 {
        for j in 0..4 {
            aug[i][j] = a[i][j];
        }
        aug[i][4 + i] = 1.0;
    }

    for col in 0..4 {
        let mut pivot_row = col;
        let mut pivot_abs = aug[col][col].abs();
        for (r, row) in aug.iter().enumerate().skip(col + 1) {
            let v = row[col].abs();
            if v > pivot_abs {
                pivot_abs = v;
                pivot_row = r;
            }
        }
        if pivot_abs <= 1e-12 {
            return None;
        }
        if pivot_row != col {
            aug.swap(pivot_row, col);
        }

        let piv = aug[col][col];
        for c in 0..8 {
            aug[col][c] /= piv;
        }

        for r in 0..4 {
            if r == col {
                continue;
            }
            let factor = aug[r][col];
            if factor == 0.0 {
                continue;
            }
            for c in 0..8 {
                aug[r][c] -= factor * aug[col][c];
            }
        }
    }

    let mut inv = [[0.0f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            inv[i][j] = aug[i][4 + j];
        }
    }
    Some(inv)
}

fn identity_4() -> [[f64; 4]; 4] {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

fn entropy(fractions: &BTreeMap<String, f64>) -> f64 {
    let mut h = 0.0f64;
    for p in fractions.values() {
        if *p > 0.0 {
            h -= p * p.ln();
        }
    }
    h
}
