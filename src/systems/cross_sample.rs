use std::collections::BTreeMap;
use std::io::BufRead;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::io::canonicalize_f64;
use crate::util::select::median_in_place;

const EPS: f64 = 1e-9;
const MAD_SCALE: f64 = 1.4826;
const DSP_AXIS_ORDER: [&str; 7] = ["RSS", "SII", "OSL", "TSM", "PCP", "ASM", "MSM"];
const REGIME_ORDER: [&str; 5] = [
    "BaselineCompensated",
    "ImmuneSuppressiveEmbedded",
    "AdaptiveHighStress",
    "HighStress_FailedCompensation",
    "SystemicCollapseRisk",
];

#[derive(Debug, Clone)]
pub struct SampleSystemsInput {
    pub sample: String,
    pub metrics_rows: Vec<MetricsRow>,
    pub axis_raw_values: BTreeMap<String, Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct MetricsRow {
    pub cell_id: String,
    pub stress_vector: f64,
    pub cpi: f64,
    pub afs: f64,
    pub imsc: f64,
    pub regime_class: String,
    pub basin_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DspAxisDelta {
    pub axis: String,
    pub value: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DspPairVector {
    pub sample_a: String,
    pub sample_b: String,
    pub vector: Vec<DspAxisDelta>,
    pub norm: f64,
    pub norm_normalized: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrossSampleAnalysis {
    pub samples: Vec<String>,
    pub regime_js_matrix: Vec<Vec<f64>>,
    pub basin_overlap_matrix: Vec<Vec<f64>>,
    pub dsp_vectors: BTreeMap<String, DspPairVector>,
    pub sdi_matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct DspTsvRow {
    pub sample_a: String,
    pub sample_b: String,
    pub axis: String,
    pub value: f64,
}

#[derive(Debug, Clone, Copy)]
struct RobustScale {
    median: f64,
    mad: f64,
}

#[derive(Debug, Clone)]
struct BasinCentroid {
    id: String,
    centroid: [f64; 4],
}

#[derive(Debug, Clone)]
struct PairDraft {
    i: usize,
    j: usize,
    js: f64,
    overlap: f64,
    norm: f64,
    dsp: Vec<f64>,
}

pub fn load_sample_systems_input(
    sample: &str,
    sample_root: &Path,
) -> Result<SampleSystemsInput, String> {
    let integration_dir = sample_root.join("kira-organelle").join("integration");
    let metrics_path = integration_dir.join("metrics.tsv");
    let normalized_path = integration_dir.join("metrics_normalized.tsv");

    let metrics_rows = parse_metrics_tsv(&metrics_path)?;
    let axis_raw_values = parse_metrics_normalized_tsv(&normalized_path)?;

    Ok(SampleSystemsInput {
        sample: sample.to_string(),
        metrics_rows,
        axis_raw_values,
    })
}

pub fn compute_cross_sample_analysis(samples: &[SampleSystemsInput]) -> CrossSampleAnalysis {
    let n_samples = samples.len();
    let sample_names = samples.iter().map(|s| s.sample.clone()).collect::<Vec<_>>();
    if n_samples == 0 {
        return CrossSampleAnalysis {
            samples: sample_names,
            regime_js_matrix: Vec::new(),
            basin_overlap_matrix: Vec::new(),
            dsp_vectors: BTreeMap::new(),
            sdi_matrix: Vec::new(),
        };
    }

    let feature_stats = compute_feature_harmonization(samples);
    let axis_stats = compute_axis_harmonization(samples);

    let mut regime_probs = Vec::<Vec<f64>>::with_capacity(n_samples);
    let mut basin_centroids = Vec::<Vec<BasinCentroid>>::with_capacity(n_samples);
    let mut axis_medians = Vec::<[f64; DSP_AXIS_ORDER.len()]>::with_capacity(n_samples);

    for sample in samples {
        regime_probs.push(regime_probabilities(sample));
        basin_centroids.push(sample_basin_centroids(sample, &feature_stats));
        axis_medians.push(sample_axis_medians(sample, &axis_stats));
    }

    let mut regime_js_matrix = vec![vec![0.0; n_samples]; n_samples];
    let mut basin_overlap_matrix = vec![vec![0.0; n_samples]; n_samples];
    let mut sdi_matrix = vec![vec![0.0; n_samples]; n_samples];
    for i in 0..n_samples {
        basin_overlap_matrix[i][i] = 1.0;
    }
    let mut pair_drafts = Vec::<PairDraft>::new();

    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            let js = canonicalize_f64(js_divergence(&regime_probs[i], &regime_probs[j]));
            let overlap = canonicalize_f64(basin_overlap_score(
                &basin_centroids[i],
                &basin_centroids[j],
            ));

            let mut dsp = [0.0f64; DSP_AXIS_ORDER.len()];
            let mut norm_sq = 0.0f64;
            for (axis_idx, _) in DSP_AXIS_ORDER.iter().enumerate() {
                let v = axis_medians[i][axis_idx] - axis_medians[j][axis_idx];
                dsp[axis_idx] = canonicalize_f64(v);
                norm_sq += v * v;
            }
            let norm = canonicalize_f64(norm_sq.sqrt());

            regime_js_matrix[i][j] = js;
            regime_js_matrix[j][i] = js;
            basin_overlap_matrix[i][j] = overlap;
            basin_overlap_matrix[j][i] = overlap;
            pair_drafts.push(PairDraft {
                i,
                j,
                js,
                overlap,
                norm,
                dsp: dsp.to_vec(),
            });
        }
    }

    let mut pair_norms = pair_drafts.iter().map(|p| p.norm).collect::<Vec<_>>();
    let norm_stats = robust_scale(&mut pair_norms);
    let norm_mad = norm_stats.mad;

    let mut dsp_vectors = BTreeMap::<String, DspPairVector>::new();
    for pair in pair_drafts {
        let norm_normalized = canonicalize_f64(pair.norm / (norm_mad + EPS));
        let sdi =
            canonicalize_f64(0.4 * pair.js + 0.3 * (1.0 - pair.overlap) + 0.3 * norm_normalized)
                .max(0.0);
        sdi_matrix[pair.i][pair.j] = sdi;
        sdi_matrix[pair.j][pair.i] = sdi;

        let mut vector = Vec::<DspAxisDelta>::with_capacity(DSP_AXIS_ORDER.len());
        for (axis_idx, axis) in DSP_AXIS_ORDER.iter().enumerate() {
            vector.push(DspAxisDelta {
                axis: (*axis).to_string(),
                value: pair.dsp[axis_idx],
            });
        }
        let sample_a = sample_names[pair.i].clone();
        let sample_b = sample_names[pair.j].clone();
        let key = format!("{sample_a}__vs__{sample_b}");
        dsp_vectors.insert(
            key,
            DspPairVector {
                sample_a,
                sample_b,
                vector,
                norm: pair.norm,
                norm_normalized,
            },
        );
    }

    CrossSampleAnalysis {
        samples: sample_names,
        regime_js_matrix,
        basin_overlap_matrix,
        dsp_vectors,
        sdi_matrix,
    }
}

pub fn render_cross_sample_dsp_tsv(analysis: &CrossSampleAnalysis) -> String {
    let mut out = String::new();
    out.push_str("sample_a\tsample_b\taxis\tvalue\n");
    for pair in analysis.dsp_vectors.values() {
        for axis in &pair.vector {
            out.push_str(&pair.sample_a);
            out.push('\t');
            out.push_str(&pair.sample_b);
            out.push('\t');
            out.push_str(&axis.axis);
            out.push('\t');
            out.push_str(&format!("{:.6}", canonicalize_f64(axis.value)));
            out.push('\n');
        }
        out.push_str(&pair.sample_a);
        out.push('\t');
        out.push_str(&pair.sample_b);
        out.push('\t');
        out.push_str("DSP_NORM");
        out.push('\t');
        out.push_str(&format!("{:.6}", canonicalize_f64(pair.norm)));
        out.push('\n');
        out.push_str(&pair.sample_a);
        out.push('\t');
        out.push_str(&pair.sample_b);
        out.push('\t');
        out.push_str("DSP_NORM_NORMALIZED");
        out.push('\t');
        out.push_str(&format!("{:.6}", canonicalize_f64(pair.norm_normalized)));
        out.push('\n');
    }
    out
}

pub fn cross_sample_dsp_rows(analysis: &CrossSampleAnalysis) -> Vec<DspTsvRow> {
    let mut rows = Vec::<DspTsvRow>::new();
    for pair in analysis.dsp_vectors.values() {
        for axis in &pair.vector {
            rows.push(DspTsvRow {
                sample_a: pair.sample_a.clone(),
                sample_b: pair.sample_b.clone(),
                axis: axis.axis.clone(),
                value: axis.value,
            });
        }
    }
    rows
}

fn parse_metrics_tsv(path: &Path) -> Result<Vec<MetricsRow>, String> {
    let file =
        std::fs::File::open(path).map_err(|e| format!("failed opening {}: {e}", path.display()))?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();
    let Some(header_line) = lines.next() else {
        return Ok(Vec::new());
    };
    let header_line = header_line.map_err(|e| format!("failed reading {}: {e}", path.display()))?;
    let header = header_line.split('\t').collect::<Vec<_>>();
    let idx = header
        .iter()
        .enumerate()
        .map(|(i, c)| ((*c).to_string(), i))
        .collect::<BTreeMap<_, _>>();

    let idx_cell = idx.get("cell_id").copied().unwrap_or(0);
    let idx_stress = idx.get("StressVector").copied();
    let idx_cpi = idx.get("CPI").copied();
    let idx_afs = idx.get("AFS").copied();
    let idx_imsc = idx.get("IMSC").copied();
    let idx_regime = idx.get("RegimeClass").copied();
    let idx_basin = idx.get("BasinId").copied();

    let mut out = Vec::<MetricsRow>::new();
    for line in lines {
        let line = line.map_err(|e| format!("failed reading {}: {e}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let cell_id = fields
            .get(idx_cell)
            .copied()
            .unwrap_or_default()
            .to_string();
        let stress_vector = idx_stress
            .and_then(|i| fields.get(i).copied())
            .and_then(parse_finite_f64)
            .unwrap_or(0.0);
        let cpi = idx_cpi
            .and_then(|i| fields.get(i).copied())
            .and_then(parse_finite_f64)
            .unwrap_or(0.0);
        let afs = idx_afs
            .and_then(|i| fields.get(i).copied())
            .and_then(parse_finite_f64)
            .unwrap_or(0.0);
        let imsc = idx_imsc
            .and_then(|i| fields.get(i).copied())
            .and_then(parse_finite_f64)
            .unwrap_or(0.0);
        let regime_class = idx_regime
            .and_then(|i| fields.get(i).copied())
            .filter(|v| !v.is_empty())
            .unwrap_or("BaselineCompensated")
            .to_string();
        let basin_id = idx_basin
            .and_then(|i| fields.get(i).copied())
            .filter(|v| !v.is_empty())
            .unwrap_or("NA")
            .to_string();

        out.push(MetricsRow {
            cell_id,
            stress_vector,
            cpi,
            afs,
            imsc,
            regime_class,
            basin_id,
        });
    }
    out.sort_by(|a, b| a.cell_id.cmp(&b.cell_id));
    Ok(out)
}

fn parse_metrics_normalized_tsv(path: &Path) -> Result<BTreeMap<String, Vec<f64>>, String> {
    if !path.is_file() {
        return Ok(BTreeMap::new());
    }
    let file =
        std::fs::File::open(path).map_err(|e| format!("failed opening {}: {e}", path.display()))?;
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();
    let Some(header_line) = lines.next() else {
        return Ok(BTreeMap::new());
    };
    let header_line = header_line.map_err(|e| format!("failed reading {}: {e}", path.display()))?;
    let header = header_line.split('\t').collect::<Vec<_>>();
    let idx = header
        .iter()
        .enumerate()
        .map(|(i, c)| ((*c).to_string(), i))
        .collect::<BTreeMap<_, _>>();
    let idx_metric = idx.get("metric").copied().ok_or_else(|| {
        format!(
            "{} missing required 'metric' column in metrics_normalized.tsv",
            path.display()
        )
    })?;
    let idx_raw = idx.get("raw_value").copied().ok_or_else(|| {
        format!(
            "{} missing required 'raw_value' column in metrics_normalized.tsv",
            path.display()
        )
    })?;

    let mut out = BTreeMap::<String, Vec<f64>>::new();
    for axis in DSP_AXIS_ORDER {
        out.insert(axis.to_string(), Vec::new());
    }

    for line in lines {
        let line = line.map_err(|e| format!("failed reading {}: {e}", path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let metric = fields.get(idx_metric).copied().unwrap_or_default();
        if !matches!(
            metric,
            "RSS" | "SII" | "OSL" | "TSM" | "PCP" | "ASM" | "MSM"
        ) {
            continue;
        }
        let Some(raw) = fields.get(idx_raw).copied().and_then(parse_finite_f64) else {
            continue;
        };
        out.entry(metric.to_string()).or_default().push(raw);
    }
    Ok(out)
}

fn compute_feature_harmonization(samples: &[SampleSystemsInput]) -> [RobustScale; 4] {
    let mut pooled = [
        Vec::<f64>::new(),
        Vec::<f64>::new(),
        Vec::<f64>::new(),
        Vec::<f64>::new(),
    ];
    for sample in samples {
        for row in &sample.metrics_rows {
            pooled[0].push(finite_or_zero(row.stress_vector));
            pooled[1].push(finite_or_zero(row.cpi));
            pooled[2].push(finite_or_zero(row.afs));
            pooled[3].push(finite_or_zero(row.imsc));
        }
    }
    [
        robust_scale(&mut pooled[0]),
        robust_scale(&mut pooled[1]),
        robust_scale(&mut pooled[2]),
        robust_scale(&mut pooled[3]),
    ]
}

fn compute_axis_harmonization(
    samples: &[SampleSystemsInput],
) -> [RobustScale; DSP_AXIS_ORDER.len()] {
    let mut pooled = [
        Vec::<f64>::new(),
        Vec::<f64>::new(),
        Vec::<f64>::new(),
        Vec::<f64>::new(),
        Vec::<f64>::new(),
        Vec::<f64>::new(),
        Vec::<f64>::new(),
    ];
    for sample in samples {
        for (axis_idx, axis) in DSP_AXIS_ORDER.iter().enumerate() {
            if let Some(vals) = sample.axis_raw_values.get(*axis) {
                pooled[axis_idx].extend(vals.iter().copied().filter(|v| v.is_finite()));
            }
        }
    }
    [
        robust_scale(&mut pooled[0]),
        robust_scale(&mut pooled[1]),
        robust_scale(&mut pooled[2]),
        robust_scale(&mut pooled[3]),
        robust_scale(&mut pooled[4]),
        robust_scale(&mut pooled[5]),
        robust_scale(&mut pooled[6]),
    ]
}

fn regime_probabilities(sample: &SampleSystemsInput) -> Vec<f64> {
    let mut counts = BTreeMap::<String, u64>::new();
    for row in &sample.metrics_rows {
        *counts.entry(row.regime_class.clone()).or_insert(0) += 1;
    }
    let total = sample.metrics_rows.len() as f64;
    if total == 0.0 {
        return vec![0.0; REGIME_ORDER.len()];
    }
    REGIME_ORDER
        .iter()
        .map(|k| canonicalize_f64(counts.get(*k).copied().unwrap_or(0) as f64 / total))
        .collect()
}

fn sample_basin_centroids(
    sample: &SampleSystemsInput,
    feature_stats: &[RobustScale; 4],
) -> Vec<BasinCentroid> {
    let mut accum = BTreeMap::<String, ([f64; 4], usize)>::new();
    for row in &sample.metrics_rows {
        if row.basin_id.is_empty() || row.basin_id == "NA" {
            continue;
        }
        let f = [
            harmonized_z(finite_or_zero(row.stress_vector), feature_stats[0]),
            harmonized_z(finite_or_zero(row.cpi), feature_stats[1]),
            harmonized_z(finite_or_zero(row.afs), feature_stats[2]),
            harmonized_z(finite_or_zero(row.imsc), feature_stats[3]),
        ];
        let entry = accum
            .entry(row.basin_id.clone())
            .or_insert(([0.0, 0.0, 0.0, 0.0], 0));
        for (idx, value) in f.iter().enumerate() {
            entry.0[idx] += *value;
        }
        entry.1 += 1;
    }

    let mut out = Vec::<BasinCentroid>::new();
    for (id, (sum, n)) in accum {
        if n == 0 {
            continue;
        }
        out.push(BasinCentroid {
            id,
            centroid: [
                canonicalize_f64(sum[0] / n as f64),
                canonicalize_f64(sum[1] / n as f64),
                canonicalize_f64(sum[2] / n as f64),
                canonicalize_f64(sum[3] / n as f64),
            ],
        });
    }
    out.sort_by(|a, b| a.id.cmp(&b.id));
    out
}

fn sample_axis_medians(
    sample: &SampleSystemsInput,
    axis_stats: &[RobustScale; DSP_AXIS_ORDER.len()],
) -> [f64; DSP_AXIS_ORDER.len()] {
    let mut out = [0.0f64; DSP_AXIS_ORDER.len()];
    for (axis_idx, axis) in DSP_AXIS_ORDER.iter().enumerate() {
        let Some(values) = sample.axis_raw_values.get(*axis) else {
            continue;
        };
        if values.is_empty() {
            continue;
        }
        let mut z = Vec::<f64>::with_capacity(values.len());
        for value in values {
            if !value.is_finite() {
                continue;
            }
            z.push(harmonized_z(*value, axis_stats[axis_idx]));
        }
        if z.is_empty() {
            continue;
        }
        out[axis_idx] = canonicalize_f64(median_in_place(&mut z).unwrap_or(0.0));
    }
    out
}

fn basin_overlap_score(a: &[BasinCentroid], b: &[BasinCentroid]) -> f64 {
    if a.is_empty() && b.is_empty() {
        return 1.0;
    }
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }

    let mut edges = Vec::<(f64, usize, usize)>::new();
    let mut dist_values = Vec::<f64>::new();
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            let d = euclidean4(ai.centroid, bj.centroid);
            edges.push((d, i, j));
            dist_values.push(d);
        }
    }
    let threshold = median_in_place(&mut dist_values).unwrap_or(0.0);
    edges.sort_by(|x, y| {
        x.0.partial_cmp(&y.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a[x.1].id.cmp(&a[y.1].id))
            .then(b[x.2].id.cmp(&b[y.2].id))
    });

    let mut used_a = vec![false; a.len()];
    let mut used_b = vec![false; b.len()];
    let mut aligned = 0usize;
    for (d, i, j) in edges {
        if d >= threshold || used_a[i] || used_b[j] {
            continue;
        }
        used_a[i] = true;
        used_b[j] = true;
        aligned += 1;
    }
    let denom = a.len().min(b.len());
    if denom == 0 {
        0.0
    } else {
        canonicalize_f64(aligned as f64 / denom as f64)
    }
}

fn js_divergence(p: &[f64], q: &[f64]) -> f64 {
    let n = p.len().min(q.len());
    if n == 0 {
        return 0.0;
    }
    let mut m = vec![0.0f64; n];
    for i in 0..n {
        m[i] = 0.5 * (p[i] + q[i]);
    }
    canonicalize_f64(0.5 * kl_divergence(p, &m) + 0.5 * kl_divergence(q, &m))
}

fn kl_divergence(p: &[f64], q: &[f64]) -> f64 {
    let n = p.len().min(q.len());
    let mut out = 0.0f64;
    for i in 0..n {
        let pi = p[i];
        let qi = q[i];
        if pi <= 0.0 || qi <= 0.0 {
            continue;
        }
        out += pi * (pi / qi).ln();
    }
    canonicalize_f64(out)
}

fn robust_scale(values: &mut [f64]) -> RobustScale {
    if values.is_empty() {
        return RobustScale {
            median: 0.0,
            mad: 0.0,
        };
    }
    let median = median_in_place(values).unwrap_or(0.0);
    let mut dev = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .map(|v| (v - median).abs())
        .collect::<Vec<_>>();
    let mad = median_in_place(&mut dev).unwrap_or(0.0);
    RobustScale {
        median: canonicalize_f64(median),
        mad: canonicalize_f64(mad),
    }
}

fn harmonized_z(value: f64, stats: RobustScale) -> f64 {
    if !value.is_finite() {
        return 0.0;
    }
    if stats.mad <= 0.0 {
        return 0.0;
    }
    canonicalize_f64((value - stats.median) / (MAD_SCALE * stats.mad + EPS))
}

fn parse_finite_f64(raw: &str) -> Option<f64> {
    let parsed = raw.parse::<f64>().ok()?;
    if parsed.is_finite() {
        Some(parsed)
    } else {
        None
    }
}

fn finite_or_zero(v: f64) -> f64 {
    if v.is_finite() { v } else { 0.0 }
}

fn euclidean4(a: [f64; 4], b: [f64; 4]) -> f64 {
    let d0 = a[0] - b[0];
    let d1 = a[1] - b[1];
    let d2 = a[2] - b[2];
    let d3 = a[3] - b[3];
    canonicalize_f64((d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn js_zero_for_identical_distribution() {
        let p = vec![0.2, 0.3, 0.5];
        let q = vec![0.2, 0.3, 0.5];
        let js = js_divergence(&p, &q);
        assert_eq!(js, 0.0);
    }

    #[test]
    fn basin_overlap_is_one_if_both_empty() {
        let score = basin_overlap_score(&[], &[]);
        assert_eq!(score, 1.0);
    }
}
