use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::cells::types::CellState;
use crate::cells::types::CellsState;
use crate::model::organelle::OrganelleId;

const ADAPTIVE_THRESHOLD: f64 = 0.33;
const RESISTANT_THRESHOLD: f64 = 0.66;
const WEIGHT_EPSILON: f64 = 1e-6;

const METRIC_MITO: &str = "mitochondrial_stress_adaptation_score";
const METRIC_TRANSLATION: &str = "translation_commitment_score";
const METRIC_SPLICE: &str = "splice_irreversibility_index";

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FiiWeights {
    pub mitochondrial: f64,
    pub translation: f64,
    pub splice: f64,
}

impl Default for FiiWeights {
    fn default() -> Self {
        Self {
            mitochondrial: 1.0 / 3.0,
            translation: 1.0 / 3.0,
            splice: 1.0 / 3.0,
        }
    }
}

impl FiiWeights {
    pub fn parse(raw: &str) -> Result<Self, String> {
        let mut values = BTreeMap::new();
        for pair in raw.split(',') {
            let mut it = pair.splitn(2, ':');
            let key = it
                .next()
                .map(str::trim)
                .ok_or_else(|| "invalid --fii-weights format".to_string())?;
            let value = it
                .next()
                .map(str::trim)
                .ok_or_else(|| "invalid --fii-weights format".to_string())?;
            let parsed = value
                .parse::<f64>()
                .map_err(|_| format!("invalid numeric value for '{key}': {value}"))?;
            values.insert(key.to_ascii_lowercase(), parsed);
        }

        let mitochondrial = *values.get("mito").unwrap_or(&f64::NAN);
        let translation = *values.get("translation").unwrap_or(&f64::NAN);
        let splice = *values.get("splice").unwrap_or(&f64::NAN);
        let weights = Self {
            mitochondrial,
            translation,
            splice,
        };
        weights.validate()?;
        Ok(weights)
    }

    pub fn validate(&self) -> Result<(), String> {
        for (name, value) in [
            ("mito", self.mitochondrial),
            ("translation", self.translation),
            ("splice", self.splice),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(format!("invalid --fii-weights value for {name}: {value}"));
            }
        }
        let sum = self.mitochondrial + self.translation + self.splice;
        if (sum - 1.0).abs() > WEIGHT_EPSILON {
            return Err(format!(
                "invalid --fii-weights sum: expected 1.0, got {sum:.8}"
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FiiRegime {
    Adaptive,
    Transition,
    Resistant,
}

impl FiiRegime {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Adaptive => "Adaptive",
            Self::Transition => "Transition",
            Self::Resistant => "Resistant",
        }
    }
}

#[derive(Debug, Clone)]
pub struct FiiCellRow {
    pub cell_id: String,
    pub mitochondrial_stress_adaptation_score: Option<f64>,
    pub translation_commitment_score: Option<f64>,
    pub splice_irreversibility_index: Option<f64>,
    pub functional_irreversibility_index: f64,
    pub fii_regime: FiiRegime,
    pub fii_low_confidence: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FiiDistributionSummary {
    pub mean: f64,
    pub median: f64,
    pub p25: f64,
    pub p75: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FiiSummary {
    pub weights: FiiWeights,
    pub distribution: FiiDistributionSummary,
    pub regime_fractions: BTreeMap<String, f64>,
    pub low_confidence_fraction: f64,
}

#[derive(Debug, Clone)]
pub struct FiiComputation {
    pub rows: Vec<FiiCellRow>,
    pub summary: FiiSummary,
}

pub fn compute_fii(cells: &CellsState, weights: FiiWeights) -> Result<FiiComputation, String> {
    weights.validate()?;
    let mut rows = Vec::with_capacity(cells.cells.len());
    let mut values = Vec::with_capacity(cells.cells.len());
    let mut low_conf_count = 0usize;
    let mut regime_counts = BTreeMap::<String, usize>::new();

    for cell in &cells.cells {
        let mito = read_cell_metric(cell, OrganelleId::Mitochondria, METRIC_MITO);
        let translation = read_cell_metric(cell, OrganelleId::Ribosome, METRIC_TRANSLATION);
        let splice = read_cell_metric(cell, OrganelleId::Spliceosome, METRIC_SPLICE);

        let mito = mito.map(clamp01);
        let translation = translation.map(clamp01);
        let splice = splice.map(clamp01);

        let mut weighted_sum = 0.0;
        let mut weight_sum = 0.0;
        if let Some(v) = mito {
            weighted_sum += weights.mitochondrial * v;
            weight_sum += weights.mitochondrial;
        }
        if let Some(v) = translation {
            weighted_sum += weights.translation * v;
            weight_sum += weights.translation;
        }
        if let Some(v) = splice {
            weighted_sum += weights.splice * v;
            weight_sum += weights.splice;
        }

        // Missing inputs are allowed but marked low-confidence; normalize over available inputs.
        let fii = if weight_sum > 0.0 {
            clamp01(weighted_sum / weight_sum)
        } else {
            0.0
        };
        let low_conf = mito.is_none() || translation.is_none() || splice.is_none();
        if low_conf {
            low_conf_count += 1;
        }

        let regime = classify_regime(fii);
        *regime_counts
            .entry(regime.as_str().to_ascii_lowercase())
            .or_insert(0) += 1;

        rows.push(FiiCellRow {
            cell_id: cell.id.clone(),
            mitochondrial_stress_adaptation_score: mito,
            translation_commitment_score: translation,
            splice_irreversibility_index: splice,
            functional_irreversibility_index: fii,
            fii_regime: regime,
            fii_low_confidence: low_conf,
        });
        values.push(fii);
    }

    let n = rows.len() as f64;
    let distribution = FiiDistributionSummary {
        mean: if n > 0.0 {
            values.iter().sum::<f64>() / n
        } else {
            0.0
        },
        median: percentile(&values, 0.50),
        p25: percentile(&values, 0.25),
        p75: percentile(&values, 0.75),
    };

    let mut regime_fractions = BTreeMap::new();
    for key in ["adaptive", "transition", "resistant"] {
        let count = *regime_counts.get(key).unwrap_or(&0) as f64;
        regime_fractions.insert(key.to_string(), if n > 0.0 { count / n } else { 0.0 });
    }

    Ok(FiiComputation {
        rows,
        summary: FiiSummary {
            weights,
            distribution,
            regime_fractions,
            low_confidence_fraction: if n > 0.0 {
                low_conf_count as f64 / n
            } else {
                0.0
            },
        },
    })
}

pub fn render_fii_tsv(rows: &[FiiCellRow]) -> String {
    let mut out = String::from(
        "cell_id\tmitochondrial_stress_adaptation_score\ttranslation_commitment_score\tsplice_irreversibility_index\tfunctional_irreversibility_index\tfii_regime\tfii_low_confidence\n",
    );
    for row in rows {
        out.push_str(&row.cell_id);
        out.push('\t');
        write_opt_float(&mut out, row.mitochondrial_stress_adaptation_score);
        out.push('\t');
        write_opt_float(&mut out, row.translation_commitment_score);
        out.push('\t');
        write_opt_float(&mut out, row.splice_irreversibility_index);
        out.push('\t');
        out.push_str(&format!("{:.6}", row.functional_irreversibility_index));
        out.push('\t');
        out.push_str(row.fii_regime.as_str());
        out.push('\t');
        out.push_str(if row.fii_low_confidence {
            "true"
        } else {
            "false"
        });
        out.push('\n');
    }
    out
}

pub fn classify_regime(fii: f64) -> FiiRegime {
    if fii < ADAPTIVE_THRESHOLD {
        FiiRegime::Adaptive
    } else if fii < RESISTANT_THRESHOLD {
        FiiRegime::Transition
    } else {
        FiiRegime::Resistant
    }
}

fn read_cell_metric(cell: &CellState, organelle: OrganelleId, metric_name: &str) -> Option<f64> {
    let org = cell.per_organelle.get(&organelle)?;
    org.axes.get(metric_name).copied()
}

fn clamp01(v: f64) -> f64 {
    if !v.is_finite() {
        return 0.0;
    }
    v.clamp(0.0, 1.0)
}

fn percentile(values: &[f64], q: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let rank = ((sorted.len() - 1) as f64 * q).round() as usize;
    sorted[rank]
}

fn write_opt_float(out: &mut String, value: Option<f64>) {
    if let Some(v) = value {
        out.push_str(&format!("{:.6}", v));
    }
}

#[cfg(test)]
mod tests {
    use super::{FiiRegime, FiiWeights, classify_regime};

    #[test]
    fn parse_weights() {
        let w = FiiWeights::parse("mito:0.34,translation:0.33,splice:0.33").expect("parse");
        assert!((w.mitochondrial - 0.34).abs() < 1e-9);
        assert!((w.translation - 0.33).abs() < 1e-9);
        assert!((w.splice - 0.33).abs() < 1e-9);
    }

    #[test]
    fn reject_bad_weight_sum() {
        let err = FiiWeights::parse("mito:0.5,translation:0.5,splice:0.5").expect_err("err");
        assert!(err.contains("sum"));
    }

    #[test]
    fn regime_boundaries() {
        assert_eq!(classify_regime(0.0), FiiRegime::Adaptive);
        assert_eq!(classify_regime(0.329999), FiiRegime::Adaptive);
        assert_eq!(classify_regime(0.33), FiiRegime::Transition);
        assert_eq!(classify_regime(0.659999), FiiRegime::Transition);
        assert_eq!(classify_regime(0.66), FiiRegime::Resistant);
        assert_eq!(classify_regime(1.0), FiiRegime::Resistant);
    }
}
