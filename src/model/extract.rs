use std::collections::BTreeMap;

use serde_json::Value;

use crate::contracts::types::Issue;
use crate::model::organelle::OrganelleId;
use crate::warn_missing_issue;

use crate::state::{AxisSummary, OrganelleState, RegimeDistribution};

pub fn extract_organelle_state(
    tool: &str,
    summary: &Value,
    issues: &mut Vec<Issue>,
) -> Option<OrganelleState> {
    let organelle = match OrganelleId::from_tool_name(tool) {
        Some(v) => v,
        None => {
            issues.push(warn_missing_issue(
                Some(tool),
                "UNKNOWN_TOOL_MAPPING",
                "tool is discovered but has no organelle mapping".to_string(),
                None,
            ));
            return None;
        }
    };

    let axes = extract_axes(tool, summary, issues);
    let regimes = extract_regimes(tool, summary, issues);
    let qc = extract_qc(tool, summary, issues);

    Some(OrganelleState {
        organelle,
        tool: tool.to_string(),
        axes,
        regimes,
        qc,
    })
}

fn extract_axes(tool: &str, summary: &Value, issues: &mut Vec<Issue>) -> Vec<AxisSummary> {
    let Some(distributions) = summary.get("distributions") else {
        // Backward-compatible mito summary shape: { "axes_median": { ... } }
        if let Some(axes_median) = summary.get("axes_median").and_then(Value::as_object) {
            let mut axes = Vec::new();
            for (axis_name, median_val) in axes_median {
                if let Some(median) = median_val.as_f64() {
                    axes.push(AxisSummary {
                        name: axis_name.clone(),
                        median,
                        p90: None,
                        p99: None,
                    });
                }
            }
            axes.sort_by(|a, b| a.name.cmp(&b.name));
            return axes;
        }
        issues.push(warn_missing_issue(
            Some(tool),
            "MISSING_DISTRIBUTIONS",
            "summary.json has no distributions object".to_string(),
            None,
        ));
        return Vec::new();
    };

    let Some(obj) = distributions.as_object() else {
        issues.push(warn_missing_issue(
            Some(tool),
            "MALFORMED_DISTRIBUTIONS",
            "distributions exists but is not an object".to_string(),
            None,
        ));
        return Vec::new();
    };

    let mut axes = Vec::new();
    for (axis_name, axis_val) in obj {
        if let Some(axis_obj) = axis_val.as_object() {
            let Some(median) = axis_obj.get("median").and_then(Value::as_f64) else {
                issues.push(warn_missing_issue(
                    Some(tool),
                    "MISSING_AXIS_MEDIAN",
                    format!("distribution axis '{}' has no numeric median", axis_name),
                    None,
                ));
                continue;
            };

            let p90 = axis_obj.get("p90").and_then(Value::as_f64);
            let p99 = axis_obj.get("p99").and_then(Value::as_f64);

            axes.push(AxisSummary {
                name: axis_name.clone(),
                median,
                p90,
                p99,
            });
            continue;
        }

        // Newer summary shape used by nuclearqc:
        // distributions: { axes: [{name, median, p90, p99}, ...], composites: [...] }
        if let Some(axis_arr) = axis_val.as_array() {
            for axis_item in axis_arr {
                let Some(axis_obj) = axis_item.as_object() else {
                    continue;
                };
                let Some(name) = axis_obj.get("name").and_then(Value::as_str) else {
                    continue;
                };
                let Some(median) = axis_obj.get("median").and_then(Value::as_f64) else {
                    continue;
                };
                let p90 = axis_obj.get("p90").and_then(Value::as_f64);
                let p99 = axis_obj.get("p99").and_then(Value::as_f64);
                axes.push(AxisSummary {
                    name: name.to_string(),
                    median,
                    p90,
                    p99,
                });
            }
            continue;
        }

        issues.push(warn_missing_issue(
            Some(tool),
            "MALFORMED_AXIS",
            format!("distribution axis '{}' is not an object", axis_name),
            None,
        ));
    }

    axes.sort_by(|a, b| a.name.cmp(&b.name));
    axes
}

fn extract_regimes(
    tool: &str,
    summary: &Value,
    issues: &mut Vec<Issue>,
) -> Option<RegimeDistribution> {
    let regimes = summary.get("regimes")?;

    let Some(regimes_obj) = regimes.as_object() else {
        issues.push(warn_missing_issue(
            Some(tool),
            "MALFORMED_REGIMES",
            "regimes exists but is not an object".to_string(),
            None,
        ));
        return None;
    };

    // Preferred shape: regimes.counts + regimes.fractions.
    let mut counts = BTreeMap::new();
    let mut fractions = BTreeMap::new();

    if let Some(counts_obj) = regimes_obj.get("counts").and_then(Value::as_object) {
        for (k, v) in counts_obj {
            if let Some(n) = as_count(v) {
                counts.insert(k.clone(), n);
            }
        }
    } else if let Some(regime_counts_obj) = summary.get("regime_counts").and_then(Value::as_object)
    {
        // Compatibility shape used by nuclearqc/autolys.
        for (k, v) in regime_counts_obj {
            if let Some(n) = as_count(v) {
                counts.insert(k.clone(), n);
            }
        }
    }

    if let Some(fractions_obj) = regimes_obj.get("fractions").and_then(Value::as_object) {
        for (k, v) in fractions_obj {
            if let Some(n) = v.as_f64() {
                fractions.insert(k.clone(), n);
            }
        }
    } else {
        // Compatibility shape used by nuclearqc/autolys: regimes itself is fraction map.
        for (k, v) in regimes_obj {
            if let Some(n) = v.as_f64() {
                fractions.insert(k.clone(), n);
            }
        }
    }

    if counts.is_empty() {
        issues.push(warn_missing_issue(
            Some(tool),
            "MISSING_REGIME_COUNTS",
            "regimes.counts missing or malformed".to_string(),
            None,
        ));
        return None;
    }
    if fractions.is_empty() {
        issues.push(warn_missing_issue(
            Some(tool),
            "MISSING_REGIME_FRACTIONS",
            "regimes.fractions missing or malformed".to_string(),
            None,
        ));
        return None;
    }

    Some(RegimeDistribution { counts, fractions })
}

fn as_count(v: &Value) -> Option<u64> {
    if let Some(n) = v.as_u64() {
        return Some(n);
    }
    let f = v.as_f64()?;
    if !f.is_finite() || f < 0.0 {
        return None;
    }
    Some(f.round() as u64)
}

fn extract_qc(tool: &str, summary: &Value, issues: &mut Vec<Issue>) -> BTreeMap<String, f64> {
    let Some(qc_value) = summary.get("qc") else {
        return BTreeMap::new();
    };
    let Some(qc_obj) = qc_value.as_object() else {
        issues.push(warn_missing_issue(
            Some(tool),
            "MALFORMED_QC",
            "qc exists but is not an object".to_string(),
            None,
        ));
        return BTreeMap::new();
    };

    let mut qc = BTreeMap::new();
    for (k, v) in qc_obj {
        if !k.ends_with("_fraction") {
            continue;
        }
        match v.as_f64() {
            Some(n) => {
                qc.insert(k.clone(), n);
            }
            None => {
                issues.push(warn_missing_issue(
                    Some(tool),
                    "MALFORMED_QC_VALUE",
                    format!("qc field '{}' is not numeric", k),
                    None,
                ));
            }
        }
    }
    qc
}
