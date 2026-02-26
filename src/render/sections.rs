use std::collections::BTreeMap;

use crate::cells::types::CELLS_SCHEMA_V1;
use crate::compare::types::ComparisonState;
use crate::render::format::{escape_html, fmt_delta, fmt_f64, fmt_opt_f64, organelle_label};
use crate::state::AggregatedState;

pub fn render_header(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
    out: &mut String,
) {
    let mode = if comparison.is_some() {
        "comparison"
    } else {
        "single"
    };
    out.push_str("<section><h2>Header</h2>");
    out.push_str("<p><strong>kira-organelle report</strong></p>");
    out.push_str("<table><tbody>");
    row(out, "Timestamp", &escape_html(&state.timestamp));
    row(out, "Mode", mode);
    row(out, "Input A", &escape_html(&state.inputs.a));
    row(
        out,
        "Input B",
        &escape_html(state.inputs.b.as_deref().unwrap_or("-")),
    );
    out.push_str("</tbody></table></section>");
}

pub fn render_global_summary(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
    out: &mut String,
) {
    out.push_str("<section><h2>Global Summary</h2>");

    match comparison {
        None => {
            out.push_str("<table><thead><tr><th>organelle</th><th>axes</th><th>dominant_regime</th><th>key_qc_flags</th></tr></thead><tbody>");
            for org in &state.organelle_states {
                let org_name = serde_json::to_value(org.organelle)
                    .ok()
                    .and_then(|v| v.as_str().map(ToOwned::to_owned))
                    .unwrap_or_else(|| "unknown".to_string());

                let dominant = org
                    .regimes
                    .as_ref()
                    .and_then(|r| {
                        r.fractions
                            .iter()
                            .max_by(|a, b| a.1.total_cmp(b.1))
                            .map(|(k, _)| k.as_str())
                    })
                    .unwrap_or("-");

                let flags = org
                    .qc
                    .iter()
                    .filter(|(_, v)| **v > 0.0)
                    .map(|(k, v)| format!("{}={}", k, fmt_f64(*v)))
                    .collect::<Vec<_>>();
                let flags_text = if flags.is_empty() {
                    "-".to_string()
                } else {
                    flags.join(", ")
                };

                out.push_str("<tr>");
                cell(out, organelle_label(&org_name));
                cell(out, &org.axes.len().to_string());
                cell(out, &escape_html(dominant));
                cell(out, &escape_html(&flags_text));
                out.push_str("</tr>");
            }
            out.push_str("</tbody></table>");
        }
        Some(cmp) => {
            out.push_str("<table><thead><tr><th>organelle</th><th>changed_axes</th><th>largest_abs_delta_axis</th><th>direction</th></tr></thead><tbody>");
            for org in &cmp.organelle_deltas {
                let org_name = serde_json::to_value(org.organelle)
                    .ok()
                    .and_then(|v| v.as_str().map(ToOwned::to_owned))
                    .unwrap_or_else(|| "unknown".to_string());

                let mut largest_axis = "-".to_string();
                let mut largest_abs = f64::NEG_INFINITY;
                let mut sum_delta = 0.0;
                for axis in &org.axes {
                    let abs = axis.delta.median.abs();
                    if abs > largest_abs {
                        largest_abs = abs;
                        largest_axis = axis.name.clone();
                    }
                    sum_delta += axis.delta.median;
                }
                let direction = if org.axes.is_empty() {
                    "flat"
                } else if sum_delta > 0.0 {
                    "net up"
                } else if sum_delta < 0.0 {
                    "net down"
                } else {
                    "flat"
                };

                out.push_str("<tr>");
                cell(out, organelle_label(&org_name));
                cell(out, &org.axes.len().to_string());
                cell(out, &escape_html(&largest_axis));
                cell(out, direction);
                out.push_str("</tr>");
            }
            out.push_str("</tbody></table>");
        }
    }

    out.push_str("</section>");
}

pub fn render_organelle_sections(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
    out: &mut String,
) {
    out.push_str("<section><h2>Organelle Sections</h2>");

    match comparison {
        None => {
            for org in &state.organelle_states {
                let org_name = serde_json::to_value(org.organelle)
                    .ok()
                    .and_then(|v| v.as_str().map(ToOwned::to_owned))
                    .unwrap_or_else(|| "unknown".to_string());
                out.push_str("<h3>");
                out.push_str(organelle_label(&org_name));
                out.push_str("</h3>");

                out.push_str("<table><thead><tr><th>axis</th><th>median</th><th>p90</th><th>p99</th></tr></thead><tbody>");
                for axis in &org.axes {
                    out.push_str("<tr>");
                    cell(out, &escape_html(&axis.name));
                    cell(out, &fmt_f64(axis.median));
                    cell(out, &fmt_opt_f64(axis.p90));
                    cell(out, &fmt_opt_f64(axis.p99));
                    out.push_str("</tr>");
                }
                out.push_str("</tbody></table>");

                if let Some(reg) = &org.regimes {
                    out.push_str(
                        "<table><thead><tr><th>regime</th><th>fraction</th></tr></thead><tbody>",
                    );
                    for (k, v) in &reg.fractions {
                        out.push_str("<tr>");
                        cell(out, &escape_html(k));
                        cell(out, &fmt_f64(*v));
                        out.push_str("</tr>");
                    }
                    out.push_str("</tbody></table>");
                }
            }
        }
        Some(cmp) => {
            for org in &cmp.organelle_deltas {
                let org_name = serde_json::to_value(org.organelle)
                    .ok()
                    .and_then(|v| v.as_str().map(ToOwned::to_owned))
                    .unwrap_or_else(|| "unknown".to_string());
                out.push_str("<h3>");
                out.push_str(organelle_label(&org_name));
                out.push_str("</h3>");

                out.push_str("<table><thead><tr><th>axis</th><th>median_a</th><th>median_b</th><th>delta</th></tr></thead><tbody>");
                for axis in &org.axes {
                    out.push_str("<tr>");
                    cell(out, &escape_html(&axis.name));
                    cell(out, &fmt_f64(axis.a.median));
                    cell(out, &fmt_f64(axis.b.median));
                    cell(out, &fmt_delta(axis.delta.median));
                    out.push_str("</tr>");
                }
                out.push_str("</tbody></table>");

                if let Some(reg) = &org.regimes {
                    out.push_str("<table><thead><tr><th>regime</th><th>fraction_a</th><th>fraction_b</th><th>delta</th></tr></thead><tbody>");
                    for regime in regime_union_keys(&reg.a.fractions, &reg.b.fractions) {
                        let a = reg.a.fractions.get(&regime).copied().unwrap_or(0.0);
                        let b = reg.b.fractions.get(&regime).copied().unwrap_or(0.0);
                        let d = reg.delta.get(&regime).copied().unwrap_or(0.0);
                        out.push_str("<tr>");
                        cell(out, &escape_html(&regime));
                        cell(out, &fmt_f64(a));
                        cell(out, &fmt_f64(b));
                        cell(out, &fmt_delta(d));
                        out.push_str("</tr>");
                    }
                    out.push_str("</tbody></table>");
                }
            }
        }
    }

    out.push_str("</section>");
}

pub fn render_comparison_sections(comparison: Option<&ComparisonState>, out: &mut String) {
    if let Some(cmp) = comparison {
        out.push_str("<section><h2>Comparison Sections</h2>");
        out.push_str("<h3>Cell-level Delta Summary</h3>");
        out.push_str("<table><thead><tr><th>organelle</th><th>axis</th><th>median_delta</th><th>p90_abs_delta</th></tr></thead><tbody>");

        let mut axes = cmp.cell_deltas.axes.clone();
        axes.sort_by(|a, b| {
            b.median_delta
                .abs()
                .total_cmp(&a.median_delta.abs())
                .then_with(|| a.axis.cmp(&b.axis))
        });

        for row in &axes {
            let org_name = serde_json::to_value(row.organelle)
                .ok()
                .and_then(|v| v.as_str().map(ToOwned::to_owned))
                .unwrap_or_else(|| "unknown".to_string());
            out.push_str("<tr>");
            cell(out, organelle_label(&org_name));
            cell(out, &escape_html(&row.axis));
            cell(out, &fmt_delta(row.median_delta));
            cell(out, &fmt_f64(row.p90_abs_delta));
            out.push_str("</tr>");
        }
        out.push_str("</tbody></table></section>");
    }
}

pub fn render_issues(issues: &[crate::contracts::types::Issue], out: &mut String) {
    out.push_str("<section><h2>Issues &amp; Warnings</h2>");

    if issues.is_empty() {
        out.push_str("<p>No issues detected</p>");
        out.push_str("</section>");
        return;
    }

    out.push_str("<table><thead><tr><th>severity</th><th>tool</th><th>code</th><th>message</th></tr></thead><tbody>");
    for issue in issues {
        out.push_str("<tr>");
        cell(
            out,
            match issue.severity {
                crate::contracts::types::Severity::Warn => "warn",
                crate::contracts::types::Severity::Error => "error",
            },
        );
        cell(out, &escape_html(issue.tool.as_deref().unwrap_or("-")));
        cell(out, &escape_html(&issue.code));
        cell(out, &escape_html(&issue.message));
        out.push_str("</tr>");
    }
    out.push_str("</tbody></table></section>");
}

pub fn render_footer(
    version: &str,
    has_comparison: bool,
    comparison_schema: Option<&str>,
    out: &mut String,
) {
    out.push_str("<section><h2>Footer / Metadata</h2>");
    out.push_str("<table><tbody>");
    row(out, "kira-organelle version", version);
    row(out, "state schema", crate::state::STATE_SCHEMA_V1);
    row(out, "cells schema", CELLS_SCHEMA_V1);
    row(
        out,
        "comparison schema",
        if has_comparison {
            comparison_schema.unwrap_or("-")
        } else {
            "-"
        },
    );
    row(out, "note", "Deterministic static report build");
    out.push_str("</tbody></table></section>");
}

fn row(out: &mut String, key: &str, value: &str) {
    out.push_str("<tr>");
    out.push_str("<th>");
    out.push_str(&escape_html(key));
    out.push_str("</th>");
    out.push_str("<td>");
    out.push_str(value);
    out.push_str("</td>");
    out.push_str("</tr>");
}

fn cell(out: &mut String, value: &str) {
    out.push_str("<td>");
    out.push_str(value);
    out.push_str("</td>");
}

fn regime_union_keys(a: &BTreeMap<String, f64>, b: &BTreeMap<String, f64>) -> Vec<String> {
    let mut keys = a.keys().cloned().collect::<Vec<_>>();
    for key in b.keys() {
        if !a.contains_key(key) {
            keys.push(key.clone());
        }
    }
    keys.sort();
    keys
}
