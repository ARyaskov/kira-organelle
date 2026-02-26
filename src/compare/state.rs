use std::collections::{BTreeMap, BTreeSet};

use crate::contracts::types::Issue;
use crate::state::{AggregatedState, AxisSummary, OrganelleState};
use crate::warn_missing_issue;

use super::types::{AxisDelta, AxisPoint, OrganelleDelta, RegimeDelta, RegimeSide};

pub fn compute_organelle_deltas(
    a: &AggregatedState,
    b: &AggregatedState,
    issues: &mut Vec<Issue>,
) -> Vec<OrganelleDelta> {
    let a_map = a
        .organelle_states
        .iter()
        .map(|s| (s.organelle, s))
        .collect::<BTreeMap<_, _>>();
    let b_map = b
        .organelle_states
        .iter()
        .map(|s| (s.organelle, s))
        .collect::<BTreeMap<_, _>>();

    let mut organelles = BTreeSet::new();
    for key in a_map.keys() {
        organelles.insert(*key);
    }
    for key in b_map.keys() {
        organelles.insert(*key);
    }

    let mut out = Vec::new();
    for organelle in organelles {
        let a_state = a_map.get(&organelle).copied();
        let b_state = b_map.get(&organelle).copied();

        if a_state.is_none() || b_state.is_none() {
            issues.push(warn_missing_issue(
                None,
                "MISSING_ORGANELLE_SIDE",
                format!("organelle {:?} is missing in one side", organelle),
                None,
            ));
        }

        let axes = match (a_state, b_state) {
            (Some(a_s), Some(b_s)) => compute_axis_deltas(a_s, b_s, issues),
            _ => Vec::new(),
        };

        let regimes = compute_regime_delta(a_state, b_state);

        out.push(OrganelleDelta {
            organelle,
            axes,
            regimes,
        });
    }

    out
}

fn compute_axis_deltas(
    a: &OrganelleState,
    b: &OrganelleState,
    issues: &mut Vec<Issue>,
) -> Vec<AxisDelta> {
    let a_axes = a
        .axes
        .iter()
        .map(|axis| (axis.name.clone(), axis))
        .collect::<BTreeMap<_, _>>();
    let b_axes = b
        .axes
        .iter()
        .map(|axis| (axis.name.clone(), axis))
        .collect::<BTreeMap<_, _>>();

    for key in a_axes.keys() {
        if !b_axes.contains_key(key) {
            issues.push(warn_missing_issue(
                None,
                "MISSING_AXIS_IN_B",
                format!("axis '{}' missing in B", key),
                None,
            ));
        }
    }
    for key in b_axes.keys() {
        if !a_axes.contains_key(key) {
            issues.push(warn_missing_issue(
                None,
                "MISSING_AXIS_IN_A",
                format!("axis '{}' missing in A", key),
                None,
            ));
        }
    }

    let mut out = Vec::new();
    for (name, a_axis) in &a_axes {
        let Some(b_axis) = b_axes.get(name) else {
            continue;
        };
        out.push(AxisDelta {
            name: name.clone(),
            a: axis_to_point(a_axis),
            b: axis_to_point(b_axis),
            delta: AxisPoint {
                median: b_axis.median - a_axis.median,
                p90: opt_delta(a_axis.p90, b_axis.p90),
                p99: opt_delta(a_axis.p99, b_axis.p99),
            },
        });
    }

    out
}

fn compute_regime_delta(
    a: Option<&OrganelleState>,
    b: Option<&OrganelleState>,
) -> Option<RegimeDelta> {
    let a_fractions = a
        .and_then(|s| s.regimes.as_ref())
        .map(|r| r.fractions.clone())
        .unwrap_or_default();
    let b_fractions = b
        .and_then(|s| s.regimes.as_ref())
        .map(|r| r.fractions.clone())
        .unwrap_or_default();

    if a_fractions.is_empty() && b_fractions.is_empty() {
        return None;
    }

    let mut keys = BTreeSet::new();
    for k in a_fractions.keys() {
        keys.insert(k.clone());
    }
    for k in b_fractions.keys() {
        keys.insert(k.clone());
    }

    let mut delta = BTreeMap::new();
    for k in keys {
        let a_v = a_fractions.get(&k).copied().unwrap_or(0.0);
        let b_v = b_fractions.get(&k).copied().unwrap_or(0.0);
        delta.insert(k, b_v - a_v);
    }

    Some(RegimeDelta {
        a: RegimeSide {
            fractions: a_fractions,
        },
        b: RegimeSide {
            fractions: b_fractions,
        },
        delta,
    })
}

fn axis_to_point(axis: &AxisSummary) -> AxisPoint {
    AxisPoint {
        median: axis.median,
        p90: axis.p90,
        p99: axis.p99,
    }
}

fn opt_delta(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(a_v), Some(b_v)) => Some(b_v - a_v),
        _ => None,
    }
}
