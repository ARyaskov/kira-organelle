use std::collections::{BTreeMap, BTreeSet};

use crate::compare::types::ComparisonState;
use crate::interpret::types::{
    EvidenceMetric, EvidenceRef, EvidenceSource, InterpretationSignal, SignalDirection,
    SignalScope, SignalSeverity,
};
use crate::model::organelle::OrganelleId;
use crate::state::AggregatedState;

pub const THRESH_DECAY_ABS: f64 = 0.60;
pub const THRESH_BIOENERGETICS_LOW: f64 = 0.40;
pub const THRESH_DELTA_UP: f64 = 0.10;
pub const THRESH_DELTA_DOWN: f64 = -0.10;
pub const THRESH_STRESS_ABS: f64 = 0.60;

pub fn evaluate_signals(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
) -> Vec<InterpretationSignal> {
    let mut out = Vec::new();
    if let Some(s) = energetic_stress_shift(state, comparison) {
        out.push(s);
    }
    if let Some(s) = proteostasis_overload(state, comparison) {
        out.push(s);
    }
    if let Some(s) = secretory_pressure_increase(state, comparison) {
        out.push(s);
    }
    if let Some(s) = global_stress_adaptation(state, comparison) {
        out.push(s);
    }
    out.extend(stress_localization_flags(state));
    out
}

fn stress_localization_flags(state: &AggregatedState) -> Vec<InterpretationSignal> {
    let Some(sli) = state.stress_localization.as_ref() else {
        return Vec::new();
    };

    let mut out = Vec::new();
    let organelles = vec![
        OrganelleId::Spliceosome,
        OrganelleId::Proteostasis,
        OrganelleId::Ribosome,
        OrganelleId::Nucleus,
        OrganelleId::Secretion,
        OrganelleId::Energetics,
    ];

    let (id, severity) = if sli.sli >= 0.75 {
        ("HIGH_LOCALIZATION", SignalSeverity::Low)
    } else if sli.sli >= 0.50 {
        ("MIXED_LOCALIZATION", SignalSeverity::Medium)
    } else {
        ("LOW_LOCALIZATION", SignalSeverity::High)
    };

    out.push(InterpretationSignal {
        id: id.to_string(),
        scope: SignalScope::Global,
        organelles: organelles.clone(),
        severity,
        direction: SignalDirection::Stable,
        evidence: vec![
            evidence_median(OrganelleId::Spliceosome, "sli", sli.sli),
            evidence_median(
                OrganelleId::Secretion,
                "rho_internal_external",
                sli.rho_internal_external,
            ),
            evidence_median(
                OrganelleId::Energetics,
                "rho_internal_fragility",
                sli.rho_internal_fragility,
            ),
        ],
    });

    if sli.gating.p_ext_top20 > 0.3 {
        out.push(InterpretationSignal {
            id: "STRESS_LEAKAGE_EXTERNAL".to_string(),
            scope: SignalScope::CrossOrganelle,
            organelles: vec![
                OrganelleId::Spliceosome,
                OrganelleId::Secretion,
                OrganelleId::Nucleus,
            ],
            severity: if sli.gating.p_ext_top20 >= 0.5 {
                SignalSeverity::High
            } else {
                SignalSeverity::Medium
            },
            direction: SignalDirection::Up,
            evidence: vec![evidence_median(
                OrganelleId::Secretion,
                "p_ext_top20",
                sli.gating.p_ext_top20,
            )],
        });
    }

    if sli.gating.p_frag_top20 > 0.3 {
        out.push(InterpretationSignal {
            id: "STRESS_ENERGETIC_COUPLING".to_string(),
            scope: SignalScope::CrossOrganelle,
            organelles: vec![
                OrganelleId::Spliceosome,
                OrganelleId::Energetics,
                OrganelleId::Proteostasis,
            ],
            severity: if sli.gating.p_frag_top20 >= 0.5 {
                SignalSeverity::High
            } else {
                SignalSeverity::Medium
            },
            direction: SignalDirection::Up,
            evidence: vec![evidence_median(
                OrganelleId::Energetics,
                "p_frag_top20",
                sli.gating.p_frag_top20,
            )],
        });
    }

    out
}

fn energetic_stress_shift(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
) -> Option<InterpretationSignal> {
    match comparison {
        Some(cmp) => {
            let decay = find_delta(cmp, OrganelleId::Mitochondria, "decay_score")?;
            let bio = find_delta(cmp, OrganelleId::Mitochondria, "bioenergetics")?;
            let mut evidence = Vec::new();
            let mut triggered = false;

            if decay >= THRESH_DELTA_UP {
                triggered = true;
                evidence.push(evidence_delta(
                    OrganelleId::Mitochondria,
                    "decay_score",
                    decay,
                ));
            }
            if bio <= THRESH_DELTA_DOWN {
                triggered = true;
                evidence.push(evidence_delta(
                    OrganelleId::Mitochondria,
                    "bioenergetics",
                    bio,
                ));
            }
            if !triggered {
                return None;
            }

            sort_evidence(&mut evidence);
            Some(InterpretationSignal {
                id: "energetic_stress_shift".to_string(),
                scope: SignalScope::Organelle,
                organelles: vec![OrganelleId::Mitochondria],
                severity: severity_from_max_abs(&evidence),
                direction: if evidence.len() > 1 {
                    SignalDirection::Mixed
                } else if evidence[0].axis == "bioenergetics" {
                    SignalDirection::Down
                } else {
                    SignalDirection::Up
                },
                evidence,
            })
        }
        None => {
            let decay = find_median(state, OrganelleId::Mitochondria, "decay_score")?;
            let bio = find_median(state, OrganelleId::Mitochondria, "bioenergetics")?;
            let mut evidence = Vec::new();
            if decay >= THRESH_DECAY_ABS {
                evidence.push(evidence_median(
                    OrganelleId::Mitochondria,
                    "decay_score",
                    decay,
                ));
            }
            if bio <= THRESH_BIOENERGETICS_LOW {
                evidence.push(evidence_median(
                    OrganelleId::Mitochondria,
                    "bioenergetics",
                    bio,
                ));
            }
            if evidence.is_empty() {
                return None;
            }
            sort_evidence(&mut evidence);

            Some(InterpretationSignal {
                id: "energetic_stress_shift".to_string(),
                scope: SignalScope::Organelle,
                organelles: vec![OrganelleId::Mitochondria],
                severity: severity_from_max_abs(&evidence),
                direction: SignalDirection::Mixed,
                evidence,
            })
        }
    }
}

fn proteostasis_overload(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
) -> Option<InterpretationSignal> {
    let (p, a, source, metric) = match comparison {
        Some(cmp) => (
            find_delta(cmp, OrganelleId::Proteostasis, "proteostasis_load")?,
            find_delta(cmp, OrganelleId::Autophagy, "stress_index")?,
            EvidenceSource::Comparison,
            EvidenceMetric::Delta,
        ),
        None => (
            find_median(state, OrganelleId::Proteostasis, "proteostasis_load")?,
            find_median(state, OrganelleId::Autophagy, "stress_index")?,
            EvidenceSource::State,
            EvidenceMetric::Median,
        ),
    };

    let pass = if comparison.is_some() {
        p >= THRESH_DELTA_UP && a >= THRESH_DELTA_UP
    } else {
        p >= THRESH_STRESS_ABS && a >= THRESH_STRESS_ABS
    };
    if !pass {
        return None;
    }

    let mut evidence = vec![
        EvidenceRef {
            source,
            organelle: OrganelleId::Proteostasis,
            axis: "proteostasis_load".to_string(),
            metric,
            value: p,
        },
        EvidenceRef {
            source,
            organelle: OrganelleId::Autophagy,
            axis: "stress_index".to_string(),
            metric,
            value: a,
        },
    ];
    sort_evidence(&mut evidence);

    Some(InterpretationSignal {
        id: "proteostasis_overload".to_string(),
        scope: SignalScope::CrossOrganelle,
        organelles: vec![OrganelleId::Proteostasis, OrganelleId::Autophagy],
        severity: severity_from_max_abs(&evidence),
        direction: SignalDirection::Up,
        evidence,
    })
}

fn secretory_pressure_increase(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
) -> Option<InterpretationSignal> {
    let (er, sec, source, metric) = match comparison {
        Some(cmp) => (
            find_delta(cmp, OrganelleId::Secretion, "er_golgi_pressure")?,
            find_delta(cmp, OrganelleId::Secretion, "stress_secretion_index")?,
            EvidenceSource::Comparison,
            EvidenceMetric::Delta,
        ),
        None => (
            find_median(state, OrganelleId::Secretion, "er_golgi_pressure")?,
            find_median(state, OrganelleId::Secretion, "stress_secretion_index")?,
            EvidenceSource::State,
            EvidenceMetric::Median,
        ),
    };

    let mut evidence = Vec::new();
    let thr = if comparison.is_some() {
        THRESH_DELTA_UP
    } else {
        THRESH_STRESS_ABS
    };

    if er >= thr {
        evidence.push(EvidenceRef {
            source,
            organelle: OrganelleId::Secretion,
            axis: "er_golgi_pressure".to_string(),
            metric,
            value: er,
        });
    }
    if sec >= thr {
        evidence.push(EvidenceRef {
            source,
            organelle: OrganelleId::Secretion,
            axis: "stress_secretion_index".to_string(),
            metric,
            value: sec,
        });
    }
    if evidence.is_empty() {
        return None;
    }

    sort_evidence(&mut evidence);
    Some(InterpretationSignal {
        id: "secretory_pressure_increase".to_string(),
        scope: SignalScope::Organelle,
        organelles: vec![OrganelleId::Secretion],
        severity: severity_from_max_abs(&evidence),
        direction: SignalDirection::Up,
        evidence,
    })
}

fn global_stress_adaptation(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
) -> Option<InterpretationSignal> {
    let mut stressed = BTreeSet::new();
    let mut evidence = Vec::new();

    match comparison {
        Some(cmp) => {
            for org in &cmp.organelle_deltas {
                for axis in &org.axes {
                    if axis.name.contains("stress") && axis.delta.median >= THRESH_DELTA_UP {
                        stressed.insert(org.organelle);
                        evidence.push(evidence_delta(org.organelle, &axis.name, axis.delta.median));
                    }
                }
            }
        }
        None => {
            for org in &state.organelle_states {
                for axis in &org.axes {
                    if axis.name.contains("stress") && axis.median >= THRESH_STRESS_ABS {
                        stressed.insert(org.organelle);
                        evidence.push(evidence_median(org.organelle, &axis.name, axis.median));
                    }
                }
            }
        }
    }

    if stressed.len() < 3 {
        return None;
    }

    sort_evidence(&mut evidence);
    let mut organelles = stressed.into_iter().collect::<Vec<_>>();
    organelles.sort();

    Some(InterpretationSignal {
        id: "global_stress_adaptation".to_string(),
        scope: SignalScope::Global,
        organelles,
        severity: severity_from_max_abs(&evidence),
        direction: SignalDirection::Up,
        evidence,
    })
}

fn find_median(state: &AggregatedState, organelle: OrganelleId, axis: &str) -> Option<f64> {
    state
        .organelle_states
        .iter()
        .find(|o| o.organelle == organelle)
        .and_then(|o| o.axes.iter().find(|a| a.name == axis))
        .map(|a| a.median)
}

fn find_delta(cmp: &ComparisonState, organelle: OrganelleId, axis: &str) -> Option<f64> {
    cmp.organelle_deltas
        .iter()
        .find(|o| o.organelle == organelle)
        .and_then(|o| o.axes.iter().find(|a| a.name == axis))
        .map(|a| a.delta.median)
}

fn evidence_median(org: OrganelleId, axis: &str, value: f64) -> EvidenceRef {
    EvidenceRef {
        source: EvidenceSource::State,
        organelle: org,
        axis: axis.to_string(),
        metric: EvidenceMetric::Median,
        value,
    }
}

fn evidence_delta(org: OrganelleId, axis: &str, value: f64) -> EvidenceRef {
    EvidenceRef {
        source: EvidenceSource::Comparison,
        organelle: org,
        axis: axis.to_string(),
        metric: EvidenceMetric::Delta,
        value,
    }
}

fn sort_evidence(evidence: &mut [EvidenceRef]) {
    evidence.sort_by(|a, b| {
        a.organelle
            .cmp(&b.organelle)
            .then_with(|| a.axis.cmp(&b.axis))
            .then_with(|| a.metric.cmp(&b.metric))
            .then_with(|| a.source.cmp(&b.source))
            .then_with(|| a.value.total_cmp(&b.value))
    });
}

fn severity_from_max_abs(evidence: &[EvidenceRef]) -> SignalSeverity {
    let max = evidence
        .iter()
        .map(|e| e.value.abs())
        .fold(0.0f64, f64::max);

    if max >= 0.8 {
        SignalSeverity::High
    } else if max >= 0.4 {
        SignalSeverity::Medium
    } else {
        SignalSeverity::Low
    }
}

pub fn confidence_score(state: &AggregatedState, signals: &[InterpretationSignal]) -> f64 {
    let low_conf = estimate_low_conf_fraction(state);
    let coverage = state.organelle_states.len() as f64 / 7.0;
    let signal_factor = (signals.len() as f64 / 4.0).min(1.0);

    let score = (1.0 - low_conf).clamp(0.0, 1.0) * 0.5
        + coverage.clamp(0.0, 1.0) * 0.3
        + signal_factor * 0.2;
    score.clamp(0.0, 1.0)
}

fn estimate_low_conf_fraction(state: &AggregatedState) -> f64 {
    let mut vals = Vec::new();
    for org in &state.organelle_states {
        for (k, v) in &org.qc {
            if k.contains("low_confidence_fraction") {
                vals.push(*v);
            }
        }
    }
    if vals.is_empty() {
        return 0.0;
    }
    vals.iter().sum::<f64>() / vals.len() as f64
}

pub fn dedup_and_sort_signals(signals: &mut Vec<InterpretationSignal>) {
    let mut map = BTreeMap::new();
    for signal in signals.drain(..) {
        map.insert(signal.id.clone(), signal);
    }

    let mut out = map.into_values().collect::<Vec<_>>();
    out.sort_by(|a, b| {
        a.scope
            .cmp(&b.scope)
            .then_with(|| a.severity.cmp(&b.severity))
            .then_with(|| a.id.cmp(&b.id))
    });
    *signals = out;
}
