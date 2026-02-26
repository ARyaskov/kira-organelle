use crate::compare::types::ComparisonState;
use crate::contracts::types::Issue;
use crate::state::AggregatedState;

use super::rules::{confidence_score, dedup_and_sort_signals, evaluate_signals};
use super::types::{INTERPRETATION_SCHEMA_V1, InterpretationState};

pub fn build_interpretation(
    state: &AggregatedState,
    comparison: Option<&ComparisonState>,
    issues: Vec<Issue>,
) -> InterpretationState {
    let mut signals = evaluate_signals(state, comparison);
    dedup_and_sort_signals(&mut signals);

    let confidence = confidence_score(state, &signals);

    InterpretationState {
        schema: INTERPRETATION_SCHEMA_V1.to_string(),
        mode: if comparison.is_some() {
            "comparison".to_string()
        } else {
            "single".to_string()
        },
        signals,
        confidence,
        issues,
        timestamp: crate::io::deterministic_rfc3339_utc(),
    }
}
