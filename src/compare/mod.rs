pub mod cells;
pub mod state;
pub mod types;

use std::fs;
use std::path::Path;

use crate::cells::types::CellsState;
use crate::contracts::types::Issue;
use crate::state::AggregatedState;
use crate::warn_missing_issue;

use types::{COMPARISON_SCHEMA_V1, ComparisonState};

pub fn build_comparison(
    a_state: &AggregatedState,
    a_cells: &CellsState,
    b_state: &AggregatedState,
    b_cells: &CellsState,
    input_a: &Path,
    input_b: &Path,
    mut issues: Vec<Issue>,
) -> ComparisonState {
    if a_cells.cell_key != b_cells.cell_key {
        issues.push(warn_missing_issue(
            None,
            "CELL_KEY_MISMATCH",
            "cells.json uses different cell_key between A and B".to_string(),
            None,
        ));
    }

    let organelle_deltas = state::compute_organelle_deltas(a_state, b_state, &mut issues);
    let cell_deltas = cells::compute_cell_delta_summary(a_cells, b_cells);

    ComparisonState {
        schema: COMPARISON_SCHEMA_V1.to_string(),
        inputs_a: input_a.to_string_lossy().to_string(),
        inputs_b: input_b.to_string_lossy().to_string(),
        organelle_deltas,
        cell_deltas,
        issues,
        timestamp: crate::io::deterministic_rfc3339_utc(),
    }
}

pub fn load_side_from_output(
    input_root: &Path,
    issues: &mut Vec<Issue>,
) -> Option<(AggregatedState, CellsState)> {
    let out_dir = input_root.join("kira-organelle");
    let state_path = out_dir.join("state.json");
    let cells_path = out_dir.join("cells.json");

    let state_raw = match fs::read_to_string(&state_path) {
        Ok(v) => v,
        Err(_) => {
            issues.push(warn_missing_issue(
                None,
                "MISSING_SIDE_STATE",
                "no precomputed state.json for side; will use raw fallback".to_string(),
                Some(state_path.to_string_lossy().to_string()),
            ));
            return None;
        }
    };
    let cells_raw = match fs::read_to_string(&cells_path) {
        Ok(v) => v,
        Err(_) => {
            issues.push(warn_missing_issue(
                None,
                "MISSING_SIDE_CELLS",
                "no precomputed cells.json for side; will use raw fallback".to_string(),
                Some(cells_path.to_string_lossy().to_string()),
            ));
            return None;
        }
    };

    let state: AggregatedState = match serde_json::from_str(&state_raw) {
        Ok(v) => v,
        Err(e) => {
            issues.push(warn_missing_issue(
                None,
                "INVALID_SIDE_STATE",
                format!("failed parsing side state.json: {e}"),
                Some(state_path.to_string_lossy().to_string()),
            ));
            return None;
        }
    };

    let cells: CellsState = match serde_json::from_str(&cells_raw) {
        Ok(v) => v,
        Err(e) => {
            issues.push(warn_missing_issue(
                None,
                "INVALID_SIDE_CELLS",
                format!("failed parsing side cells.json: {e}"),
                Some(cells_path.to_string_lossy().to_string()),
            ));
            return None;
        }
    };

    Some((state, cells))
}
