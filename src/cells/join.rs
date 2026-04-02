use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use crate::contracts::types::{Issue, ToolContracts};
use crate::model::organelle::OrganelleId;
use crate::warn_missing_issue;

use super::read::{ToolCellData, read_tool_primary_metrics};
use super::types::{CellKey, CellState, CellsState, OrganelleAxes};

pub fn build_cells_state(
    input_root: &Path,
    tools: &[ToolContracts],
    base_issues: &[Issue],
) -> CellsState {
    let mut issues = base_issues.to_vec();

    let mut tool_rows = Vec::new();
    for tool in tools {
        let Some(organelle) = OrganelleId::from_tool_name(&tool.name) else {
            issues.push(warn_missing_issue(
                Some(&tool.name),
                "UNKNOWN_TOOL_MAPPING",
                "tool is discovered but has no organelle mapping".to_string(),
                None,
            ));
            continue;
        };

        match read_tool_primary_metrics(
            input_root,
            &tool.name,
            organelle,
            tool.primary_metrics_path.as_deref(),
            &mut issues,
        ) {
            Some(data) => tool_rows.push(data),
            None => continue,
        }
    }

    if tool_rows.is_empty() {
        return CellsState::empty(issues);
    }

    let cell_key = detect_global_cell_key(&tool_rows);
    let mut global_cells: BTreeMap<String, CellState> = BTreeMap::new();
    let mut ingestion_diagnostics = Vec::new();

    let mut axes_by_organelle: BTreeMap<OrganelleId, BTreeSet<String>> = BTreeMap::new();

    for data in &tool_rows {
        ingestion_diagnostics.push(data.ingestion_diagnostics.clone());
        let axis_set = axes_by_organelle.entry(data.organelle).or_default();
        for axis in &data.axes_union {
            axis_set.insert(axis.clone());
        }

        for row in &data.rows {
            let entry = global_cells
                .entry(row.id.clone())
                .or_insert_with(|| CellState {
                    id: row.id.clone(),
                    per_organelle: BTreeMap::new(),
                });
            entry
                .per_organelle
                .insert(data.organelle, row.organelle_state.clone());
        }
    }

    let organelle_axes = axes_by_organelle
        .into_iter()
        .map(|(organelle, axes)| OrganelleAxes {
            organelle,
            axes: axes.into_iter().collect(),
        })
        .collect::<Vec<_>>();

    let cells = global_cells.into_values().collect::<Vec<_>>();
    ingestion_diagnostics.sort_by(|a, b| a.tool.cmp(&b.tool));

    CellsState {
        schema: super::types::CELLS_SCHEMA_V1.to_string(),
        cell_key,
        n_cells: cells.len(),
        organelle_axes,
        ingestion_diagnostics,
        normalization_context: None,
        cells,
        issues,
    }
}

fn detect_global_cell_key(data: &[ToolCellData]) -> CellKey {
    if data.iter().any(|d| d.cell_key == CellKey::Barcode) {
        CellKey::Barcode
    } else {
        CellKey::Sample
    }
}
