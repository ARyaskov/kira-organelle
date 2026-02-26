use std::collections::BTreeMap;

use crate::cells::types::CellsState;
use crate::util::select::{median_in_place, percentile_nearest_rank_in_place};

use super::types::{CellAxisDeltaSummary, CellDeltaSummary};

pub fn compute_cell_delta_summary(a: &CellsState, b: &CellsState) -> CellDeltaSummary {
    let a_map = a
        .cells
        .iter()
        .map(|cell| (cell.id.as_str(), cell))
        .collect::<BTreeMap<_, _>>();
    let b_map = b
        .cells
        .iter()
        .map(|cell| (cell.id.as_str(), cell))
        .collect::<BTreeMap<_, _>>();

    let common_ids = a_map
        .keys()
        .filter(|id| b_map.contains_key(**id))
        .copied()
        .collect::<Vec<_>>();

    let mut deltas: BTreeMap<(crate::model::organelle::OrganelleId, String), Vec<f64>> =
        BTreeMap::new();

    for id in &common_ids {
        let a_cell = a_map.get(id).expect("exists");
        let b_cell = b_map.get(id).expect("exists");

        for (organelle, a_org) in &a_cell.per_organelle {
            let Some(b_org) = b_cell.per_organelle.get(organelle) else {
                continue;
            };

            for (axis, a_val) in &a_org.axes {
                let Some(b_val) = b_org.axes.get(axis) else {
                    continue;
                };
                deltas
                    .entry((*organelle, axis.clone()))
                    .or_default()
                    .push(*b_val - *a_val);
            }
        }
    }

    let axes = deltas
        .into_iter()
        .filter_map(|((organelle, axis), values)| {
            if values.is_empty() {
                return None;
            }
            let mut median_buf = values.clone();
            let median_delta = median_in_place(&mut median_buf)?;
            let mut abs_buf = values;
            for v in &mut abs_buf {
                *v = v.abs();
            }
            let p90_abs_delta = percentile_nearest_rank_in_place(&mut abs_buf, 0.90)?;

            Some(CellAxisDeltaSummary {
                organelle,
                axis,
                median_delta,
                p90_abs_delta,
            })
        })
        .collect::<Vec<_>>();

    CellDeltaSummary {
        n_common_cells: common_ids.len() as u64,
        axes,
    }
}
