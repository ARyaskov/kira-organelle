use std::collections::BTreeMap;

use kira_organelle::cells::types::{
    CellKey, CellState, CellsState, OrganelleAxes, OrganelleCellState,
};
use kira_organelle::model::organelle::OrganelleId;
use kira_organelle::stress_localization::compute;

fn build_synthetic_cells(n: usize, coupled: bool, zero_variance: bool) -> CellsState {
    let mut cells = Vec::new();
    for i in 0..n {
        let base = if zero_variance {
            0.5
        } else {
            i as f64 / n as f64
        };
        let ext = if coupled {
            base
        } else if zero_variance {
            0.5
        } else {
            1.0 - base
        };
        let frag = if coupled {
            base
        } else if zero_variance {
            0.5
        } else {
            1.0 - base
        };

        let mut per_organelle = BTreeMap::new();
        per_organelle.insert(
            OrganelleId::Spliceosome,
            OrganelleCellState {
                axes: BTreeMap::from([("splice_stress_index".to_string(), base)]),
                regime: None,
                confidence: Some(1.0),
                flags: vec![],
            },
        );
        per_organelle.insert(
            OrganelleId::Proteostasis,
            OrganelleCellState {
                axes: BTreeMap::from([("proteo_stress_index".to_string(), base)]),
                regime: None,
                confidence: Some(1.0),
                flags: vec![],
            },
        );
        per_organelle.insert(
            OrganelleId::Ribosome,
            OrganelleCellState {
                axes: BTreeMap::from([("translation_stress_index".to_string(), base)]),
                regime: None,
                confidence: Some(1.0),
                flags: vec![],
            },
        );
        per_organelle.insert(
            OrganelleId::Nucleus,
            OrganelleCellState {
                axes: BTreeMap::from([("nuclear_stress_index".to_string(), base)]),
                regime: None,
                confidence: Some(1.0),
                flags: vec![],
            },
        );
        per_organelle.insert(
            OrganelleId::Secretion,
            OrganelleCellState {
                axes: BTreeMap::from([
                    ("stress_secretion_index".to_string(), ext),
                    ("er_golgi_pressure".to_string(), ext),
                ]),
                regime: None,
                confidence: Some(1.0),
                flags: vec![],
            },
        );
        per_organelle.insert(
            OrganelleId::Energetics,
            OrganelleCellState {
                axes: BTreeMap::from([
                    ("barrier_min".to_string(), 1.0 - frag),
                    ("unstable_trap_flag".to_string(), frag),
                    ("trap_depth".to_string(), frag),
                ]),
                regime: None,
                confidence: Some(1.0),
                flags: vec![],
            },
        );

        cells.push(CellState {
            id: format!("cell_{i}"),
            per_organelle,
        });
    }

    CellsState {
        schema: "kira-organelle-cells-v1".to_string(),
        cell_key: CellKey::Barcode,
        n_cells: cells.len(),
        organelle_axes: vec![
            OrganelleAxes {
                organelle: OrganelleId::Spliceosome,
                axes: vec!["splice_stress_index".to_string()],
            },
            OrganelleAxes {
                organelle: OrganelleId::Proteostasis,
                axes: vec!["proteo_stress_index".to_string()],
            },
            OrganelleAxes {
                organelle: OrganelleId::Ribosome,
                axes: vec!["translation_stress_index".to_string()],
            },
            OrganelleAxes {
                organelle: OrganelleId::Nucleus,
                axes: vec!["nuclear_stress_index".to_string()],
            },
            OrganelleAxes {
                organelle: OrganelleId::Secretion,
                axes: vec![
                    "stress_secretion_index".to_string(),
                    "er_golgi_pressure".to_string(),
                ],
            },
            OrganelleAxes {
                organelle: OrganelleId::Energetics,
                axes: vec![
                    "barrier_min".to_string(),
                    "unstable_trap_flag".to_string(),
                    "trap_depth".to_string(),
                ],
            },
        ],
        cells,
        issues: vec![],
    }
}

#[test]
fn synthetic_decoupled_has_high_sli() {
    let cells = build_synthetic_cells(200, false, false);
    let mut issues = Vec::new();
    let sli = compute(&cells, &mut issues).expect("sli");
    assert!(issues.is_empty());
    assert!(sli.sli > 0.75, "sli={}", sli.sli);
}

#[test]
fn synthetic_coupled_has_low_sli() {
    let cells = build_synthetic_cells(200, true, false);
    let mut issues = Vec::new();
    let sli = compute(&cells, &mut issues).expect("sli");
    assert!(issues.is_empty());
    assert!(sli.sli < 0.50, "sli={}", sli.sli);
}

#[test]
fn zero_variance_rhos_are_zero() {
    let cells = build_synthetic_cells(200, true, true);
    let mut issues = Vec::new();
    let sli = compute(&cells, &mut issues).expect("sli");
    assert!(issues.is_empty());
    assert_eq!(sli.rho_internal_external, 0.0);
    assert_eq!(sli.rho_internal_fragility, 0.0);
    assert!((sli.sli - 1.0).abs() < 1e-9);
}
