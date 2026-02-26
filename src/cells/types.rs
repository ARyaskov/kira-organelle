use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

use crate::contracts::types::{Issue, Severity};
use crate::model::organelle::OrganelleId;

pub const CELLS_SCHEMA_V1: &str = "kira-organelle-cells-v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CellKey {
    Barcode,
    Sample,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganelleCellState {
    pub axes: BTreeMap<String, f64>,
    pub regime: Option<String>,
    pub confidence: Option<f64>,
    pub flags: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellState {
    pub id: String,
    pub per_organelle: BTreeMap<OrganelleId, OrganelleCellState>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganelleAxes {
    pub organelle: OrganelleId,
    pub axes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellsState {
    pub schema: String,
    pub cell_key: CellKey,
    pub n_cells: usize,
    pub organelle_axes: Vec<OrganelleAxes>,
    pub cells: Vec<CellState>,
    pub issues: Vec<Issue>,
}

impl CellsState {
    pub fn empty(issues: Vec<Issue>) -> Self {
        Self {
            schema: CELLS_SCHEMA_V1.to_string(),
            cell_key: CellKey::Sample,
            n_cells: 0,
            organelle_axes: Vec::new(),
            cells: Vec::new(),
            issues,
        }
    }

    pub fn to_json_value(&self) -> Value {
        let mut root = Map::new();
        root.insert("schema".to_string(), Value::String(self.schema.clone()));
        root.insert(
            "cell_key".to_string(),
            serde_json::to_value(self.cell_key).unwrap_or(Value::Null),
        );
        root.insert("n_cells".to_string(), Value::from(self.n_cells as u64));

        let organelle_axes = self
            .organelle_axes
            .iter()
            .map(|oa| {
                let mut item = Map::new();
                item.insert(
                    "organelle".to_string(),
                    serde_json::to_value(oa.organelle).unwrap_or(Value::Null),
                );
                item.insert(
                    "axes".to_string(),
                    Value::Array(oa.axes.iter().cloned().map(Value::String).collect()),
                );
                Value::Object(item)
            })
            .collect::<Vec<_>>();
        root.insert("organelle_axes".to_string(), Value::Array(organelle_axes));

        let cells = self
            .cells
            .iter()
            .map(|cell| {
                let mut item = Map::new();
                item.insert("id".to_string(), Value::String(cell.id.clone()));

                let mut per_organelle = Map::new();
                for (org, org_state) in &cell.per_organelle {
                    let mut org_item = Map::new();
                    org_item.insert(
                        "axes".to_string(),
                        Value::Object(
                            org_state
                                .axes
                                .iter()
                                .map(|(k, v)| (k.clone(), Value::from(*v)))
                                .collect(),
                        ),
                    );
                    match &org_state.regime {
                        Some(v) => {
                            org_item.insert("regime".to_string(), Value::String(v.clone()));
                        }
                        None => {
                            org_item.insert("regime".to_string(), Value::Null);
                        }
                    }
                    match org_state.confidence {
                        Some(v) => {
                            org_item.insert("confidence".to_string(), Value::from(v));
                        }
                        None => {
                            org_item.insert("confidence".to_string(), Value::Null);
                        }
                    }
                    org_item.insert(
                        "flags".to_string(),
                        Value::Array(org_state.flags.iter().cloned().map(Value::String).collect()),
                    );

                    let key = serde_json::to_value(org)
                        .ok()
                        .and_then(|v| v.as_str().map(ToOwned::to_owned))
                        .unwrap_or_else(|| "unknown".to_string());
                    per_organelle.insert(key, Value::Object(org_item));
                }
                item.insert("per_organelle".to_string(), Value::Object(per_organelle));
                Value::Object(item)
            })
            .collect::<Vec<_>>();
        root.insert("cells".to_string(), Value::Array(cells));

        let issues = self
            .issues
            .iter()
            .map(|issue| {
                let mut item = Map::new();
                item.insert(
                    "severity".to_string(),
                    Value::String(
                        match issue.severity {
                            Severity::Warn => "warn",
                            Severity::Error => "error",
                        }
                        .to_string(),
                    ),
                );
                match &issue.tool {
                    Some(v) => {
                        item.insert("tool".to_string(), Value::String(v.clone()));
                    }
                    None => {
                        item.insert("tool".to_string(), Value::Null);
                    }
                }
                item.insert("code".to_string(), Value::String(issue.code.clone()));
                item.insert("message".to_string(), Value::String(issue.message.clone()));
                match &issue.path {
                    Some(v) => {
                        item.insert("path".to_string(), Value::String(v.clone()));
                    }
                    None => {
                        item.insert("path".to_string(), Value::Null);
                    }
                }
                Value::Object(item)
            })
            .collect::<Vec<_>>();
        root.insert("issues".to_string(), Value::Array(issues));

        Value::Object(root)
    }
}
