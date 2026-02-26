use std::collections::BTreeMap;

use serde_json::{Map, Value};

use crate::contracts::types::{Issue, Severity};
use crate::model::organelle::OrganelleId;

pub const COMPARISON_SCHEMA_V1: &str = "kira-organelle-comparison-v1";

#[derive(Debug, Clone)]
pub struct AxisPoint {
    pub median: f64,
    pub p90: Option<f64>,
    pub p99: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct AxisDelta {
    pub name: String,
    pub a: AxisPoint,
    pub b: AxisPoint,
    pub delta: AxisPoint,
}

#[derive(Debug, Clone)]
pub struct RegimeSide {
    pub fractions: BTreeMap<String, f64>,
}

#[derive(Debug, Clone)]
pub struct RegimeDelta {
    pub a: RegimeSide,
    pub b: RegimeSide,
    pub delta: BTreeMap<String, f64>,
}

#[derive(Debug, Clone)]
pub struct OrganelleDelta {
    pub organelle: OrganelleId,
    pub axes: Vec<AxisDelta>,
    pub regimes: Option<RegimeDelta>,
}

#[derive(Debug, Clone)]
pub struct CellAxisDeltaSummary {
    pub organelle: OrganelleId,
    pub axis: String,
    pub median_delta: f64,
    pub p90_abs_delta: f64,
}

#[derive(Debug, Clone)]
pub struct CellDeltaSummary {
    pub n_common_cells: u64,
    pub axes: Vec<CellAxisDeltaSummary>,
}

#[derive(Debug, Clone)]
pub struct ComparisonState {
    pub schema: String,
    pub inputs_a: String,
    pub inputs_b: String,
    pub organelle_deltas: Vec<OrganelleDelta>,
    pub cell_deltas: CellDeltaSummary,
    pub issues: Vec<Issue>,
    pub timestamp: String,
}

impl ComparisonState {
    pub fn to_json_value(&self) -> Value {
        let mut root = Map::new();
        root.insert("schema".to_string(), Value::String(self.schema.clone()));

        let mut inputs = Map::new();
        inputs.insert("a".to_string(), Value::String(self.inputs_a.clone()));
        inputs.insert("b".to_string(), Value::String(self.inputs_b.clone()));
        root.insert("inputs".to_string(), Value::Object(inputs));

        let organelle_deltas = self
            .organelle_deltas
            .iter()
            .map(|od| {
                let mut item = Map::new();
                item.insert(
                    "organelle".to_string(),
                    serde_json::to_value(od.organelle).unwrap_or(Value::Null),
                );

                let axes = od
                    .axes
                    .iter()
                    .map(|axis| {
                        let mut axis_item = Map::new();
                        axis_item.insert("name".to_string(), Value::String(axis.name.clone()));
                        axis_item.insert("a".to_string(), axis_point_to_value(&axis.a));
                        axis_item.insert("b".to_string(), axis_point_to_value(&axis.b));
                        axis_item.insert("delta".to_string(), axis_point_to_value(&axis.delta));
                        Value::Object(axis_item)
                    })
                    .collect::<Vec<_>>();
                item.insert("axes".to_string(), Value::Array(axes));

                match &od.regimes {
                    Some(reg) => {
                        let mut reg_item = Map::new();
                        reg_item.insert(
                            "a".to_string(),
                            Value::Object({
                                let mut x = Map::new();
                                x.insert(
                                    "fractions".to_string(),
                                    Value::Object(
                                        reg.a
                                            .fractions
                                            .iter()
                                            .map(|(k, v)| (k.clone(), Value::from(*v)))
                                            .collect(),
                                    ),
                                );
                                x
                            }),
                        );
                        reg_item.insert(
                            "b".to_string(),
                            Value::Object({
                                let mut x = Map::new();
                                x.insert(
                                    "fractions".to_string(),
                                    Value::Object(
                                        reg.b
                                            .fractions
                                            .iter()
                                            .map(|(k, v)| (k.clone(), Value::from(*v)))
                                            .collect(),
                                    ),
                                );
                                x
                            }),
                        );
                        reg_item.insert(
                            "delta".to_string(),
                            Value::Object(
                                reg.delta
                                    .iter()
                                    .map(|(k, v)| (k.clone(), Value::from(*v)))
                                    .collect(),
                            ),
                        );
                        item.insert("regimes".to_string(), Value::Object(reg_item));
                    }
                    None => {
                        item.insert("regimes".to_string(), Value::Null);
                    }
                }

                Value::Object(item)
            })
            .collect::<Vec<_>>();
        root.insert(
            "organelle_deltas".to_string(),
            Value::Array(organelle_deltas),
        );

        let mut cell_deltas = Map::new();
        cell_deltas.insert(
            "n_common_cells".to_string(),
            Value::from(self.cell_deltas.n_common_cells),
        );
        cell_deltas.insert(
            "axes".to_string(),
            Value::Array(
                self.cell_deltas
                    .axes
                    .iter()
                    .map(|x| {
                        let mut item = Map::new();
                        item.insert(
                            "organelle".to_string(),
                            serde_json::to_value(x.organelle).unwrap_or(Value::Null),
                        );
                        item.insert("axis".to_string(), Value::String(x.axis.clone()));
                        item.insert("median_delta".to_string(), Value::from(x.median_delta));
                        item.insert("p90_abs_delta".to_string(), Value::from(x.p90_abs_delta));
                        Value::Object(item)
                    })
                    .collect(),
            ),
        );
        root.insert("cell_deltas".to_string(), Value::Object(cell_deltas));

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

        root.insert(
            "timestamp".to_string(),
            Value::String(self.timestamp.clone()),
        );

        Value::Object(root)
    }
}

fn axis_point_to_value(p: &AxisPoint) -> Value {
    let mut item = Map::new();
    item.insert("median".to_string(), Value::from(p.median));
    match p.p90 {
        Some(v) => {
            item.insert("p90".to_string(), Value::from(v));
        }
        None => {
            item.insert("p90".to_string(), Value::Null);
        }
    }
    match p.p99 {
        Some(v) => {
            item.insert("p99".to_string(), Value::from(v));
        }
        None => {
            item.insert("p99".to_string(), Value::Null);
        }
    }
    Value::Object(item)
}
