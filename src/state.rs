use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

use crate::contracts::types::{Issue, Severity};
use crate::fii::FiiSummary;
use crate::model::organelle::OrganelleId;

pub const STATE_SCHEMA_V1: &str = "kira-organelle-state-v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolMeta {
    pub name: String,
    pub version: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Inputs {
    pub a: String,
    pub b: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AxisSummary {
    pub name: String,
    pub median: f64,
    pub p90: Option<f64>,
    pub p99: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegimeDistribution {
    pub counts: BTreeMap<String, u64>,
    pub fractions: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganelleState {
    pub organelle: OrganelleId,
    pub tool: String,
    pub axes: Vec<AxisSummary>,
    pub regimes: Option<RegimeDistribution>,
    pub qc: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AggregatedState {
    pub schema: String,
    pub tool: ToolMeta,
    pub inputs: Inputs,
    pub organelle_states: Vec<OrganelleState>,
    #[serde(default)]
    pub stress_localization: Option<StressLocalization>,
    #[serde(default)]
    pub functional_irreversibility: Option<FiiSummary>,
    pub issues: Vec<Issue>,
    pub timestamp: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressLocalizationGating {
    pub p_ext_top20: f64,
    pub p_frag_top20: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressLocalization {
    pub sli: f64,
    pub rho_internal_external: f64,
    pub rho_internal_fragility: f64,
    pub gating: StressLocalizationGating,
}

impl AggregatedState {
    pub fn to_json_value(&self) -> Value {
        let mut root = Map::new();
        root.insert("schema".to_string(), Value::String(self.schema.clone()));

        let mut tool = Map::new();
        tool.insert("name".to_string(), Value::String(self.tool.name.clone()));
        tool.insert(
            "version".to_string(),
            Value::String(self.tool.version.clone()),
        );
        root.insert("tool".to_string(), Value::Object(tool));

        let mut inputs = Map::new();
        inputs.insert("a".to_string(), Value::String(self.inputs.a.clone()));
        match &self.inputs.b {
            Some(v) => {
                inputs.insert("b".to_string(), Value::String(v.clone()));
            }
            None => {
                inputs.insert("b".to_string(), Value::Null);
            }
        }
        root.insert("inputs".to_string(), Value::Object(inputs));

        let organelle_states = self
            .organelle_states
            .iter()
            .map(organelle_state_to_value)
            .collect::<Vec<_>>();
        root.insert(
            "organelle_states".to_string(),
            Value::Array(organelle_states),
        );

        match &self.stress_localization {
            Some(v) => {
                let mut stress = Map::new();
                stress.insert("sli".to_string(), Value::from(v.sli));
                stress.insert(
                    "rho_internal_external".to_string(),
                    Value::from(v.rho_internal_external),
                );
                stress.insert(
                    "rho_internal_fragility".to_string(),
                    Value::from(v.rho_internal_fragility),
                );
                let mut gating = Map::new();
                gating.insert("p_ext_top20".to_string(), Value::from(v.gating.p_ext_top20));
                gating.insert(
                    "p_frag_top20".to_string(),
                    Value::from(v.gating.p_frag_top20),
                );
                stress.insert("gating".to_string(), Value::Object(gating));
                root.insert("stress_localization".to_string(), Value::Object(stress));
            }
            None => {
                root.insert("stress_localization".to_string(), Value::Null);
            }
        }

        match &self.functional_irreversibility {
            Some(v) => {
                root.insert(
                    "functional_irreversibility".to_string(),
                    serde_json::to_value(v).unwrap_or(Value::Null),
                );
            }
            None => {
                root.insert("functional_irreversibility".to_string(), Value::Null);
            }
        }

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

fn organelle_state_to_value(state: &OrganelleState) -> Value {
    let mut item = Map::new();
    item.insert(
        "organelle".to_string(),
        serde_json::to_value(state.organelle).unwrap_or(Value::Null),
    );
    item.insert("tool".to_string(), Value::String(state.tool.clone()));

    let axes = state
        .axes
        .iter()
        .map(|axis| {
            let mut obj = Map::new();
            obj.insert("name".to_string(), Value::String(axis.name.clone()));
            obj.insert("median".to_string(), Value::from(axis.median));
            match axis.p90 {
                Some(v) => {
                    obj.insert("p90".to_string(), Value::from(v));
                }
                None => {
                    obj.insert("p90".to_string(), Value::Null);
                }
            }
            match axis.p99 {
                Some(v) => {
                    obj.insert("p99".to_string(), Value::from(v));
                }
                None => {
                    obj.insert("p99".to_string(), Value::Null);
                }
            }
            Value::Object(obj)
        })
        .collect::<Vec<_>>();
    item.insert("axes".to_string(), Value::Array(axes));

    match &state.regimes {
        Some(regime) => {
            let mut regimes = Map::new();
            regimes.insert(
                "counts".to_string(),
                Value::Object(
                    regime
                        .counts
                        .iter()
                        .map(|(k, v)| (k.clone(), Value::from(*v)))
                        .collect(),
                ),
            );
            regimes.insert(
                "fractions".to_string(),
                Value::Object(
                    regime
                        .fractions
                        .iter()
                        .map(|(k, v)| (k.clone(), Value::from(*v)))
                        .collect(),
                ),
            );
            item.insert("regimes".to_string(), Value::Object(regimes));
        }
        None => {
            item.insert("regimes".to_string(), Value::Null);
        }
    }

    item.insert(
        "qc".to_string(),
        Value::Object(
            state
                .qc
                .iter()
                .map(|(k, v)| (k.clone(), Value::from(*v)))
                .collect(),
        ),
    );

    Value::Object(item)
}
