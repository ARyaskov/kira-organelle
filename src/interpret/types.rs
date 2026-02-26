use serde_json::{Map, Value};

use crate::contracts::types::{Issue, Severity};
use crate::model::organelle::OrganelleId;

pub const INTERPRETATION_SCHEMA_V1: &str = "kira-organelle-interpretation-v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SignalScope {
    Organelle,
    CrossOrganelle,
    Global,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SignalSeverity {
    Low,
    Medium,
    High,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SignalDirection {
    Up,
    Down,
    Mixed,
    Stable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EvidenceSource {
    State,
    Comparison,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EvidenceMetric {
    Median,
    Delta,
}

#[derive(Debug, Clone)]
pub struct EvidenceRef {
    pub source: EvidenceSource,
    pub organelle: OrganelleId,
    pub axis: String,
    pub metric: EvidenceMetric,
    pub value: f64,
}

#[derive(Debug, Clone)]
pub struct InterpretationSignal {
    pub id: String,
    pub scope: SignalScope,
    pub organelles: Vec<OrganelleId>,
    pub severity: SignalSeverity,
    pub direction: SignalDirection,
    pub evidence: Vec<EvidenceRef>,
}

#[derive(Debug, Clone)]
pub struct InterpretationState {
    pub schema: String,
    pub mode: String,
    pub signals: Vec<InterpretationSignal>,
    pub confidence: f64,
    pub issues: Vec<Issue>,
    pub timestamp: String,
}

impl InterpretationState {
    pub fn to_json_value(&self) -> Value {
        let mut root = Map::new();
        root.insert("schema".to_string(), Value::String(self.schema.clone()));
        root.insert("mode".to_string(), Value::String(self.mode.clone()));

        let signals = self
            .signals
            .iter()
            .map(|s| {
                let mut item = Map::new();
                item.insert("id".to_string(), Value::String(s.id.clone()));
                item.insert(
                    "scope".to_string(),
                    Value::String(scope_str(s.scope).to_string()),
                );
                item.insert(
                    "organelles".to_string(),
                    Value::Array(
                        s.organelles
                            .iter()
                            .map(|o| serde_json::to_value(o).unwrap_or(Value::Null))
                            .collect(),
                    ),
                );
                item.insert(
                    "severity".to_string(),
                    Value::String(severity_str(s.severity).to_string()),
                );
                item.insert(
                    "direction".to_string(),
                    Value::String(direction_str(s.direction).to_string()),
                );

                let evidence = s
                    .evidence
                    .iter()
                    .map(|e| {
                        let mut ev = Map::new();
                        ev.insert(
                            "source".to_string(),
                            Value::String(source_str(e.source).to_string()),
                        );
                        ev.insert(
                            "organelle".to_string(),
                            serde_json::to_value(e.organelle).unwrap_or(Value::Null),
                        );
                        ev.insert("axis".to_string(), Value::String(e.axis.clone()));
                        ev.insert(
                            "metric".to_string(),
                            Value::String(metric_str(e.metric).to_string()),
                        );
                        ev.insert("value".to_string(), Value::from(e.value));
                        Value::Object(ev)
                    })
                    .collect::<Vec<_>>();
                item.insert("evidence".to_string(), Value::Array(evidence));
                Value::Object(item)
            })
            .collect::<Vec<_>>();
        root.insert("signals".to_string(), Value::Array(signals));

        root.insert("confidence".to_string(), Value::from(self.confidence));

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

pub fn scope_str(v: SignalScope) -> &'static str {
    match v {
        SignalScope::Organelle => "organelle",
        SignalScope::CrossOrganelle => "cross-organelle",
        SignalScope::Global => "global",
    }
}

pub fn severity_str(v: SignalSeverity) -> &'static str {
    match v {
        SignalSeverity::Low => "low",
        SignalSeverity::Medium => "medium",
        SignalSeverity::High => "high",
    }
}

pub fn direction_str(v: SignalDirection) -> &'static str {
    match v {
        SignalDirection::Up => "up",
        SignalDirection::Down => "down",
        SignalDirection::Mixed => "mixed",
        SignalDirection::Stable => "stable",
    }
}

pub fn source_str(v: EvidenceSource) -> &'static str {
    match v {
        EvidenceSource::State => "state",
        EvidenceSource::Comparison => "comparison",
    }
}

pub fn metric_str(v: EvidenceMetric) -> &'static str {
    match v {
        EvidenceMetric::Median => "median",
        EvidenceMetric::Delta => "delta",
    }
}
