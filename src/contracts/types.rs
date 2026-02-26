use std::path::PathBuf;

use serde::{Deserialize, Serialize};
use serde_json::Value;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    Warn,
    Error,
}

impl Severity {
    pub fn is_error(&self) -> bool {
        matches!(self, Severity::Error)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Issue {
    pub severity: Severity,
    pub tool: Option<String>,
    pub code: String,
    pub message: String,
    pub path: Option<String>,
}

#[derive(Debug, Clone)]
pub struct ToolContracts {
    pub name: String,
    pub dir: PathBuf,
    pub summary_path: PathBuf,
    pub summary_json: Value,
    pub pipeline_step_path: Option<PathBuf>,
    pub primary_metrics_path: Option<String>,
}
