use std::fs;
use std::path::Path;

use serde_json::Value;

use crate::{strict_missing_issue, warn_missing_issue};

use super::types::{Issue, ToolContracts};

pub fn read_tool_contracts(
    tool: &str,
    tool_dir: &Path,
    strict: bool,
    issues: &mut Vec<Issue>,
) -> Result<Option<ToolContracts>, String> {
    let summary_path = tool_dir.join("summary.json");
    if !summary_path.is_file() {
        let issue = if strict {
            strict_missing_issue(
                Some(tool),
                "MISSING_SUMMARY",
                "required summary.json missing".to_string(),
                Some(summary_path.to_string_lossy().to_string()),
            )
        } else {
            warn_missing_issue(
                Some(tool),
                "MISSING_SUMMARY",
                "required summary.json missing".to_string(),
                Some(summary_path.to_string_lossy().to_string()),
            )
        };
        issues.push(issue);

        if strict {
            return Err(format!("strict mode: missing summary.json for {tool}"));
        }
        return Ok(None);
    }

    let summary_json = read_json_value(&summary_path).inspect_err(|e| {
        issues.push(if strict {
            strict_missing_issue(
                Some(tool),
                "SUMMARY_PARSE_ERROR",
                e.clone(),
                Some(summary_path.to_string_lossy().to_string()),
            )
        } else {
            warn_missing_issue(
                Some(tool),
                "SUMMARY_PARSE_ERROR",
                e.clone(),
                Some(summary_path.to_string_lossy().to_string()),
            )
        });
    })?;

    let pipeline_step_path = tool_dir.join("pipeline_step.json");
    let (pipeline_step_path_opt, primary_metrics_path) = if pipeline_step_path.is_file() {
        match read_json_value(&pipeline_step_path) {
            Ok(v) => {
                let primary = v
                    .get("artifacts")
                    .and_then(|a| a.get("primary_metrics"))
                    .and_then(Value::as_str)
                    .map(ToOwned::to_owned);
                let _summary_artifact = v
                    .get("artifacts")
                    .and_then(|a| a.get("summary"))
                    .and_then(Value::as_str)
                    .map(ToOwned::to_owned);
                (Some(pipeline_step_path.clone()), primary)
            }
            Err(e) => {
                issues.push(if strict {
                    strict_missing_issue(
                        Some(tool),
                        "PIPELINE_STEP_PARSE_ERROR",
                        e,
                        Some(pipeline_step_path.to_string_lossy().to_string()),
                    )
                } else {
                    warn_missing_issue(
                        Some(tool),
                        "PIPELINE_STEP_PARSE_ERROR",
                        "pipeline_step.json exists but could not be parsed".to_string(),
                        Some(pipeline_step_path.to_string_lossy().to_string()),
                    )
                });
                (Some(pipeline_step_path.clone()), None)
            }
        }
    } else {
        (None, None)
    };

    Ok(Some(ToolContracts {
        name: tool.to_string(),
        dir: tool_dir.to_path_buf(),
        summary_path,
        summary_json,
        pipeline_step_path: pipeline_step_path_opt,
        primary_metrics_path,
    }))
}

fn read_json_value(path: &Path) -> Result<Value, String> {
    let raw = fs::read_to_string(path)
        .map_err(|e| format!("failed reading {}: {e}", path.to_string_lossy()))?;
    serde_json::from_str(&raw)
        .map_err(|e| format!("failed parsing {} as json: {e}", path.to_string_lossy()))
}
