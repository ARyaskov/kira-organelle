pub mod read;
pub mod types;

use std::path::Path;

use crate::{strict_missing_issue, warn_missing_issue};
use read::read_tool_contracts;
use types::{Issue, ToolContracts};

pub const EXPECTED_TOOLS: [&str; 8] = [
    "kira-mitoqc",
    "kira-nuclearqc",
    "kira-spliceqc",
    "kira-proteoqc",
    "kira-autolys",
    "kira-secretion",
    "kira-energetics",
    "kira-microenvironment",
];

const OPTIONAL_TOOLS: [&str; 1] = ["kira-riboqc"];

pub fn discover_and_read(
    input: &Path,
    strict: bool,
    issues: &mut Vec<Issue>,
) -> Result<Vec<ToolContracts>, String> {
    let mut found = Vec::new();

    for tool in EXPECTED_TOOLS {
        let tool_dir = input.join(tool);
        if !tool_dir.is_dir() {
            let issue = if strict {
                strict_missing_issue(
                    Some(tool),
                    "MISSING_TOOL_DIR",
                    format!("missing expected tool directory: {}", tool_dir.display()),
                    Some(tool_dir.to_string_lossy().to_string()),
                )
            } else {
                warn_missing_issue(
                    Some(tool),
                    "MISSING_TOOL_DIR",
                    format!("missing expected tool directory: {}", tool_dir.display()),
                    Some(tool_dir.to_string_lossy().to_string()),
                )
            };
            issues.push(issue);
            continue;
        }

        if let Some(tool_contracts) = read_tool_contracts(tool, &tool_dir, strict, issues)? {
            found.push(tool_contracts);
        }
    }

    // Optional tool contracts are loaded when present but do not emit missing-dir issues.
    for tool in OPTIONAL_TOOLS {
        let tool_dir = input.join(tool);
        if !tool_dir.is_dir() {
            continue;
        }
        if let Some(tool_contracts) = read_tool_contracts(tool, &tool_dir, strict, issues)? {
            found.push(tool_contracts);
        }
    }

    if strict && issues.iter().any(|i| i.severity.is_error()) {
        return Err("strict mode validation failed".to_string());
    }

    Ok(found)
}
