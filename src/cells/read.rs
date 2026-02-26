use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::contracts::types::Issue;
use crate::model::organelle::OrganelleId;
use crate::util::tsv::TsvReader;
use crate::warn_missing_issue;

use super::types::{CellKey, OrganelleCellState};

#[derive(Debug, Clone)]
pub struct ToolCellRow {
    pub id: String,
    pub organelle_state: OrganelleCellState,
}

#[derive(Debug, Clone)]
pub struct ToolCellData {
    pub tool: String,
    pub organelle: OrganelleId,
    pub cell_key: CellKey,
    pub axes_union: Vec<String>,
    pub rows: Vec<ToolCellRow>,
}

pub fn read_tool_primary_metrics(
    input_root: &Path,
    tool: &str,
    organelle: OrganelleId,
    primary_metrics_path: Option<&str>,
    issues: &mut Vec<Issue>,
) -> Option<ToolCellData> {
    let Some(primary_metrics_path) = primary_metrics_path else {
        issues.push(warn_missing_issue(
            Some(tool),
            "MISSING_PRIMARY_METRICS",
            "pipeline_step.json does not declare artifacts.primary_metrics".to_string(),
            None,
        ));
        return None;
    };

    let path = resolve_primary_metrics_path(input_root, tool, primary_metrics_path);
    let mut reader = match TsvReader::open(&path) {
        Ok(v) => v,
        Err(e) => {
            issues.push(warn_missing_issue(
                Some(tool),
                "PRIMARY_METRICS_READ_ERROR",
                format!("failed to open primary metrics TSV: {e}"),
                Some(path.to_string_lossy().to_string()),
            ));
            return None;
        }
    };

    let mut header_ranges = Vec::new();
    match reader.read_record(&mut header_ranges) {
        Ok(true) => {}
        Ok(false) => {
            issues.push(warn_missing_issue(
                Some(tool),
                "EMPTY_PRIMARY_METRICS",
                "primary metrics TSV is empty".to_string(),
                Some(path.to_string_lossy().to_string()),
            ));
            return None;
        }
        Err(e) => {
            issues.push(warn_missing_issue(
                Some(tool),
                "PRIMARY_METRICS_READ_ERROR",
                format!("failed reading header: {e}"),
                Some(path.to_string_lossy().to_string()),
            ));
            return None;
        }
    }

    let header_len = TsvReader::fields_len(&header_ranges);
    let mut headers = Vec::with_capacity(header_len);
    for idx in 0..header_len {
        headers.push(
            reader
                .field(&header_ranges, idx)
                .unwrap_or_default()
                .to_string(),
        );
    }
    let header_index = build_header_index(&headers);

    let cell_key = if header_index.contains_key("barcode") {
        CellKey::Barcode
    } else if header_index.contains_key("sample") {
        CellKey::Sample
    } else {
        issues.push(warn_missing_issue(
            Some(tool),
            "MISSING_ID_COLUMN",
            "primary metrics TSV has neither 'barcode' nor 'sample' column".to_string(),
            Some(path.to_string_lossy().to_string()),
        ));
        return None;
    };

    let id_idx = match cell_key {
        CellKey::Barcode => *header_index.get("barcode").unwrap_or(&0),
        CellKey::Sample => *header_index.get("sample").unwrap_or(&0),
    };

    let regime_idx = header_index.get("regime").copied();
    let confidence_idx = header_index.get("confidence").copied();
    let flags_idx = header_index
        .get("flags")
        .copied()
        .or_else(|| header_index.get("flag").copied());

    let mut axis_columns = Vec::with_capacity(header_len);
    for (idx, col) in headers.iter().enumerate() {
        let lower = col.to_ascii_lowercase();
        if idx == id_idx
            || regime_idx == Some(idx)
            || confidence_idx == Some(idx)
            || flags_idx == Some(idx)
            || lower.is_empty()
        {
            continue;
        }
        axis_columns.push((idx, col.clone()));
    }

    let mut axis_union = BTreeMap::new();
    let mut rows = Vec::new();
    let mut fields = Vec::with_capacity(header_len);
    let mut line_no = 1usize;

    loop {
        line_no += 1;
        let read_ok = match reader.read_record(&mut fields) {
            Ok(v) => v,
            Err(e) => {
                issues.push(warn_missing_issue(
                    Some(tool),
                    "PRIMARY_METRICS_ROW_READ_ERROR",
                    format!("failed reading row {}: {e}", line_no),
                    Some(path.to_string_lossy().to_string()),
                ));
                continue;
            }
        };
        if !read_ok {
            break;
        }

        if row_is_empty(&reader, &fields) {
            continue;
        }

        let id = reader
            .field(&fields, id_idx)
            .map(str::trim)
            .unwrap_or_default();
        if id.is_empty() {
            issues.push(warn_missing_issue(
                Some(tool),
                "EMPTY_CELL_ID",
                format!("row {} has empty id", line_no),
                Some(path.to_string_lossy().to_string()),
            ));
            continue;
        }

        let mut axes = BTreeMap::new();
        for (idx, axis_name) in &axis_columns {
            let raw = reader
                .field(&fields, *idx)
                .map(str::trim)
                .unwrap_or_default();
            if raw.is_empty() {
                continue;
            }
            match raw.parse::<f64>() {
                Ok(v) => {
                    axes.insert(axis_name.clone(), v);
                    axis_union.insert(axis_name.clone(), ());
                }
                Err(_) => {
                    issues.push(warn_missing_issue(
                        Some(tool),
                        "INVALID_NUMERIC_AXIS_VALUE",
                        format!(
                            "row {} column '{}' is not a valid number",
                            line_no, axis_name
                        ),
                        Some(path.to_string_lossy().to_string()),
                    ));
                }
            }
        }

        let regime = regime_idx
            .and_then(|idx| reader.field(&fields, idx))
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(ToOwned::to_owned);

        let confidence = confidence_idx
            .and_then(|idx| reader.field(&fields, idx))
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .and_then(|s| match s.parse::<f64>() {
                Ok(v) => Some(v),
                Err(_) => {
                    issues.push(warn_missing_issue(
                        Some(tool),
                        "INVALID_CONFIDENCE",
                        format!("row {} confidence is not a valid number", line_no),
                        Some(path.to_string_lossy().to_string()),
                    ));
                    None
                }
            });

        let flags = flags_idx
            .and_then(|idx| reader.field(&fields, idx))
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|raw| {
                raw.split([';', ','])
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(ToOwned::to_owned)
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();

        rows.push(ToolCellRow {
            id: id.to_string(),
            organelle_state: OrganelleCellState {
                axes,
                regime,
                confidence,
                flags,
            },
        });
    }

    Some(ToolCellData {
        tool: tool.to_string(),
        organelle,
        cell_key,
        axes_union: axis_union.into_keys().collect(),
        rows,
    })
}

fn resolve_primary_metrics_path(input_root: &Path, tool: &str, raw: &str) -> PathBuf {
    let raw_path = Path::new(raw);
    if raw_path.is_absolute() {
        return raw_path.to_path_buf();
    }
    input_root.join(tool).join(raw_path)
}

fn build_header_index(headers: &[String]) -> BTreeMap<String, usize> {
    let mut index = BTreeMap::new();
    for (i, h) in headers.iter().enumerate() {
        index.insert(h.to_ascii_lowercase(), i);
    }
    index
}

fn row_is_empty(reader: &TsvReader, fields: &[std::ops::Range<usize>]) -> bool {
    for i in 0..fields.len() {
        if !reader
            .field(fields, i)
            .unwrap_or_default()
            .trim()
            .is_empty()
        {
            return false;
        }
    }
    true
}
