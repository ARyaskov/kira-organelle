use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use crate::cells::types::CellsState;
use crate::model::organelle::OrganelleId;
use crate::run::experiments::sanitize_name;
use crate::state::AggregatedState;

const AXES_HEADER: [&str; 6] = [
    "mito",
    "proteostasis",
    "splice",
    "secretion",
    "energetics",
    "autophagy",
];
const CPI_LAMBDA: f64 = 0.2;

#[derive(Debug, Clone)]
struct IntegrationRow {
    entity_id: String,
    timepoint: usize,
    axes: [f64; 6],
    source: Option<String>,
}

#[derive(Debug, Clone)]
struct ExpressionRow {
    entity_id: String,
    gene: String,
    value: f64,
}

pub fn write_single_integration_artifacts(
    state: &AggregatedState,
    cells: &CellsState,
    out_dir: &Path,
) -> Result<(), String> {
    let entity = derive_entity_id_from_path(&state.inputs.a);
    let row = IntegrationRow {
        entity_id: entity.clone(),
        timepoint: 0,
        axes: extract_canonical_axes(state),
        source: None,
    };
    let expression = aggregate_expression_proxy(cells, &entity);
    write_integration_bundle(&[row], &expression, out_dir)
}

pub fn write_multi_integration_artifacts(
    experiments: &[(String, PathBuf)],
    out_dir: &Path,
    manifest: Option<&Path>,
) -> Result<(), String> {
    if experiments.is_empty() {
        return Ok(());
    }

    let entity = derive_entity_id_from_pathbuf(out_dir);
    let ordered = order_experiments(experiments, manifest)?;
    let mut rows = Vec::with_capacity(ordered.len());
    let mut expression = Vec::new();

    for (idx, exp) in ordered.iter().enumerate() {
        let exp_name = &exp.name;
        let exp_out = &exp.out_dir;
        let state_path = exp_out.join("kira-organelle").join("state.json");
        let state_raw = std::fs::read(&state_path)
            .map_err(|e| format!("failed reading {}: {e}", state_path.display()))?;
        let state: AggregatedState = serde_json::from_slice(&state_raw)
            .map_err(|e| format!("failed parsing {}: {e}", state_path.display()))?;
        rows.push(IntegrationRow {
            entity_id: entity.clone(),
            timepoint: exp.timepoint.unwrap_or(idx),
            axes: extract_canonical_axes(&state),
            source: Some(exp_name.clone()),
        });

        let cells_path = exp_out.join("kira-organelle").join("cells.json");
        let cells_raw = match std::fs::read(&cells_path) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let Ok(cells) = serde_json::from_slice::<CellsState>(&cells_raw) else {
            continue;
        };
        expression.extend(aggregate_expression_proxy(&cells, exp_name));
    }

    write_integration_bundle(&rows, &expression, out_dir)
}

#[derive(Debug, Clone)]
struct OrderedExperiment {
    name: String,
    out_dir: PathBuf,
    order_rank: usize,
    timepoint: Option<usize>,
}

#[derive(Debug, Clone)]
struct ManifestOrder {
    order_rank: usize,
    timepoint: Option<usize>,
}

fn order_experiments(
    experiments: &[(String, PathBuf)],
    manifest: Option<&Path>,
) -> Result<Vec<OrderedExperiment>, String> {
    let Some(manifest_path) = manifest else {
        let out = experiments
            .iter()
            .enumerate()
            .map(|(idx, (name, out_dir))| OrderedExperiment {
                name: name.clone(),
                out_dir: out_dir.clone(),
                order_rank: idx,
                timepoint: None,
            })
            .collect::<Vec<_>>();
        return Ok(out);
    };

    let default_base = usize::MAX / 4;
    let mut out = experiments
        .iter()
        .enumerate()
        .map(|(idx, (name, out_dir))| OrderedExperiment {
            name: name.clone(),
            out_dir: out_dir.clone(),
            order_rank: default_base + idx,
            timepoint: None,
        })
        .collect::<Vec<_>>();

    let order_map = parse_manifest_order(manifest_path)?;
    for item in &mut out {
        if let Some(mapped) = order_map.get(&sanitize_name(&item.name)) {
            item.order_rank = mapped.order_rank;
            item.timepoint = mapped.timepoint;
        }
    }
    out.sort_by(|a, b| a.order_rank.cmp(&b.order_rank).then(a.name.cmp(&b.name)));
    Ok(out)
}

fn parse_manifest_order(path: &Path) -> Result<BTreeMap<String, ManifestOrder>, String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("failed reading manifest {}: {e}", path.display()))?;
    let ext = path
        .extension()
        .and_then(|v| v.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    if ext == "json" {
        return parse_manifest_json(&raw);
    }
    let delim = if ext == "csv" { ',' } else { '\t' };
    parse_manifest_delimited(&raw, delim)
}

fn parse_manifest_json(raw: &str) -> Result<BTreeMap<String, ManifestOrder>, String> {
    let value: Value =
        serde_json::from_str(raw).map_err(|e| format!("failed parsing manifest json: {e}"))?;
    let arr = value
        .as_array()
        .ok_or_else(|| "manifest json must be an array".to_string())?;
    let mut out = BTreeMap::new();
    for (idx, item) in arr.iter().enumerate() {
        let obj = match item.as_object() {
            Some(v) => v,
            None => continue,
        };
        let order_rank = obj
            .get("order")
            .or_else(|| obj.get("order_rank"))
            .and_then(Value::as_u64)
            .map(|v| v as usize)
            .unwrap_or(idx);
        let timepoint = obj
            .get("timepoint")
            .and_then(parse_manifest_timepoint_value)
            .or(Some(order_rank));
        let label = obj
            .get("label")
            .and_then(Value::as_str)
            .map(ToOwned::to_owned)
            .or_else(|| {
                obj.get("path")
                    .and_then(Value::as_str)
                    .and_then(path_tail_as_name)
            });
        let Some(label) = label else { continue };
        out.insert(
            sanitize_name(&label),
            ManifestOrder {
                order_rank,
                timepoint,
            },
        );
    }
    Ok(out)
}

fn parse_manifest_delimited(
    raw: &str,
    delim: char,
) -> Result<BTreeMap<String, ManifestOrder>, String> {
    let mut lines = raw.lines();
    let header = lines
        .next()
        .ok_or_else(|| "manifest is empty; expected header row".to_string())?;
    let cols = header.split(delim).map(str::trim).collect::<Vec<_>>();
    let mut idx_path = None;
    let mut idx_label = None;
    let mut idx_order = None;
    let mut idx_timepoint = None;
    for (idx, col) in cols.iter().enumerate() {
        match col.to_ascii_lowercase().as_str() {
            "path" | "sample" | "dir" | "sample_dir" => idx_path = Some(idx),
            "label" | "sample_label" => idx_label = Some(idx),
            "order" | "order_rank" => idx_order = Some(idx),
            "timepoint" => idx_timepoint = Some(idx),
            _ => {}
        }
    }
    let mut out = BTreeMap::new();
    for (line_idx, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split(delim).map(str::trim).collect::<Vec<_>>();
        let order_rank = idx_order
            .and_then(|i| fields.get(i).copied())
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(line_idx);
        let timepoint = idx_timepoint
            .and_then(|i| fields.get(i).copied())
            .and_then(parse_manifest_timepoint_str)
            .or(Some(order_rank));

        let label = idx_label
            .and_then(|i| fields.get(i).copied())
            .filter(|v| !v.is_empty())
            .map(ToOwned::to_owned)
            .or_else(|| {
                idx_path
                    .and_then(|i| fields.get(i).copied())
                    .filter(|v| !v.is_empty())
                    .and_then(path_tail_as_name)
            });
        let Some(label) = label else { continue };
        out.insert(
            sanitize_name(&label),
            ManifestOrder {
                order_rank,
                timepoint,
            },
        );
    }
    Ok(out)
}

fn parse_manifest_timepoint_value(value: &Value) -> Option<usize> {
    if let Some(v) = value.as_u64() {
        return Some(v as usize);
    }
    if let Some(v) = value.as_str() {
        return parse_manifest_timepoint_str(v);
    }
    None
}

fn parse_manifest_timepoint_str(v: &str) -> Option<usize> {
    v.trim().parse::<usize>().ok()
}

fn path_tail_as_name(path: &str) -> Option<String> {
    let p = Path::new(path);
    p.file_name()
        .and_then(|v| v.to_str())
        .filter(|v| !v.is_empty())
        .map(ToOwned::to_owned)
}

fn write_integration_bundle(
    rows: &[IntegrationRow],
    expression_rows: &[ExpressionRow],
    out_dir: &Path,
) -> Result<(), String> {
    if rows.is_empty() {
        return Ok(());
    }

    let integration_dir = out_dir.join("integration");
    std::fs::create_dir_all(&integration_dir)
        .map_err(|e| format!("failed creating {}: {e}", integration_dir.display()))?;

    let timeseries_path = integration_dir.join("timeseries.tsv");
    let timeseries_tsv = render_timeseries_tsv(rows);
    crate::io::write_bytes_atomic(&timeseries_path, timeseries_tsv.as_bytes())
        .map_err(|e| format!("failed writing {}: {e}", timeseries_path.display()))?;

    let expression_path = integration_dir.join("expression_aggregated.tsv");
    let expression_tsv = render_expression_tsv(expression_rows);
    crate::io::write_bytes_atomic(&expression_path, expression_tsv.as_bytes())
        .map_err(|e| format!("failed writing {}: {e}", expression_path.display()))?;

    let summary_path = integration_dir.join("summary.json");
    let summary = build_summary_json(rows, !expression_rows.is_empty());
    crate::io::write_json_atomic(&summary_path, &summary)
        .map_err(|e| format!("failed writing {}: {e}", summary_path.display()))?;

    Ok(())
}

fn render_timeseries_tsv(rows: &[IntegrationRow]) -> String {
    let mut out = String::new();
    out.push_str("entity_id\ttimepoint");
    for axis in AXES_HEADER {
        out.push('\t');
        out.push_str(axis);
    }
    out.push('\n');

    for row in rows {
        out.push_str(&row.entity_id);
        out.push('\t');
        out.push_str(&row.timepoint.to_string());
        for value in row.axes {
            out.push('\t');
            out.push_str(&format!("{value:.6}"));
        }
        out.push('\n');
    }

    out
}

fn render_expression_tsv(rows: &[ExpressionRow]) -> String {
    let mut out = String::new();
    out.push_str("entity_id\tgene\tvalue\n");
    for row in rows {
        out.push_str(&row.entity_id);
        out.push('\t');
        out.push_str(&row.gene);
        out.push('\t');
        out.push_str(&format!("{:.6}", row.value));
        out.push('\n');
    }
    out
}

fn build_summary_json(rows: &[IntegrationRow], expression_available: bool) -> Value {
    let latest = rows.last().expect("checked by caller");
    let collapse_curve = rows.iter().map(|r| cpi(r.axes)).collect::<Vec<_>>();
    let collapse_peak = collapse_curve
        .iter()
        .copied()
        .reduce(f64::max)
        .unwrap_or(0.0)
        .clamp(0.0, 1.0);

    let mut axes_obj = serde_json::Map::new();
    for (idx, name) in AXES_HEADER.iter().enumerate() {
        axes_obj.insert((*name).to_string(), json!(latest.axes[idx]));
    }

    let mut root = serde_json::Map::new();
    root.insert(
        "schema".to_string(),
        Value::String("kira-organelle-integration-v1".to_string()),
    );
    root.insert(
        "entity".to_string(),
        Value::String(latest.entity_id.clone()),
    );
    root.insert("timepoint".to_string(), json!(latest.timepoint));
    root.insert("axes".to_string(), Value::Object(axes_obj));
    root.insert(
        "collapse_proximity_baseline".to_string(),
        json!(cpi(latest.axes)),
    );
    root.insert("collapse_proximity_peak".to_string(), json!(collapse_peak));
    root.insert(
        "expression_available".to_string(),
        Value::Bool(expression_available),
    );
    root.insert(
        "expression_source".to_string(),
        Value::String("organelle_axes_proxy".to_string()),
    );
    if let Some(source) = &latest.source {
        root.insert(
            "source_experiment".to_string(),
            Value::String(source.clone()),
        );
    }
    Value::Object(root)
}

fn cpi(axes: [f64; 6]) -> f64 {
    let mean_axes = axes.iter().sum::<f64>() / 6.0;
    let entropy = normalized_entropy(&axes);
    (mean_axes + CPI_LAMBDA * entropy).clamp(0.0, 1.0)
}

fn normalized_entropy(axes: &[f64; 6]) -> f64 {
    let sum = axes.iter().sum::<f64>();
    if sum <= 0.0 {
        return 0.0;
    }

    let mut h = 0.0;
    for value in axes {
        if *value <= 0.0 {
            continue;
        }
        let p = *value / sum;
        h -= p * p.ln();
    }
    let max_h = (axes.len() as f64).ln();
    if max_h <= 0.0 {
        return 0.0;
    }
    (h / max_h).clamp(0.0, 1.0)
}

fn aggregate_expression_proxy(cells: &CellsState, entity_id: &str) -> Vec<ExpressionRow> {
    let mut sums = BTreeMap::<String, (f64, usize)>::new();
    for cell in &cells.cells {
        for (org, org_state) in &cell.per_organelle {
            let org_name = organelle_to_str(*org);
            for (axis, value) in &org_state.axes {
                if !value.is_finite() {
                    continue;
                }
                let key = format!("{org_name}::{axis}");
                let entry = sums.entry(key).or_insert((0.0, 0));
                entry.0 += *value;
                entry.1 += 1;
            }
        }
    }

    sums.into_iter()
        .filter_map(|(gene, (sum, n))| {
            if n == 0 {
                return None;
            }
            Some(ExpressionRow {
                entity_id: entity_id.to_string(),
                gene,
                value: (sum / n as f64).clamp(0.0, 1.0),
            })
        })
        .collect()
}

fn extract_canonical_axes(state: &AggregatedState) -> [f64; 6] {
    [
        mean_axis_for(state, OrganelleId::Mitochondria),
        mean_axis_for(state, OrganelleId::Proteostasis),
        mean_axis_for(state, OrganelleId::Spliceosome),
        mean_axis_for(state, OrganelleId::Secretion),
        mean_axis_for(state, OrganelleId::Energetics),
        mean_axis_for(state, OrganelleId::Autophagy),
    ]
}

fn mean_axis_for(state: &AggregatedState, organelle: OrganelleId) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for item in &state.organelle_states {
        if item.organelle != organelle {
            continue;
        }
        for axis in &item.axes {
            if !axis.median.is_finite() {
                continue;
            }
            sum += axis.median;
            n += 1;
        }
    }
    if n == 0 {
        0.0
    } else {
        (sum / n as f64).clamp(0.0, 1.0)
    }
}

fn organelle_to_str(id: OrganelleId) -> &'static str {
    match id {
        OrganelleId::Mitochondria => "mitochondria",
        OrganelleId::Nucleus => "nucleus",
        OrganelleId::Spliceosome => "spliceosome",
        OrganelleId::Ribosome => "ribosome",
        OrganelleId::Proteostasis => "proteostasis",
        OrganelleId::Autophagy => "autophagy",
        OrganelleId::Secretion => "secretion",
        OrganelleId::Energetics => "energetics",
        OrganelleId::Microenvironment => "microenvironment",
    }
}

fn derive_entity_id_from_path(path: &str) -> String {
    derive_entity_id_from_pathbuf(Path::new(path))
}

fn derive_entity_id_from_pathbuf(path: &Path) -> String {
    path.file_name()
        .and_then(|v| v.to_str())
        .filter(|v| !v.is_empty())
        .map(ToOwned::to_owned)
        .unwrap_or_else(|| "sample_0".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn manifest_tsv_order_is_applied() {
        let raw = "path\torder\ttimepoint\nB\t1\t10\nA\t0\t5\n";
        let parsed = parse_manifest_delimited(raw, '\t').expect("manifest parse");
        let a = parsed.get("A").expect("A");
        let b = parsed.get("B").expect("B");
        assert_eq!(a.order_rank, 0);
        assert_eq!(a.timepoint, Some(5));
        assert_eq!(b.order_rank, 1);
        assert_eq!(b.timepoint, Some(10));
    }

    #[test]
    fn order_falls_back_to_input_for_missing_manifest_entries() {
        let dir = tempfile::tempdir().expect("tempdir");
        let manifest = dir.path().join("manifest.tsv");
        std::fs::write(&manifest, "label\torder\nB\t0\n").expect("manifest write");

        let exps = vec![
            ("A".to_string(), dir.path().join("out-a")),
            ("B".to_string(), dir.path().join("out-b")),
        ];
        let ordered = order_experiments(&exps, Some(&manifest)).expect("order");
        assert_eq!(ordered[0].name, "B");
        assert_eq!(ordered[0].timepoint, Some(0));
        assert_eq!(ordered[1].name, "A");
        assert_eq!(ordered[1].timepoint, None);
    }
}
