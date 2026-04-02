use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use crate::cells::types::CellsState;
use crate::model::organelle::OrganelleId;
use crate::normalize::global_robust::GlobalNormalizationSummary;
use crate::registry::output_schema::{OutputSchemaAvailability, metrics_tsv_columns};
use crate::run::experiments::sanitize_name;
use crate::state::AggregatedState;
use crate::systems::composites::{
    SystemsMetricRow, SystemsNormalizedRow, SystemsStageProfile, SystemsStateModel,
    compute_systems_state_model_profiled,
};
use crate::systems::cross_sample::{
    CrossSampleAnalysis, compute_cross_sample_analysis, load_sample_systems_input,
    render_cross_sample_dsp_tsv,
};
use crate::systems::formal_model::build_systems_model_export;
use crate::version;

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
    export_systems_model: Option<&Path>,
    profile_out: Option<&mut SystemsStageProfile>,
) -> Result<(), String> {
    let entity = derive_entity_id_from_path(&state.inputs.a);
    let row = IntegrationRow {
        entity_id: entity.clone(),
        timepoint: 0,
        axes: extract_canonical_axes(state),
        source: None,
    };
    let expression = aggregate_expression_proxy(cells, &entity);
    let (systems_state, systems_profile) = compute_systems_state_model_profiled(cells);
    if let Some(out) = profile_out {
        *out = systems_profile;
    }
    if let (Some(export_path), Some(systems)) = (export_systems_model, systems_state.as_ref()) {
        let systems_export = build_systems_model_export(systems);
        let value = serde_json::to_value(systems_export).unwrap_or(Value::Null);
        crate::io::write_json_atomic(export_path, &value)
            .map_err(|e| format!("failed writing {}: {e}", export_path.display()))?;
    }
    write_integration_bundle(
        &[row],
        &expression,
        out_dir,
        systems_state.as_ref(),
        Some(cells),
        state.global_normalization.as_ref(),
        None,
    )
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
    let mut sample_inputs = Vec::with_capacity(ordered.len());

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
        if let Ok(sample_input) = load_sample_systems_input(exp_name, exp_out) {
            sample_inputs.push(sample_input);
        }
    }

    let cross_sample = if sample_inputs.is_empty() {
        None
    } else {
        Some(compute_cross_sample_analysis(&sample_inputs))
    };
    if let Some(analysis) = cross_sample.as_ref() {
        let integration_dir = out_dir.join("integration");
        std::fs::create_dir_all(&integration_dir)
            .map_err(|e| format!("failed creating {}: {e}", integration_dir.display()))?;
        let dsp_path = integration_dir.join("cross_sample_dsp.tsv");
        let tsv = render_cross_sample_dsp_tsv(analysis);
        crate::io::write_bytes_atomic(&dsp_path, tsv.as_bytes())
            .map_err(|e| format!("failed writing {}: {e}", dsp_path.display()))?;
    }

    write_integration_bundle(
        &rows,
        &expression,
        out_dir,
        None,
        None,
        None,
        cross_sample.as_ref(),
    )
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
    systems_state: Option<&SystemsStateModel>,
    cells: Option<&CellsState>,
    global_normalization: Option<&GlobalNormalizationSummary>,
    cross_sample_analysis: Option<&CrossSampleAnalysis>,
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

    if let Some(systems) = systems_state {
        let metrics_path = integration_dir.join("metrics.tsv");
        let metrics_tsv = render_systems_metrics_tsv(&systems.metrics_rows);
        crate::io::write_bytes_atomic(&metrics_path, metrics_tsv.as_bytes())
            .map_err(|e| format!("failed writing {}: {e}", metrics_path.display()))?;

        let normalized_path = integration_dir.join("metrics_normalized.tsv");
        let normalized_tsv = render_systems_normalized_tsv(&systems.normalized_rows);
        crate::io::write_bytes_atomic(&normalized_path, normalized_tsv.as_bytes())
            .map_err(|e| format!("failed writing {}: {e}", normalized_path.display()))?;

        let systems_export = build_systems_model_export(systems);
        let systems_model_value = serde_json::to_value(&systems_export).unwrap_or(Value::Null);
        let systems_model_path = integration_dir.join("systems_model.json");
        crate::io::write_json_atomic(&systems_model_path, &systems_model_value)
            .map_err(|e| format!("failed writing {}: {e}", systems_model_path.display()))?;

        let landscape_path = integration_dir.join("landscape.html");
        let landscape_html = build_landscape_html(systems, &systems_export);
        crate::io::write_bytes_atomic(&landscape_path, landscape_html.as_bytes())
            .map_err(|e| format!("failed writing {}: {e}", landscape_path.display()))?;
    }

    let summary_path = integration_dir.join("summary.json");
    let summary = build_summary_json(
        rows,
        !expression_rows.is_empty(),
        systems_state,
        cells,
        global_normalization,
        cross_sample_analysis,
    );
    crate::io::write_json_atomic(&summary_path, &summary)
        .map_err(|e| format!("failed writing {}: {e}", summary_path.display()))?;

    Ok(())
}

fn build_landscape_html(
    systems: &SystemsStateModel,
    systems_export: &crate::systems::formal_model::SystemsModelExport,
) -> String {
    let mut points = Vec::with_capacity(systems.metrics_rows.len());
    for row in &systems.metrics_rows {
        points.push(json!({
            "id": row.cell_id,
            "x": row.potential,
            "y": row.stability_gradient,
            "regime": row.regime_class,
            "transition": row.transition_candidate,
            "cpi": row.cpi
        }));
    }

    let graph_nodes = systems_export
        .state_graph
        .nodes
        .iter()
        .map(|n| json!({"id": n.id, "kind": n.kind}))
        .collect::<Vec<_>>();
    let graph_edges = systems_export
        .state_graph
        .edges
        .iter()
        .map(|e| json!({"from": e.from, "to": e.to, "weight": e.weight}))
        .collect::<Vec<_>>();

    let matrix_regimes = systems_export
        .state_transition_matrix
        .regimes
        .iter()
        .cloned()
        .collect::<Vec<_>>();
    let matrix_values = systems_export
        .state_transition_matrix
        .matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|v| canonicalize_json_f64(*v))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut membership = BTreeMap::<&str, (&str, bool, bool)>::new();
    for row in &systems.metrics_rows {
        membership.insert(
            row.cell_id.as_str(),
            (
                row.regime_class.as_str(),
                row.transition_candidate,
                row.rare_state,
            ),
        );
    }
    let groups = [
        "AllCells",
        "TransitionCandidates",
        "SystemicCollapseRisk",
        "RareStates",
    ];
    let mut metric_acc = BTreeMap::<&str, [f64; 4]>::new();
    let mut metric_cnt = BTreeMap::<&str, [usize; 4]>::new();
    for row in &systems.normalized_rows {
        if !row.z_global.is_finite() {
            continue;
        }
        let Some((regime, transition, rare)) = membership.get(row.cell_id.as_str()).copied() else {
            continue;
        };
        let sum = metric_acc.entry(row.metric.as_str()).or_insert([0.0; 4]);
        let cnt = metric_cnt.entry(row.metric.as_str()).or_insert([0usize; 4]);
        sum[0] += row.z_global;
        cnt[0] += 1;
        if transition {
            sum[1] += row.z_global;
            cnt[1] += 1;
        }
        if regime == "SystemicCollapseRisk" {
            sum[2] += row.z_global;
            cnt[2] += 1;
        }
        if rare {
            sum[3] += row.z_global;
            cnt[3] += 1;
        }
    }
    let mut metric_rank = metric_acc
        .iter()
        .map(|(metric, sum)| {
            let cnt = metric_cnt.get(metric).copied().unwrap_or([0usize; 4]);
            let means = [
                if cnt[0] > 0 {
                    sum[0] / cnt[0] as f64
                } else {
                    0.0
                },
                if cnt[1] > 0 {
                    sum[1] / cnt[1] as f64
                } else {
                    0.0
                },
                if cnt[2] > 0 {
                    sum[2] / cnt[2] as f64
                } else {
                    0.0
                },
                if cnt[3] > 0 {
                    sum[3] / cnt[3] as f64
                } else {
                    0.0
                },
            ];
            let rank_score = means[1].abs() + means[2].abs() + means[3].abs();
            ((*metric).to_string(), means, rank_score)
        })
        .collect::<Vec<_>>();
    metric_rank.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.0.cmp(&b.0)));
    let top = metric_rank.into_iter().take(16).collect::<Vec<_>>();
    let heat_metrics = top.iter().map(|x| x.0.clone()).collect::<Vec<_>>();
    let heat_values = top
        .iter()
        .map(|x| {
            x.1.iter()
                .map(|v| canonicalize_json_f64(*v))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let payload = json!({
        "points": points,
        "graph": {
            "nodes": graph_nodes,
            "edges": graph_edges
        },
        "transition_matrix": {
            "regimes": matrix_regimes,
            "values": matrix_values
        },
        "normalized_heatmap": {
            "groups": groups,
            "metrics": heat_metrics,
            "values": heat_values
        }
    });
    let payload_json = serde_json::to_string(&payload).unwrap_or_else(|_| "{}".to_string());

    let template = r###"<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>kira-organelle systems report</title><style>
:root{
 --ink:#12202f; --muted:#5b6a79; --line:#d9e1ea; --panel:#ffffff; --canvas:#f8fafc;
 --accent:#0e7490; --accent-2:#155e75; --warn:#b91c1c;
 --reg-base:#4c78a8; --reg-immune:#72b7b2; --reg-adapt:#f58518; --reg-failed:#e45756; --reg-collapse:#b279a2;
 --math-font:"STIX Two Math","Cambria Math","Latin Modern Math","Times New Roman",serif;
}
*{box-sizing:border-box}
body{
 margin:0;padding:24px;
 font-family:"Avenir Next","Segoe UI",Roboto,Helvetica,Arial,sans-serif;
 color:var(--ink);
 background:radial-gradient(1200px 700px at 100% -10%, #d8ecff 0%, #edf4fb 40%, #f4f7fb 68%, #f8fafc 100%);
}
.wrap{max-width:1260px;margin:0 auto;}
.hero{
 border:1px solid var(--line); border-radius:18px; padding:18px 20px; background:linear-gradient(145deg,#ffffff 0%,#f4f8fc 100%);
 box-shadow:0 20px 55px rgba(15,23,42,0.08);
}
h1{margin:0 0 8px; font-size:33px; letter-spacing:-0.02em;}
h2{margin:0 0 10px; font-size:25px; letter-spacing:-0.01em;}
h3{margin:0; font-size:14px; text-transform:uppercase; color:var(--muted); letter-spacing:0.06em;}
p,li{line-height:1.55; margin:0;}
small{color:var(--muted);}
code{background:#eef4fa;padding:2px 6px;border-radius:6px}
.lead{color:#223447;font-size:15px}
.cards{display:grid;gap:12px;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));margin-top:14px}
.card{background:#ffffff;border:1px solid var(--line);border-radius:12px;padding:10px 12px;box-shadow:0 6px 18px rgba(15,23,42,0.06)}
.card .v{font-size:24px;font-weight:700;margin-top:4px}
.tabs{position:sticky;top:8px;z-index:20;background:rgba(248,250,252,0.9);backdrop-filter:blur(6px);padding:10px 0 2px;display:flex;flex-wrap:wrap;gap:8px}
.tab{border:1px solid #b7c9da;border-radius:999px;background:#f4f8fc;color:#214159;padding:8px 12px;font-weight:600;cursor:pointer}
.tab.active{background:linear-gradient(180deg,#0e7490,#155e75);color:#fff;border-color:#0f5f75}
.panel{
 margin-top:14px; border:1px solid var(--line); border-radius:16px; background:var(--panel);
 box-shadow:0 15px 38px rgba(15,23,42,0.08); padding:16px;
}
.panel.hidden{display:none}
.formula{
 display:inline-block; font-family:var(--math-font); font-size:1.08em; font-style:italic;
 background:linear-gradient(180deg,#f9fbff 0%, #edf4fb 100%); border:1px solid #ced8e4; border-radius:8px; padding:4px 8px; margin:0 2px;
}
.split{display:grid;grid-template-columns:2fr 1fr;gap:12px}
.controls{display:flex;align-items:center;gap:10px;flex-wrap:wrap;color:#284155;font-size:13px;margin:8px 0 6px}
input[type="range"]{accent-color:var(--accent)}
svg,canvas{width:100%;display:block;border:1px solid var(--line);border-radius:10px;background:var(--canvas)}
.legend{display:flex;gap:10px;flex-wrap:wrap;font-size:12px;margin:10px 0 0}
.chip{display:inline-flex;align-items:center;gap:6px;padding:4px 7px;border:1px solid #d4deea;border-radius:999px;background:#f8fbff}
.dot{width:10px;height:10px;border-radius:10px;display:inline-block}
table{border-collapse:collapse;width:100%;font-size:12px;margin-top:8px}
th,td{border:1px solid var(--line);padding:5px 7px;text-align:right}
th:first-child,td:first-child{text-align:left}
.note{margin-top:8px;font-size:12px;color:#3f5468}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
@media (max-width:980px){.split{grid-template-columns:1fr} body{padding:14px} h1{font-size:28px}}
</style></head><body><div class="wrap">
<div class="hero">
<h1>Systems Landscape And Transition Geometry</h1>
<p class="lead">Integrated deterministic report from <code>integration/metrics.tsv</code>, <code>integration/metrics_normalized.tsv</code>, and <code>integration/systems_model.json</code>. The report is designed for visual triage first, then mathematical inspection.</p>
<div class="cards" id="summaryCards"></div>
</div>
<div class="tabs" id="tabs">
 <button class="tab active" data-target="panel1">1) Landscape</button>
 <button class="tab" data-target="panel2">2) State Graph</button>
 <button class="tab" data-target="panel3">3) Regime Transition Matrix</button>
 <button class="tab" data-target="panel4">4) Normalized Signals</button>
</div>
<section id="panel1" class="panel">
 <h2>1) Landscape Scatter</h2>
 <p>Each point is one cell in latent system coordinates: <span class="formula">x = Potential</span>, <span class="formula">y = StabilityGradient</span>. Radius maps to collapse pressure via <span class="formula">r = 2 + 6·σ(CPI)</span>, where <span class="formula">σ(z) = 1 / (1 + e<sup>-z</sup>)</span>. Red ring denotes <code>TransitionCandidate = true</code>.</p>
 <canvas id="scatter" width="1120" height="560"></canvas>
 <div class="legend" id="scatterLegend"></div>
 <p class="note">Interpretation shortcut: right-shift means higher drive (Potential), up-shift means stronger local stability barrier (StabilityGradient).</p>
</section>
<section id="panel2" class="panel hidden">
 <h2>2) State Graph</h2>
 <p>Two-level view: (A) clear regime-to-regime flow (always visible), (B) optional basin layer for fine structure. Directed edge width scales with transition weight.</p>
 <div class="controls">
  <label>Basin-edge percentile:
   <input id="edgePercentile" type="range" min="35" max="90" value="60" step="1">
  </label>
  <label><input id="showBasins" type="checkbox" checked> Show basin layer</label>
  <span id="edgePercentileLabel" class="mono"></span>
 </div>
 <div class="split">
  <svg id="graph" viewBox="0 0 1120 700"></svg>
  <div>
   <h3>Top Directed Edges</h3>
   <table id="edgeTable"></table>
   <p class="note">If this panel still looks dense, increase threshold to focus on strongest deterministic transitions only.</p>
  </div>
 </div>
 <div class="legend" id="graphLegend"></div>
</section>
<section id="panel3" class="panel hidden">
 <h2>3) Regime Transition Matrix</h2>
 <p>Transition probability is estimated by local neighborhood dynamics: <span class="formula">P<sub>ij</sub> = N(i → j) / Σ<sub>j</sub> N(i → j)</span>.</p>
 <canvas id="matrix" width="1120" height="430"></canvas>
</section>
<section id="panel4" class="panel hidden">
 <h2>4) Normalized Signal Heatmap</h2>
 <p>Top normalized metrics (z-space) across subsets <span class="mono">AllCells</span>, <span class="mono">TransitionCandidates</span>, <span class="mono">SystemicCollapseRisk</span>, <span class="mono">RareStates</span>. Robust normalization form: <span class="formula">z = (x - median) / MAD</span>.</p>
 <canvas id="heat" width="1120" height="560"></canvas>
 <table id="heatTable"></table>
</section>
<script>const DATA=__PAYLOAD__;
const COLORS={BaselineCompensated:getComputedStyle(document.documentElement).getPropertyValue('--reg-base').trim(),ImmuneSuppressiveEmbedded:getComputedStyle(document.documentElement).getPropertyValue('--reg-immune').trim(),AdaptiveHighStress:getComputedStyle(document.documentElement).getPropertyValue('--reg-adapt').trim(),HighStress_FailedCompensation:getComputedStyle(document.documentElement).getPropertyValue('--reg-failed').trim(),SystemicCollapseRisk:getComputedStyle(document.documentElement).getPropertyValue('--reg-collapse').trim()};
const REGIME_ORDER=['BaselineCompensated','ImmuneSuppressiveEmbedded','AdaptiveHighStress','HighStress_FailedCompensation','SystemicCollapseRisk'];
const REGIME_X=new Map(REGIME_ORDER.map((r,i)=>[r,i]));
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const scale=(v,vmin,vmax,a,b)=>vmax<=vmin?(a+b)/2:a+(v-vmin)*(b-a)/(vmax-vmin);
const qtile=(arr,q)=>{if(!arr.length)return 0;const a=arr.slice().sort((x,y)=>x-y);return a[Math.floor((a.length-1)*q)];};
function initTabs(){
 const tabs=[...document.querySelectorAll('.tab')];
 tabs.forEach(t=>t.addEventListener('click',()=>{
   tabs.forEach(x=>x.classList.remove('active'));
   document.querySelectorAll('.panel').forEach(p=>p.classList.add('hidden'));
   t.classList.add('active');
   document.getElementById(t.dataset.target).classList.remove('hidden');
 }));
}
function summaryCards(){
 const pts=DATA.points||[];
 if(!pts.length) return;
 const n=pts.length;
 const transitions=pts.filter(p=>p.transition).length;
 const collapse=pts.filter(p=>p.regime==='SystemicCollapseRisk').length;
 const xs=pts.map(p=>p.x).sort((a,b)=>a-b), ys=pts.map(p=>p.y).sort((a,b)=>a-b);
 const med=(a)=>a.length?a[Math.floor(a.length/2)]:0;
 const cards=[
  ['Cells',String(n)],
  ['Transition Candidates',`${transitions} (${(100*transitions/n).toFixed(1)}%)`],
  ['Systemic Collapse Risk',`${collapse} (${(100*collapse/n).toFixed(1)}%)`],
  ['Median Potential',med(xs).toFixed(3)],
  ['Median StabilityGradient',med(ys).toFixed(3)]
 ];
 const root=document.getElementById('summaryCards');
 root.innerHTML=cards.map(([k,v])=>`<div class="card"><small>${k}</small><div class="v">${v}</div></div>`).join('');
}
function drawScatter(){
 const cv=document.getElementById('scatter'); const ctx=cv.getContext('2d');
 const pts=DATA.points||[]; if(!pts.length) return;
 const W=cv.width,H=cv.height,pad=56;
 const xs=pts.map(p=>p.x), ys=pts.map(p=>p.y);
 const xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys), ymax=Math.max(...ys);
 const step=6;
 ctx.clearRect(0,0,W,H);
 for(let x=pad;x<=W-pad;x+=Math.max(60,(W-2*pad)/step)){ctx.strokeStyle='#e5edf5';ctx.beginPath();ctx.moveTo(x,pad);ctx.lineTo(x,H-pad);ctx.stroke();}
 for(let y=pad;y<=H-pad;y+=Math.max(50,(H-2*pad)/step)){ctx.strokeStyle='#e5edf5';ctx.beginPath();ctx.moveTo(pad,y);ctx.lineTo(W-pad,y);ctx.stroke();}
 ctx.strokeStyle='#5a6c80'; ctx.lineWidth=1.2; ctx.strokeRect(pad,pad,W-2*pad,H-2*pad);
 for(const p of pts){
  const x=scale(p.x,xmin,xmax,pad,W-pad);
  const y=scale(p.y,ymin,ymax,H-pad,pad);
  const r=2+6*(1/(1+Math.exp(-p.cpi)));
  ctx.beginPath(); ctx.arc(x,y,r,0,Math.PI*2); ctx.fillStyle=COLORS[p.regime]||'#8b97a6'; ctx.globalAlpha=0.58; ctx.fill();
  if(p.transition){ ctx.beginPath(); ctx.arc(x,y,r+1.2,0,Math.PI*2); ctx.strokeStyle='#be123c'; ctx.globalAlpha=0.86; ctx.lineWidth=1.3; ctx.stroke(); }
 }
 ctx.globalAlpha=1; ctx.fillStyle='#1c2d40'; ctx.font='13px "Avenir Next",sans-serif';
 ctx.fillText('Potential',W/2-25,H-18);
 ctx.save(); ctx.translate(20,H/2+38); ctx.rotate(-Math.PI/2); ctx.fillText('StabilityGradient',0,0); ctx.restore();
 ctx.fillStyle='#3e556b'; ctx.font='11px ui-monospace, SFMono-Regular, Menlo, monospace';
 ctx.fillText(`x:[${xmin.toFixed(2)}, ${xmax.toFixed(2)}]`,pad+4,H-pad+18);
 ctx.fillText(`y:[${ymin.toFixed(2)}, ${ymax.toFixed(2)}]`,pad+170,H-pad+18);
 const lg=document.getElementById('scatterLegend'); lg.innerHTML='';
 for(const [k,v] of Object.entries(COLORS)){const d=document.createElement('span'); d.className='chip'; d.innerHTML=`<span class="dot" style="background:${v}"></span>${k}`; lg.appendChild(d);}
 const t=document.createElement('span'); t.className='chip'; t.innerHTML='<span class="dot" style="background:#fff;border:2px solid #be123c"></span>TransitionCandidate'; lg.appendChild(t);
}
function drawGraph(){
 const svg=document.getElementById('graph');
 const nodes=(DATA.graph&&DATA.graph.nodes)||[];
 const edges=(DATA.graph&&DATA.graph.edges)||[];
 const mat=(DATA.transition_matrix)||{};
 const regimes=mat.regimes||REGIME_ORDER.slice();
 const tvals=mat.values||[];
 while(svg.firstChild) svg.removeChild(svg.firstChild);
 if(!nodes.length && !regimes.length) return;
 const W=1120,H=700;
 const showBasins=document.getElementById('showBasins').checked;
 const edgePct=Number(document.getElementById('edgePercentile').value||60)/100;
 const thr=qtile(edges.map(e=>e.weight),edgePct);
 let drawEdges=edges.filter(e=>e.weight>=thr);
 const minEdges=Math.max(12, regimes.length*3);
 if(drawEdges.length<minEdges){
  drawEdges=edges.slice().sort((a,b)=>b.weight-a.weight).slice(0,minEdges);
 }
 const effThr=drawEdges.length?Math.min(...drawEdges.map(e=>e.weight)):thr;
 document.getElementById('edgePercentileLabel').textContent=`q=${edgePct.toFixed(2)} threshold=${thr.toFixed(3)} effective_min=${effThr.toFixed(3)} edges=${drawEdges.length}`;
 const nodesById=new Map(nodes.map(n=>[n.id,n]));
 const regimeNodes=(nodes.filter(n=>n.kind==='RegimeClass').length?nodes.filter(n=>n.kind==='RegimeClass'):regimes.map(id=>({id,kind:'RegimeClass'})));
 const basinNodes=nodes.filter(n=>n.kind!=='RegimeClass');
  const pos=new Map();
 const xMin=120, xMax=W-120;
 const regimeY=170;
 for(const rn of regimeNodes){
  const idx=REGIME_X.has(rn.id)?REGIME_X.get(rn.id):Math.floor(REGIME_ORDER.length/2);
  const x=scale(idx,0,REGIME_ORDER.length-1,xMin,xMax);
  pos.set(rn.id,{x,y:regimeY});
 }
 const influence=new Map(); // weighted by all edges for stable positioning
 for(const b of basinNodes){ influence.set(b.id,new Array(REGIME_ORDER.length).fill(0)); }
 for(const e of edges){
  const a=nodesById.get(e.from), b=nodesById.get(e.to);
  if(!a||!b) continue;
  if(a.kind!=='RegimeClass'&&b.kind==='RegimeClass'&&REGIME_X.has(b.id)){influence.get(a.id)[REGIME_X.get(b.id)] += e.weight;}
  if(b.kind!=='RegimeClass'&&a.kind==='RegimeClass'&&REGIME_X.has(a.id)){influence.get(b.id)[REGIME_X.get(a.id)] += e.weight;}
 }
 const lanes=Array.from({length:REGIME_ORDER.length},()=>[]);
 for(const b of basinNodes){
  const s=influence.get(b.id) || [];
  let best=0; for(let i=1;i<s.length;i++) if((s[i]||0)>(s[best]||0)) best=i;
  lanes[best].push(b.id);
 }
 for(let lane=0;lane<lanes.length;lane++){
  const ids=lanes[lane];
  ids.sort();
  for(let i=0;i<ids.length;i++){
   const count=ids.length;
   const center=H*0.56;
   const spread=Math.min(320,Math.max(120,count*12));
   const y=center+(count===1?0:(i/(count-1)-0.5)*spread);
   const hash=[...ids[i]].reduce((a,c)=>((a*131 + c.charCodeAt(0))>>>0),7);
   const jitter=((hash%9)-4)*7;
   const x=scale(lane,0,REGIME_ORDER.length-1,xMin+24,xMax-24)+jitter;
   pos.set(ids[i],{x,y});
  }
 }
 const defs=document.createElementNS('http://www.w3.org/2000/svg','defs');
 const marker=document.createElementNS('http://www.w3.org/2000/svg','marker');
 marker.setAttribute('id','arrow'); marker.setAttribute('viewBox','0 0 10 10'); marker.setAttribute('refX','8'); marker.setAttribute('refY','5');
 marker.setAttribute('markerWidth','6'); marker.setAttribute('markerHeight','6'); marker.setAttribute('orient','auto-start-reverse');
 const path=document.createElementNS('http://www.w3.org/2000/svg','path'); path.setAttribute('d','M0,0 L10,5 L0,10 z'); path.setAttribute('fill','#6b7d91');
 marker.appendChild(path); defs.appendChild(marker); svg.appendChild(defs);
 const regimeFlow=[];
 for(let i=0;i<regimes.length;i++){
  for(let j=0;j<regimes.length;j++){
   const v=(tvals[i]&&tvals[i][j])||0;
   if(i===j || v<=0) continue;
   regimeFlow.push({from:regimes[i],to:regimes[j],weight:v});
  }
 }
 const regimeThr=qtile(regimeFlow.map(e=>e.weight),0.35);
 const drawRegimeEdges=regimeFlow.filter(e=>e.weight>=regimeThr);
 const maxRw=Math.max(0.0001,...drawRegimeEdges.map(e=>e.weight),0.0001);
 for(const e of drawRegimeEdges){
  const a=pos.get(e.from), b=pos.get(e.to); if(!a||!b) continue;
  const color=COLORS[e.from]||'#546a80';
  const p=document.createElementNS('http://www.w3.org/2000/svg','path');
  const dx=b.x-a.x, bend=clamp(Math.abs(dx)*0.14,16,52);
  const c1x=a.x+dx*0.33, c2x=a.x+dx*0.66;
  p.setAttribute('d',`M ${a.x} ${a.y+12} C ${c1x} ${a.y+bend}, ${c2x} ${b.y+bend}, ${b.x} ${b.y+12}`);
  p.setAttribute('fill','none');
  p.setAttribute('stroke',color);
  p.setAttribute('stroke-opacity',String(clamp(0.25 + e.weight/maxRw,0.3,0.92)));
  p.setAttribute('stroke-width',String(1.8 + 7.5*e.weight/maxRw));
  p.setAttribute('marker-end','url(#arrow)');
  svg.appendChild(p);
 }
 const maxW=Math.max(0.0001,...drawEdges.map(e=>e.weight),0.0001);
 if(showBasins){
 for(const e of drawEdges){
  const a=pos.get(e.from), b=pos.get(e.to); if(!a||!b) continue;
  const src=nodesById.get(e.from), dst=nodesById.get(e.to);
  const color=(src&&src.kind==='RegimeClass'&&COLORS[src.id])?COLORS[src.id]:((dst&&dst.kind==='RegimeClass'&&COLORS[dst.id])?COLORS[dst.id]:'#6b7d91');
  const p=document.createElementNS('http://www.w3.org/2000/svg','path');
  const dx=b.x-a.x, bend=clamp(Math.abs(dx)*0.22,14,70)*(a.y<b.y?1:-1);
  const c1x=a.x+dx*0.33, c2x=a.x+dx*0.66;
  p.setAttribute('d',`M ${a.x} ${a.y} C ${c1x} ${a.y+bend}, ${c2x} ${b.y-bend}, ${b.x} ${b.y}`);
  p.setAttribute('fill','none');
  p.setAttribute('stroke',color);
  p.setAttribute('stroke-opacity',String(clamp(0.2 + e.weight/maxW,0.26,0.9)));
  p.setAttribute('stroke-width',String(0.8 + 4.2*e.weight/maxW));
  p.setAttribute('marker-end','url(#arrow)');
  svg.appendChild(p);
 }
 }
 for(const n of nodes){
  const p0=pos.get(n.id); if(!p0) continue;
  const isReg=n.kind==='RegimeClass';
  if(!showBasins && !isReg) continue;
  const c=document.createElementNS('http://www.w3.org/2000/svg','circle');
  c.setAttribute('cx',String(p0.x)); c.setAttribute('cy',String(p0.y));
  c.setAttribute('r',isReg?'13':'4.8');
  c.setAttribute('fill',isReg?(COLORS[n.id]||'#334155'):'#7f8ea0');
  c.setAttribute('stroke',isReg?'#fff':'#e7edf5'); c.setAttribute('stroke-width',isReg?'2':'1');
  svg.appendChild(c);
  if(isReg){
    const t=document.createElementNS('http://www.w3.org/2000/svg','text');
    t.setAttribute('x',String(p0.x)); t.setAttribute('y',String(p0.y-20));
    t.setAttribute('text-anchor','middle'); t.setAttribute('font-size','12');
    t.setAttribute('font-family','Avenir Next, Segoe UI, sans-serif'); t.setAttribute('fill','#1f3347');
    t.textContent=n.id; svg.appendChild(t);
  }
 }
 const edgeTable=document.getElementById('edgeTable');
 const topReg=drawRegimeEdges.slice().sort((a,b)=>b.weight-a.weight).slice(0,8);
 const topBas=drawEdges.slice().sort((a,b)=>b.weight-a.weight).slice(0,8);
 const rows=topReg.map(e=>`<tr><td>${e.from}</td><td>${e.to}</td><td>${Number(e.weight).toFixed(3)}</td></tr>`).join('');
 const rows2=showBasins?topBas.map(e=>`<tr><td>${e.from}</td><td>${e.to}</td><td>${Number(e.weight).toFixed(3)}</td></tr>`).join(''):'';
 edgeTable.innerHTML='<thead><tr><th>from</th><th>to</th><th>w</th></tr></thead><tbody>'+rows+(showBasins?rows2:'')+'</tbody>';
 const gl=document.getElementById('graphLegend');
 gl.innerHTML='';
 for(const r of REGIME_ORDER){
  if(!nodesById.has(r)) continue;
  const el=document.createElement('span'); el.className='chip'; el.innerHTML=`<span class="dot" style="background:${COLORS[r]}"></span>${r}`; gl.appendChild(el);
 }
 const basin=document.createElement('span'); basin.className='chip'; basin.innerHTML='<span class="dot" style="background:#7f8ea0"></span>Basin node';
 gl.appendChild(basin);
}
function drawMatrix(){
 const cv=document.getElementById('matrix'); const ctx=cv.getContext('2d'); const mat=DATA.transition_matrix; if(!mat) return;
 const regs=mat.regimes||[]; const vals=mat.values||[]; const n=regs.length; const left=220, top=50, cell=56;
 ctx.clearRect(0,0,cv.width,cv.height);
 let vmax=0; for(let i=0;i<n;i++) for(let j=0;j<n;j++) vmax=Math.max(vmax,Math.abs((vals[i]&&vals[i][j])||0)); if(vmax===0)vmax=1;
 for(let i=0;i<n;i++){
  for(let j=0;j<n;j++){
   const v=(vals[i]&&vals[i][j])||0;
   const t=clamp(v/vmax,0,1);
   const r=Math.round(232-120*t), g=Math.round(241-130*t), b=Math.round(248-20*t);
   ctx.fillStyle=`rgb(${r},${g},${b})`;
   ctx.fillRect(left+j*cell,top+i*cell,cell,cell);
   ctx.strokeStyle='#d8e2ee'; ctx.strokeRect(left+j*cell,top+i*cell,cell,cell);
   ctx.fillStyle=t>0.58?'#f8fbff':'#13263a'; ctx.font='11px ui-monospace, SFMono-Regular, Menlo, monospace';
   ctx.fillText(v.toFixed(2),left+j*cell+10,top+i*cell+31);
  }
 }
 ctx.fillStyle='#21384f'; ctx.font='12px "Avenir Next",sans-serif';
 for(let i=0;i<n;i++){
  ctx.fillText(regs[i],16,top+i*cell+32);
  ctx.save(); ctx.translate(left+i*cell+32,35); ctx.rotate(-0.68); ctx.fillText(regs[i],0,0); ctx.restore();
 }
 ctx.fillText('from-state',16,20); ctx.fillText('to-state',left,20);
}
function drawHeat(){
 const cv=document.getElementById('heat'); const ctx=cv.getContext('2d'); const h=DATA.normalized_heatmap; if(!h) return;
 const rows=h.metrics||[]; const cols=h.groups||[]; const vals=h.values||[];
 const left=250, top=36, cw=195, ch=28;
 ctx.clearRect(0,0,cv.width,cv.height);
 let vmax=0; for(const r of vals) for(const v of r) vmax=Math.max(vmax,Math.abs(v||0)); if(vmax===0) vmax=1;
 for(let i=0;i<rows.length;i++){
  for(let j=0;j<cols.length;j++){
   const v=(vals[i]&&vals[i][j])||0;
   const t=clamp((v+vmax)/(2*vmax),0,1);
   const rr=Math.round(44+197*t), gg=Math.round(80+150*(1-Math.abs(t-0.5)*2)), bb=Math.round(214-174*t);
   ctx.fillStyle=`rgb(${rr},${gg},${bb})`;
   ctx.fillRect(left+j*cw,top+i*ch,cw,ch);
   ctx.strokeStyle='#e0e8f1'; ctx.strokeRect(left+j*cw,top+i*ch,cw,ch);
   ctx.fillStyle=Math.abs(v)>vmax*0.55?'#ffffff':'#13263a';
   ctx.font='11px ui-monospace, SFMono-Regular, Menlo, monospace';
   ctx.fillText(v.toFixed(2),left+j*cw+8,top+i*ch+18);
  }
  ctx.fillStyle='#233a52'; ctx.font='12px "Avenir Next",sans-serif'; ctx.fillText(rows[i],10,top+i*ch+18);
 }
 for(let j=0;j<cols.length;j++){ctx.fillStyle='#1f3246';ctx.font='12px "Avenir Next",sans-serif';ctx.fillText(cols[j],left+j*cw+8,20);}
 const tb=document.getElementById('heatTable');
 tb.innerHTML='<thead><tr><th>Metric</th>'+cols.map(c=>`<th>${c}</th>`).join('')+'</tr></thead><tbody>' +
  rows.map((r,i)=>'<tr><td>'+r+'</td>'+cols.map((_,j)=>'<td>'+(((vals[i]||[])[j]||0).toFixed(3))+'</td>').join('')+'</tr>').join('') +
  '</tbody>';
}
initTabs();
summaryCards();
drawScatter();
drawMatrix();
drawHeat();
drawGraph();
document.getElementById('edgePercentile').addEventListener('input',drawGraph);
document.getElementById('showBasins').addEventListener('change',drawGraph);
</script></div></body></html>"###;
    template.replace("__PAYLOAD__", &payload_json)
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
            out.push_str(&fmt6(value));
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
        out.push_str(&fmt6(row.value));
        out.push('\n');
    }
    out
}

fn render_systems_metrics_tsv(rows: &[SystemsMetricRow]) -> String {
    let include_fragility = rows.iter().any(|r| r.dominant_axis != "NA");
    let availability = OutputSchemaAvailability {
        include_aggregate: true,
        include_fragility,
        include_landscape: true,
        include_therapeutic: true,
        include_dynamic: true,
    };
    let cols = metrics_tsv_columns(availability);
    let mut out = String::new();
    out.push_str(&cols.join("\t"));
    out.push('\n');
    for row in rows {
        out.push_str(&row.cell_id);
        out.push('\t');
        out.push_str(&row.cluster);
        out.push('\t');
        out.push_str(&fmt6(row.stress_vector));
        out.push('\t');
        out.push_str(&fmt6(row.compensation_deficit));
        out.push('\t');
        out.push_str(&fmt6(row.cpi));
        out.push('\t');
        out.push_str(&fmt6(row.afs));
        out.push('\t');
        out.push_str(&fmt6(row.imsc));
        out.push('\t');
        out.push_str(&row.regime_class);
        if availability.include_aggregate {
            out.push('\t');
            out.push_str(if row.rare_state { "true" } else { "false" });
            out.push('\t');
            if row.mahalanobis_distance.is_finite() {
                out.push_str(&fmt6(row.mahalanobis_distance));
            } else {
                out.push_str("NaN");
            }
        }
        if availability.include_fragility {
            out.push('\t');
            out.push_str(&row.dominant_axis);
            out.push('\t');
            out.push_str(&fmt6(row.dominant_axis_sensitivity));
        }
        if availability.include_landscape {
            out.push('\t');
            out.push_str(&fmt6(row.potential));
            out.push('\t');
            out.push_str(&fmt6(row.stability_gradient));
            out.push('\t');
            out.push_str(&fmt6(row.tpi_landscape));
            out.push('\t');
            out.push_str(&row.basin_id);
            out.push('\t');
            out.push_str(if row.transition_candidate {
                "true"
            } else {
                "false"
            });
        }
        if availability.include_therapeutic {
            out.push('\t');
            out.push_str(&row.dominant_vulnerable_axis);
        }
        if availability.include_dynamic {
            out.push('\t');
            out.push_str(&fmt6(row.lsi));
            out.push('\t');
            out.push_str(&fmt6(row.max_trajectory));
            out.push('\t');
            out.push_str(&fmt6(row.bee));
        }
        out.push('\n');
    }
    out
}

fn render_systems_normalized_tsv(rows: &[SystemsNormalizedRow]) -> String {
    let mut out = String::new();
    out.push_str("cell_id\tcluster\tmetric\traw_value\tz_global\n");
    for row in rows {
        out.push_str(&row.cell_id);
        out.push('\t');
        out.push_str(&row.cluster);
        out.push('\t');
        out.push_str(&row.metric);
        out.push('\t');
        out.push_str(&fmt6(row.raw_value));
        out.push('\t');
        out.push_str(&fmt6(row.z_global));
        out.push('\n');
    }
    out
}

fn build_summary_json(
    rows: &[IntegrationRow],
    expression_available: bool,
    systems_state: Option<&SystemsStateModel>,
    cells: Option<&CellsState>,
    global_normalization: Option<&GlobalNormalizationSummary>,
    cross_sample_analysis: Option<&CrossSampleAnalysis>,
) -> Value {
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
        "model".to_string(),
        json!({
            "engine": version::ENGINE_NAME,
            "model_version": version::MODEL_VERSION,
            "normalization": version::NORMALIZATION_VERSION,
            "composites_version": version::COMPOSITES_VERSION,
            "fragility_version": version::FRAGILITY_VERSION,
        }),
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
    if let Some(systems) = systems_state {
        root.insert(
            "systems_state_model".to_string(),
            systems_state_model_to_json(systems),
        );
    }
    if let Some(cells) = cells {
        let diag = cells
            .ingestion_diagnostics
            .iter()
            .map(|d| {
                json!({
                    "tool": d.tool,
                    "source_path": d.source_path,
                    "recognized_columns": d.recognized_columns,
                    "missing_expected_columns": d.missing_expected_columns,
                    "unknown_columns": d.unknown_columns,
                })
            })
            .collect::<Vec<_>>();
        root.insert("ingestion_diagnostics".to_string(), Value::Array(diag));
    }
    if let Some(global_normalization) = global_normalization {
        root.insert(
            "global_normalization".to_string(),
            serde_json::to_value(global_normalization).unwrap_or(Value::Null),
        );
    }
    if let Some(cross_sample_analysis) = cross_sample_analysis {
        root.insert(
            "cross_sample_analysis".to_string(),
            serde_json::to_value(cross_sample_analysis).unwrap_or(Value::Null),
        );
    }
    Value::Object(root)
}

fn systems_state_model_to_json(model: &SystemsStateModel) -> Value {
    let mut global_stats = serde_json::Map::new();
    global_stats.insert("n_cells".to_string(), json!(model.metrics_rows.len()));
    global_stats.insert(
        "cpi_p50".to_string(),
        json!(canonicalize_json_f64(model.cpi_p50)),
    );
    global_stats.insert(
        "cpi_p90".to_string(),
        json!(canonicalize_json_f64(model.cpi_p90)),
    );
    global_stats.insert(
        "stress_p50".to_string(),
        json!(canonicalize_json_f64(model.stress_p50)),
    );
    global_stats.insert(
        "stress_p90".to_string(),
        json!(canonicalize_json_f64(model.stress_p90)),
    );
    global_stats.insert(
        "system_entropy".to_string(),
        json!(canonicalize_json_f64(model.system_entropy)),
    );
    global_stats.insert(
        "normalized_entropy".to_string(),
        json!(canonicalize_json_f64(model.normalized_entropy)),
    );
    global_stats.insert(
        "rare_fraction".to_string(),
        json!(canonicalize_json_f64(model.rare_fraction_global)),
    );

    let mut cluster_stats = serde_json::Map::new();
    let mut sorted = model.cluster_stats.clone();
    sorted.sort_by(|a, b| cluster_key_cmp(&a.cluster, &b.cluster));
    for c in sorted {
        let median_cpi = c.metrics.get("CPI").map(|m| m.median).unwrap_or(0.0);
        let median_stress = c
            .metrics
            .get("StressVector")
            .map(|m| m.median)
            .unwrap_or(0.0);
        cluster_stats.insert(
            c.cluster.clone(),
            json!({
                "n_cells": c.n_cells,
                "median_cpi": canonicalize_json_f64(median_cpi),
                "median_stress": canonicalize_json_f64(median_stress),
                "heterogeneity_index": canonicalize_json_f64(c.heterogeneity_index),
                "tail_fraction_cpi": canonicalize_json_f64(c.tail_fraction_cpi),
                "tail_fraction_stress": canonicalize_json_f64(c.tail_fraction_stress),
                "rare_fraction": canonicalize_json_f64(c.rare_fraction),
                "regime_fraction": c.regime_fraction,
                "metrics": c.metrics,
            }),
        );
    }

    let mut frag_cluster = serde_json::Map::new();
    let mut frag_clusters = model
        .fragility
        .cluster_axis_ranking
        .iter()
        .collect::<Vec<_>>();
    frag_clusters.sort_by(|(a, _), (b, _)| cluster_key_cmp(a, b));
    for (cluster_id, ranking) in frag_clusters {
        frag_cluster.insert(cluster_id.clone(), json!(ranking));
    }

    json!({
        "global_stats": global_stats,
        "cluster_stats": cluster_stats,
        "regime_distribution": {
            "counts": model.global_regime_counts,
            "fractions": model.global_regime_fraction,
        },
        "global_cpi_p50": canonicalize_json_f64(model.cpi_p50),
        "global_cpi_p90": canonicalize_json_f64(model.cpi_p90),
        "global_stress_p50": canonicalize_json_f64(model.stress_p50),
        "global_stress_p90": canonicalize_json_f64(model.stress_p90),
        "rare_fraction_global": canonicalize_json_f64(model.rare_fraction_global),
        "fragility": {
            "delta": model.fragility.delta,
            "global_axis_ranking": model.fragility.global_axis_ranking,
            "cluster_axis_ranking": frag_cluster,
        },
        "therapeutic_projection": model.therapeutic_projection,
        "dynamic_stability": model.dynamic_stability,
        "validation": model.validation,
        "compatibility": {
            "missing_metrics": model.compatibility_missing_metrics,
            "degraded_axes": model.compatibility_degraded_axes,
        },
        "landscape": model.landscape,
        "missing_metric_counts": model.missing_metric_counts,
    })
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

fn fmt6(v: f64) -> String {
    if !v.is_finite() {
        return "NaN".to_string();
    }
    format!("{:.6}", canonicalize_json_f64(v))
}

fn canonicalize_json_f64(v: f64) -> f64 {
    crate::io::canonicalize_f64(v)
}

fn cluster_key_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    match (a.parse::<i64>(), b.parse::<i64>()) {
        (Ok(ia), Ok(ib)) => ia.cmp(&ib).then(a.cmp(b)),
        _ => a.cmp(b),
    }
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
