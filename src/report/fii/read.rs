use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::report::fii::binning::{BinSpec, Heatmap2D, Histogram1D, clamp01};
use crate::util::tsv::TsvReader;

const FII_BINS: usize = 50;
const COUPLING_BINS: usize = 50;

#[derive(Debug, Clone)]
pub struct SampleInput {
    pub root: PathBuf,
    pub label: String,
    pub order: usize,
    pub timepoint: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SampleSummary {
    pub mean: f64,
    pub median: f64,
    pub p25: f64,
    pub p75: f64,
    pub adaptive: f64,
    pub transition: f64,
    pub resistant: f64,
    pub low_confidence_fraction: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct HeatmapPayload {
    pub counts: Vec<u64>,
    pub mean_x: Vec<f64>,
    pub mean_y: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CouplingPayload {
    pub component: String,
    pub bins_x: usize,
    pub bins_y: usize,
    pub all: HeatmapPayload,
    pub adaptive: HeatmapPayload,
    pub transition: HeatmapPayload,
    pub resistant: HeatmapPayload,
    pub all_high_conf: HeatmapPayload,
    pub adaptive_high_conf: HeatmapPayload,
    pub transition_high_conf: HeatmapPayload,
    pub resistant_high_conf: HeatmapPayload,
}

#[derive(Debug, Clone, Serialize)]
pub struct SampleReportData {
    pub id: String,
    pub label: String,
    pub path: String,
    pub order: usize,
    pub timepoint: Option<String>,
    pub n_cells: u64,
    pub fii_histogram: Vec<u64>,
    pub fii_bins: usize,
    pub fii_min: f64,
    pub fii_max: f64,
    pub fii_bin_width: f64,
    pub summary: SampleSummary,
    pub coupling: Vec<CouplingPayload>,
}

#[derive(Debug, Deserialize)]
struct ManifestJsonItem {
    path: String,
    label: Option<String>,
    order: Option<usize>,
    timepoint: Option<String>,
}

#[derive(Debug, Clone)]
struct RowValues {
    fii: f64,
    mito: Option<f64>,
    translation: Option<f64>,
    splice: Option<f64>,
    regime: String,
    low_confidence: bool,
}

pub fn resolve_inputs(
    inputs: &[PathBuf],
    manifest: Option<&Path>,
) -> Result<Vec<SampleInput>, String> {
    if let Some(manifest_path) = manifest {
        return read_manifest(manifest_path);
    }

    let mut out = Vec::new();
    for (idx, root) in inputs.iter().enumerate() {
        let label = root
            .file_name()
            .and_then(|v| v.to_str())
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| root.display().to_string());
        out.push(SampleInput {
            root: root.clone(),
            label,
            order: idx,
            timepoint: None,
        });
    }
    Ok(out)
}

fn read_manifest(path: &Path) -> Result<Vec<SampleInput>, String> {
    let raw = std::fs::read_to_string(path).map_err(|e| {
        format!(
            "MISSING or READ failure for manifest {}: {e}",
            path.display()
        )
    })?;
    let ext = path
        .extension()
        .and_then(|v| v.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();

    if ext == "json" {
        let parsed: Vec<ManifestJsonItem> =
            serde_json::from_str(&raw).map_err(|e| format!("PARSE manifest JSON failed: {e}"))?;
        let mut out = parsed
            .into_iter()
            .enumerate()
            .map(|(idx, it)| SampleInput {
                root: PathBuf::from(&it.path),
                label: it.label.unwrap_or(it.path.clone()),
                order: it.order.unwrap_or(idx),
                timepoint: it.timepoint,
            })
            .collect::<Vec<_>>();
        out.sort_by_key(|s| s.order);
        return Ok(out);
    }

    let delim = if ext == "csv" { ',' } else { '\t' };
    let mut lines = raw.lines();
    let header = lines
        .next()
        .ok_or_else(|| "MISSING manifest header row".to_string())?;
    let cols = header.split(delim).map(str::trim).collect::<Vec<_>>();
    let mut idx_path = None;
    let mut idx_label = None;
    let mut idx_order = None;
    let mut idx_time = None;
    for (idx, col) in cols.iter().enumerate() {
        match col.to_ascii_lowercase().as_str() {
            "path" | "sample" | "dir" | "sample_dir" => idx_path = Some(idx),
            "label" | "sample_label" => idx_label = Some(idx),
            "order" | "order_rank" => idx_order = Some(idx),
            "timepoint" => idx_time = Some(idx),
            _ => {}
        }
    }
    let idx_path = idx_path.ok_or_else(|| "MISSING manifest 'path' column".to_string())?;
    let mut out = Vec::new();
    for (line_idx, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split(delim).map(str::trim).collect::<Vec<_>>();
        let path_val = fields
            .get(idx_path)
            .copied()
            .ok_or_else(|| format!("INVALID manifest row {} path", line_idx + 2))?;
        let label = idx_label
            .and_then(|i| fields.get(i).copied())
            .filter(|v| !v.is_empty())
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| path_val.to_string());
        let order = idx_order
            .and_then(|i| fields.get(i).copied())
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(line_idx);
        let timepoint = idx_time
            .and_then(|i| fields.get(i).copied())
            .filter(|v| !v.is_empty())
            .map(ToOwned::to_owned);
        out.push(SampleInput {
            root: PathBuf::from(path_val),
            label,
            order,
            timepoint,
        });
    }
    out.sort_by_key(|s| s.order);
    Ok(out)
}

pub fn read_sample(sample: &SampleInput) -> Result<SampleReportData, String> {
    let tsv_path = sample.root.join("functional_irreversibility_index.tsv");
    if !tsv_path.is_file() {
        return Err(format!(
            "MISSING functional_irreversibility_index.tsv in {}",
            sample.root.display()
        ));
    }

    let state_path = if sample.root.join("summary.json").is_file() {
        sample.root.join("summary.json")
    } else {
        sample.root.join("state.json")
    };
    if !state_path.is_file() {
        return Err(format!(
            "MISSING summary.json/state.json in {}",
            sample.root.display()
        ));
    }

    let rows = read_fii_rows(&tsv_path)?;
    let summary = read_summary(&state_path).unwrap_or_else(|| summary_from_rows(&rows));

    let mut hist = Histogram1D::new(BinSpec {
        bins: FII_BINS,
        min: 0.0,
        max: 1.0,
    });

    let mut all_maps = init_coupling_maps();
    for row in &rows {
        hist.add(row.fii);
        update_coupling_maps(&mut all_maps, row);
    }

    let coupling = all_maps
        .into_iter()
        .map(|(name, maps)| build_coupling_payload(name, maps))
        .collect::<Vec<_>>();

    Ok(SampleReportData {
        id: sanitize_id(&sample.label),
        label: sample.label.clone(),
        path: sample.root.display().to_string(),
        order: sample.order,
        timepoint: sample.timepoint.clone(),
        n_cells: rows.len() as u64,
        fii_histogram: hist.counts,
        fii_bins: FII_BINS,
        fii_min: 0.0,
        fii_max: 1.0,
        fii_bin_width: 1.0 / FII_BINS as f64,
        summary,
        coupling,
    })
}

fn read_fii_rows(path: &Path) -> Result<Vec<RowValues>, String> {
    let mut reader =
        TsvReader::open(path).map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let mut header = Vec::new();
    let has_header = reader
        .read_record(&mut header)
        .map_err(|e| format!("READ header failed for {}: {e}", path.display()))?;
    if !has_header {
        return Err(format!("MISSING header in {}", path.display()));
    }

    let cols = (0..TsvReader::fields_len(&header))
        .map(|i| {
            reader
                .field(&header, i)
                .unwrap_or_default()
                .trim()
                .to_string()
        })
        .collect::<Vec<_>>();
    let idx = col_map(&cols);

    let required = [
        "cell_id",
        "functional_irreversibility_index",
        "fii_regime",
        "fii_low_confidence",
    ];
    for name in required {
        if !idx.contains_key(name) {
            return Err(format!(
                "MISSING required column '{}' in {}",
                name,
                path.display()
            ));
        }
    }

    let mut rows = Vec::new();
    let mut fields = Vec::new();
    while reader
        .read_record(&mut fields)
        .map_err(|e| format!("READ row failed for {}: {e}", path.display()))?
    {
        if fields.is_empty() {
            continue;
        }
        let fii = parse_req_float(
            reader
                .field(
                    &fields,
                    *idx.get("functional_irreversibility_index")
                        .expect("checked"),
                )
                .unwrap_or_default(),
            "functional_irreversibility_index",
            path,
        )?;
        let regime = reader
            .field(&fields, *idx.get("fii_regime").expect("checked"))
            .unwrap_or_default()
            .trim()
            .to_string();
        let low_confidence = parse_bool(
            reader
                .field(&fields, *idx.get("fii_low_confidence").expect("checked"))
                .unwrap_or_default(),
        );
        let mito = idx
            .get("mitochondrial_stress_adaptation_score")
            .and_then(|i| reader.field(&fields, *i))
            .and_then(parse_opt_float);
        let translation = idx
            .get("translation_commitment_score")
            .and_then(|i| reader.field(&fields, *i))
            .and_then(parse_opt_float);
        let splice = idx
            .get("splice_irreversibility_index")
            .and_then(|i| reader.field(&fields, *i))
            .and_then(parse_opt_float);
        rows.push(RowValues {
            fii: clamp01(fii),
            mito,
            translation,
            splice,
            regime,
            low_confidence,
        });
    }
    Ok(rows)
}

fn parse_req_float(raw: &str, col: &str, path: &Path) -> Result<f64, String> {
    raw.trim()
        .parse::<f64>()
        .map_err(|_| format!("PARSE invalid '{}' in {}", col, path.display()))
}

fn parse_opt_float(raw: &str) -> Option<f64> {
    let r = raw.trim();
    if r.is_empty() {
        None
    } else {
        r.parse::<f64>().ok().map(clamp01)
    }
}

fn parse_bool(raw: &str) -> bool {
    matches!(
        raw.trim().to_ascii_lowercase().as_str(),
        "true" | "1" | "yes"
    )
}

fn read_summary(path: &Path) -> Option<SampleSummary> {
    let raw = std::fs::read_to_string(path).ok()?;
    let v: Value = serde_json::from_str(&raw).ok()?;
    let fi = v.get("functional_irreversibility")?;
    let dist = fi.get("distribution")?;
    let regimes = fi.get("regime_fractions")?;
    Some(SampleSummary {
        mean: dist.get("mean")?.as_f64()?,
        median: dist.get("median")?.as_f64()?,
        p25: dist.get("p25")?.as_f64()?,
        p75: dist.get("p75")?.as_f64()?,
        adaptive: regimes.get("adaptive")?.as_f64()?,
        transition: regimes.get("transition")?.as_f64()?,
        resistant: regimes.get("resistant")?.as_f64()?,
        low_confidence_fraction: fi.get("low_confidence_fraction")?.as_f64()?,
    })
}

fn summary_from_rows(rows: &[RowValues]) -> SampleSummary {
    let mut values = rows.iter().map(|r| r.fii).collect::<Vec<_>>();
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len() as f64;
    let mean = if n > 0.0 {
        values.iter().sum::<f64>() / n
    } else {
        0.0
    };
    let p = |q: f64| -> f64 {
        if values.is_empty() {
            0.0
        } else {
            let idx = ((values.len() - 1) as f64 * q).round() as usize;
            values[idx]
        }
    };
    let mut regimes = BTreeMap::<&str, u64>::new();
    let mut low = 0u64;
    for row in rows {
        let key = match row.regime.as_str() {
            "Adaptive" => "adaptive",
            "Transition" => "transition",
            "Resistant" => "resistant",
            _ => "transition",
        };
        *regimes.entry(key).or_insert(0) += 1;
        if row.low_confidence {
            low += 1;
        }
    }

    SampleSummary {
        mean,
        median: p(0.5),
        p25: p(0.25),
        p75: p(0.75),
        adaptive: fraction(*regimes.get("adaptive").unwrap_or(&0), rows.len()),
        transition: fraction(*regimes.get("transition").unwrap_or(&0), rows.len()),
        resistant: fraction(*regimes.get("resistant").unwrap_or(&0), rows.len()),
        low_confidence_fraction: fraction(low, rows.len()),
    }
}

fn fraction(count: u64, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64
    }
}

fn col_map(cols: &[String]) -> BTreeMap<String, usize> {
    let mut map = BTreeMap::new();
    for (idx, col) in cols.iter().enumerate() {
        map.insert(col.to_string(), idx);
    }
    map
}

fn sanitize_id(label: &str) -> String {
    label
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect()
}

struct ComponentMaps {
    all: Heatmap2D,
    adaptive: Heatmap2D,
    transition: Heatmap2D,
    resistant: Heatmap2D,
    all_high_conf: Heatmap2D,
    adaptive_high_conf: Heatmap2D,
    transition_high_conf: Heatmap2D,
    resistant_high_conf: Heatmap2D,
}

fn init_coupling_maps() -> BTreeMap<String, ComponentMaps> {
    let mut out = BTreeMap::new();
    let x = BinSpec {
        bins: COUPLING_BINS,
        min: 0.0,
        max: 1.0,
    };
    let y = BinSpec {
        bins: COUPLING_BINS,
        min: 0.0,
        max: 1.0,
    };
    for name in ["mitochondrial", "translation", "splice"] {
        let maps = ComponentMaps {
            all: Heatmap2D::new(x, y),
            adaptive: Heatmap2D::new(x, y),
            transition: Heatmap2D::new(x, y),
            resistant: Heatmap2D::new(x, y),
            all_high_conf: Heatmap2D::new(x, y),
            adaptive_high_conf: Heatmap2D::new(x, y),
            transition_high_conf: Heatmap2D::new(x, y),
            resistant_high_conf: Heatmap2D::new(x, y),
        };
        out.insert(name.to_string(), maps);
    }
    out
}

fn update_coupling_maps(maps: &mut BTreeMap<String, ComponentMaps>, row: &RowValues) {
    for (name, comp_value) in [
        ("mitochondrial", row.mito),
        ("translation", row.translation),
        ("splice", row.splice),
    ] {
        let Some(x) = comp_value else { continue };
        let m = maps.get_mut(name).expect("known component");
        m.all.add(x, row.fii);
        match row.regime.as_str() {
            "Adaptive" => m.adaptive.add(x, row.fii),
            "Transition" => m.transition.add(x, row.fii),
            "Resistant" => m.resistant.add(x, row.fii),
            _ => {}
        }
        if !row.low_confidence {
            m.all_high_conf.add(x, row.fii);
            match row.regime.as_str() {
                "Adaptive" => m.adaptive_high_conf.add(x, row.fii),
                "Transition" => m.transition_high_conf.add(x, row.fii),
                "Resistant" => m.resistant_high_conf.add(x, row.fii),
                _ => {}
            }
        }
    }
}

fn build_coupling_payload(component: String, maps: ComponentMaps) -> CouplingPayload {
    CouplingPayload {
        component,
        bins_x: COUPLING_BINS,
        bins_y: COUPLING_BINS,
        all: to_payload(maps.all),
        adaptive: to_payload(maps.adaptive),
        transition: to_payload(maps.transition),
        resistant: to_payload(maps.resistant),
        all_high_conf: to_payload(maps.all_high_conf),
        adaptive_high_conf: to_payload(maps.adaptive_high_conf),
        transition_high_conf: to_payload(maps.transition_high_conf),
        resistant_high_conf: to_payload(maps.resistant_high_conf),
    }
}

fn to_payload(map: Heatmap2D) -> HeatmapPayload {
    let mut mean_x = vec![0.0; map.counts.len()];
    let mut mean_y = vec![0.0; map.counts.len()];
    for i in 0..map.counts.len() {
        if map.counts[i] > 0 {
            let denom = map.counts[i] as f64;
            mean_x[i] = map.sum_x[i] / denom;
            mean_y[i] = map.sum_y[i] / denom;
        }
    }
    HeatmapPayload {
        counts: map.counts,
        mean_x,
        mean_y,
    }
}
