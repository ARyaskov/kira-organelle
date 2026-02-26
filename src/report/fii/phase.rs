use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::cli::ReportPhaseArgs;
use crate::io;
use crate::util::tsv::TsvReader;

use super::read::{SampleInput, resolve_inputs};

const PROJ_FII_VEL: &str = "fii_vs_velocity";
const PROJ_FII_ACC: &str = "fii_vs_acceleration";
const PROJ_VEL_ACC: &str = "velocity_vs_acceleration";

#[derive(Debug, Clone, Serialize)]
struct PhasePoint {
    sample_label: String,
    order_rank: usize,
    fii_mean: f64,
    fii_median: f64,
    fii_velocity: Option<f64>,
    fii_acceleration: Option<f64>,
    resistant_fraction: f64,
    low_confidence_fraction: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct AxisRange {
    min: f64,
    max: f64,
    bins: usize,
}

#[derive(Debug, Clone)]
struct CellPoint {
    fii: f64,
    regime: String,
    low_conf: bool,
}

#[derive(Debug, Clone, Serialize)]
struct BinRow {
    projection: String,
    sample_label: String,
    bin_x: usize,
    bin_y: usize,
    count: u64,
    frac_adaptive: f64,
    frac_transition: f64,
    frac_resistant: f64,
    frac_low_confidence: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PhaseConfig {
    velocity_threshold: f64,
    resistant_onset_fii: f64,
    acceleration_threshold: f64,
    fii_q: f64,
    vel_q: f64,
    resistant_fraction_jump: f64,
    volatility_abs_velocity_threshold: f64,
}

impl Default for PhaseConfig {
    fn default() -> Self {
        Self {
            velocity_threshold: 0.10,
            resistant_onset_fii: 0.66,
            acceleration_threshold: 0.08,
            fii_q: 0.50,
            vel_q: 0.12,
            resistant_fraction_jump: 0.20,
            volatility_abs_velocity_threshold: 0.10,
        }
    }
}

#[derive(Debug, Clone)]
struct FlagRow {
    sample_label: String,
    order_rank: usize,
    flag: String,
    metric: String,
    value: f64,
    threshold: f64,
    note: String,
}

pub fn run_report_phase(args: &ReportPhaseArgs) -> Result<(), String> {
    if args.inputs.is_empty() && args.manifest.is_none() {
        return Err("INVALID report-phase args: --inputs or --manifest required".to_string());
    }
    if let Some(manifest) = args.manifest.as_deref() {
        validate_manifest_ordering(manifest)?;
    }

    let mut specs = resolve_inputs(&args.inputs, args.manifest.as_deref())?;
    specs.sort_by(|a, b| a.order.cmp(&b.order).then(a.label.cmp(&b.label)));
    validate_sample_ranks(&specs)?;

    let mut per_sample_cells = BTreeMap::<String, Vec<CellPoint>>::new();
    let mut phase_points = Vec::<PhasePoint>::new();
    for spec in &specs {
        let cells = read_fii_cells(spec)?;
        let point = summarize_phase_point(spec, &cells);
        per_sample_cells.insert(spec.label.clone(), cells);
        phase_points.push(point);
    }
    attach_velocity_acceleration(&mut phase_points);

    let vel_range = derive_centered_range(
        phase_points
            .iter()
            .filter_map(|p| p.fii_velocity)
            .collect::<Vec<_>>(),
        40,
    );
    let acc_range = derive_centered_range(
        phase_points
            .iter()
            .filter_map(|p| p.fii_acceleration)
            .collect::<Vec<_>>(),
        40,
    );
    let fii_range = AxisRange {
        min: 0.0,
        max: 1.0,
        bins: 40,
    };

    let bins = build_phase_bins(
        &phase_points,
        &per_sample_cells,
        fii_range,
        vel_range,
        acc_range,
    );
    let cfg = read_config(args.config.as_deref())?;
    let flags = build_flags(&phase_points, &cfg);

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;

    let points_tsv = render_points_tsv(&phase_points);
    io::write_bytes_atomic(
        &args.out.join("phase_portrait_points.tsv"),
        points_tsv.as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for phase_portrait_points.tsv: {e}"))?;

    let bins_tsv = render_bins_tsv(&bins);
    io::write_bytes_atomic(
        &args.out.join("phase_portrait_bins.tsv"),
        bins_tsv.as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for phase_portrait_bins.tsv: {e}"))?;

    let summary_json = build_summary_json(&phase_points, &bins, fii_range, vel_range, acc_range);
    io::write_json_atomic(&args.out.join("phase_portrait_summary.json"), &summary_json)
        .map_err(|e| format!("WRITE failed for phase_portrait_summary.json: {e}"))?;

    let flags_tsv = render_flags_tsv(&flags);
    io::write_bytes_atomic(
        &args.out.join("early_warning_flags.tsv"),
        flags_tsv.as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for early_warning_flags.tsv: {e}"))?;

    Ok(())
}

fn validate_sample_ranks(specs: &[SampleInput]) -> Result<(), String> {
    let mut seen = BTreeSet::new();
    for s in specs {
        if !seen.insert(s.order) {
            return Err(format!(
                "INVALID manifest ordering: duplicated order_rank {}",
                s.order
            ));
        }
    }
    for (idx, s) in specs.iter().enumerate() {
        if s.order != idx {
            return Err(format!(
                "INVALID manifest ordering: missing/non-monotonic rank at {} (got {})",
                idx, s.order
            ));
        }
    }
    Ok(())
}

fn validate_manifest_ordering(path: &Path) -> Result<(), String> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ manifest failed {}: {e}", path.display()))?;
    let ext = path
        .extension()
        .and_then(|v| v.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    let mut seq = Vec::new();
    if ext == "json" {
        let v: Value =
            serde_json::from_str(&raw).map_err(|e| format!("PARSE manifest JSON failed: {e}"))?;
        let arr = v
            .as_array()
            .ok_or_else(|| "INVALID manifest JSON: expected array".to_string())?;
        for item in arr {
            let rank = item
                .get("order_rank")
                .or_else(|| item.get("order"))
                .and_then(Value::as_u64)
                .ok_or_else(|| "INVALID manifest JSON: missing order_rank/order".to_string())?
                as usize;
            seq.push(rank);
        }
    } else {
        let delim = if ext == "csv" { ',' } else { '\t' };
        let mut lines = raw.lines();
        let header = lines
            .next()
            .ok_or_else(|| "MISSING manifest header row".to_string())?;
        let cols = header.split(delim).map(str::trim).collect::<Vec<_>>();
        let rank_idx = cols
            .iter()
            .position(|c| matches!((*c).to_ascii_lowercase().as_str(), "order_rank" | "order"))
            .ok_or_else(|| "MISSING manifest order_rank/order column".to_string())?;
        for line in lines {
            if line.trim().is_empty() {
                continue;
            }
            let fields = line.split(delim).map(str::trim).collect::<Vec<_>>();
            let rank = fields
                .get(rank_idx)
                .ok_or_else(|| "INVALID manifest row: missing order".to_string())?
                .parse::<usize>()
                .map_err(|_| "INVALID manifest row: non-numeric order".to_string())?;
            seq.push(rank);
        }
    }
    if seq.is_empty() {
        return Err("MISSING manifest samples".to_string());
    }
    for i in 1..seq.len() {
        if seq[i] <= seq[i - 1] {
            return Err("INVALID manifest ordering: non-monotonic order_rank".to_string());
        }
    }
    Ok(())
}

fn read_fii_cells(spec: &SampleInput) -> Result<Vec<CellPoint>, String> {
    let path = spec.root.join("functional_irreversibility_index.tsv");
    if !path.is_file() {
        return Err(format!(
            "MISSING functional_irreversibility_index.tsv in {}",
            spec.root.display()
        ));
    }
    let mut reader =
        TsvReader::open(&path).map_err(|e| format!("READ failed for {}: {e}", path.display()))?;
    let mut header = Vec::new();
    if !reader
        .read_record(&mut header)
        .map_err(|e| format!("READ header failed for {}: {e}", path.display()))?
    {
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
    let idx = cols
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i))
        .collect::<BTreeMap<_, _>>();
    let fii_idx = *idx.get("functional_irreversibility_index").ok_or_else(|| {
        format!(
            "MISSING column functional_irreversibility_index in {}",
            path.display()
        )
    })?;
    let regime_idx = *idx
        .get("fii_regime")
        .ok_or_else(|| format!("MISSING column fii_regime in {}", path.display()))?;
    let low_idx = *idx
        .get("fii_low_confidence")
        .ok_or_else(|| format!("MISSING column fii_low_confidence in {}", path.display()))?;

    let mut out = Vec::new();
    let mut fields = Vec::new();
    while reader
        .read_record(&mut fields)
        .map_err(|e| format!("READ row failed for {}: {e}", path.display()))?
    {
        let fii = reader
            .field(&fields, fii_idx)
            .unwrap_or_default()
            .trim()
            .parse::<f64>()
            .ok()
            .and_then(|v| {
                if v.is_finite() {
                    Some(v.clamp(0.0, 1.0))
                } else {
                    None
                }
            });
        let Some(fii) = fii else { continue };
        let regime = reader
            .field(&fields, regime_idx)
            .unwrap_or_default()
            .trim()
            .to_string();
        let low_conf = matches!(
            reader
                .field(&fields, low_idx)
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase()
                .as_str(),
            "true" | "1" | "yes"
        );
        out.push(CellPoint {
            fii,
            regime,
            low_conf,
        });
    }
    Ok(out)
}

fn summarize_phase_point(spec: &SampleInput, cells: &[CellPoint]) -> PhasePoint {
    if cells.is_empty() {
        return PhasePoint {
            sample_label: spec.label.clone(),
            order_rank: spec.order,
            fii_mean: 0.0,
            fii_median: 0.0,
            fii_velocity: None,
            fii_acceleration: None,
            resistant_fraction: 0.0,
            low_confidence_fraction: 1.0,
        };
    }
    let mut fii = cells.iter().map(|c| c.fii).collect::<Vec<_>>();
    fii.sort_by(|a, b| a.total_cmp(b));
    let mean = fii.iter().sum::<f64>() / fii.len() as f64;
    let median = percentile(&fii, 0.5);
    let resistant =
        cells.iter().filter(|c| c.regime == "Resistant").count() as f64 / cells.len() as f64;
    let low = cells.iter().filter(|c| c.low_conf).count() as f64 / cells.len() as f64;
    PhasePoint {
        sample_label: spec.label.clone(),
        order_rank: spec.order,
        fii_mean: mean,
        fii_median: median,
        fii_velocity: None,
        fii_acceleration: None,
        resistant_fraction: resistant,
        low_confidence_fraction: low,
    }
}

fn attach_velocity_acceleration(points: &mut [PhasePoint]) {
    for i in 1..points.len() {
        points[i].fii_velocity = Some(points[i].fii_mean - points[i - 1].fii_mean);
    }
    for i in 2..points.len() {
        let v_i = points[i].fii_velocity.unwrap_or(0.0);
        let v_prev = points[i - 1].fii_velocity.unwrap_or(0.0);
        points[i].fii_acceleration = Some(v_i - v_prev);
    }
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * q).round() as usize;
    sorted[idx]
}

fn derive_centered_range(values: Vec<f64>, bins: usize) -> AxisRange {
    if values.is_empty() {
        return AxisRange {
            min: -0.5,
            max: 0.5,
            bins,
        };
    }
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    for v in values {
        min_v = min_v.min(v);
        max_v = max_v.max(v);
    }
    let span = (max_v - min_v).abs();
    let pad = if span > 0.0 { span * 0.05 } else { 0.05 };
    AxisRange {
        min: min_v - pad,
        max: max_v + pad,
        bins,
    }
}

#[derive(Default, Clone)]
struct BinAgg {
    count: u64,
    adaptive: u64,
    transition: u64,
    resistant: u64,
    low_conf: u64,
}

fn build_phase_bins(
    points: &[PhasePoint],
    per_sample_cells: &BTreeMap<String, Vec<CellPoint>>,
    fii_range: AxisRange,
    vel_range: AxisRange,
    acc_range: AxisRange,
) -> Vec<BinRow> {
    let mut map = BTreeMap::<(String, String, usize, usize), BinAgg>::new();

    for point in points {
        let cells = per_sample_cells
            .get(&point.sample_label)
            .cloned()
            .unwrap_or_default();
        for cell in cells {
            if let Some(v) = point.fii_velocity {
                add_bin(
                    &mut map,
                    PROJ_FII_VEL,
                    &point.sample_label,
                    index_for(fii_range, cell.fii),
                    index_for(vel_range, v),
                    &cell,
                );
                add_bin(
                    &mut map,
                    PROJ_FII_VEL,
                    "ALL",
                    index_for(fii_range, cell.fii),
                    index_for(vel_range, v),
                    &cell,
                );
            }
            if let Some(a) = point.fii_acceleration {
                add_bin(
                    &mut map,
                    PROJ_FII_ACC,
                    &point.sample_label,
                    index_for(fii_range, cell.fii),
                    index_for(acc_range, a),
                    &cell,
                );
                add_bin(
                    &mut map,
                    PROJ_FII_ACC,
                    "ALL",
                    index_for(fii_range, cell.fii),
                    index_for(acc_range, a),
                    &cell,
                );
            }
            if let (Some(v), Some(a)) = (point.fii_velocity, point.fii_acceleration) {
                add_bin(
                    &mut map,
                    PROJ_VEL_ACC,
                    &point.sample_label,
                    index_for(vel_range, v),
                    index_for(acc_range, a),
                    &cell,
                );
                add_bin(
                    &mut map,
                    PROJ_VEL_ACC,
                    "ALL",
                    index_for(vel_range, v),
                    index_for(acc_range, a),
                    &cell,
                );
            }
        }
    }

    map.into_iter()
        .map(|((projection, sample_label, bin_x, bin_y), agg)| BinRow {
            projection,
            sample_label,
            bin_x,
            bin_y,
            count: agg.count,
            frac_adaptive: agg.adaptive as f64 / agg.count as f64,
            frac_transition: agg.transition as f64 / agg.count as f64,
            frac_resistant: agg.resistant as f64 / agg.count as f64,
            frac_low_confidence: agg.low_conf as f64 / agg.count as f64,
        })
        .collect()
}

fn add_bin(
    map: &mut BTreeMap<(String, String, usize, usize), BinAgg>,
    projection: &str,
    sample_label: &str,
    bin_x: usize,
    bin_y: usize,
    cell: &CellPoint,
) {
    let key = (
        projection.to_string(),
        sample_label.to_string(),
        bin_x,
        bin_y,
    );
    let agg = map.entry(key).or_default();
    agg.count += 1;
    match cell.regime.as_str() {
        "Adaptive" => agg.adaptive += 1,
        "Transition" => agg.transition += 1,
        "Resistant" => agg.resistant += 1,
        _ => {}
    }
    if cell.low_conf {
        agg.low_conf += 1;
    }
}

fn index_for(range: AxisRange, value: f64) -> usize {
    if range.max <= range.min {
        return 0;
    }
    let v = value.clamp(range.min, range.max);
    let ratio = (v - range.min) / (range.max - range.min);
    let mut idx = (ratio * range.bins as f64).floor() as usize;
    if idx >= range.bins {
        idx = range.bins - 1;
    }
    idx
}

fn build_summary_json(
    points: &[PhasePoint],
    bins: &[BinRow],
    fii_range: AxisRange,
    vel_range: AxisRange,
    acc_range: AxisRange,
) -> Value {
    let by_proj = |name: &str| {
        bins.iter()
            .filter(|b| b.projection == name && b.sample_label == "ALL")
            .map(|b| {
                json!({
                    "bin_x": b.bin_x,
                    "bin_y": b.bin_y,
                    "count": b.count,
                    "frac_adaptive": b.frac_adaptive,
                    "frac_transition": b.frac_transition,
                    "frac_resistant": b.frac_resistant,
                    "frac_low_confidence": b.frac_low_confidence
                })
            })
            .collect::<Vec<_>>()
    };
    json!({
        "meta": {
            "version": 1,
            "bins": {
                "fii": fii_range.bins,
                "velocity": vel_range.bins,
                "acceleration": acc_range.bins
            },
            "ranges": {
                "fii": [fii_range.min, fii_range.max],
                "velocity": [vel_range.min, vel_range.max],
                "acceleration": [acc_range.min, acc_range.max]
            }
        },
        "samples": points.iter().map(|p| {
            json!({
                "label": p.sample_label,
                "rank": p.order_rank,
                "fii_mean": p.fii_mean,
                "fii_median": p.fii_median,
                "velocity": p.fii_velocity,
                "acceleration": p.fii_acceleration,
                "resistant_fraction": p.resistant_fraction,
                "low_confidence_fraction": p.low_confidence_fraction
            })
        }).collect::<Vec<_>>(),
        "bins": {
            PROJ_FII_VEL: by_proj(PROJ_FII_VEL),
            PROJ_FII_ACC: by_proj(PROJ_FII_ACC),
            PROJ_VEL_ACC: by_proj(PROJ_VEL_ACC)
        }
    })
}

fn read_config(path: Option<&Path>) -> Result<PhaseConfig, String> {
    let mut cfg = PhaseConfig::default();
    let Some(path) = path else { return Ok(cfg) };
    let raw = std::fs::read_to_string(path)
        .map_err(|e| format!("READ failed for config {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for config {}: {e}", path.display()))?;
    let get = |k: &str| v.get(k).and_then(Value::as_f64);
    if let Some(x) = get("velocity_threshold") {
        cfg.velocity_threshold = x;
    }
    if let Some(x) = get("resistant_onset_fii") {
        cfg.resistant_onset_fii = x;
    }
    if let Some(x) = get("acceleration_threshold") {
        cfg.acceleration_threshold = x;
    }
    if let Some(x) = get("fii_q") {
        cfg.fii_q = x;
    }
    if let Some(x) = get("vel_q") {
        cfg.vel_q = x;
    }
    if let Some(x) = get("resistant_fraction_jump") {
        cfg.resistant_fraction_jump = x;
    }
    if let Some(x) = get("volatility_abs_velocity_threshold") {
        cfg.volatility_abs_velocity_threshold = x;
    }
    Ok(cfg)
}

fn build_flags(points: &[PhasePoint], cfg: &PhaseConfig) -> Vec<FlagRow> {
    let mut out = Vec::new();
    for i in 1..points.len() {
        let p = &points[i];
        if let Some(v) = p.fii_velocity
            && v >= cfg.velocity_threshold
            && p.fii_mean < cfg.resistant_onset_fii
            && p.fii_acceleration.unwrap_or(0.0) >= 0.0
        {
            out.push(FlagRow {
                sample_label: p.sample_label.clone(),
                order_rank: p.order_rank,
                flag: "VELOCITY_SPIKE".to_string(),
                metric: "fii_velocity".to_string(),
                value: v,
                threshold: cfg.velocity_threshold,
                note: "rising velocity before resistant onset".to_string(),
            });
        }
        if let Some(a) = p.fii_acceleration
            && a >= cfg.acceleration_threshold
        {
            out.push(FlagRow {
                sample_label: p.sample_label.clone(),
                order_rank: p.order_rank,
                flag: "ACCELERATION_PEAK".to_string(),
                metric: "fii_acceleration".to_string(),
                value: a,
                threshold: cfg.acceleration_threshold,
                note: "acceleration exceeds threshold".to_string(),
            });
        }
        if let Some(v) = p.fii_velocity
            && p.fii_mean >= cfg.fii_q
            && v >= cfg.vel_q
        {
            let prev = &points[i - 1];
            let entered_first_time =
                prev.fii_mean < cfg.fii_q || prev.fii_velocity.unwrap_or(0.0) < cfg.vel_q;
            let resistant_jump = p.resistant_fraction - prev.resistant_fraction;
            if entered_first_time || resistant_jump >= cfg.resistant_fraction_jump {
                out.push(FlagRow {
                    sample_label: p.sample_label.clone(),
                    order_rank: p.order_rank,
                    flag: "PHASE_DRIFT_TO_RESISTANT".to_string(),
                    metric: "danger_quadrant".to_string(),
                    value: resistant_jump.max(0.0),
                    threshold: cfg.resistant_fraction_jump,
                    note: "entered danger quadrant or resistant fraction jumped".to_string(),
                });
            }
        }
    }
    for i in 2..points.len() {
        let prev_v = points[i - 1].fii_velocity.unwrap_or(0.0);
        let curr_v = points[i].fii_velocity.unwrap_or(0.0);
        if prev_v.abs() >= cfg.volatility_abs_velocity_threshold
            && curr_v.abs() >= cfg.volatility_abs_velocity_threshold
            && prev_v.signum() != curr_v.signum()
        {
            out.push(FlagRow {
                sample_label: points[i].sample_label.clone(),
                order_rank: points[i].order_rank,
                flag: "TRANSITION_VOLATILITY".to_string(),
                metric: "fii_velocity_sign_flip".to_string(),
                value: curr_v,
                threshold: cfg.volatility_abs_velocity_threshold,
                note: "high-magnitude sign flip across adjacent steps".to_string(),
            });
        }
    }
    out.sort_by(|a, b| {
        a.order_rank
            .cmp(&b.order_rank)
            .then(a.flag.cmp(&b.flag))
            .then(a.sample_label.cmp(&b.sample_label))
    });
    out
}

fn render_points_tsv(points: &[PhasePoint]) -> String {
    let mut out = String::from(
        "sample_label\torder_rank\tfii_mean\tfii_median\tfii_velocity\tfii_acceleration\tresistant_fraction\tlow_confidence_fraction\n",
    );
    for p in points {
        out.push_str(&format!(
            "{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\n",
            p.sample_label,
            p.order_rank,
            p.fii_mean,
            p.fii_median,
            opt6(p.fii_velocity),
            opt6(p.fii_acceleration),
            p.resistant_fraction,
            p.low_confidence_fraction
        ));
    }
    out
}

fn render_bins_tsv(rows: &[BinRow]) -> String {
    let mut sorted = rows.to_vec();
    sorted.sort_by(|a, b| {
        a.projection
            .cmp(&b.projection)
            .then(a.sample_label.cmp(&b.sample_label))
            .then(a.bin_x.cmp(&b.bin_x))
            .then(a.bin_y.cmp(&b.bin_y))
    });
    let mut out = String::from(
        "projection\tsample_label\tbin_x\tbin_y\tcount\tfrac_adaptive\tfrac_transition\tfrac_resistant\tfrac_low_confidence\n",
    );
    for r in sorted {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\n",
            r.projection,
            r.sample_label,
            r.bin_x,
            r.bin_y,
            r.count,
            r.frac_adaptive,
            r.frac_transition,
            r.frac_resistant,
            r.frac_low_confidence
        ));
    }
    out
}

fn render_flags_tsv(rows: &[FlagRow]) -> String {
    let mut out = String::from("sample_label\torder_rank\tflag\tmetric\tvalue\tthreshold\tnote\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\n",
            r.sample_label, r.order_rank, r.flag, r.metric, r.value, r.threshold, r.note
        ));
    }
    out
}

fn opt6(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.6}")).unwrap_or_default()
}
