use std::fmt::Write as _;

use serde::Serialize;
use serde_json::{Value, json};

use crate::cli::ComputeCocsArgs;
use crate::io;

#[derive(Debug, Clone)]
pub struct SamplePoint {
    pub sample_label: String,
    pub order_rank: usize,
    pub nucleus: Option<f64>,
    pub splice: Option<f64>,
    pub proteostasis: Option<f64>,
    pub mito: Option<f64>,
    pub tme: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CocsRow {
    pub sample_label: String,
    pub order_rank: usize,
    pub cocs_global: f64,
    pub cocs_core: f64,
    pub cocs_extended: f64,
}

pub fn run_compute_cocs(args: &ComputeCocsArgs) -> Result<(), String> {
    let raw = std::fs::read_to_string(&args.input)
        .map_err(|e| format!("READ failed for {}: {e}", args.input.display()))?;
    let root: Value = serde_json::from_str(&raw)
        .map_err(|e| format!("PARSE failed for {}: {e}", args.input.display()))?;
    let mut points = read_points(&root)?;
    points.sort_by(|a, b| {
        a.order_rank
            .cmp(&b.order_rank)
            .then(a.sample_label.cmp(&b.sample_label))
    });
    validate_order(&points)?;

    let rows = compute_rows(&points);

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;
    io::write_bytes_atomic(
        &args.out.join("cross_organelle_coupling.tsv"),
        render_tsv(&rows).as_bytes(),
    )
    .map_err(|e| format!("WRITE failed for cross_organelle_coupling.tsv: {e}"))?;
    let json_value = json!({
        "samples": rows.iter().map(|r| {
            json!({
                "sample_label": r.sample_label,
                "order_rank": r.order_rank,
                "COCS_global": r.cocs_global,
                "COCS_core": r.cocs_core,
                "COCS_extended": r.cocs_extended
            })
        }).collect::<Vec<_>>()
    });
    io::write_json_atomic(&args.out.join("cross_organelle_coupling.json"), &json_value)
        .map_err(|e| format!("WRITE failed for cross_organelle_coupling.json: {e}"))?;
    Ok(())
}

fn read_points(root: &Value) -> Result<Vec<SamplePoint>, String> {
    let arr = root
        .get("samples")
        .and_then(Value::as_array)
        .ok_or_else(|| "MISSING samples array".to_string())?;
    let mut out = Vec::new();
    for item in arr {
        let label = item
            .get("sample_label")
            .or_else(|| item.get("label"))
            .and_then(Value::as_str)
            .ok_or_else(|| "MISSING sample label".to_string())?
            .to_string();
        let rank =
            item.get("order_rank")
                .or_else(|| item.get("rank"))
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("MISSING order_rank for {}", label))? as usize;
        out.push(SamplePoint {
            sample_label: label,
            order_rank: rank,
            nucleus: pick_component(item, &["FII_nucleus", "fii_nucleus", "nucleus_fii"]),
            splice: pick_component(item, &["FII_splice", "fii_splice", "splice_fii"]),
            proteostasis: pick_component(
                item,
                &["FII_proteostasis", "fii_proteostasis", "proteostasis_fii"],
            ),
            mito: pick_component(item, &["FII_mito", "fii_mito", "mito_fii"]),
            tme: pick_component(item, &["FII_tme", "fii_tme", "tme_fii"]),
        });
    }
    Ok(out)
}

fn pick_component(v: &Value, keys: &[&str]) -> Option<f64> {
    for key in keys {
        let Some(x) = v.get(*key).and_then(Value::as_f64) else {
            continue;
        };
        if x.is_finite() {
            return Some(x);
        }
    }
    None
}

fn validate_order(points: &[SamplePoint]) -> Result<(), String> {
    if points.len() < 2 {
        return Err("MISSING samples: at least 2 ordered samples required for COCS".to_string());
    }
    for i in 1..points.len() {
        if points[i].order_rank <= points[i - 1].order_rank {
            return Err("INVALID order_rank sequence for COCS".to_string());
        }
    }
    Ok(())
}

pub fn compute_rows(points: &[SamplePoint]) -> Vec<CocsRow> {
    let mut rows = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let (core_pairs, extended_pairs) = pairwise_corrs_split(points, i);
        let global: Vec<f64> = core_pairs
            .iter()
            .chain(extended_pairs.iter())
            .copied()
            .collect();
        rows.push(CocsRow {
            sample_label: points[i].sample_label.clone(),
            order_rank: points[i].order_rank,
            cocs_global: mean_or_zero(&global),
            cocs_core: mean_or_zero(&core_pairs),
            cocs_extended: mean_or_zero(&extended_pairs),
        });
    }
    rows
}

const PAIRS: [(&str, &str); 10] = [
    ("nucleus", "splice"),
    ("nucleus", "proteostasis"),
    ("nucleus", "mito"),
    ("nucleus", "tme"),
    ("splice", "proteostasis"),
    ("splice", "mito"),
    ("splice", "tme"),
    ("proteostasis", "mito"),
    ("proteostasis", "tme"),
    ("mito", "tme"),
];

fn is_core_pair(a: &str, b: &str) -> bool {
    matches!(
        (a, b),
        ("nucleus", "splice") | ("nucleus", "proteostasis") | ("splice", "proteostasis")
    )
}

fn pairwise_corrs_split(points: &[SamplePoint], end_idx: usize) -> (Vec<f64>, Vec<f64>) {
    let mut core = Vec::with_capacity(3);
    let mut extended = Vec::with_capacity(7);
    for (a, b) in PAIRS {
        let (xs, ys) = delta_vectors(points, end_idx, a, b);
        if xs.len() < 2 {
            continue;
        }
        let value = spearman_abs(&xs, &ys);
        if is_core_pair(a, b) {
            core.push(value);
        } else {
            extended.push(value);
        }
    }
    (core, extended)
}

fn delta_vectors(points: &[SamplePoint], end_idx: usize, a: &str, b: &str) -> (Vec<f64>, Vec<f64>) {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for i in 1..=end_idx {
        let prev = &points[i - 1];
        let curr = &points[i];
        let dt = (curr.order_rank as f64 - prev.order_rank as f64).max(1.0);
        let da = component(curr, a)
            .zip(component(prev, a))
            .map(|(c, p)| (c - p) / dt);
        let db = component(curr, b)
            .zip(component(prev, b))
            .map(|(c, p)| (c - p) / dt);
        let (Some(x), Some(y)) = (da, db) else {
            continue;
        };
        if x.is_finite() && y.is_finite() {
            xs.push(x);
            ys.push(y);
        }
    }
    (xs, ys)
}

fn component(p: &SamplePoint, key: &str) -> Option<f64> {
    match key {
        "nucleus" => p.nucleus,
        "splice" => p.splice,
        "proteostasis" => p.proteostasis,
        "mito" => p.mito,
        "tme" => p.tme,
        _ => None,
    }
}

fn mean_or_zero(v: &[f64]) -> f64 {
    if v.is_empty() {
        0.0
    } else {
        (v.iter().sum::<f64>() / v.len() as f64).clamp(0.0, 1.0)
    }
}

fn spearman_abs(xs: &[f64], ys: &[f64]) -> f64 {
    let rx = rank_average_ties(xs);
    let ry = rank_average_ties(ys);
    pearson_abs(&rx, &ry)
}

fn rank_average_ties(v: &[f64]) -> Vec<f64> {
    let mut idx = (0..v.len()).collect::<Vec<_>>();
    idx.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
    let mut ranks = vec![0.0; v.len()];
    let mut pos = 0usize;
    while pos < idx.len() {
        let start = pos;
        let val = v[idx[pos]];
        while pos < idx.len() && (v[idx[pos]] - val).abs() <= f64::EPSILON {
            pos += 1;
        }
        let avg_rank = (start + pos - 1) as f64 / 2.0 + 1.0;
        for k in start..pos {
            ranks[idx[k]] = avg_rank;
        }
    }
    ranks
}

fn pearson_abs(xs: &[f64], ys: &[f64]) -> f64 {
    if xs.len() != ys.len() || xs.len() < 2 {
        return 0.0;
    }
    let mx = xs.iter().sum::<f64>() / xs.len() as f64;
    let my = ys.iter().sum::<f64>() / ys.len() as f64;
    let mut num = 0.0;
    let mut dx2 = 0.0;
    let mut dy2 = 0.0;
    for i in 0..xs.len() {
        let dx = xs[i] - mx;
        let dy = ys[i] - my;
        num += dx * dy;
        dx2 += dx * dx;
        dy2 += dy * dy;
    }
    let den = (dx2 * dy2).sqrt();
    if den <= 1e-12 {
        0.0
    } else {
        (num / den).abs().clamp(0.0, 1.0)
    }
}

pub fn render_tsv(rows: &[CocsRow]) -> String {
    let mut out = String::with_capacity(rows.len() * 64 + 96);
    out.push_str("sample_label\torder_rank\tCOCS_global\tCOCS_core\tCOCS_extended\n");
    for r in rows {
        let _ = writeln!(
            &mut out,
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}",
            r.sample_label, r.order_rank, r.cocs_global, r.cocs_core, r.cocs_extended
        );
    }
    out
}
