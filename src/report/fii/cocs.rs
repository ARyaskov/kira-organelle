use serde::Serialize;
use serde_json::{Value, json};

use crate::cli::ComputeCocsArgs;
use crate::io;

#[derive(Debug, Clone)]
struct SamplePoint {
    sample_label: String,
    order_rank: usize,
    nucleus: Option<f64>,
    splice: Option<f64>,
    proteostasis: Option<f64>,
    mito: Option<f64>,
    tme: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
struct CocsRow {
    sample_label: String,
    order_rank: usize,
    cocs_global: f64,
    cocs_core: f64,
    cocs_extended: f64,
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

fn compute_rows(points: &[SamplePoint]) -> Vec<CocsRow> {
    let mut rows = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let global_pairs = pairwise_corrs(points, i, true);
        let core_pairs = pairwise_corrs(points, i, false);
        let extended = mean_or_zero(&global_pairs);
        rows.push(CocsRow {
            sample_label: points[i].sample_label.clone(),
            order_rank: points[i].order_rank,
            cocs_global: extended,
            cocs_core: mean_or_zero(&core_pairs),
            cocs_extended: extended,
        });
    }
    rows
}

fn pairwise_corrs(points: &[SamplePoint], end_idx: usize, include_extended: bool) -> Vec<f64> {
    let mut corr = Vec::new();
    let pairs: [(&str, &str); 10] = [
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
    for (a, b) in pairs {
        let core = matches!(
            (a, b),
            ("nucleus", "splice") | ("nucleus", "proteostasis") | ("splice", "proteostasis")
        );
        if !include_extended && !core {
            continue;
        }
        let (xs, ys) = delta_vectors(points, end_idx, a, b);
        if xs.len() < 2 {
            continue;
        }
        corr.push(spearman_abs(&xs, &ys));
    }
    corr
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

fn render_tsv(rows: &[CocsRow]) -> String {
    let mut out = String::from("sample_label\torder_rank\tCOCS_global\tCOCS_core\tCOCS_extended\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\n",
            r.sample_label, r.order_rank, r.cocs_global, r.cocs_core, r.cocs_extended
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(label: &str, rank: usize, n: f64, s: f64, pr: f64, m: f64, t: f64) -> SamplePoint {
        SamplePoint {
            sample_label: label.to_string(),
            order_rank: rank,
            nucleus: Some(n),
            splice: Some(s),
            proteostasis: Some(pr),
            mito: Some(m),
            tme: Some(t),
        }
    }

    #[test]
    fn coupled_vs_uncoupled() {
        let uncoupled = vec![
            p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
            p("s1", 1, 0.2, 0.05, 0.13, 0.11, 0.07),
            p("s2", 2, 0.25, 0.2, 0.11, 0.18, 0.09),
            p("s3", 3, 0.3, 0.1, 0.2, 0.15, 0.14),
        ];
        let coupled = vec![
            p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
            p("s1", 1, 0.2, 0.2, 0.2, 0.2, 0.2),
            p("s2", 2, 0.3, 0.3, 0.3, 0.3, 0.3),
            p("s3", 3, 0.45, 0.45, 0.45, 0.45, 0.45),
        ];
        let ru = compute_rows(&uncoupled);
        let rc = compute_rows(&coupled);
        assert!(rc.last().expect("c").cocs_global > ru.last().expect("u").cocs_global);
    }

    #[test]
    fn monotonic_higher_than_noisy() {
        let mono = vec![
            p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
            p("s1", 1, 0.14, 0.13, 0.15, 0.145, 0.135),
            p("s2", 2, 0.23, 0.21, 0.25, 0.235, 0.225),
            p("s3", 3, 0.26, 0.24, 0.29, 0.27, 0.255),
            p("s4", 4, 0.38, 0.35, 0.42, 0.39, 0.37),
        ];
        let noisy = vec![
            p("s0", 0, 0.10, 0.10, 0.20, 0.20, 0.20),
            p("s1", 1, 0.25, 0.13, 0.12, 0.31, 0.16),
            p("s2", 2, 0.13, 0.33, 0.26, 0.26, 0.23),
            p("s3", 3, 0.35, 0.16, 0.28, 0.10, 0.41),
            p("s4", 4, 0.17, 0.21, 0.07, 0.34, 0.32),
            p("s5", 5, 0.36, 0.10, 0.16, 0.30, 0.19),
        ];
        let rm = compute_rows(&mono);
        let rn = compute_rows(&noisy);
        assert!(rm.last().expect("m").cocs_global >= rn.last().expect("n").cocs_global);
    }

    #[test]
    fn missing_organelle_values_are_robust() {
        let mut pts = vec![
            p("s0", 0, 0.1, 0.1, 0.1, 0.1, 0.1),
            p("s1", 1, 0.2, 0.2, 0.2, 0.2, 0.2),
            p("s2", 2, 0.3, 0.3, 0.3, 0.3, 0.3),
        ];
        pts[1].tme = None;
        pts[2].mito = None;
        let rows = compute_rows(&pts);
        let last = rows.last().expect("last");
        assert!(last.cocs_global.is_finite());
        assert!((0.0..=1.0).contains(&last.cocs_global));
    }

    #[test]
    fn deterministic() {
        let pts = vec![
            p("s0", 0, 0.1, 0.2, 0.3, 0.4, 0.5),
            p("s1", 1, 0.2, 0.3, 0.4, 0.5, 0.6),
            p("s2", 2, 0.4, 0.5, 0.6, 0.7, 0.8),
        ];
        let a = compute_rows(&pts);
        let b = compute_rows(&pts);
        assert_eq!(render_tsv(&a), render_tsv(&b));
    }
}
