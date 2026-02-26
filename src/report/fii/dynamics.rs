use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::Serialize;
use serde_json::{Value, json};

use crate::cli::ComputeStateDynamicsArgs;
use crate::io;
use crate::util::tsv::TsvReader;

use super::read::{SampleInput, resolve_inputs};

#[derive(Debug, Clone)]
struct SampleDynamics {
    label: String,
    order_rank: usize,
    fii_mean: f64,
    fii_median: f64,
    mito_mean: Option<f64>,
    translation_mean: Option<f64>,
    splice_mean: Option<f64>,
    low_confidence: bool,
}

#[derive(Debug, Clone, Serialize)]
struct DynamicsPoint {
    sample: String,
    value: f64,
}

#[derive(Debug, Clone, Serialize)]
struct DynamicsSummary {
    velocity: Vec<DynamicsPoint>,
    acceleration: Vec<DynamicsPoint>,
}

pub fn run_compute_state_dynamics(args: &ComputeStateDynamicsArgs) -> Result<(), String> {
    if args.inputs.is_empty() && args.manifest.is_none() {
        return Err(
            "INVALID compute-state-dynamics args: --inputs or --manifest required".to_string(),
        );
    }

    let mut specs = resolve_inputs(&args.inputs, args.manifest.as_deref())?;
    validate_ordering(&mut specs)?;

    let mut rows = Vec::with_capacity(specs.len());
    for spec in &specs {
        rows.push(read_sample_dynamics(spec)?);
    }
    rows.sort_by_key(|r| r.order_rank);

    std::fs::create_dir_all(&args.out)
        .map_err(|e| format!("WRITE failed creating {}: {e}", args.out.display()))?;

    let velocity_tsv = render_velocity_tsv(&rows);
    let velocity_path = args.out.join("fii_state_velocity.tsv");
    io::write_bytes_atomic(&velocity_path, velocity_tsv.as_bytes())
        .map_err(|e| format!("WRITE failed for {}: {e}", velocity_path.display()))?;

    let acceleration_tsv = render_acceleration_tsv(&rows);
    let acceleration_path = args.out.join("fii_state_acceleration.tsv");
    io::write_bytes_atomic(&acceleration_path, acceleration_tsv.as_bytes())
        .map_err(|e| format!("WRITE failed for {}: {e}", acceleration_path.display()))?;

    let dynamics = build_dynamics_summary(&rows);
    let summary_path = args.out.join("summary.json");
    let mut root = if summary_path.is_file() {
        let raw = std::fs::read_to_string(&summary_path)
            .map_err(|e| format!("READ failed for {}: {e}", summary_path.display()))?;
        serde_json::from_str::<Value>(&raw)
            .map_err(|e| format!("PARSE failed for {}: {e}", summary_path.display()))?
    } else {
        json!({})
    };
    let obj = root
        .as_object_mut()
        .ok_or_else(|| format!("INVALID summary root in {}", summary_path.display()))?;
    obj.insert(
        "fii_dynamics".to_string(),
        serde_json::to_value(dynamics).map_err(|e| format!("PARSE fii_dynamics JSON: {e}"))?,
    );
    io::write_json_atomic(&summary_path, &root)
        .map_err(|e| format!("WRITE failed for {}: {e}", summary_path.display()))?;

    Ok(())
}

fn validate_ordering(specs: &mut [SampleInput]) -> Result<(), String> {
    specs.sort_by_key(|s| s.order);
    let mut seen = BTreeMap::<usize, String>::new();
    for spec in specs {
        if let Some(prev) = seen.insert(spec.order, spec.label.clone()) {
            return Err(format!(
                "INVALID manifest ordering: duplicated order_rank {} for '{}' and '{}'",
                spec.order, prev, spec.label
            ));
        }
    }
    Ok(())
}

fn read_sample_dynamics(spec: &SampleInput) -> Result<SampleDynamics, String> {
    let tsv_path = spec.root.join("functional_irreversibility_index.tsv");
    if !tsv_path.is_file() {
        return Err(format!(
            "MISSING functional_irreversibility_index.tsv in {}",
            spec.root.display()
        ));
    }

    let mut reader = TsvReader::open(&tsv_path)
        .map_err(|e| format!("READ failed for {}: {e}", tsv_path.display()))?;
    let mut header = Vec::new();
    if !reader
        .read_record(&mut header)
        .map_err(|e| format!("READ header failed for {}: {e}", tsv_path.display()))?
    {
        return Err(format!("MISSING header in {}", tsv_path.display()));
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
    let idx = build_idx(&cols);

    let req = ["functional_irreversibility_index"];
    for col in req {
        if !idx.contains_key(col) {
            return Err(format!(
                "MISSING required column '{}' in {}",
                col,
                tsv_path.display()
            ));
        }
    }

    let mut fields = Vec::new();
    let mut fii_values = Vec::new();
    let mut mito_sum = 0.0;
    let mut mito_n = 0usize;
    let mut translation_sum = 0.0;
    let mut translation_n = 0usize;
    let mut splice_sum = 0.0;
    let mut splice_n = 0usize;
    let mut low_conf = false;

    while reader
        .read_record(&mut fields)
        .map_err(|e| format!("READ row failed for {}: {e}", tsv_path.display()))?
    {
        let fii = reader
            .field(
                &fields,
                *idx.get("functional_irreversibility_index")
                    .expect("checked"),
            )
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
        if let Some(v) = fii {
            fii_values.push(v);
        }

        if let Some(v) = parse_opt_at(
            &reader,
            &fields,
            idx.get("mitochondrial_stress_adaptation_score").copied(),
        ) {
            mito_sum += v;
            mito_n += 1;
        }
        if let Some(v) = parse_opt_at(
            &reader,
            &fields,
            idx.get("translation_commitment_score").copied(),
        ) {
            translation_sum += v;
            translation_n += 1;
        }
        if let Some(v) = parse_opt_at(
            &reader,
            &fields,
            idx.get("splice_irreversibility_index").copied(),
        ) {
            splice_sum += v;
            splice_n += 1;
        }
    }

    if fii_values.is_empty() {
        return Ok(SampleDynamics {
            label: spec.label.clone(),
            order_rank: spec.order,
            fii_mean: 0.0,
            fii_median: 0.0,
            mito_mean: None,
            translation_mean: None,
            splice_mean: None,
            low_confidence: true,
        });
    }

    fii_values.sort_by(|a, b| a.total_cmp(b));
    let mean = fii_values.iter().sum::<f64>() / fii_values.len() as f64;
    let median = percentile(&fii_values, 0.5);
    if mito_n == 0 || translation_n == 0 || splice_n == 0 {
        low_conf = true;
    }

    Ok(SampleDynamics {
        label: spec.label.clone(),
        order_rank: spec.order,
        fii_mean: mean,
        fii_median: median,
        mito_mean: if mito_n > 0 {
            Some(mito_sum / mito_n as f64)
        } else {
            None
        },
        translation_mean: if translation_n > 0 {
            Some(translation_sum / translation_n as f64)
        } else {
            None
        },
        splice_mean: if splice_n > 0 {
            Some(splice_sum / splice_n as f64)
        } else {
            None
        },
        low_confidence: low_conf,
    })
}

fn parse_opt_at(
    reader: &TsvReader,
    fields: &[std::ops::Range<usize>],
    idx: Option<usize>,
) -> Option<f64> {
    let idx = idx?;
    let raw = reader.field(fields, idx)?.trim();
    if raw.is_empty() {
        return None;
    }
    raw.parse::<f64>().ok().and_then(|v| {
        if v.is_finite() {
            Some(v.clamp(0.0, 1.0))
        } else {
            None
        }
    })
}

fn build_idx(cols: &[String]) -> BTreeMap<String, usize> {
    let mut out = BTreeMap::new();
    for (idx, col) in cols.iter().enumerate() {
        out.insert(col.clone(), idx);
    }
    out
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * q).round() as usize;
    sorted[idx]
}

fn render_velocity_tsv(rows: &[SampleDynamics]) -> String {
    let mut out = String::from(
        "sample_label\torder_rank\tfii_mean\tfii_median\tfii_velocity\tmito_velocity\ttranslation_velocity\tsplice_velocity\tlow_confidence\n",
    );
    for (i, row) in rows.iter().enumerate() {
        let prev = if i > 0 { Some(&rows[i - 1]) } else { None };
        let fii_vel = prev.map(|p| row.fii_mean - p.fii_mean).unwrap_or(0.0);
        let mito_vel = delta_opt(row.mito_mean, prev.and_then(|p| p.mito_mean));
        let translation_vel =
            delta_opt(row.translation_mean, prev.and_then(|p| p.translation_mean));
        let splice_vel = delta_opt(row.splice_mean, prev.and_then(|p| p.splice_mean));
        out.push_str(&format!(
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\n",
            row.label,
            row.order_rank,
            row.fii_mean,
            row.fii_median,
            fii_vel,
            fmt_opt(mito_vel),
            fmt_opt(translation_vel),
            fmt_opt(splice_vel),
            if row.low_confidence { "true" } else { "false" }
        ));
    }
    out
}

fn render_acceleration_tsv(rows: &[SampleDynamics]) -> String {
    let mut out = String::from(
        "sample_label\torder_rank\tfii_acceleration\tmito_acceleration\ttranslation_acceleration\tsplice_acceleration\tlow_confidence\n",
    );
    for i in 0..rows.len() {
        let row = &rows[i];
        let fii_acc = if i >= 2 {
            let v_i = rows[i].fii_mean - rows[i - 1].fii_mean;
            let v_prev = rows[i - 1].fii_mean - rows[i - 2].fii_mean;
            v_i - v_prev
        } else {
            0.0
        };
        let mito_acc = if i >= 2 {
            second_diff_opt(
                rows[i].mito_mean,
                rows[i - 1].mito_mean,
                rows[i - 2].mito_mean,
            )
        } else {
            None
        };
        let translation_acc = if i >= 2 {
            second_diff_opt(
                rows[i].translation_mean,
                rows[i - 1].translation_mean,
                rows[i - 2].translation_mean,
            )
        } else {
            None
        };
        let splice_acc = if i >= 2 {
            second_diff_opt(
                rows[i].splice_mean,
                rows[i - 1].splice_mean,
                rows[i - 2].splice_mean,
            )
        } else {
            None
        };
        out.push_str(&format!(
            "{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\n",
            row.label,
            row.order_rank,
            fii_acc,
            fmt_opt(mito_acc),
            fmt_opt(translation_acc),
            fmt_opt(splice_acc),
            if row.low_confidence { "true" } else { "false" }
        ));
    }
    out
}

fn build_dynamics_summary(rows: &[SampleDynamics]) -> DynamicsSummary {
    let mut velocity = Vec::new();
    let mut acceleration = Vec::new();
    for i in 1..rows.len() {
        velocity.push(DynamicsPoint {
            sample: rows[i].label.clone(),
            value: rows[i].fii_mean - rows[i - 1].fii_mean,
        });
    }
    for i in 2..rows.len() {
        let v_i = rows[i].fii_mean - rows[i - 1].fii_mean;
        let v_prev = rows[i - 1].fii_mean - rows[i - 2].fii_mean;
        acceleration.push(DynamicsPoint {
            sample: rows[i].label.clone(),
            value: v_i - v_prev,
        });
    }
    DynamicsSummary {
        velocity,
        acceleration,
    }
}

fn delta_opt(curr: Option<f64>, prev: Option<f64>) -> Option<f64> {
    Some(curr? - prev?)
}

fn second_diff_opt(curr: Option<f64>, prev: Option<f64>, prev2: Option<f64>) -> Option<f64> {
    let v_i = curr? - prev?;
    let v_prev = prev? - prev2?;
    Some(v_i - v_prev)
}

fn fmt_opt(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.6}")).unwrap_or_default()
}

#[allow(dead_code)]
fn _keep_path(_: &Path, _: &PathBuf) {}
