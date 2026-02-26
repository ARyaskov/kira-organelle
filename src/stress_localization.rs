use crate::cells::types::CellsState;
use crate::contracts::types::{Issue, Severity};
use crate::model::organelle::OrganelleId;
use crate::state::{AggregatedState, StressLocalization, StressLocalizationGating};

const EPS: f64 = 1e-9;
const MIN_VALID_CELLS: usize = 100;

pub fn compute_and_attach(state: &mut AggregatedState, cells: &CellsState) {
    let result = compute(cells, &mut state.issues);
    state.stress_localization = result;
}

pub fn compute(cells: &CellsState, issues: &mut Vec<Issue>) -> Option<StressLocalization> {
    if cells.n_cells == 0 {
        issues.push(issue_warn(
            "SLI_NO_CELLS",
            "cannot compute SLI: cells.json has no cells".to_string(),
        ));
        return None;
    }

    let requirements = [
        RequiredMetric::new(
            OrganelleId::Spliceosome,
            "splice_stress_index",
            &[
                "splice_stress_index",
                "splicing_stress_index",
                "stress_splicing_index",
                "sis",
            ],
        ),
        RequiredMetric::new(
            OrganelleId::Proteostasis,
            "proteo_stress_index",
            &[
                "proteo_stress_index",
                "proteostasis_stress_index",
                "stress_proteostasis_index",
                "proteostasis_load",
                "stress_index",
            ],
        ),
        RequiredMetric::new(
            OrganelleId::Nucleus,
            "nuclear_stress_index",
            &[
                "nuclear_stress_index",
                "nuclear_stress",
                "a4_trs",
                "stress_index",
            ],
        ),
        RequiredMetric::new(
            OrganelleId::Secretion,
            "stress_secretion_index",
            &["stress_secretion_index"],
        ),
        RequiredMetric::new(
            OrganelleId::Secretion,
            "er_golgi_pressure",
            &["er_golgi_pressure"],
        ),
    ];

    let mut missing = Vec::new();
    for req in requirements {
        if !cells_have_axis(cells, req.organelle, req.aliases) {
            missing.push(format!("{:?}:{}", req.organelle, req.canonical));
        }
    }
    if !missing.is_empty() {
        issues.push(issue_error(
            "SLI_MISSING_INPUT_AXES",
            format!(
                "cannot compute SLI: missing required axes [{}]",
                missing.join(", ")
            ),
        ));
        return None;
    }

    let mut internal_raw = Vec::new();
    let mut external_raw = Vec::new();
    let mut fragility_raw = Vec::new();

    for cell in &cells.cells {
        let per = &cell.per_organelle;
        if has_low_confidence_flag(per) {
            continue;
        }

        let v_splice = axis_value(
            per,
            OrganelleId::Spliceosome,
            &[
                "splice_stress_index",
                "splicing_stress_index",
                "stress_splicing_index",
                "sis",
            ],
        );
        let v_proteo = axis_value(
            per,
            OrganelleId::Proteostasis,
            &[
                "proteo_stress_index",
                "proteostasis_stress_index",
                "stress_proteostasis_index",
                "proteostasis_load",
                "stress_index",
            ],
        );
        let v_trans = axis_value(
            per,
            OrganelleId::Ribosome,
            &[
                "translation_stress_index",
                "ribo_stress_index",
                "translation_load",
            ],
        );
        let v_nuclear = axis_value(
            per,
            OrganelleId::Nucleus,
            &[
                "nuclear_stress_index",
                "nuclear_stress",
                "a4_trs",
                "stress_index",
            ],
        );
        let v_sec = axis_value(per, OrganelleId::Secretion, &["stress_secretion_index"]);
        let v_er = axis_value(per, OrganelleId::Secretion, &["er_golgi_pressure"]);
        let v_barrier = axis_value(
            per,
            OrganelleId::Energetics,
            &["barrier_min", "mean_barrier_min", "median_barrier_min"],
        );
        let v_unstable = axis_value(
            per,
            OrganelleId::Energetics,
            &["unstable_trap_flag", "unstable_trap"],
        );
        let v_trap = axis_value(
            per,
            OrganelleId::Energetics,
            &["trap_depth", "deep_trap_fraction", "unstable_trap_fraction"],
        );
        let v_autophagy = axis_value(
            per,
            OrganelleId::Autophagy,
            &["stress_autophagy_index", "autophagy_flux_proxy"],
        );

        if let (Some(splice), Some(proteo), Some(nuclear), Some(sec), Some(er)) =
            (v_splice, v_proteo, v_nuclear, v_sec, v_er)
        {
            let mut internal = vec![splice, proteo, nuclear];
            if let Some(trans) = v_trans {
                internal.push(trans);
            }

            let mut fragility = Vec::new();
            if let Some(barrier) = v_barrier {
                fragility.push(1.0 - barrier);
            }
            if let Some(unstable) = v_unstable {
                fragility.push(if unstable >= 0.5 { 1.0 } else { 0.0 });
            }
            if let Some(trap) = v_trap {
                fragility.push(trap);
            }
            if fragility.is_empty() {
                if let Some(auto) = v_autophagy {
                    fragility.push(auto);
                } else {
                    fragility.push(sec);
                }
            }

            internal_raw.push(median_dynamic(&internal));
            external_raw.push(median_of2(sec, er));
            fragility_raw.push(median_dynamic(&fragility));
        }
    }

    if internal_raw.len() < MIN_VALID_CELLS {
        issues.push(issue_warn(
            "SLI_INSUFFICIENT_CELLS",
            format!(
                "cannot interpret SLI: valid QC-passing cells < {} (got {})",
                MIN_VALID_CELLS,
                internal_raw.len()
            ),
        ));
        return None;
    }

    let i_axis = robust_z(&internal_raw);
    let e_axis = robust_z(&external_raw);
    let f_axis = robust_z(&fragility_raw);

    let rho_ie = spearman_rho(&i_axis, &e_axis);
    let rho_if = spearman_rho(&i_axis, &f_axis);
    let coupled = (0.5 * rho_ie.max(0.0) + 0.4 * rho_if.max(0.0)).clamp(0.0, 1.0);
    let sli = 1.0 - coupled;

    let i_top20_thr = percentile_nearest_rank(&i_axis, 0.80);
    let e_thr = percentile_nearest_rank(&e_axis, 0.75);
    let f_thr = percentile_nearest_rank(&f_axis, 0.75);

    let mut n_top = 0usize;
    let mut n_ext = 0usize;
    let mut n_frag = 0usize;
    for idx in 0..i_axis.len() {
        if i_axis[idx] >= i_top20_thr {
            n_top += 1;
            if e_axis[idx] >= e_thr {
                n_ext += 1;
            }
            if f_axis[idx] >= f_thr {
                n_frag += 1;
            }
        }
    }

    let p_ext = if n_top == 0 {
        0.0
    } else {
        n_ext as f64 / n_top as f64
    };
    let p_frag = if n_top == 0 {
        0.0
    } else {
        n_frag as f64 / n_top as f64
    };

    Some(StressLocalization {
        sli: sli.clamp(0.0, 1.0),
        rho_internal_external: rho_ie,
        rho_internal_fragility: rho_if,
        gating: StressLocalizationGating {
            p_ext_top20: p_ext.clamp(0.0, 1.0),
            p_frag_top20: p_frag.clamp(0.0, 1.0),
        },
    })
}

#[derive(Clone, Copy)]
struct RequiredMetric<'a> {
    organelle: OrganelleId,
    canonical: &'a str,
    aliases: &'a [&'a str],
}

impl<'a> RequiredMetric<'a> {
    const fn new(organelle: OrganelleId, canonical: &'a str, aliases: &'a [&'a str]) -> Self {
        Self {
            organelle,
            canonical,
            aliases,
        }
    }
}

fn cells_have_axis(cells: &CellsState, organelle: OrganelleId, aliases: &[&str]) -> bool {
    cells
        .organelle_axes
        .iter()
        .find(|x| x.organelle == organelle)
        .map(|x| {
            x.axes
                .iter()
                .any(|axis| aliases.iter().any(|a| axis.eq_ignore_ascii_case(a)))
        })
        .unwrap_or(false)
}

fn axis_value(
    per: &std::collections::BTreeMap<OrganelleId, crate::cells::types::OrganelleCellState>,
    organelle: OrganelleId,
    aliases: &[&str],
) -> Option<f64> {
    let state = per.get(&organelle)?;
    for alias in aliases {
        if let Some(v) = state.axes.get(*alias) {
            return Some(*v);
        }
        if let Some((_, v)) = state
            .axes
            .iter()
            .find(|(k, _)| k.eq_ignore_ascii_case(alias))
        {
            return Some(*v);
        }
    }
    None
}

fn has_low_confidence_flag(
    per: &std::collections::BTreeMap<OrganelleId, crate::cells::types::OrganelleCellState>,
) -> bool {
    per.values().any(|v| {
        v.flags
            .iter()
            .any(|f| f.to_ascii_uppercase().contains("LOW_CONFIDENCE"))
    })
}

fn robust_z(x: &[f64]) -> Vec<f64> {
    let med = median(x);
    let dev = x.iter().map(|v| (v - med).abs()).collect::<Vec<_>>();
    let mad = median(&dev);
    x.iter()
        .map(|v| {
            let z = (v - med) / (mad + EPS);
            let zc = z.clamp(-3.0, 3.0);
            (zc + 3.0) / 6.0
        })
        .collect()
}

fn spearman_rho(a: &[f64], b: &[f64]) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return 0.0;
    }
    let ra = ranks_with_ties(a);
    let rb = ranks_with_ties(b);
    pearson(&ra, &rb)
}

fn ranks_with_ties(v: &[f64]) -> Vec<f64> {
    let mut indexed = v
        .iter()
        .enumerate()
        .map(|(i, x)| (i, *x))
        .collect::<Vec<_>>();
    indexed.sort_by(|a, b| a.1.total_cmp(&b.1));

    let mut out = vec![0.0; v.len()];
    let mut i = 0usize;
    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && indexed[j].1 == indexed[i].1 {
            j += 1;
        }
        let rank_start = i as f64 + 1.0;
        let rank_end = j as f64;
        let avg = (rank_start + rank_end) / 2.0;
        for k in i..j {
            out[indexed[k].0] = avg;
        }
        i = j;
    }
    out
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let mean_a = a.iter().sum::<f64>() / n;
    let mean_b = b.iter().sum::<f64>() / n;
    let mut cov = 0.0;
    let mut var_a = 0.0;
    let mut var_b = 0.0;
    for i in 0..a.len() {
        let da = a[i] - mean_a;
        let db = b[i] - mean_b;
        cov += da * db;
        var_a += da * da;
        var_b += db * db;
    }
    if var_a <= EPS || var_b <= EPS {
        return 0.0;
    }
    cov / (var_a.sqrt() * var_b.sqrt())
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut data = values.to_vec();
    data.sort_by(|a, b| a.total_cmp(b));
    let n = data.len();
    if n % 2 == 1 {
        data[n / 2]
    } else {
        (data[n / 2 - 1] + data[n / 2]) / 2.0
    }
}

fn percentile_nearest_rank(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut data = values.to_vec();
    data.sort_by(|a, b| a.total_cmp(b));
    let n = data.len();
    let rank = ((n as f64) * p).ceil() as usize;
    let idx = rank.saturating_sub(1).min(n - 1);
    data[idx]
}

fn median_of2(a: f64, b: f64) -> f64 {
    (a + b) / 2.0
}

fn median_dynamic(values: &[f64]) -> f64 {
    median(values)
}

fn issue_error(code: &str, message: String) -> Issue {
    Issue {
        severity: Severity::Error,
        tool: Some("kira-organelle".to_string()),
        code: code.to_string(),
        message,
        path: None,
    }
}

fn issue_warn(code: &str, message: String) -> Issue {
    Issue {
        severity: Severity::Warn,
        tool: Some("kira-organelle".to_string()),
        code: code.to_string(),
        message,
        path: None,
    }
}
