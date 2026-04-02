use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::io::canonicalize_f64;
use crate::systems::composites::SystemsStateModel;
use crate::util::select::median_in_place;
use crate::version;

const ODE_LAMBDA: f64 = 0.5;
const ODE_BETA: f64 = 0.1;
const REGIME_ORDER: [&str; 5] = [
    "BaselineCompensated",
    "ImmuneSuppressiveEmbedded",
    "AdaptiveHighStress",
    "HighStress_FailedCompensation",
    "SystemicCollapseRisk",
];

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemsModelExport {
    pub model_version: String,
    pub normalization: String,
    pub composites_version: String,
    pub state_graph: DeterministicStateGraph,
    pub ode_proxy: OdeProxy,
    pub state_transition_matrix: StateTransitionMatrix,
    pub basin_connectivity: BTreeMap<String, Vec<BasinConnectivityEdge>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeterministicStateGraph {
    pub nodes: Vec<GraphNode>,
    pub edges: Vec<GraphEdge>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphNode {
    pub id: String,
    pub kind: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphEdge {
    pub from: String,
    pub to: String,
    pub weight: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OdeProxy {
    pub lambda: f64,
    pub beta: f64,
    pub mean_vector: [f64; 4],
    pub flow_matrix: [[f64; 4]; 4],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateTransitionMatrix {
    pub regimes: Vec<String>,
    pub matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BasinConnectivityEdge {
    pub neighbor: String,
    pub distance: f64,
}

#[derive(Debug, Clone)]
struct BasinCentroid {
    id: String,
    centroid: [f64; 4],
}

pub fn build_systems_model_export(systems: &SystemsStateModel) -> SystemsModelExport {
    let features = systems
        .metrics_rows
        .iter()
        .map(|r| [r.stress_vector, r.cpi, r.afs, r.imsc])
        .collect::<Vec<_>>();

    let basin_centroids = compute_basin_centroids(systems);
    let state_graph = build_state_graph(systems, &features, &basin_centroids);
    let ode_proxy = build_ode_proxy(&features);
    let state_transition_matrix = build_state_transition_matrix(systems, &features);
    let basin_connectivity = build_basin_connectivity(&basin_centroids);

    SystemsModelExport {
        model_version: version::MODEL_VERSION.to_string(),
        normalization: version::NORMALIZATION_VERSION.to_string(),
        composites_version: version::COMPOSITES_VERSION.to_string(),
        state_graph,
        ode_proxy,
        state_transition_matrix,
        basin_connectivity,
    }
}

fn compute_basin_centroids(systems: &SystemsStateModel) -> Vec<BasinCentroid> {
    let mut accum = BTreeMap::<String, ([f64; 4], usize)>::new();
    for row in &systems.metrics_rows {
        if row.basin_id == "NA" {
            continue;
        }
        let entry = accum
            .entry(row.basin_id.clone())
            .or_insert(([0.0, 0.0, 0.0, 0.0], 0));
        entry.0[0] += row.stress_vector;
        entry.0[1] += row.cpi;
        entry.0[2] += row.afs;
        entry.0[3] += row.imsc;
        entry.1 += 1;
    }

    let mut out = accum
        .into_iter()
        .filter_map(|(id, (sum, n))| {
            if n == 0 {
                None
            } else {
                Some(BasinCentroid {
                    id,
                    centroid: [
                        canonicalize_f64(sum[0] / n as f64),
                        canonicalize_f64(sum[1] / n as f64),
                        canonicalize_f64(sum[2] / n as f64),
                        canonicalize_f64(sum[3] / n as f64),
                    ],
                })
            }
        })
        .collect::<Vec<_>>();
    out.sort_by(|a, b| a.id.cmp(&b.id));
    out
}

fn build_state_graph(
    systems: &SystemsStateModel,
    features: &[[f64; 4]],
    basin_centroids: &[BasinCentroid],
) -> DeterministicStateGraph {
    let mut nodes = Vec::<GraphNode>::new();
    for basin in basin_centroids {
        nodes.push(GraphNode {
            id: basin.id.clone(),
            kind: "Basin".to_string(),
        });
    }
    for regime in REGIME_ORDER {
        nodes.push(GraphNode {
            id: regime.to_string(),
            kind: "RegimeClass".to_string(),
        });
    }
    nodes.sort_by(|a, b| a.id.cmp(&b.id).then(a.kind.cmp(&b.kind)));

    let mut edges = Vec::<GraphEdge>::new();
    let mut basin_distances = Vec::<f64>::new();
    for i in 0..basin_centroids.len() {
        for j in (i + 1)..basin_centroids.len() {
            basin_distances.push(euclidean4(
                basin_centroids[i].centroid,
                basin_centroids[j].centroid,
            ));
        }
    }
    let threshold = median_in_place(&mut basin_distances).unwrap_or(0.0);
    if threshold > 0.0 {
        for i in 0..basin_centroids.len() {
            for j in (i + 1)..basin_centroids.len() {
                let d = euclidean4(basin_centroids[i].centroid, basin_centroids[j].centroid);
                if d < threshold {
                    let w = canonicalize_f64(1.0 / (1.0 + d));
                    edges.push(GraphEdge {
                        from: basin_centroids[i].id.clone(),
                        to: basin_centroids[j].id.clone(),
                        weight: w,
                    });
                    edges.push(GraphEdge {
                        from: basin_centroids[j].id.clone(),
                        to: basin_centroids[i].id.clone(),
                        weight: w,
                    });
                }
            }
        }
    }

    let transition_matrix = build_state_transition_matrix(systems, features);
    for (i, from) in transition_matrix.regimes.iter().enumerate() {
        for (j, to) in transition_matrix.regimes.iter().enumerate() {
            if i == j {
                continue;
            }
            let w = transition_matrix.matrix[i][j];
            if w <= 0.0 {
                continue;
            }
            edges.push(GraphEdge {
                from: from.clone(),
                to: to.clone(),
                weight: canonicalize_f64(w),
            });
        }
    }

    edges.sort_by(|a, b| a.from.cmp(&b.from).then(a.to.cmp(&b.to)));
    DeterministicStateGraph { nodes, edges }
}

fn build_ode_proxy(features: &[[f64; 4]]) -> OdeProxy {
    if features.is_empty() {
        return OdeProxy {
            lambda: ODE_LAMBDA,
            beta: ODE_BETA,
            mean_vector: [0.0; 4],
            flow_matrix: [[0.0; 4]; 4],
        };
    }

    let n = features.len() as f64;
    let mut mean = [0.0f64; 4];
    for f in features {
        for i in 0..4 {
            mean[i] += f[i];
        }
    }
    for v in &mut mean {
        *v = canonicalize_f64(*v / n);
    }

    let mut cov = [[0.0f64; 4]; 4];
    for f in features {
        let d0 = f[0] - mean[0];
        let d1 = f[1] - mean[1];
        let d2 = f[2] - mean[2];
        let d3 = f[3] - mean[3];
        let d = [d0, d1, d2, d3];
        for i in 0..4 {
            for j in 0..4 {
                cov[i][j] += d[i] * d[j];
            }
        }
    }
    let denom = (features.len().saturating_sub(1)).max(1) as f64;
    for row in &mut cov {
        for v in row {
            *v = canonicalize_f64(*v / denom);
        }
    }

    let mut flow = [[0.0f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut v = ODE_BETA * cov[i][j];
            if i == j {
                v -= ODE_LAMBDA;
            }
            flow[i][j] = canonicalize_f64(v);
        }
    }

    OdeProxy {
        lambda: ODE_LAMBDA,
        beta: ODE_BETA,
        mean_vector: mean,
        flow_matrix: flow,
    }
}

fn build_state_transition_matrix(
    systems: &SystemsStateModel,
    features: &[[f64; 4]],
) -> StateTransitionMatrix {
    let n = systems.metrics_rows.len();
    let regimes = REGIME_ORDER
        .iter()
        .map(|v| (*v).to_string())
        .collect::<Vec<_>>();
    let mut matrix = vec![vec![0.0f64; regimes.len()]; regimes.len()];
    if n < 2 {
        for (i, row) in matrix.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        return StateTransitionMatrix { regimes, matrix };
    }

    let mut order = (0..n).collect::<Vec<_>>();
    order.sort_by(|a, b| {
        feature_cmp(features[*a], features[*b])
            .then(
                systems.metrics_rows[*a]
                    .cell_id
                    .cmp(&systems.metrics_rows[*b].cell_id),
            )
            .then(a.cmp(b))
    });
    let mut pos = vec![0usize; n];
    for (rank, idx) in order.iter().copied().enumerate() {
        pos[idx] = rank;
    }

    let mut counts = vec![vec![0u64; regimes.len()]; regimes.len()];
    let mut row_totals = vec![0u64; regimes.len()];
    for idx in 0..n {
        let nn = knn_for_idx(idx, features, &order, &pos, 1);
        if nn.is_empty() {
            continue;
        }
        let src = regime_index(&systems.metrics_rows[idx].regime_class);
        let dst = regime_index(&systems.metrics_rows[nn[0].0].regime_class);
        counts[src][dst] += 1;
        row_totals[src] += 1;
    }

    for i in 0..regimes.len() {
        if row_totals[i] == 0 {
            matrix[i][i] = 1.0;
            continue;
        }
        for j in 0..regimes.len() {
            matrix[i][j] = canonicalize_f64(counts[i][j] as f64 / row_totals[i] as f64);
        }
    }

    StateTransitionMatrix { regimes, matrix }
}

fn build_basin_connectivity(
    basin_centroids: &[BasinCentroid],
) -> BTreeMap<String, Vec<BasinConnectivityEdge>> {
    let mut out = BTreeMap::<String, Vec<BasinConnectivityEdge>>::new();
    if basin_centroids.is_empty() {
        return out;
    }
    for (i, basin) in basin_centroids.iter().enumerate() {
        let mut best: Option<(String, f64)> = None;
        for (j, other) in basin_centroids.iter().enumerate() {
            if i == j {
                continue;
            }
            let d = euclidean4(basin.centroid, other.centroid);
            let replace = match &best {
                None => true,
                Some((best_id, best_d)) => d < *best_d || (d == *best_d && other.id < *best_id),
            };
            if replace {
                best = Some((other.id.clone(), d));
            }
        }
        let neighbors = match best {
            Some((neighbor, distance)) => vec![BasinConnectivityEdge {
                neighbor,
                distance: canonicalize_f64(distance),
            }],
            None => Vec::new(),
        };
        out.insert(basin.id.clone(), neighbors);
    }
    out
}

fn regime_index(name: &str) -> usize {
    REGIME_ORDER
        .iter()
        .position(|x| *x == name)
        .unwrap_or_default()
}

fn feature_cmp(a: [f64; 4], b: [f64; 4]) -> std::cmp::Ordering {
    a[0].partial_cmp(&b[0])
        .unwrap_or(std::cmp::Ordering::Equal)
        .then(a[1].partial_cmp(&b[1]).unwrap_or(std::cmp::Ordering::Equal))
        .then(a[2].partial_cmp(&b[2]).unwrap_or(std::cmp::Ordering::Equal))
        .then(a[3].partial_cmp(&b[3]).unwrap_or(std::cmp::Ordering::Equal))
}

fn euclidean4(a: [f64; 4], b: [f64; 4]) -> f64 {
    let d0 = a[0] - b[0];
    let d1 = a[1] - b[1];
    let d2 = a[2] - b[2];
    let d3 = a[3] - b[3];
    canonicalize_f64((d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3).sqrt())
}

fn knn_for_idx(
    idx: usize,
    features: &[[f64; 4]],
    order: &[usize],
    pos: &[usize],
    k: usize,
) -> Vec<(usize, f64)> {
    if k == 0 || order.len() <= 1 {
        return Vec::new();
    }
    let n = order.len();
    let p = pos[idx];
    let target = features[idx];
    let mut left = p as isize - 1;
    let mut right = p + 1;
    let mut cand = Vec::<(f64, usize)>::with_capacity(k + 4);

    while left >= 0 || right < n {
        let lb_left = if left >= 0 {
            (target[0] - features[order[left as usize]][0]).abs()
        } else {
            f64::INFINITY
        };
        let lb_right = if right < n {
            (target[0] - features[order[right]][0]).abs()
        } else {
            f64::INFINITY
        };
        let worst = if cand.len() < k {
            f64::INFINITY
        } else {
            cand.last().map(|x| x.0).unwrap_or(f64::INFINITY)
        };
        if lb_left > worst && lb_right > worst && cand.len() >= k {
            break;
        }
        let take_left = lb_left <= lb_right;
        let next_rank = if take_left {
            let r = left as usize;
            left -= 1;
            r
        } else {
            let r = right;
            right += 1;
            r
        };
        let j = order[next_rank];
        if j == idx {
            continue;
        }
        let dist = euclidean4(target, features[j]);
        cand.push((dist, j));
        cand.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.1.cmp(&b.1))
        });
        if cand.len() > k {
            cand.truncate(k);
        }
    }
    cand.into_iter().map(|(d, j)| (j, d)).collect()
}
