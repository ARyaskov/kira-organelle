# Kira Organelle Systems State Metrics

## Systems State Integrator Outputs

`kira-organelle` emits cross-axis systems metrics to `integration/metrics.tsv`:

- `StressVector`: Euclidean multi-axis stress magnitude across nuclear, splicing, mitochondrial, translation, proteostasis, autolysosomal, and microenvironmental stress contributors.
- `CompensationDeficit`: mean positive imbalance between stress drivers and compensatory capacity across translation/proteostasis, mito/autophagy, nuclear/splicing, and hypoxia/energetics gaps.
- `CPI`: collapse proximity index combining stress magnitude, compensation deficit, and proteostasis collapse pressure (`PCP`).
- `AFS`: adaptation-vs-failure score computed as adaptive reserve minus compensation deficit.
- `IMSC`: immune-metabolic stress coupling from immune suppression/metabolic stress (`IMSI`) and microenvironment suppression (`MSM`).
- `RegimeClass`: deterministic regime label derived from threshold rules (no ML).
- `RareState`: boolean rare-system-state flag per cluster from Mahalanobis tailing.
- `MahalanobisDistance`: distance in 4D systems space `[StressVector, CPI, AFS, IMSC]`.
- `DominantAxis`: per-cell axis with largest absolute finite-difference CPI sensitivity.
- `DominantAxisSensitivity`: magnitude of that dominant axis sensitivity.
- `Potential`: deterministic cellular potential energy (`max(0, CPI - AFS)`).
- `StabilityGradient`: local potential slope from deterministic kNN in systems space.
- `TPI_landscape`: transition-proximity index from gradient and basin-threshold proximity.
- `BasinId`: deterministic low-potential basin assignment (`Basin_N` or `NA`).
- `TransitionCandidate`: true when `TPI_landscape >= median + 2*MAD`.
- `DominantVulnerableAxis`: per-cell axis with highest deterministic AVS (Axis Vulnerability Score).
- `LSI`: Lyapunov-like local stability index in systems space, clamped to `[0,1]`.
- `MaxTrajectory`: max deterministic perturbation trajectory magnitude across major stress axes.
- `BEE`: basin escape energy proxy (`max(0, T_basin - Potential)` for basin-assigned cells).

## Global Re-normalization

Before cross-axis integration, each detected metric is robustly re-normalized across merged per-cell data:

`Z_global = (x - median(x)) / (1.4826 * MAD(x) + eps)`

- If `MAD == 0`, `Z_global` is set to `0`.
- Missing metrics are treated as unavailable and excluded from means/gaps.
- Raw and normalized values are exported in `integration/metrics_normalized.tsv`.

## Robust Z-score Rules (Stage 1)

- `eps = 1e-9` for denominator stability.
- Finite values are used for `median`/`MAD`; non-finite values are ignored.
- If `n_valid < MIN_VALID_CELLS` (`20`):
  - metric is marked `reliable=false`
  - z-scores are stored as `NaN` for all cells for that metric
- If `MAD == 0`:
  - metric is marked `scale_degenerate=true`
  - z-scores are `0` for finite inputs and `NaN` for missing inputs
- Missing input values always remain `NaN` in normalized storage.

## Regime Classification Rules

Rules are deterministic and applied in this order:

1. `CPI >= 3.0` -> `Systemic Collapse Risk`
2. `StressVector >= 3.0 && AFS < 0` -> `High Stress / Failed Compensation`
3. `StressVector >= 3.0 && AFS >= 0` -> `Adaptive High-Stress State`
4. `IMSC >= 2.0` -> `Immune-Suppressive Embedded State`
5. Else -> `Baseline / Compensated`

## Summary JSON

`integration/summary.json` contains:

- `systems_state_model.global_stats` (`n_cells`, CPI/stress p50+p90, entropy, global rare fraction)
- `systems_state_model.cluster_stats` (per-cluster robust summaries, tail fractions, heterogeneity index, rare fraction)
- `systems_state_model.regime_distribution` (global regime counts/fractions)
- `systems_state_model.missing_metric_counts` (composite-input missingness diagnostics)

## Topology Analytics (Stage 3)

- Cluster robust summaries are computed for `StressVector`, `CPI`, `CompensationDeficit`, `AFS`, `IMSC`: `median`, `p10`, `p90`, `mad`.
- Tail thresholds are global robust cutoffs:
  - `T_CPI = median(CPI) + 2*MAD(CPI)`
  - `T_Stress = median(StressVector) + 2*MAD(StressVector)`
- Per cluster:
  - `tail_fraction_cpi = frac(CPI >= T_CPI)`
  - `tail_fraction_stress = frac(StressVector >= T_Stress)`
  - `heterogeneity_index = mean(var(StressVector), var(CPI), var(AFS)) + tail_fraction_cpi`
- Rare-state detection:
  - Mahalanobis distance per cluster with covariance regularization (`lambda = 1e-3 * trace(Sigma)/4`) when needed.
  - `RareState = D >= median(D) + 2*MAD(D)` (cluster-local threshold), skipped for clusters `< 10` cells.
- System entropy:
  - `H = -sum_i p_i ln(p_i)` over regime fractions
  - `H_norm = H / ln(5)` for the 5 deterministic regime classes

## Fragility Analytics (Stage 4)

- Sensitivity axes are evaluated in stable metric-id order and require metric presence plus Stage-1 reliability.
- CPI local sensitivity for axis `X` uses centered finite differences in Z-space:
  - `Sens(X,c) = (CPI(X+delta) - CPI(X-delta)) / (2*delta)`, with `delta=0.25`
  - `AbsSens(X,c) = |Sens(X,c)|`
- Missing per-cell axis values yield zero sensitivity for that cell/axis.
- `DominantAxis` is `argmax AbsSens` per cell, with deterministic tie-break on lowest `MetricId`.
- `systems_state_model.fragility` in `summary.json` exposes:
  - `delta`
  - `global_axis_ranking` (top 10 by median `AbsSens`)
  - `cluster_axis_ranking` (top 5 per cluster by median `AbsSens`)

## Landscape Analytics (Stage 6)

- Systems feature vector per cell: `[StressVector, CPI, AFS, IMSC]`.
- Potential:
  - `Potential = max(0, CPI - AFS)`
- Local gradient:
  - deterministic `k=10` nearest-neighbor search in feature space
  - `Gradient = mean(|Potential(c)-Potential(nn)|) / mean(dist(c,nn))`
- Basin detection:
  - low-energy cells satisfy `Potential <= median(Potential)`
  - adjacency on deterministic kNN edges where `dist < eps`
  - `eps = median(kNN edge distances)`
  - connected components with size `< 5` are ignored
- Basin metrics:
  - `depth = T_basin - median(Potential in basin)`
  - `width = mean internal kNN-edge distance in basin`
  - `stability = depth / (1 + width)`
- Transition proximity:
  - `TPI_landscape = Gradient * 1/(1+|Potential - T_basin|)`
- `TransitionCandidate` threshold uses `median + 2*MAD`

## Cross-Sample Alignment (Stage 7)

Available in multi-sample integration mode.

- Harmonized Z-space:
  - pooled robust normalization across all samples (no per-sample batch correction)
  - for each metric `M`: `Z_harmonized = (x - median_all) / (1.4826*MAD_all + eps)`
- Regime distribution comparison:
  - pairwise Jensen-Shannon divergence matrix over deterministic regime fractions
- Basin alignment:
  - basin centroids in harmonized systems space
  - deterministic overlap score per sample pair
- Differential Stress Projection (DSP):
  - axis-wise median shifts between samples in harmonized stress-axis space (`RSS,SII,OSL,TSM,PCP,ASM,MSM`)
  - vector norm reported per pair
- Systems Divergence Index (SDI):
- `SDI = 0.4*JS + 0.3*(1-overlap) + 0.3*normalized(||DSP||)`
  - DSP norm is normalized by robust spread (`MAD`) across pairwise DSP norms

## Therapeutic Projection (Stage 8)

- Axis Vulnerability Score (AVS), per cell and axis:
  - `AVS_X = pos(Z(X)) * AbsSens(X) * (1/(1+CompCapacity_X))`
  - if compensatory axis is unavailable, the denominator term is omitted
- Target Opportunity Index (TOI), per cluster and axis:
  - `TOI_X = median(AVS_X)`
  - `TOI_norm_X = TOI_X / (1 + median_cluster(TOI_*))`
- Compensation Risk Index (CRI), per cluster and axis:
  - deterministic proxy using positive sensitivity-correlation to compensatory axes
  - no ML or stochastic coupling model
- Therapeutic Priority Score (TPS):
  - `TPS_X = TOI_norm_X * (1 - CRI_X)`
- Combination Synergy Potential (CSP), per cluster and axis-pair:
  - `CSP_XY = (TOI_X + TOI_Y) * |corr_sens(X,Y)|`
  - top 5 combinations are reported per cluster

`systems_state_model.therapeutic_projection` in `integration/summary.json` includes:
- `cluster_rankings`
- `global_axis_ranking`
- `top_combinations`

## Dynamic Stability Modeling (Stage 9)

- Systems vector: `[StressVector, CPI, AFS, IMSC]`.
- LSI per cell:
  - `LSI = 1 / (1 + LocalDispersion * (1 + LocalPotentialVariance))`
  - computed from deterministic `k=10` nearest neighbors in systems space.
- Attractor robustness (per basin):
  - `ARS = Depth / (1 + Width)`, and `ARS_norm = ARS / max(ARS)`.
- Perturbation trajectory projection:
  - stress-axis perturbation in Z-space with fixed `Δ=0.5`
  - recompute local systems composites deterministically
  - store per-cell `MaxTrajectory` only.
- System Criticality Index:
  - `SCI_system = 0.4*norm(mean_gradient) + 0.3*norm(rare_fraction) + 0.3*norm(system_entropy)`
  - normalization is deterministic median-based scaling to `[0,1]`.

`systems_state_model.dynamic_stability` includes:
- `mean_lsi`
- `global_sci_system`
- `basin_robustness` (`basin_id`, `ars`, `ars_norm`)

## Validation and Biological Grounding (Stage 11)

`systems_state_model.validation` includes deterministic checks:

- `invariants`: counts and examples for CPI non-negativity, basin-escape consistency, and collapse-regime threshold consistency.
- `robustness_index`: envelope stability from deterministic parameter perturbations (`k=10` vs `k=8`, and fragility delta `0.25` vs `0.2`).
- `effective_degrees_of_freedom`: covariance-rank estimate over active axes.
- `condition_number`: covariance conditioning diagnostic.
- `biological_constraints`: configured directional correlation checks (`OSL~PCP`, `HSI~ISS`, `TSM~CCI`) with violation list.
- `sensitivity_stability_index`: ranking stability after deterministic top-5% CPI removal.
- `overparameterization_flag`: complexity-vs-DoF warning.

## Caveats

- Integration is purely transcriptional and metric-based.
- This is not a predictive clinical tool.
- Intended for reproducible hypothesis generation and cross-axis state interpretation.
