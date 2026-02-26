# Pipeline

## High-Level Pipeline Diagram
```text
Raw data
  ↓
kira-mitoqc
  ↓
kira-nuclearqc
  ↓
kira-spliceqc
  ↓
kira-proteoqc
  ↓
kira-autolys
  ↓
kira-secretion
  ↓
kira-organelle
```

## Pipeline Stages (External Tools)
1. `kira-mitoqc`
- Focus: mitochondrial quality and stress-related metrics.
- Contract artifacts: writes `summary.json`; may write `pipeline_step.json` with `artifacts.primary_metrics`.

2. `kira-nuclearqc`
- Focus: nucleus-associated QC metrics.
- Contract artifacts: writes `summary.json`; may write `pipeline_step.json` with `artifacts.primary_metrics`.

3. `kira-spliceqc`
- Focus: splicing-related QC metrics.
- Contract artifacts: writes `summary.json`; may write `pipeline_step.json` with `artifacts.primary_metrics`.

4. `kira-proteoqc`
- Focus: proteostasis/protein quality metrics.
- Contract artifacts: writes `summary.json`; may write `pipeline_step.json` with `artifacts.primary_metrics`.

5. `kira-autolys`
- Focus: autophagy/autolysosomal stress metrics.
- Contract artifacts: writes `summary.json`; may write `pipeline_step.json` with `artifacts.primary_metrics`.

6. `kira-secretion`
- Focus: secretory pathway stress/pressure metrics.
- Contract artifacts: writes `summary.json`; may write `pipeline_step.json` with `artifacts.primary_metrics`.

Optional supported tool in aggregation mode:
- `kira-riboqc` (loaded when present; not required in the default `run` chain)

## kira-organelle Internal Stages

### Stage 0 — Contract Discovery & Validation
- Inputs: per-tool directories under input root.
- Reads: `summary.json`, optional `pipeline_step.json`.
- Outputs: validation issues collected in-memory and propagated.
- Determinism: fixed expected tool order, stable issue serialization order.

### Stage 1 — Organelle State Normalization
- Inputs: tool `summary.json` values.
- Output: `state.json` (`kira-organelle-state-v1`).
- Determinism: organelle order by enum, axes sorted, map fields serialized stably.

### Stage 2 — Cell-Level Ingestion & Join
- Inputs: `pipeline_step.json` → `artifacts.primary_metrics`, then `primary_metrics.tsv`.
- Output: `cells.json` (`kira-organelle-cells-v1`).
- Determinism: streaming TSV read, lexicographic cell IDs, stable organelle/axis ordering.

### Stage 2.5 — Functional Irreversibility Index (Optional)
- Inputs (per-cell proxies):
  - `mitochondrial_stress_adaptation_score` (`kira-mitoqc`)
  - `translation_commitment_score` (`kira-riboqc`)
  - `splice_irreversibility_index` (`kira-spliceqc`)
- Enabled by CLI `--fii-weights mito:<w>,translation:<w>,splice:<w>` where weights sum to `1.0`.
- Computation:
  - clamp every input to `[0,1]`;
  - weighted sum (normalized over available inputs when some are missing);
  - clamp final FII to `[0,1]`.
- Regimes (fixed thresholds):
  - `Adaptive` if `FII < 0.33`
  - `Transition` if `0.33 <= FII < 0.66`
  - `Resistant` if `FII >= 0.66`
- Outputs:
  - `functional_irreversibility_index.tsv` (per-cell)
  - `state.json.functional_irreversibility` summary block
- Determinism: fixed thresholds, explicit arithmetic, stable row order.
- Interpretation scope: proxy-only signal, not a direct measure of biochemical irreversibility or temporal direction.

### Stage 3 — Comparison & Delta Model
- Inputs: A/B `state.json` + `cells.json` (or raw fallback compute path).
- Output: `comparison.json` (`kira-organelle-comparison-v1`).
- Determinism: explicit A→B deltas, fixed ordering, stable quantile selection.

### Stage 4 — HTML Rendering
- Inputs: in-memory state (+ comparison when present).
- Output: `report.html` (single self-contained static page).
- Determinism: fixed section order, fixed numeric formatting, no JS, no network.

### Stage 5 — Orchestration (`run`)
- Inputs: raw dataset path and run options.
- Output layout: per-tool directories plus `kira-organelle/` aggregate outputs.
- Determinism: fixed tool execution order; deterministic dry-run execution plan.

### Stage 6 — Performance & Hardening
- Inputs: same contracts as prior stages.
- Outputs: unchanged schemas and artifact formats.
- Determinism guarantees preserved while improving TSV ingestion, quantile selection, atomic IO, and diagnostics.

### Stage 7 — Interpretation Layer
- Inputs: `state.json` and optional `comparison.json`.
- Output: `interpretation.json` (`kira-organelle-interpretation-v1`).
- Determinism: fixed rule order, threshold constants, stable signal/evidence ordering.

### Stage 8 — FII Landscape Report (Post-Processing)
- Command: `kira-organelle report-fii --inputs <dir[,dir...]> --out <out_dir> [--manifest <path>] [--title <str>]`
- Inputs per sample:
  - `functional_irreversibility_index.tsv`
  - `summary.json` or `state.json` with `functional_irreversibility` section
- Outputs:
  - `fii_landscape.html` (single self-contained offline HTML)
  - `fii_landscape_index.json`
- Performance model:
  - streaming TSV read
  - stable 1D bins for FII distribution
  - stable 2D bins for coupling (component vs FII)
  - raw cell tables are not embedded by default (pre-binned payload only)
- Interpretation scope:
  - comparative visualization of proxy signals
  - not causal modeling and not direct biochemical irreversibility measurement

### Stage 9 — State Dynamics (Post-Processing)
- Command: `kira-organelle compute-state-dynamics --inputs <dir[,dir...]> --out <out_dir> [--manifest <path>]`
- Inputs per ordered sample:
  - `functional_irreversibility_index.tsv`
  - `summary.json` (or compatible sample directory with FII TSV)
- Derived quantities:
  - `fii_mean[i]`, `fii_median[i]`
  - velocity: `fii_mean[i] - fii_mean[i-1]` for `i >= 1`
  - acceleration: `velocity[i] - velocity[i-1]` for `i >= 2`
- Outputs:
  - `fii_state_velocity.tsv`
  - `fii_state_acceleration.tsv`
  - additive `summary.json.fii_dynamics`
- Determinism: stable sample ordering, deterministic parsing, fixed numeric formatting.
- Limitation:
  - State Velocity and State Acceleration are comparative metrics across ordered samples.
  - They do not represent true temporal derivatives, RNA velocity, or single-cell trajectories.

### Stage 10 — Phase Portrait + Early-Warning Flags (Post-Processing)
- Command: `kira-organelle report-phase --inputs <dir[,dir...]> [--manifest <path>] [--config <phase_config.json>] --out <out_dir>`
- Inputs:
  - `functional_irreversibility_index.tsv` per sample
  - explicit ordering from `--manifest` or CLI order
- Derived phase-space projections:
  - `FII` vs `Velocity`
  - `FII` vs `Acceleration`
  - `Velocity` vs `Acceleration`
- Outputs:
  - `phase_portrait_points.tsv`
  - `phase_portrait_bins.tsv`
  - `phase_portrait_summary.json`
  - `early_warning_flags.tsv`
- Flag system (deterministic thresholds):
  - `VELOCITY_SPIKE`
  - `ACCELERATION_PEAK`
  - `PHASE_DRIFT_TO_RESISTANT`
  - `TRANSITION_VOLATILITY`
- Limitation:
  - phase portrait is a comparative state-space summary across ordered conditions.
  - no causal inference, no pseudotime, no RNA velocity, no single-cell trajectory reconstruction.

### Stage 11 — Decision Layer (Post-Processing)
- Command: `kira-organelle decision --inputs <phase_portrait_summary.json> [--config <decision_config.json>] --out <out_dir>`
- Inputs:
  - `phase_portrait_summary.json`
  - optional `early_warning_flags.tsv` in the same directory
- Deterministic tier priority:
  1. `FIXED_RESISTANT`
  2. `PRE_RESISTANT`
  3. `TRANSITION_RISK`
  4. `STABLE_ADAPTIVE`
- Outputs:
  - `sample_decision.tsv`
  - `sample_decision.json`
- Confidence:
  - deterministic `[0..1]` score from normalized FII / velocity / acceleration + flag fraction
  - confidence does not override tier priority
- Limitation:
  - decision tiers are functional risk summaries inferred from analytics.
  - they are not clinical predictions or treatment recommendations.

### Stage 12 — Decision Timeline / Stability (Post-Processing)
- Command: `kira-organelle report-decision-timeline --inputs <decision_output_dir> --out <out_dir>`
- Inputs:
  - `sample_decision.json` (preferred) or `sample_decision.tsv`
- Outputs:
  - `decision_timeline.tsv`
  - `decision_stability.tsv`
  - `decision_stability.json`
- Deterministic metrics:
  - longest persistence run per tier
  - flip count between adjacent samples
  - volatility = `flip_count / (n_samples - 1)`
  - early flip index into `{PRE_RESISTANT, FIXED_RESISTANT}`
  - confidence-weighted stability score over tier-rank jumps
- Limitation:
  - stability/volatility are consistency descriptors of rule-based classifications
  - no forecasting, no causal inference, no ML

### Stage 13 — Irreversibility Localization Index (ILI) (Post-Processing)
- Command: `kira-organelle compute-ili --inputs <phase_portrait_summary.json> [--config <ili_config.json>] --out <out_dir>`
- Inputs:
  - ordered `samples` from `phase_portrait_summary.json`
  - optional component fields per sample:
    - `FII_nucleus`
    - `FII_splice`
    - `FII_proteostasis`
    - `FII_mito`
    - `FII_tme`
- Transition logic for each non-baseline sample:
  - compute positive `delta_component / delta_order`
  - choose argmax component if signal is above threshold and not tied within epsilon
  - otherwise assign `UNRESOLVED`
- Output categories:
  - `NUCLEUS_DRIVEN`
  - `SPLICE_DRIVEN`
  - `PROTEOSTASIS_DRIVEN`
  - `METABOLIC_DRIVEN`
  - `TME_DRIVEN`
  - `UNRESOLVED`
- Outputs:
  - `irreversibility_localization.tsv`
  - `irreversibility_localization.json`
- Limitation:
  - ILI labels the leading proxy signal among available components in ordered transitions.
  - it does not establish biological causality.

### Stage 14 — Commitment Asymmetry Index (CAI) (Post-Processing)
- Command: `kira-organelle compute-cai --inputs <dir[,dir...]> [--manifest <path>] [--config <cai_config.json>] [--cai-weights <spec>] --out <out_dir>`
- Inputs per sample:
  - `functional_irreversibility_index.tsv`
- Distribution components (deterministic):
  - quantile-based skewness component (Bowley-type, right-tail only)
  - tail-heaviness component from upper-quantile spread relative to IQR
  - tail mass above hard threshold (`FII >= 0.66`)
- Aggregation:
  - `CAI = clamp(w1*skew + w2*tail_heaviness + w3*tail_mass_norm, 0, 1)`
  - default weights: `0.33 / 0.33 / 0.34`
- Outputs:
  - `commitment_asymmetry.tsv`
  - `commitment_asymmetry.json`
- Integration scope:
  - CAI is attached as context in decision-layer/report overlays.
  - CAI does not drive decision-tier assignment.
- Limitation:
  - CAI describes distribution shape (diffuse vs subpopulation-driven fixation pattern).
  - it is not a direct severity score and does not establish causality.

### Stage 15 — Plasticity Reserve Index (PRI) (Post-Processing)
- Command: `kira-organelle compute-pri --inputs <dir[,dir...]> [--manifest <path>] [--config <pri_config.json>] [--pri-weights <spec>] --out <out_dir>`
- Inputs per sample:
  - nuclear plasticity proxy: `fraction_PlasticAdaptive` (or `1 - fraction_CommittedState`)
  - splice integrity proxy: `low_splice_noise_fraction`
  - translation selectivity proxy: `translation_selectivity_index`
- Deterministic aggregation:
  - components are clamped to `[0,1]`
  - `PRI = clamp(w_n*NPC + w_s*SIC + w_t*TSC, 0, 1)`
  - default weights: `w_n=0.4`, `w_s=0.3`, `w_t=0.3`
- Outputs:
  - `plasticity_reserve.tsv`
  - `plasticity_reserve.json`
- Integration scope:
  - PRI is attached as context in decision timeline/report overlays.
  - PRI does not override decision-tier assignment.
- Limitation:
  - PRI is a proxy for adaptive reserve; it is not causal inference or prediction.

### Stage 16 — Cross-Organelle Coupling Strength (COCS) (Post-Processing)
- Command: `kira-organelle compute-cocs --inputs <phase_portrait_summary.json> --out <out_dir>`
- Inputs:
  - ordered samples with organelle-specific irreversibility fields:
    - `FII_nucleus`
    - `FII_splice`
    - `FII_proteostasis`
    - `FII_mito`
    - `FII_tme`
- Deterministic computation:
  - per transition deltas `ΔFII_organelle`
  - pairwise robust coupling via `abs(Spearman(ΔA, ΔB))`
  - aggregate means:
    - `COCS_global`
    - `COCS_core` (nucleus/splice/proteostasis subset)
    - `COCS_extended` (all available pairs)
- Outputs:
  - `cross_organelle_coupling.tsv`
  - `cross_organelle_coupling.json`
- Integration scope:
  - COCS is attached as context in decision timeline/report overlays.
  - COCS does not override decision-tier assignment.
- Limitation:
  - COCS is a structural synchrony proxy, not causal directionality.

### Stage 17 — Decision Concordance Index (DCI) (Post-Processing)
- Command: `kira-organelle compute-dci --inputs <phase_portrait_summary.json> [--config <dci_config.json>] --out <out_dir>`
- Inputs:
  - ordered samples with organelle-level irreversibility fields (`FII_nucleus`, `FII_splice`, `FII_proteostasis`, `FII_mito`, `FII_tme`)
- Deterministic vote normalization:
  - map each organelle value to ordinal vote:
    - `0` adaptive-like (`< low_threshold`)
    - `1` transitional (`[low_threshold, high_threshold)`)
    - `2` fixed-like (`>= high_threshold`)
  - defaults: `low_threshold=0.33`, `high_threshold=0.66`
- Pairwise concordance:
  - `agreement(A,B)=1-|vote_A-vote_B|/2`
  - `DCI = mean(agreement across valid organelle pairs)`
  - `DCI_core` uses nucleus/splice/proteostasis subset
  - `DCI_extended` mirrors full-set aggregation
- Outputs:
  - `decision_concordance.tsv`
  - `decision_concordance.json`
- Integration scope:
  - DCI is attached as context in decision timeline/report overlays.
  - DCI does not override decision-tier assignment.
- Limitation:
  - DCI measures agreement structure only; not causality, directionality, or severity.

## Artifact Contracts

### `summary.json` (per-tool)
- Required for each discovered tool directory.
- Parsed schema-agnostically.
- Used for organelle-level axis/regime/QC extraction.

### `primary_metrics.tsv`
- Located via `pipeline_step.json` `artifacts.primary_metrics`.
- Parsed as TSV in streaming mode.
- Used for cell/sample matrix and cell-level comparisons.

### `pipeline_step.json`
- Optional for tools, required to locate `primary_metrics.tsv` path in contract mode.
- `kira-organelle` writes its own `pipeline_step.json` with produced artifacts.

`kira-organelle` ignores non-contract files for computation.

## Execution Modes

### `aggregate`
- Aggregates existing tool outputs under `--input`.
- Optional `--input-b` computes comparison.
- Writes `state.json`, `cells.json`, `integration/timeseries.tsv`, `integration/summary.json`, `integration/expression_aggregated.tsv`, `report.html`, `interpretation.json`, optional `comparison.json`, and `pipeline_step.json`.

### `run`
- Orchestrates full tool chain in fixed order, then runs internal `aggregate` on produced outputs.

### `--validate-only`
- `aggregate` parses and validates contracts but writes no artifacts.

### `--dry-run`
- `run` prints a deterministic execution plan and does not execute tool steps.

## Failure Handling
- `--strict`:
  - `aggregate`: missing/invalid required contracts fail execution.
  - `run`: aborts on first failed tool step.
- Non-strict mode:
  - collects issues (`warn`) and continues where possible.
- Exit codes are mapped by error taxonomy (`Io`, `Parse`, `ContractViolation`, `ToolFailure`, `Internal`) to stable process exit statuses.

## Design Principles
- Determinism first.
- Explicit file contracts.
- Rule-based interpretation only.
- Clear separation of computation (state/cells/comparison) from interpretation (signals).
