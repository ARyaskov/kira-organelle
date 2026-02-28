# kira-organelle

## Project Overview
`kira-organelle` is the aggregation and orchestration layer for the Kira organelle QC stack. It consumes the per-tool contract artifacts emitted by `kira-mitoqc`, `kira-nuclearqc`, `kira-spliceqc`, `kira-riboqc`, `kira-proteoqc`, `kira-autolys`, and `kira-secretion`, then produces unified machine-readable outputs and a static report.

The project is deterministic and rule-based by design: identical inputs produce byte-identical outputs. Interpretation is explicit threshold logic over computed metrics and deltas. There is no ML model, no learned scoring, and no black-box inference.

## Installation
Rust requirement: `1.95+`.

Build from source:

```bash
cargo build --release
```

Optional install:

```bash
cargo install --path .
```

## Quick Start

### A) Aggregate existing pipeline outputs

```bash
kira-organelle aggregate \
  --input ./out \
  --out ./out/kira-organelle
```

Produces (depending on mode/data availability):
- `state.json`
- `cells.json`
- `integration/timeseries.tsv`
- `integration/summary.json`
- `integration/expression_aggregated.tsv`
- `functional_irreversibility_index.tsv` (optional; when `--fii-weights` is set)
- `report.html`
- `interpretation.json`
- `comparison.json` (only when `--input-b` is provided)
- `pipeline_step.json`

### B) Full pipeline run (end-to-end)

```bash
kira-organelle run \
  --input ./data \
  --integration-manifest ./manifest.tsv \
  --out ./out
```

`run` executes tools in fixed order:
1. `kira-mitoqc`
2. `kira-nuclearqc`
3. `kira-spliceqc`
4. `kira-proteoqc`
5. `kira-autolys`
6. `kira-secretion`

`kira-riboqc` remains supported in `aggregate` mode and is loaded when present.

Then it invokes internal `aggregate` to generate final outputs under `./out/kira-organelle`.
When `--integration-manifest` is provided, `integration/timeseries.tsv` follows manifest `order` / `timepoint`.

Shared cache behavior:
- default cache path: `<input>/kira-organelle.bin`
- propagated to all tool invocations in `run` mode as `--cache <path>`
- disabled with `--no-cache`

Binary lookup behavior for each external tool:
- `KIRA_ORGANELLE_BIN_<TOOL>` override (if set)
- executable in the same directory as `kira-organelle`
- executable in `PATH`

### C) Comparison (A vs B)

```bash
kira-organelle aggregate \
  --input ./out_before \
  --input-b ./out_after \
  --out ./out_diff
```

Comparison mode adds `comparison.json` and comparison sections in `report.html`.

### D) FII Landscape HTML

```bash
kira-organelle report-fii \
  --inputs ./out/baseline,./out/c2,./out/eot \
  --out ./out/fii-report \
  --title "FII Landscape"
```

Outputs:
- `fii_landscape.html` (single self-contained offline file)
- `fii_landscape_index.json` (machine-readable sample index)

The HTML now includes a **Decision Overlay** section:
- decision-tier markers aligned with timeline and phase coordinates,
- confidence encoded by marker opacity,
- tooltip with tier/confidence/drivers and ILI/CAI/PRI/COCS/DCI (when available),
- toggles: overlay/trajectory/dim-base.

Interpretation caveat:
Decision Overlay visualizes deterministic decision-layer outputs. It does not represent biological trajectories or time-continuous dynamics.

Optional manifest:

```bash
kira-organelle report-fii \
  --manifest ./manifest.tsv \
  --out ./out/fii-report
```

Manifest accepted formats:
- `json`: array of `{path,label,order,timepoint?}`
- `tsv/csv`: columns including `path`, optional `label`, `order`, `timepoint`

### E) State Dynamics (Velocity / Acceleration)

```bash
kira-organelle compute-state-dynamics \
  --inputs ./out/baseline,./out/c2,./out/eot \
  --out ./out/fii-dynamics
```

or with explicit manifest ordering:

```bash
kira-organelle compute-state-dynamics \
  --manifest ./manifest.csv \
  --out ./out/fii-dynamics
```

Outputs:
- `fii_state_velocity.tsv`
- `fii_state_acceleration.tsv`
- `summary.json` with additive `fii_dynamics` block

### F) Phase Portrait + Early-Warning Flags

```bash
kira-organelle report-phase \
  --inputs ./out/baseline,./out/c2,./out/eot \
  --out ./out/phase-report
```

Optional thresholds config:

```bash
kira-organelle report-phase \
  --manifest ./manifest.csv \
  --config ./phase_config.json \
  --out ./out/phase-report
```

Outputs:
- `phase_portrait_points.tsv`
- `phase_portrait_bins.tsv`
- `phase_portrait_summary.json`
- `early_warning_flags.tsv`

### G) Decision Layer

```bash
kira-organelle decision \
  --inputs ./out/phase-report/phase_portrait_summary.json \
  --out ./out/decision
```

Optional threshold overrides:

```bash
kira-organelle decision \
  --inputs ./out/phase-report/phase_portrait_summary.json \
  --config ./decision_config.json \
  --out ./out/decision
```

Outputs:
- `sample_decision.tsv`
- `sample_decision.json`

### H) Decision Timeline / Stability Report

```bash
kira-organelle report-decision-timeline \
  --inputs ./out/decision \
  --out ./out/decision-timeline
```

Outputs:
- `decision_timeline.tsv`
- `decision_stability.tsv`
- `decision_stability.json`

### I) Irreversibility Localization Index (ILI)

```bash
kira-organelle compute-ili \
  --inputs ./out/phase-report/phase_portrait_summary.json \
  --out ./out/ili
```

Optional thresholds config:

```bash
kira-organelle compute-ili \
  --inputs ./out/phase-report/phase_portrait_summary.json \
  --config ./ili_config.json \
  --out ./out/ili
```

Outputs:
- `irreversibility_localization.tsv`
- `irreversibility_localization.json`

Canonical categories:
- `NUCLEUS_DRIVEN`
- `SPLICE_DRIVEN`
- `PROTEOSTASIS_DRIVEN`
- `METABOLIC_DRIVEN`
- `TME_DRIVEN`
- `UNRESOLVED`

Interpretation caveat:
ILI identifies the leading irreversibility signal in ordered-sample transitions. It is not a causal attribution.

### J) Commitment Asymmetry Index (CAI)

```bash
kira-organelle compute-cai \
  --inputs ./out/baseline,./out/c2,./out/eot \
  --out ./out/cai
```

Optional inputs ordering via manifest and optional config/weights:

```bash
kira-organelle compute-cai \
  --manifest ./manifest.csv \
  --config ./cai_config.json \
  --cai-weights skewness:0.33,tail_heaviness:0.33,tail_mass:0.34 \
  --out ./out/cai
```

Outputs:
- `commitment_asymmetry.tsv`
- `commitment_asymmetry.json`

Interpretation:
- low CAI: diffuse population-wide adaptation
- high CAI: subpopulation-driven hard fixation pattern

Interpretation caveat:
CAI measures FII distribution structure (asymmetry/tail pattern), not irreversibility severity and not causality.

### K) Plasticity Reserve Index (PRI)

```bash
kira-organelle compute-pri \
  --inputs ./out/baseline,./out/c2,./out/eot \
  --out ./out/pri
```

Optional manifest and weight overrides:

```bash
kira-organelle compute-pri \
  --manifest ./manifest.csv \
  --pri-weights nuclear:0.4,splice:0.3,translation:0.3 \
  --out ./out/pri
```

Outputs:
- `plasticity_reserve.tsv`
- `plasticity_reserve.json`

Interpretation:
- high PRI: preserved adaptive maneuverability
- low PRI: reduced plastic reserve and increased risk of fixed commitment

Interpretation caveat:
PRI is a deterministic aggregate proxy of adaptive reserve. It contextualizes decision states and does not override decision tiers.

### L) Cross-Organelle Coupling Strength (COCS)

```bash
kira-organelle compute-cocs \
  --inputs ./out/phase-report/phase_portrait_summary.json \
  --out ./out/cocs
```

Outputs:
- `cross_organelle_coupling.tsv`
- `cross_organelle_coupling.json`

Interpretation:
- low COCS: modular/semi-independent organelle dynamics
- high COCS: tightly coupled system dynamics

Interpretation caveat:
COCS is a structural coupling proxy over organelle irreversibility deltas; it is not causality or pathway direction.

### M) Decision Concordance Index (DCI)

```bash
kira-organelle compute-dci \
  --inputs ./out/phase-report/phase_portrait_summary.json \
  --out ./out/dci
```

Outputs:
- `decision_concordance.tsv`
- `decision_concordance.json`

Interpretation:
- low DCI: subsystem disagreement / internal conflict
- high DCI: subsystem agreement / structurally stable system-level decision

Interpretation caveat:
DCI is an agreement metric between organelle-level decision votes; it does not encode severity or causality.

## Command Reference

### `aggregate`
- `--input <DIR>`: required input pipeline directory.
- `--input-b <DIR>`: optional second input for A→B comparison.
- `--out <DIR>`: output directory (default: `<input>/kira-organelle`).
- `--strict`: treat missing/invalid required contracts as hard failures.
- `--validate-only`: parse and validate without writing artifacts.
- `--json`: accepted flag; JSON artifacts are emitted in normal aggregate mode.
- `--fii-weights <SPEC>`: enable Functional Irreversibility Index (FII) and set weights, format `mito:0.34,translation:0.33,splice:0.33`.

### `run`
- `--input <DIR>`: required raw dataset input directory (MTX/features/barcodes(.gz) layout or BD Rhapsody export layout; forwarded to each tool).
- `--out <DIR>`: pipeline output root.
- `--threads <N>`: pass thread count through to tool invocations.
- `--no-cache`: disable shared cache path propagation.
- `--strict`: abort on first tool failure.
- `--dry-run`: print deterministic execution plan only.
- `--fii-weights <SPEC>`: forward FII weighting into final aggregation step.

## Output Artifacts

| File | Description |
| --- | --- |
| `state.json` | Organelle-level normalized aggregates |
| `cells.json` | Cell/sample-level organelle matrix |
| `functional_irreversibility_index.tsv` | Per-cell composite proxy for adaptive/transition/resistant fixation (optional) |
| `comparison.json` | A→B deltas (optional) |
| `interpretation.json` | Rule-based interpretation signals |
| `report.html` | Deterministic one-page static report |
| `llm_report.md` | English LLM-friendly narrative summary |
| `fii_landscape.html` | Offline interactive FII landscape report across samples (post-processing) |
| `fii_landscape_index.json` | Report index with included samples |
| `fii_state_velocity.tsv` | First-difference dynamics across ordered samples |
| `fii_state_acceleration.tsv` | Second-difference dynamics across ordered samples |
| `phase_portrait_points.tsv` | Sample-level phase-space points (`FII`, velocity, acceleration) |
| `phase_portrait_bins.tsv` | Binned phase projections with regime/low-confidence fractions |
| `phase_portrait_summary.json` | Phase portrait metadata, ranges, points, and aggregated bins |
| `early_warning_flags.tsv` | Deterministic pre-transition warning flags |
| `sample_decision.tsv` | Compact per-sample decision tier, confidence, and drivers |
| `sample_decision.json` | Machine-readable decision payload with thresholds |
| `decision_timeline.tsv` | Ordered decision trajectory across samples |
| `decision_stability.tsv` | Stability/volatility metrics summary |
| `decision_stability.json` | Machine-readable timeline + metrics |

## Determinism & Reproducibility
- Identical input artifacts produce byte-identical outputs.
- Ordering is stable (organelles, axes, cells, signals, maps).
- Quantile selection and delta summaries are deterministic.
- Report rendering is static and deterministic (no client-side computation).

## Non-Goals
- No ML-based modeling.
- No clinical recommendations.
- No interactive dashboards.

## Functional Irreversibility Index (FII)
FII is an additive proxy-only composite built from:
- `mitochondrial_stress_adaptation_score` (kira-mitoqc)
- `translation_commitment_score` (kira-riboqc)
- `splice_irreversibility_index` (kira-spliceqc)

Formula:

```text
FII = w_mito * mito + w_translation * translation + w_splice * splice
```

All inputs are clamped to `[0,1]`; output is clamped to `[0,1]`.
Regimes:
- `Adaptive`: `FII < 0.33`
- `Transition`: `0.33 <= FII < 0.66`
- `Resistant`: `FII >= 0.66`

If one or more inputs are missing for a cell, `fii_low_confidence=true` is emitted; aggregation continues deterministically.

Limitation:
Functional Irreversibility Index (FII) is a composite proxy summarizing multi-organelle functional fixation. It does not represent true biochemical irreversibility or direct temporal ordering.

Interpretation hints:
- Adaptive peak with a small Resistant tail suggests predominantly reversible/adaptive state.
- Transition broadening suggests heterogeneous intermediate fixation.
- Strong FII-component coupling heatmaps indicate which subsystem most strongly tracks irreversible shift.

The FII landscape report is comparative and descriptive, not causal inference.

State Velocity and State Acceleration are comparative metrics across ordered samples. They do not represent true temporal derivatives, RNA velocity, or single-cell trajectories.

Phase portrait is a comparative state-space summary across ordered conditions. It does not establish causality and does not represent single-cell trajectories.

Decision tiers represent functional risk states inferred from analytic summaries. They are not clinical predictions or treatment recommendations.

Decision timeline interpretation:
- low volatility + high persistence => stable regime assignment
- early flip into high-risk tiers with high confidence => early functional commitment signal
- high volatility + low confidence => unstable/noisy decision trajectory

Decision stability does not imply clinical benefit; it reflects consistency of functional classification across ordered samples.
