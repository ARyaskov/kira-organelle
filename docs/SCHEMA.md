# Output Schema Contract

## integration/metrics.tsv

Canonical column groups are defined in `src/registry/output_schema.rs`:

- `BASE_COLUMNS`
- `COMPOSITE_COLUMNS`
- `AGG_COLUMNS`
- `FRAGILITY_COLUMNS`
- `LANDSCAPE_COLUMNS`

Columns are emitted in that fixed group order. Runtime may only omit whole optional groups; it never reorders columns.

Current canonical order when all groups are available:

1. `cell_id`
2. `cluster`
3. `StressVector`
4. `CompensationDeficit`
5. `CPI`
6. `AFS`
7. `IMSC`
8. `RegimeClass`
9. `RareState`
10. `MahalanobisDistance`
11. `DominantAxis`
12. `DominantAxisSensitivity`
13. `Potential`
14. `StabilityGradient`
15. `TPI_landscape`
16. `BasinId`
17. `TransitionCandidate`
18. `DominantVulnerableAxis`
19. `LSI`
20. `MaxTrajectory`
21. `BEE`

## integration/summary.json

Top-level key order is stable:

1. `schema`
2. `model`
3. `entity`
4. `timepoint`
5. `axes`
6. `collapse_proximity_baseline`
7. `collapse_proximity_peak`
8. `expression_available`
9. `expression_source`
10. `source_experiment` (optional)
11. `systems_state_model` (optional)
12. `ingestion_diagnostics`
13. `global_normalization`
14. `cross_sample_analysis` (optional; multi-sample integration mode)

`systems_state_model` includes:

- `global_stats`
- `cluster_stats`
- `regime_distribution`
- `global_cpi_p50`, `global_cpi_p90`
- `global_stress_p50`, `global_stress_p90`
- `rare_fraction_global`
- `fragility`
- `landscape`
- `therapeutic_projection`
- `dynamic_stability`
- `validation`
- `compatibility`
- `missing_metric_counts`

`cross_sample_analysis` includes:

- `samples`: ordered sample IDs in integration input order
- `regime_js_matrix`: symmetric `N x N` matrix
- `basin_overlap_matrix`: symmetric `N x N` matrix
- `dsp_vectors`: pair-keyed (`sampleA__vs__sampleB`) deterministic DSP vectors and norms
- `sdi_matrix`: symmetric `N x N` matrix

## integration/cross_sample_dsp.tsv (optional)

Written only in multi-sample integration mode.

Columns:

1. `sample_a`
2. `sample_b`
3. `axis`
4. `value`

Rows are emitted in deterministic pair and axis order.

## Compatibility Block

`systems_state_model.compatibility`:

- `missing_metrics`: sensitivity axes with `n_valid == 0`
- `degraded_axes`: sensitivity axes with `n_valid > 0` and `reliable == false`

Older upstream tables are valid input; missing metrics degrade to `NaN` and composites fall back deterministically.
