# Formal Systems Model Export (Stage 10)

`kira-organelle` can export a deterministic formal systems model with:

- `--export-systems-model /path/to/systems_model.json`

The export is derived from already-computed systems metrics and contains no fitted or stochastic parameters.

## Export Contract

Top-level keys:

1. `model_version`
2. `normalization`
3. `composites_version`
4. `state_graph`
5. `ode_proxy`
6. `state_transition_matrix`
7. `basin_connectivity`

No timestamps are included.

## Deterministic State Graph (DSG)

- Nodes:
  - all basin IDs
  - all deterministic regime classes
- Basin edges:
  - connect basin pairs with centroid distance `< median(all basin-pair distances)`
  - edge weight: `1 / (1 + distance)`
- Regime edges:
  - weight from empirical nearest-neighbor regime transitions
- Edge ordering:
  - sorted by `(from, to)`

## ODE-like Proxy Model

For systems vector `F = [StressVector, CPI, AFS, IMSC]`:

- `mu`: global mean of `F`
- `Sigma`: global covariance of `F`
- constants: `lambda = 0.5`, `beta = 0.1`
- flow matrix: `A = -lambda*I + beta*Sigma`
- exported dynamics form: `F' = A(F - mu)`

No trajectory simulation is performed in this stage.

## State Transition Matrix (STM)

- Regime order is fixed and deterministic:
  - `BaselineCompensated`
  - `ImmuneSuppressiveEmbedded`
  - `AdaptiveHighStress`
  - `HighStress_FailedCompensation`
  - `SystemicCollapseRisk`
- Transition probabilities from nearest-neighbor cross-links:
  - `P_ij = count(Ri -> Rj) / total_outgoing(Ri)`
- Rows are stochastic; empty rows degrade to identity.

## Basin Connectivity Map

For each basin:

- nearest basin neighbor by centroid distance
- deterministic tie-break on basin ID

Output shape:

- `basin_connectivity[basin_id] = [{ "neighbor": "...", "distance": ... }]`

or empty list if no neighbor exists.
