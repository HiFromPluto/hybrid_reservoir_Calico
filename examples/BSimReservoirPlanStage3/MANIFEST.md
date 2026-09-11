# BSimReservoirPlanStage3 mechanism manifest

## Decision and scope

This package is a copy-derived, isolated Plan Stage 3 scaffold. The frozen
`BSimReservoirPlan` Java, configuration, and model files are not modified.
The center AC1 position `(500,250,5)`, formerly the repellent source, is
repurposed as a direct graded AHL source. Its repellent field is NULL: there is
no repellent source, chemotactic sign flip, or repellent toxicity. AC0 and AC2
remain attractant-source architecture at their original positions, but their
Plan Stage 3 calibration inputs are identically zero.

## Preserved mechanisms

- The continuous four-state Danino LuxI/AHL/AiiA/LuxR-AHL ODE and all kinetic
  constants are unchanged.
- Stage 2 timing, 50x25x1 chemical grid, readout grids, 18,000 s warmup,
  population clamp, growth law, and AHL field decay `0.0033 /s` are retained.
- `q = LA^2/(Kmla^2 + LA^2)` is a continuous diagnostic of the modeled
  LuxR-AHL response. `q > 0.5` is only a reporting threshold.

## Dimensional audit

`MOLECULES_PER_UM3_PER_UM = 602.2`. Intracellular AHL `y[1]` is uM and the
BSim field is molecules/um3:

1. `diffConc = y[1] * 602.2 - fieldConc` (molecules/um3).
2. `exchange_uM = diffConc / 602.2` (uM), used by the intracellular ODE.
3. Field exchange is `diffConc * CELL_WALL_DIFF * dt * cellVolume`
   (molecules).
4. Reported extracellular uM is `fieldConc / 602.2`.

The direct source follows the simple Stage10 `addQuantity` architecture:
`field.ahl.source.rate * graded_input * DT` molecules per active timestep.

## Explicit exclusions

There is no ArtificialCell, `gammaSym`, `CS_max`, refill, glucose, Monod,
nutrient field, `ALPHA_AHL`, `AHL_CONSUMPTION`, `K_QS_personal`,
`reservoir_new` Hill receiver, or luminescence. These are not implicit zeros;
they are absent from this stage.

## Population interpretation and stochasticity

For the preserved homeostatic removal law, the mean equilibrium is
`N* = K(1 + log2(ln 2))`, approximately `0.4712 K` and approximately 943 cells
for `K=2000`. This is not X-wall washout. X-boundary removals can occur only
after cell motion crosses the wall; there is no imposed flow (`FLOW_SPEED=0`).

`rng.seed` controls the explicit Plan Stage 3 `Random`, but BSim-internal
stochastic paths are not proven to share it. Runs must therefore be labeled
`stochastic_replicate`, not deterministic seed replicates.

## Output lifecycle

One `OutputOwner` owns all CSV writers. It closes each writer once after
simulation and no writes occur after closure. `window_summary.csv` has exactly
one data row per completed window. Its distribution and response statistics
pool the 16 predeclared samples spanning each complete 300 s window; they are
not endpoint snapshots. The analysis script validates expected row counts and
rectangular column counts and records them in its evidence JSON.
Calibration configs disable whole-window sample and voxel files; selected
production replicate configs retain both.

## Calibration outcome

The initial exploratory 0.25M--64M molecule/s runs exposed the endpoint-summary
defect and are retained only as pre-fix diagnostics. After correcting the
summary accumulator, the 128M/256M comparison showed that 128M could meet the
mid-input occupancy target but did not produce separation or correlation. That
comparison alone was not broad enough to establish that source-strength tuning
in general must fail.

The final attempt matched warmup exposure to measurement exactly:
`warmup.ahl.input=0.5` was applied for 75 s in every 300 s warmup cycle, not
continuously. It then reran corrected 32M, 64M, and 128M molecule/s sweeps and
tested 128M in three stochastic runs using three independently ordered,
balanced 40-window input permutations. Each permutation contains eight
instances of every input level and has near-zero correlation with window
index.

The direct source remains responsive, but the receiver does not: across the
three final runs, mean zero-lag `r(input, extracellular AHL)=0.891`, while
`r(input, mean q)=0.00294`. Mean `r(window index, mean q)=0.952`, showing that
the long activation transient remains even under exposure-matched
preconditioning. Mid-input `f(q>0.5)=0.00137`, so occupancy, separation, and
correlation all fail. These matched results support the narrower conclusion
that the preserved Danino receiver cannot resolve the specified 300 s inputs
under the tested protocol without changing kinetics or protocol.

`results/stage3_gate_evidence.json` is the final numeric source of truth.
Stage 3 is not closed and Stage 4 must not begin.
