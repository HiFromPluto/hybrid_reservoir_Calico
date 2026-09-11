# BSimReservoirPlanStage4 mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanStage3B` and adds an
input-driven death channel. Stage 3B is not modified. The archived Stage 3
Danino package is not reopened. Stage 5 is not started.

The frozen Stage 3B receiver remains:

`dR/dt = (C^2 / (K^2 + C^2) - R) / tau`

with `K=1.6 uM`, `n=2`, `tau=15 s`, analytic `t95=44.936 s`. Those values are
compile-time constants and are not retuned. AC1 remains the direct AHL source.
AC0/AC2 remain silent attractant architecture. `GROWTH_RATE` remains the
compile-time Stage 3B constant. pH taxis and pH-dependent growth modulation are
not implemented in Stage 4 so the death gate cannot displace the receiver.

## Death mechanism

A dedicated mixed-acid field is added at `(300, 375, 5)`, cross-stream from
center AC1. This is not AC0, AC1, or AC2.

The field stores acid load in molecules/um3. Physical conversion:

`acid_mM = concentration / (602.2 * 1000)`

`pH = clamp(7.1 - acid_mM / 2.0, 1, 14)`

`PH_BASE=7.1` is the Ingraham & Marr 1996 midpoint of pH 6.8-7.4.
`BUFFER_CAPACITY=2.0 mM/pH` is an assumed weakly buffered medium parameter, not
a Small/Castanie-Cornet constant.

Diffusivity is Option B buffered effective D: `200 um2/s`, the midpoint of the
100-500 um2/s buffered range in `pH_unified_mechanism.md`. This is not bulk
proton diffusivity. Decay is `0.0067 /s` so the acid field does not ratchet
across successive high-input windows.

Acid has two sources:

1. Bacterial mixed-acid production: `1e6` molecules/s/cell. This is a magnitude
   placeholder so the field is population-coupled; it is not a converted
   mmol/gDW/h literature rate. At dish scale it keeps source-off pH near 7.
2. An independently driven point source at `(300, 375, 5)` with configurable
   rate, default `2.0e11` molecules/s, pulsed on the 75/300 s Stage 3B schedule.
   The acid source is off during warmup so the clamp plateau is preserved.
   A first calibration at `6.4e11` molecules/s and `K_MAX=0.003 /s` drove
   dish-mean pH through the kill range and extinguished the population.

Kill is additive to the LacOperon homeostatic clamp. Independent Bernoulli
trials are drawn each tick. If both fire, the death is logged as input-driven
so the clamp cannot mask the channel. X-wall OOB removals are a third cause.

Acid-arm Hill, anchored on Small et al. 1994 (growth limit pH 4.5) and
Castanie-Cornet et al. 1999 (collapse pH 2.5):

`k_acid = K_MAX * (max(0, 4.5-pH)^2) / ((4.5-3.75)^2 + max(0, 4.5-pH)^2)`

Half-max is at pH 3.75, the midpoint of 2.5-4.5. `K_MAX=0.002 /s` is
calibrated, not a literature constant. An alkaline arm is present but unused
because the field only adds acid.

## Protocol

Same Stage 3B timing: 18000 s warmup, 40 x 300 s windows, 75 s pulses,
whole-window sampling every 20 s. AHL keeps its three original balanced
permutations. Acid has its own three balanced permutations generated with seed
`20260815`; these are not Stage 3B holdout sequences. Fresh Java seeds
`101/202/303` plus source-off seed `404` with `field.acid.source.rate=0`.

## Gates

All required:

1. `|r(Acid_Input, Input_Driven_Deaths)| > 0.5` in every one of >=3 runs.
2. Mean `Total_Deaths` at least 2x the source-off / clamp-only baseline.
3. Source-off control has no input-driven deaths.
4. CSV rectangularity and 40/640/640 row counts.
5. Frozen Stage 3B receiver still passes on the same protocol.

## Final outcome

Overall **FAIL**. The death channel exists and correlates; the required
2x total-death fold does not.

Production holdouts (seeds 101/202/303) plus source-off seed 404:

- `|r(Acid_Input, Input_Driven_Deaths)| = 0.867 / 0.905 / 0.895`
- Mean `Total_Deaths` fold vs silent = `0.925` (FAIL; 100.8 vs 109.0)
- Source-off input-driven deaths = 0
- CSV 40/640/640 rectangular: PASS
- Frozen receiver `r = 0.786 / 0.818 / 0.817`, mid occupancy `0.336`,
  ordered `0.254 < 0.431 < 0.561`

A stronger first calibration extinguished the population and is archived
under `results/stage4_cal1_source6.4e11/`. See `results/GATE_EVIDENCE.md`.

CSV completeness is the on-disk 40/640/640 files, not the truncated BSim
timestep stdout. Those runs were not rerun for log truncation.

The 2x `Total_Deaths` fold is kept as written. It is not replaced by
input-driven-death correlation, and the clamp is not disabled after the
fact. Stage 4 is a documented FAIL. Stage 5 remains blocked.
