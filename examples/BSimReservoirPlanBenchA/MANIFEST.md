# BSimReservoirPlanBenchA mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanStage6`. Stage 6
is not modified and remains the NARMA-10 provenance checkpoint
(**PASS**). Kinetics are not retuned. Clamp, acid source rate, AHL
source rate, `K_MAX`, `GROWTH_RATE`, and Stage 3B receiver `K=1.6`,
`n=2`, `tau=15` are not retuned. Glucose, Monod, Danino, ODE-gated
vesicle ACs, Stage 7 four-mode rank, Stage 9 product-bit, and Track B
(5 ACs, patient task, wider dish) are not imported.

Track A drives the same dish with Lorenz’63 and Mackey–Glass. New AHL
sequences and targets only. 408-feature ridge. 200-window budget.

## Three arms

Identical dish, timing, and seeds `101/202/303` except as named.

1. **Brownian null.** `BSimParticle` drift and diffusion only. No
   sensing, growth, division, death, communication, or `L`. AHL and
   acid sources still pulse (identical dish) but are not readout
   features. Ridge is spatial density `Den_*` on the 20×10 grid.
2. **Biology silent.** Full Stage 6 cells. `AHL source rate=0`,
   `acid source rate=0`, `warmup.ahl.input=0` (also forced in Java
   when `arm=silent`). Stage 6 silent CSVs are reused; silent is not
   rerun.
3. **Biology driven.** Full Stage 6 cells. AHL carries the task
   sequence. Acid held at 0.5.

## Frozen readout

Silent and driven ridge features, and only these:

- voxel `Receiver_R_*` / Mean_q (20×10)
- voxel `Lum_Mean_*` / Mean_L (20×10)
- coarse `Input_Driven_Death_*` (4×2)

Window-mean AHL, occupancy, pH, Births, Total_Deaths, Clamp_Deaths,
OOB, Population, Lum_Sum, and voxel AHL/pH/Den/occupancy/Lum_Sum stay
out of the ridge.

Optional diagnostic: driven voxel `AHL_uM_*` field-only ridge. Not a
substitute for the Brownian gate. Not a fail if biology loses to the
delay line.

## Protocol

See `PROTOCOL.md`. Frozen Stage 6 dish. Lorenz first, then
Mackey–Glass. Affine `u ∈ [0, 0.5]`. Washout 40 / train 110 / test 50.
Lambda selected on windows 128..149 only (Stage 6 val-slice leak
closed). Targets are raw `y[n+1]` (Lorenz) and `x[n+1]` (Mackey–Glass),
not affine `u`.

## Gates

All required, each task separately.

1. Driven test NRMSE lower than Brownian by a clear margin across 3
   seeds (non-overlapping mean ± s.e., and every seed).
2. Silent not ≈ driven. Silent NRMSE worse than driven.
3. Driven uses only the frozen Stage 5/6 analysis channels. Print the
   list.
4. CSV 200/3200/3200 rectangular; last sample `199;15;299.95`.
5. Stage 6 NARMA-10 remains PASS. Stage 6 Java untouched. No glucose /
   Danino / extra ACs / Track B.

If a task fails, record FAIL and stop retuning that task. Run the
other task. Do not start Track B.
