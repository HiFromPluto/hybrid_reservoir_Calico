# BSimReservoirPlanStage6 mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanStage5`. Stage 5 is
not modified and remains **PASS**. Stage 4 is not rewritten and remains
**FAIL**. The 2x `Total_Deaths` fold is not weakened. Clamp, acid source
rate, `K_MAX`, `GROWTH_RATE`, and Stage 3B receiver `K=1.6`, `n=2`,
`tau=15` are not retuned. Glucose, Monod, Danino, ODE-gated vesicle ACs,
and Stage 7 AC layout are not imported. Stage 7 is not started.

## Three arms

Identical dish, timing, and seeds `101/202/303` except as named.

1. **Brownian null.** `BSimParticle` drift and diffusion only. No sensing,
   growth, division, death, communication, or `L`. AHL and acid sources
   still pulse (identical dish) but are not readout features. Ridge is
   spatial density `Den_*` on the 20x10 grid.
2. **Biology silent.** Full Stage 5 cells. `AHL source rate=0`,
   `acid source rate=0`, `warmup.ahl.input=0` (also forced in Java when
   `arm=silent`). Measures intrinsic drift of R, L, and clamp deaths.
3. **Biology driven.** Full Stage 5 cells. AHL carries NARMA-10. Acid held
   at 0.5.

## Frozen readout

Silent and driven ridge features, and only these:

- voxel `Receiver_R_*` / Mean_q (20x10)
- voxel `Lum_Mean_*` / Mean_L (20x10)
- coarse `Input_Driven_Death_*` (4x2)

Window-mean AHL, occupancy, pH, Births, Total_Deaths, Clamp_Deaths, OOB,
Population, Lum_Sum, and voxel AHL/pH/Den/occupancy/Lum_Sum stay out of
the ridge. Total_Deaths is not a Stage 6 feature.

Optional diagnostic: driven voxel `AHL_uM_*` field-only ridge. Not a
substitute for the Brownian gate.

## Protocol extension

See `PROTOCOL.md`. `dt=0.05`, 18000 s warmup, 300 s windows, 75 s pulses,
whole-window sampling. `num.windows=200` with washout 40 / train 110 /
test 50. Ridge grid selected on training data only.

## Gates

All required. NARMA-10 is the gate; memory capacity is diagnostic.

1. Driven NARMA-10 NRMSE lower than Brownian by a clear margin across 3
   seeds (non-overlapping mean ± s.e., and every seed).
2. Silent not ≈ driven. Silent NRMSE worse than driven.
3. Driven uses only the frozen Stage 5 analysis channels. Print the list.
4. CSV 200/3200/3200 rectangular for every arm.
5. Stage 4 remains FAIL. Stage 5 remains PASS.

If this fails, stop. Do not add features, retune, or start Stage 7.
