# BSimReservoirPlanBenchA2 mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanBenchA` Java onto
the Stage 6 dish. BenchA Lorenz’63 remains **FAIL** and is not
rewritten. BenchA Mackey–Glass remains **PASS** and is not rerun.
Stage 6 is not modified and remains the NARMA-10 provenance checkpoint
(**PASS**). Kinetics are not retuned. Clamp, acid source rate, AHL
source rate, `K_MAX`, `GROWTH_RATE`, and Stage 3B receiver `K=1.6`,
`n=2`, `tau=15` are not retuned. Glucose, Monod, Danino, ODE-gated
vesicle ACs, attractant/repellent encoding, Stage 7 four-mode rank,
Stage 9 product-bit, and Track B are not imported.

Track A2 drives the same dish with subsampled Lorenz’63 (SKIP=50,
Δt_sample=1.0). New AHL sequence and target only. 408-feature ridge.
200-window budget. Drive `x`, predict `y[n+1]` at the next sample.
`y` and `z` are not injected into any field.

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
3. **Biology driven.** Full Stage 6 cells. AHL carries subsampled
   Lorenz `u`. Acid held at 0.5.

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

See `PROTOCOL.md`. Frozen Stage 6 dish. SKIP=50 is frozen. Affine
`u ∈ [0, 0.5]`. Washout 40 / train 110 / test 50. Lambda selected on
windows 128..149 only (Stage 6 val-slice leak closed). Target is raw
`y[n+1]` at the next sample, not affine `u`.

## Gates

All required.

1. Driven test NRMSE lower than Brownian by a clear margin across 3
   seeds (non-overlapping mean ± s.e., and every seed). Do not drop
   seed 101.
2. Silent not ≈ driven. Silent NRMSE worse than driven.
3. Driven uses only the frozen Stage 5/6 analysis channels. Print the
   list.
4. CSV 200/3200/3200 rectangular; last sample `199;15;299.95`.
5. BenchA Lorenz remains FAIL. BenchA Mackey–Glass remains PASS.
   Stage 6 NARMA-10 remains PASS. Stage 6 and BenchA Java untouched.
   No glucose / Danino / extra ACs / att-rep / Track B.

If this fails, record FAIL and stop. Do not retune SKIP, kinetics, or
features. Do not start Track B.
