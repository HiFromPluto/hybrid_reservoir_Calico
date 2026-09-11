# Stage 6 protocol (frozen before NRMSE)

This is a Stage 6 protocol extension of the Stage 5 dish. Mechanisms are
frozen. Stage 4 remains FAIL. Stage 5 remains PASS. Stage 7 is not started.

## Dish, timing, seeds

Copied from Stage 5 and not retuned:

- `dt=0.05`, warmup `18000 s`, windows `300 s`, pulses `75 s`, whole-window
  sampling every `20 s`
- Field grid `50x25x1`, readout `20x10` and `4x2`
- Receiver `K=1.6`, `n=2`, `tau=15`
- `ALPHA_LUX=1/1500`, `DELTA_LUX=1/1500`
- Acid source rate `2e11`, `K_MAX=0.002`, `GROWTH_RATE` compile-time constant
- Clamp `K=2000` unchanged
- Seeds `101 / 202 / 303` shared across Brownian, silent, and driven

`num.windows=200` is the only protocol extension. 40 balanced holdout windows
are too short for NARMA-10.

## Input encoding (predeclared)

- AHL carries NARMA-10. `u[n] ~ Uniform[0, 0.5]` from `random.Random(20260814)`.
  SHA-256 of the 12-decimal sequence: `d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.
  That is the standard NARMA-10 range and a subset of `[0,1]`. Uniform`[0,1]`
  is not used because the Atiya/Parlos recurrence diverges there.
- Acid is option (a): held at `0.5` every analysis window so the death channel
  is present but is not a second free input. Acid remains off during warmup,
  as in Stage 5.
- Silent arm: `field.ahl.source.rate=0`, `field.acid.source.rate=0`,
  `warmup.ahl.input=0`. Java also forces those three to zero when `arm=silent`.
- Brownian arm: same AHL/acid pulses as driven (identical dish). Particles do
  not sense the fields. Ridge uses density only.

Do not change this encoding after seeing NRMSE.

## NARMA-10 target

```
y[0] = 0
y[n+1] = 0.3 y[n] + 0.05 y[n] * sum_{i=0..9} y[n-i]
         + 1.5 u[n-9] u[n] + 0.1
```

with `y[k]=0` and `u[k]=0` for `k<0`. Window `n` state, collected after pulse
`u[n]`, predicts `y[n+1]`. Sidecar: `narma10_target.csv`.

## Washout / train / test

| Split | Windows | Count |
|---|---|---|
| Washout | 0 .. 39 | 40 |
| Train | 40 .. 149 | 110 |
| Test | 150 .. 199 | 50 |

Lambda selection uses only training data: fit on windows 40 .. 127, pick
`lambda` on windows 128 .. 149, then refit on all 110 train windows.
Test windows are never used to choose `lambda` or to standardize.

## Ridge rule (same for every arm and seed)

- Bias column of ones; intercept is not regularized
- Features standardized with training mean/std (zero-variance columns kept at 0)
- Grid: `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`
- Select the `lambda` with lowest validation NRMSE; ties take the larger `lambda`
- NRMSE = RMSE / std(target on that split, population std)
- Independent `lambda` per arm×seed, identical procedure

## Readout state

Per-window features, after whole-window sampling:

- State channels: mean of the 16 intra-window samples
- Count channel: last sample of the window (full-window accumulation)

| Arm | Features | Count |
|---|---|---|
| Brownian | `Den_*` 20×10 | 200 |
| Silent / driven | `Receiver_R_*` 20×10, `Lum_Mean_*` 20×10, `Input_Driven_Death_*` 4×2 | 408 |

Not in the ridge: window-mean AHL, occupancy, pH, Births, Total_Deaths,
Clamp_Deaths, OOB, Population, Lum_Sum, voxel AHL/pH/occupancy/Den (biology
arms), voxel Lum_Sum.

Optional diagnostic, not a gate: driven-arm ridge on voxel `AHL_uM_*` only.

## Memory capacity

Diagnostic only. Reconstruct delayed AHL `u[n-k]` for `k=1..20` with the same
ridge rule and splits. `MC = sum_k max(0, R^2_k)` on test. `k=0` is reported
as an immediate-reconstruction diagnostic and is not added to MC. NARMA-10 is
the gate.

## Gates

1. Driven NARMA-10 test NRMSE is lower than Brownian with non-overlapping
   mean ± s.e. across the 3 seeds, and driven < Brownian in every seed.
2. Silent is not ≈ driven (overlapping mean ± s.e. is FAIL). Silent NRMSE is
   higher than driven (between Brownian and driven, or worse than both).
3. Driven uses only the frozen Stage 5 analysis channels. Print the list.
4. CSV rectangularity: 200 summary rows and 3200 sample/voxel rows per arm.
5. Stage 4 remains FAIL. Stage 5 remains PASS.

If the gate fails, stop. Do not add features, retune kinetics, or start Stage 7.
