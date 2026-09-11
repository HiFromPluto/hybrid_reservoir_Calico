# Track A2 protocol (frozen before NRMSE)

This is Track A2 on the frozen Plan Stage 6 hybrid dish: subsampled
Lorenz’63 with SKIP=50. Mechanisms are not retuned. BenchA Lorenz’63
Overall: FAIL is not rewritten. BenchA Mackey–Glass Overall: PASS is
not rerun. Stage 6 remains the NARMA-10 provenance checkpoint.
Track B is not started.

Same cells, same 408-feature ridge, same 200-window budget. The only
change is the Lorenz clock: one window = 50 RK4 steps = Δt_sample=1.0.
Drive x, predict y[n+1] at the next sample. y and z are not injected
into any field. SKIP is frozen at 50. If this FAILs, stop. Do not then
try skip=25, predict z, add ACs, or add ridge features.

## Dish, timing, seeds (copied from Stage 6 / BenchA — do not change)

- `dt=0.05`, warmup `18000 s`, windows `300 s`, pulses `75 s`,
  sampling every `20 s`
- Field grid `50×25×1`, readout `20×10` and `4×2`
- Receiver `K=1.6`, `n=2`, `tau=15`
- `ALPHA_LUX=DELTA_LUX=1/1500`
- AHL source `1.28e8`, acid source `2e11`, `K_MAX=0.002`,
  `GROWTH_RATE` compile-time `4π/1800`
- Clamp `K=2000`, `INITIAL_POP=1800`
- `FLOW_SPEED=0`
- AHL AC at `(500, 250, 5)`. Acid at `(300, 375, 5)` held at `0.5`
  every analysis window, off during warmup
- Attractant AC0/AC2 stay at `u=0`
- Seeds `101 / 202 / 303` shared across Brownian, silent, and driven
- One-way `addQuantity` ACs. No vesicles, no glucose, no Danino, no
  `setGoal(glucose)`, no extra ACs

`num.windows=200` is unchanged.

## Input encoding (predeclared)

- AHL carries subsampled Lorenz `x`, affine-mapped onto `u ∈ [0, 0.5]`.
  Uniform`[0, 1]` is not used.
- Acid is held at `0.5` every analysis window (Stage 6 / BenchA copy:
  `input_acid_held05_200.txt`). Acid remains off during warmup.
- Silent arm: `field.ahl.source.rate=0`, `field.acid.source.rate=0`,
  `warmup.ahl.input=0`. Java also forces those three to zero when
  `arm=silent`. Silent CSVs do not depend on `u`; Stage 6 silent runs
  are reused and are not rerun.
- Brownian arm: same AHL/acid pulses as driven (identical dish).
  Particles do not sense the fields. Ridge uses density only.
  Brownian is rerun because the AHL file changed.

Do not change maps, ICs, transients, SKIP, or seeds after seeing NRMSE.

### Subsampled Lorenz’63

`σ=10`, `ρ=28`, `β=8/3`. RK4, `Δt_L=0.02`. IC `(x,y,z)=(1,1,1)`.
Discard 5000 transient RK4 steps. Then record 201 samples with
SKIP=50 RK4 steps between samples:

    Δt_sample = 50 × 0.02 = 1.0
    sample i is the state after 5000 + i×50 RK4 steps from the IC,
    i = 0 .. 200

BenchA’s every-step recording is not used.

Drive: `x[0..199]`. Affine to `u ∈ [0, 0.5]` using the min/max of
those 200 `x` values (not the 201st sample):

    xmin = -17.981997956481
    xmax =  16.013692471459
    u[n] = 0.5 * (x[n] − xmin) / (xmax − xmin)

Target for window `n`: `y[n+1]` (raw Lorenz `y` at the next SAMPLE,
not affine-mapped, not `y` after one RK4 step). Generator seed is
not used (deterministic ODE).

SHA-256 of the 12-decimal `u` sequence (Stage 6 recipe):
`35646e7b504940dbc859ec66f1f2c5b49017f25388c28f0d662ffc29a7f9a1e5`.

Sidecars: `input_ahl_lorenz_skip50.txt`, `input_acid_held05_200.txt`,
`lorenz_skip50_target.csv` with columns `n;u;x;y_next`.

## Washout / train / test

| Split | Windows | Count |
|---|---|---|
| Washout | 0 .. 39 | 40 |
| Train | 40 .. 149 | 110 |
| Test | 150 .. 199 | 50 |

After dropping washout, `X` has 160 rows. Fit inner train on rows
`0..87` (windows `40..127`). Pick `lambda` on rows `88..109` (windows
`128..149`) **only**. Do not use `X[88:]` — that was the Stage 6 leak
into test; it stays closed here. Then refit on all 110 train rows.
Test is never used to choose `lambda` or to standardize.

## Ridge rule (same for every arm and seed)

- Bias column of ones; intercept is not regularized
- Features standardized with training mean/std (zero-variance columns
  kept at 0)
- Grid: `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`
- Select the `lambda` with lowest validation NRMSE; ties take the
  larger `lambda`
- Independent `lambda` per arm×seed
- NRMSE = RMSE / pop-std(target on that split)

## Readout state

Per-window features, after whole-window sampling:

- State channels: mean of the 16 intra-window samples
- Count channel: last sample of the window (full-window accumulation)

| Arm | Features | Count |
|---|---|---|
| Brownian | `Den_*` 20×10 | 200 |
| Silent / driven | `Receiver_R_*` 20×10, `Lum_Mean_*` 20×10, `Input_Driven_Death_*` 4×2 | 408 |

Not in the ridge: window-mean AHL, occupancy, pH, Births, Total_Deaths,
Clamp_Deaths, OOB, Population, Lum_Sum, voxel AHL (biology ridge),
voxel pH/Den/occupancy/Lum_Sum.

Optional diagnostic, not a gate: driven-arm ridge on voxel `AHL_uM_*`
only. If field-only beats driven, record it. That is not a fail of
gate 1. Do not add features or retune to beat the delay line.

## Memory capacity

Diagnostic only. Reconstruct delayed AHL `u[n-k]` for `k=1..20` with
the same ridge rule and splits. `MC = sum_k max(0, R²_k)` on test.
`k=0` is printed, not added.

## Arms

1. **Driven.** AHL carries subsampled Lorenz `u`; acid held at 0.5.
2. **Silent.** AHL rate=0, acid rate=0, `warmup.ahl.input=0` (also
   forced in Java when `arm=silent`). Reuse Stage 6 silent CSVs.
3. **Brownian.** Same AHL/acid pulses as driven. Particles do not
   sense. Ridge is `Den` only.

## Gates (all required)

1. Driven test NRMSE < Brownian, non-overlapping mean±s.e. across 3
   seeds, AND driven < Brownian in every seed. Do not drop seed 101.
2. Silent not ≈ driven (overlapping mean±s.e. is FAIL). Silent NRMSE
   higher than driven.
3. Driven uses only the frozen 408 Stage 5/6 channels. Print the list.
4. CSV 200 / 3200 / 3200 rectangular; last sample `199;15;299.95`.
   Completeness is the CSV, not stdout.
5. BenchA Lorenz GATE_EVIDENCE still says Overall: FAIL. BenchA
   Mackey–Glass still Overall: PASS. Stage 6 NARMA-10 still PASS.
   Stage 6 and BenchA Java untouched. No glucose / Danino / extra ACs /
   att-rep / Track B.

If gate 1 fails, Overall: FAIL. Write it. Stop. Do not retune SKIP,
source rates, K, features, or the affine map. Do not switch to
`z[n+1]` or MSE.

## Runs

6 production BSim runs (driven+brownian × 3 seeds) plus the 3 copied
silent dirs.

    results/lorenz_skip50_{driven,brownian,silent}_seedXXX

Checker: `python examples/BSimReservoirPlanBenchA2/check_bencha2.py`
Evidence: `examples/BSimReservoirPlanBenchA2/results/GATE_EVIDENCE.md`

Do not start Track B.
