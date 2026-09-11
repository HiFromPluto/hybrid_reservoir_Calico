# Track A protocol (frozen before NRMSE)

This is Track A on the frozen Plan Stage 6 hybrid dish. Mechanisms are
not retuned. Stage 6 remains the NARMA-10 provenance checkpoint
(`examples/BSimReservoirPlanStage6/`, Overall: PASS). Track B (5 ACs,
patient task, wider dish) is not started.

Same cells, same 408-feature ridge, same 200-window budget. New AHL
sequences and targets only.

## Dish, timing, seeds (copied from Stage 6 — do not change)

- `dt=0.05`, warmup `18000 s`, windows `300 s`, pulses `75 s`,
  whole-window sampling every `20 s`
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

`num.windows=200` is unchanged from Stage 6.

## Input encoding (predeclared)

- AHL carries the task sequence, affine-mapped onto `u ∈ [0, 0.5]`.
  Uniform`[0, 1]` is not used.
- Acid is held at `0.5` every analysis window (Stage 6 option (a) copy:
  `input_acid_held05_200.txt`). Acid remains off during warmup.
- Silent arm: `field.ahl.source.rate=0`, `field.acid.source.rate=0`,
  `warmup.ahl.input=0`. Java also forces those three to zero when
  `arm=silent`. Silent CSVs do not depend on `u`; Stage 6 silent runs
  are reused and are not rerun.
- Brownian arm: same AHL/acid pulses as driven (identical dish).
  Particles do not sense the fields. Ridge uses density only.
  Brownian is rerun because the AHL file changed.

Do not change maps, ICs, transients, or seeds after seeing NRMSE.

### Task 1 — Lorenz’63

`σ=10`, `ρ=28`, `β=8/3`. RK4, `Δt_L=0.02`. IC `(x,y,z)=(1,1,1)`.
Discard 5000 transient steps. Then record 201 samples (indices 0..200).

Drive: `x[0..199]`. Affine to `u ∈ [0, 0.5]` using the min/max of
those 200 `x` values (not the 201st sample):

    xmin = -16.849584157354
    xmax =  14.975551014173
    u[n] = 0.5 * (x[n] − xmin) / (xmax − xmin)

Target for window `n`: `y[n+1]` (raw Lorenz `y` at the next sample,
not affine-mapped). Generator seed is not used (deterministic ODE).

SHA-256 of the 12-decimal `u` sequence (Stage 6 recipe):
`3bc69bbe32f6f2cf6a8638c4e2775133b4d27026a4223cbe783c057c19232f79`.

Sidecars: `input_ahl_lorenz200.txt`, `input_acid_held05_200.txt`,
`lorenz_target.csv` with columns `n;u;x;y_next`.

### Task 2 — Mackey–Glass

Discrete map used in RC papers (`dt=1`, `τ=17`):

    x[t+1] = x[t] + 0.2 * x[t−17] / (1 + x[t−17]**10) − 0.1 * x[t]

For `t−17 < 0`, `x = 1.2`. Run 1000 transient steps from that history,
then record 201 samples.

Drive: `x[0..199]`. Affine to `u ∈ [0, 0.5]`:

    xmin = 0.411423008740
    xmax = 1.295296067959
    u[n] = 0.5 * (x[n] − xmin) / (xmax − xmin)

Target for window `n`: `x[n+1]` (raw MG state, not affine-mapped).

SHA-256 of the 12-decimal `u` sequence:
`e06810ea70f61bbe4baf61939f93b1e1cae9d92d9eee8d757865c4935dbac780`.

Sidecars: `input_ahl_mg200.txt`, `mg_target.csv` with columns
`n;u;x;x_next`. Acid file is the Stage 6 copy already used for Lorenz.

If Lorenz FAILs, Mackey–Glass is still run. Do not retune. Do not drop
a task.

## Washout / train / test

| Split | Windows | Count |
|---|---|---|
| Washout | 0 .. 39 | 40 |
| Train | 40 .. 149 | 110 |
| Test | 150 .. 199 | 50 |

After dropping washout, `X` has 160 rows. Fit inner train on rows
`0..87` (windows `40..127`). Pick `lambda` on rows `88..109` (windows
`128..149`) **only**. Do not use `X[88:]` — that was the Stage 6 leak
into test; it is closed here. Then refit on all 110 train rows. Test
is never used to choose `lambda` or to standardize.

## Ridge rule (same for every arm, seed, and task)

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
only. If driven does not beat field-only, record it. That is not a
fail of gate 1. Do not add features or retune to beat the delay line.

## Memory capacity

Diagnostic only. Reconstruct delayed AHL `u[n-k]` for `k=1..20` with
the same ridge rule and splits. `MC = sum_k max(0, R²_k)` on test.
`k=0` is printed, not added.

## Arms

1. **Driven.** AHL carries the task sequence; acid held at 0.5.
2. **Silent.** AHL rate=0, acid rate=0, `warmup.ahl.input=0` (also
   forced in Java when `arm=silent`). Reuse Stage 6 silent CSVs.
3. **Brownian.** Same AHL/acid pulses as driven. Particles do not
   sense. Ridge is `Den` only.

## Gates (each task separately)

All required.

1. Driven test NRMSE < Brownian, non-overlapping mean±s.e. across 3
   seeds, and driven < Brownian in every seed.
2. Silent not ≈ driven (overlapping mean±s.e. is FAIL). Silent NRMSE
   higher than driven.
3. Driven uses only the frozen 408 Stage 5/6 channels. Print the list.
4. CSV 200 / 3200 / 3200 rectangular; last sample `199;15;299.95`.
   Completeness is the CSV, not stdout.
5. Stage 6 NARMA-10 `GATE_EVIDENCE.md` still says Overall: PASS.
   Stage 6 Java untouched. No glucose / Danino / extra ACs / Track B.

If a task FAILs, write that section of `GATE_EVIDENCE.md` as FAIL and
stop retuning that task. Then run the other task.

## Runs

Per task, 6 production BSim runs (driven+brownian × 3 seeds) plus the
3 copied silent dirs. Lorenz first, then Mackey–Glass. Separate
`output.dir` trees:

    results/lorenz_{driven,brownian,silent}_seedXXX
    results/mg_{driven,brownian,silent}_seedXXX

Silent copies may be shared (identical files); still point the checker
at them.

Checker: `python examples/BSimReservoirPlanBenchA/check_bencha.py`
Evidence: `examples/BSimReservoirPlanBenchA/results/GATE_EVIDENCE.md`
One evidence file with two sections (Lorenz, Mackey–Glass), each with
Overall PASS/FAIL.

Do not start Track B.
