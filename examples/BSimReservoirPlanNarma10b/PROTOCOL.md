# NARMA-10b protocol (frozen before NRMSE)

Independent NARMA-10 replication on the frozen Plan Stage 6 hybrid dish.
New bacterial RNG seeds `111 / 222 / 333`. Mechanisms are not retuned.
The NARMA `u` sequence is copied from Stage 6 and is not regenerated.

Stage 6 NARMA-10 remains the provenance checkpoint
(`examples/BSimReservoirPlanStage6/`, Overall: PASS). If this package
FAILS, Stage 6 PASS is not rewritten; both results are reported. If it
PASSes, both stand. Track B, waveform, Lorenz, vesicles, glucose, and
Danino are not started. BenchA / A2 Java is not edited.

## Dish, timing (copied from Stage 6 — do not change)

- `dt=0.05`, warmup `18000 s`, windows `300 s`, pulses `75 s`,
  whole-window sampling every `20 s`
- Field grid `50×25×1`, readout `20×10` and `4×2`
- Receiver `K=1.6`, `n=2`, `tau=15`
- `ALPHA_LUX=DELTA_LUX=1/1500`
- AHL source `1.28e8` at `(500, 250, 5)`. Acid source `2e11` at
  `(300, 375, 5)` held at `0.5` every analysis window, off during warmup
- Attractant AC0/AC2 stay at `u=0`
- `K_MAX=0.002`, `GROWTH_RATE` compile-time `4π/1800`
- Clamp `K=2000`, `INITIAL_POP=1800`, `FLOW_SPEED=0`
- One-way `addQuantity` ACs. No vesicles, no glucose, no Danino, no
  `setGoal(glucose)`, no extra ACs

`num.windows=200` is unchanged from Stage 6.

## Seeds

`111 / 222 / 333` shared across Brownian, silent, and driven.
Do **not** reuse Stage 6 seeds `101 / 202 / 303`. A same-seed rerun is
not a replication. Silent is rerun: clamp and positions depend on RNG
even with sources off. Stage 6 silent CSVs are not copied.

## Input encoding (copied, not regenerated)

- AHL carries NARMA-10. `u[n] ~ Uniform[0, 0.5]` from
  `random.Random(20260814)`. SHA-256 of the 12-decimal sequence:
  `d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.
  Uniform`[0, 1]` is not used because the Atiya/Parlos recurrence
  diverges there.
- Acid is held at `0.5` every analysis window. Acid remains off during
  warmup.
- Silent arm: `field.ahl.source.rate=0`, `field.acid.source.rate=0`,
  `warmup.ahl.input=0`. Java also forces those three to zero when
  `arm=silent`.
- Brownian arm: same AHL/acid pulses as driven (identical dish).
  Particles do not sense the fields. Ridge uses density only.

Sidecars copied from Stage 6: `input_ahl_narma200.txt`,
`input_acid_held05_200.txt`, `narma10_target.csv`.

Do not change this encoding after seeing NRMSE.

## NARMA-10 target

```
y[0] = 0
y[n+1] = 0.3 y[n] + 0.05 y[n] * sum_{i=0..9} y[n-i]
         + 1.5 u[n-9] u[n] + 0.1
```

with `y[k]=0` and `u[k]=0` for `k<0`. Window `n` state, collected after
pulse `u[n]`, predicts `y[n+1]`. The checker recomputes the SHA-256 and
the recurrence and aborts on mismatch.

## Washout / train / test

| Split | Windows | Count |
|---|---|---|
| Washout | 0 .. 39 | 40 |
| Train | 40 .. 149 | 110 |
| Test | 150 .. 199 | 50 |

After dropping washout, `X` has 160 rows. Fit inner train on rows
`0..87` (windows `40..127`). Pick `lambda` on rows `88..109` (windows
`128..149`) **only**. Do not use `X[88:]` — that was the Stage 6 leak
into test; it is closed here (BenchA rule). Then refit on all 110 train
rows. Test is never used to choose `lambda` or to standardize.

## Ridge rule (same for every arm and seed)

- Bias column of ones; intercept is not regularized
- Features standardized with training mean/std (zero-variance columns
  kept at 0)
- Grid: `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`
- Select the `lambda` with lowest validation NRMSE; ties take the
  larger `lambda`
- Independent `lambda` per arm×seed
- NRMSE = RMSE / pop-std(target on that split)
- No extra ridge features: no `Den` on biology arms, no voxel `AHL_uM`
  in the biology ridge

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
`k=0` is printed, not added. NARMA-10 is the gate.

## Arms

1. **Driven.** AHL = frozen NARMA file; acid held at 0.5;
   `warmup.ahl.input=0.5`.
2. **Silent.** AHL rate=0, acid rate=0, `warmup.ahl.input=0` (also
   forced in Java when `arm=silent`). Rerun, do not copy Stage 6 silent.
3. **Brownian.** Same AHL/acid pulses as driven. Particles do not
   sense. Ridge is `Den` only.

## Gates (same claim as Stage 6)

All required.

1. Driven test NRMSE < Brownian, non-overlapping mean±s.e. across 3
   seeds, and driven < Brownian in every seed.
2. Silent not ≈ driven (overlapping mean±s.e. is FAIL). Silent NRMSE
   higher than driven.
3. Driven uses only the frozen 408 channels. Print the list.
4. CSV 200 / 3200 / 3200 rectangular; last sample `199;15;299.95`.
   Completeness is the CSV, not stdout.
5. Stage 6 NARMA-10 `GATE_EVIDENCE.md` still says Overall: PASS and was
   not edited. Stage 6 Java untouched. Waveform / Track B / BenchA / A2
   not edited. No glucose / Danino / extra ACs.

If gate 1 fails, Overall: FAIL for **this** package. Do not retune `K`,
sources, warmup, or `u`. Do not fall back to Stage 6 seeds. Write both
NRMSE tables (this replication vs Stage 6) in `GATE_EVIDENCE.md`.
Field-only is a diagnostic, not a kill gate.

## Runs

9 production BSim runs: driven, brownian, silent × `111/222/333`.

    results/narma10b_{driven,brownian,silent}_seedXXX

Checker: `python examples/BSimReservoirPlanNarma10b/check_narma10b.py`
Evidence: `examples/BSimReservoirPlanNarma10b/results/GATE_EVIDENCE.md`

Include a side-by-side with Stage 6 means (driven `0.9306`, Brownian
`1.1625`, silent `1.1622`, field-only `1.0126`) labelled
“Stage 6, not this run.”

Do not start another task after this package.
