# Lorenz L3 protocol (frozen before new-k NRMSE)

Track E3 sample-clock screen, then at most one living scout on the
frozen HybridDish / BenchA claim dish.

This package does **not** rewrite BenchA Lorenz FAIL or BenchA2
SKIP=50 FAIL. Those `GATE_EVIDENCE.md` files are not edited. Kinetics,
K, n, tau, source rate, clamp, mortality, flow, and layout are not
retuned. Claim dish stays CENTER / FLOW=0. One AHL channel only.
y and z are not injected. Acid is not a Lorenz channel.

C1 remains DEFER. E0.3 remains NO_STORY_MOVE. Waveform1 FAIL and
Waveform2c field≥driven stay as written. Do not start C1, Waveform
confirmation, IPC BSim, AC implementation, biomarker BSim, Stage99,
or E0.3 confirmation.

## Scientific question

Is there a Lorenz sampling interval, mapped one sample per 300 s
window, where a legal linear map of the drive does not already solve
the target, but the task is not SKIP=50 destruction? If yes, does the
living 408-D readout beat persistence, the selected delay baseline,
Brownian, silent, and field on that clock?

Beating RAW delay-of-x is required for a living-layer claim on this
task. Beating field is the carrier diagnostic. Biology is not required
to beat an oracle that sees future Lorenz states.

## Frozen Lorenz integrator

σ=10, ρ=28, β=8/3, RK4, internal dt=0.02, IC (1,1,1), discard 5000
transient steps, then 201 recorded samples (200 windows). Affine
`u=0.5*(x−xmin)/(xmax−xmin)` on the 200 drive samples only,
`u∈[0,0.5]`. 75 s pulse / 300 s window, acid held 0.5. Same closed
ridge as NARMA (washout 40 / train 110 / test 50, inner val 22, lambda
grid, unreg intercept, train-only standardization, population-std
NRMSE).

Do not change σ, ρ, β, IC, or transient after seeing NRMSE.

Clock = skip k internal RK4 steps between recorded samples.
Δt_window = k × 0.02 Lorenz time units per dish window.

## Predeclared clocks

Historical anchors (do not regenerate; do not BSim):

| k | Δt | Provenance | Expected u SHA-256 |
|---|---|---|---|
| 1 | 0.02 | BenchA `input_ahl_lorenz200.txt` | `3bc69bbe32f6f2cf6a8638c4e2775133b4d27026a4223cbe783c057c19232f79` |
| 50 | 1.00 | BenchA2 `input_ahl_lorenz_skip50.txt` | `35646e7b504940dbc859ec66f1f2c5b49017f25388c28f0d662ffc29a7f9a1e5` |

New development clocks, generated and hashed before any baseline
ranking:

| k | Δt | Files |
|---|---|---|
| 5 | 0.10 | `input_ahl_lorenz_k5.txt`, `lorenz_target_k5.csv` |
| 10 | 0.20 | `input_ahl_lorenz_k10.txt`, `lorenz_target_k10.csv` |
| 20 | 0.40 | `input_ahl_lorenz_k20.txt`, `lorenz_target_k20.csv` |
| 40 | 0.80 | `input_ahl_lorenz_k40.txt`, `lorenz_target_k40.csv` |

Target CSV columns: `n;u;x;y_next;x_next`. Acid is copied from Stage 6
/ BenchA `input_acid_held05_200.txt`. Do not regenerate acid.

Hashes of the 12-decimal u sequences are written to
`results/lorenz_l3_input_hashes.csv` by `generate_lorenz_l3.py` **before**
`screen_lorenz_l3.py` ranks anything. The recorded hashes are copied
here after generation and before ranking:

| k | Role | u SHA-256 | xmin | xmax |
|---|---|---|---|---|
| 1 | anchor | `3bc69bbe32f6f2cf6a8638c4e2775133b4d27026a4223cbe783c057c19232f79` | -16.849584157354 | 14.975551014173 |
| 5 | development | `466e572660930ae36fe8f17912312569cb2c9d4af4b7feb1dfe6ba6b868c652f` | -16.849584157354 | 15.478808766204 |
| 10 | development | `69b861f256fb39527c037b994885d119ccc8ab2873a3d2bca59f7a8802c44aff` | -16.849584157354 | 15.478808766204 |
| 20 | development | `50f0ad290bc7cb4d4aeac166bdcf87a7cb325c153fc6f555d1704f0188ec7de2` | -16.021729621432 | 16.181393629053 |
| 40 | development | `4e706b8c524627cd63475348243fa59ce506125160eef0d23a9c6278353c4ba2` | -16.095795266623 | 16.155554217425 |
| 50 | anchor | `35646e7b504940dbc859ec66f1f2c5b49017f25388c28f0d662ffc29a7f9a1e5` | -17.981997956481 | 16.013692471459 |

Targets, both scored on every clock:

- AUTO_X: x[n] → x[n+1]
- CROSS_Y: x[n] → y[n+1]

## Closed ridge (unchanged)

| Split | Windows | Count |
|---|---|---|
| Washout | 0 .. 39 | 40 |
| Train | 40 .. 149 | 110 |
| Test | 150 .. 199 | 50 |

Inner train: windows 40..127. Validation: windows 128..149 only.
Lambda grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`. Ties take the
larger lambda. Unregularized intercept. Train-only standardization.
Zero-variance columns kept at 0. NRMSE = RMSE / pop-std of that
split's target. Test is never used to choose lambda, m, or
standardization.

## Module 1 — trivial-baseline screen (no Java, no BSim)

AUTO_X legal baselines:

- train-intercept
- persistence x[n]
- linear AR on `[x[n],…,x[n−m+1]]` for m in {1,3,5,7,10}
  (zero-pad before 0). Select m on validation only; ties take the
  smaller m.

CROSS_Y legal baselines:

- train-intercept
- persistence y[n] (teacher-forced; legal, usually weak)
- linear x-delay `[x[n],…,x[n−m+1]]` for the same m grid.
  Select m on validation only; ties take the smaller m.

Test-mean is an oracle reference, not ranked. Every development m is
reported. Only the validation-selected delay joins intercept and
persistence in the legal set.

Frozen labels, declared now, before looking at the new k numbers:

- **TRIVIAL:** best legal test NRMSE ≤ 0.30
- **DESTROYED:** best legal test NRMSE ≥ 0.95
- **SURVIVE:** 0.30 < best legal test NRMSE < 0.95

k=1 must come out TRIVIAL and k=50 DESTROYED if the zero-sim selected
delay test NRMSE reproduces within 1e-3:

| Clock | Task | Zero-sim selected delay test NRMSE |
|---|---|---|
| k=1 | AUTO_X | 0.0014 |
| k=1 | CROSS_Y | 0.0118 |
| k=50 | AUTO_X | 1.0047 |
| k=50 | CROSS_Y | 1.0019 |

If that sanity fails, stop and diagnose. Do not change thresholds.

If no new k SURVIVEs on either target:

```
L3_DISPOSITION: NO_INTERMEDIATE_CLOCK
```

Write `results/LORENZ_L3_CLOCK_SCREEN.md` and stop. No BSim.

If at least one (k, target) SURVIVEs, selection uses Module-1 numbers
only:

- consider SURVIVE pairs only
- choose the pair whose best legal *validation* NRMSE is largest
  (hardest remaining linear clock)
- at most one pair for living BSim
- do not select on a dish score; there is none yet
- tie-break, frozen: larger k, then CROSS_Y over AUTO_X

Write the chosen pair and the full clock table before any Java.

### Module-1 result (recorded before Java)

Anchor sanity versus zero-sim selected delay: **PASS** (all |Δ| < 1e-3).
k=1 TRIVIAL, k=50 DESTROYED.

SURVIVE pairs: (k=5, AUTO_X), (k=5, CROSS_Y), (k=10, AUTO_X).
k=10 CROSS_Y and all k≥20 pairs are DESTROYED.

Selected living scout pair, Module-1 validation only:

- **k=10, Δt=0.20, AUTO_X**
- best legal validation NRMSE = 0.8570 (LINEAR_AR m=10, λ=1e-6)
- best legal test NRMSE = 0.8233
- u SHA-256 `69b861f256fb39527c037b994885d119ccc8ab2873a3d2bca59f7a8802c44aff`

Full table: `results/LORENZ_L3_CLOCK_SCREEN.md`.
Do not switch k or target after seeing a dish NRMSE.

### Module-2 result (seed 101 only)

Occupancy **ALIVE**. Driven F408 test NRMSE `0.8256` does not beat
the legal AR `0.8233`. Living-layer **FAIL**. Brownian and silent
`1.0076`. Field `0.7973` (diagnostic; below driven). CSV 200 / 3200
/ 3200; last sample `199;15;299.95`. Seeds 202/303 were not run.
No retune. No new claim dish. Confirmation on a new Lorenz IC is a
later prompt.

## Module 2 — living scout (only the selected pair)

Package Java from the claim HybridDish / BenchA dish. Pulse protocol
stays 75/300. Do not add a continuous-hold L1 arm in this scout.

Dish, copied and frozen:

- `1000 × 500 × 10` µm, `dt=0.05`, `FLOW_SPEED=0`
- Field `50×25×1`, readout `20×10` + `4×2`
- `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL `1.28e8` at `(500, 250, 5)`. Acid `(300, 375, 5)` held 0.5 in
  analysis windows, off in warmup
- Attractants silent. Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, window `300 s`, pulse `75 s`, sample every `20 s`,
  200 windows
- Official 408. Brownian Den 20×10. Field-only driven AHL_uM 20×10
- Seeds `101 / 202 / 303`. 202/303 only after seed 101 living-layer pass

Reuse, do not rerun:

- 200-window silent seed 101 from Stage 6 / BenchA (sources off; same
  dish): `examples/BSimReservoirPlanStage6/results/stage6_silent_seed101`
- claim-dish occupancy sanity is not a Lorenz score
- do not reuse BenchA driven Lorenz voxels: those are k=1

Rerun Brownian seed 101 (new AHL file). Driven seed 101 only at first.
If 200-window silent exists for 202/303, reuse those Stage 6 silent
dirs; do not rerun silent.

Smoke: print dish, source, FLOW=0, finite AHL. Not evidence.

CSV 200 / 3200 / 3200; last sample `199;15;299.95`.

Occupancy **DEAD** if `mean_R < 0.05` or `|r(mean_R,u)| < 0.5`. If
DEAD, keep the row, do not raise rate or K, do not run 202/303.

Analysis on the selected target only:

- official F408
- field-only AHL
- FRL, R-only, L-only (descriptive)
- unmasked and occupancy-masked kinetic surrogates (frozen R/L ODEs,
  exact exponential ZOH-previous AHL, K=1.6, τ_R=15, τ_L=1500)
- the Module-1 legal baselines on this same target
- silent and Brownian

Living-layer claim on this scout (seed 101):

driven F408 test NRMSE < best legal baseline AND < Brownian AND
< silent.

If both driven and a comparator are ≥ 0.90, the margin must be larger
than 0.03. Field is diagnostic, not the living-layer definition.

If seed 101 is occupancy-DEAD or does not beat the best legal
baseline, stop. Do not run 202/303. Do not retune. Do not switch
target or k after seeing NRMSE.

If seed 101 beats the best legal baseline, Brownian, and silent: run
seeds 202 and 303 (driven+Brownian each; silent reused per seed if
200-window silent exists). No best-seed selection.

Do not freeze a new claim dish. Confirmation on a new Lorenz IC is a
later prompt.

## Commands

From the repo root:

```
python examples/BSimReservoirPlanLorenzL3/generate_lorenz_l3.py
python examples/BSimReservoirPlanLorenzL3/screen_lorenz_l3.py
```

Module 2, only after a SURVIVE pair is written into
`LORENZ_L3_CLOCK_SCREEN.md`:

```
python examples/BSimReservoirPlanLorenzL3/check_lorenz_l3.py
```

## Completion checks

- k=1 and k=50 were not rerun in BSim
- hashes recorded before baseline ranking
- TRIVIAL/DESTROYED/SURVIVE applied as frozen
- at most one living (k, target)
- no y/z/acid encoding
- no GATE_EVIDENCE edits
- 202/303 only after seed-101 living-layer pass
