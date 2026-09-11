# E5.1 — isolated one-AHL horizon-1 CRP scout (Task B only)

Frozen: 2026-08-17, before occupancy or NRMSE.

Track E5 living scout on the claim dish. It does not rewrite Track B,
E4, Narma10b, Waveform, or Lorenz. `GATE_EVIDENCE.md` files are not
edited. `K`, `n`, `tau_R`, `tau_L`, clamp, mortality, flow, and layout
are not retuned. Claim dish stays CENTER / `FLOW=0`.

Task A / `Y` classification is not run. Glucose, Danino, vesicles,
acid-as-lactate, Track B 5-AC sites, sequential batching, C1,
Lorenz 202/333, Stage99, and E4 retunes are not started.

`COHORT_STATUS` remains `SYNTHETIC_COHORT`. Not a clinical dataset.

E5.0 already authorized this scout (`PROTOCOL_FROZEN_NO_BSIM`, Task B
survived). Do not apply a new 0.35 void rule here.

## Scientific question

On isolated claim-dish runs, does the living 408-D readout forecast
next-window CRP better than Brownian, silent, and field, on the
patient-pooled NRMSE gate?

Biology is not required to beat 5-channel delay-10. That baseline sees
four channels the dish does not get. Fair input ceilings: persistence
and CRP-only delay-10, same patient-pooled metric, delays never
crossing patient bounds. Concat NRMSE is a labelled diagnostic only.

## Frozen scout (do not redraw)

From `examples/HybridDish/biomarker_protocol_audit/results/RESET_ISOLATED_DESIGN.md`.

| Arm | n | Patient indices (0-based) | SHA-256 of concatenated CSV |
|---|---|---|---|
| development | 12 | `1,4,6,7,9,13,16,18,20,22,24,25` | `83acc7102e532112c961bb45a1327504290ffddea1c06e659bfb5e68be032915` |
| confirmation | 8 | `8,12,14,15,17,19,21,23` | `c90cca8f5627c737b38024b4863415e46a384d8715e6d06a5823e8cb93ba3787` |
| seed | 111 only | — | — |

Do not overwrite the E5.0 split. Seeds 222/333 are not in this prompt.

Vendored mat:
`examples/HybridDish/external_analysis/synthetic_biomarker_BSim_inputs.mat`

Target is `Xnorm[p, t+1, 0]` (channel 0 = CRP, min-max units). Do not
use `Y`.

## Frozen chemical map (lock before occupancy or NRMSE)

One-AHL only:

```
u[n] = 0.5 * Xnorm[p, n, 0]     n = 0 .. 47
J(t) = Jmax * u[n] during the 75 s pulse, 0 after
Jmax = 1.28e8 molecules/s
AHL at (500, 250, 5)
acid at (300, 375, 5) held 0.5 in analysis windows, off in warmup
attractant / repellent silent
```

No second chemical. No glucose. No Track B sites.

Warmup (isolated-patient, not NARMA 0.5):

```
18000 s at u = 0.5 * Xnorm[p, 0, 0] = u[0]
```

so `L` starts at that patient’s baseline, not a foreign operating
point. Then 48 windows of 300 s, pulse 75 s, sample every 20 s,
`dt=0.05`.

## Claim dish (copied from Narma10b — do not change)

- `1000 × 500 × 10` µm, `FLOW_SPEED=0`, `NO_FLUX`
- Field `50×25×1`, readout `20×10` + `4×2`
- Receiver `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL `1.28e8` at `(500, 250, 5)`. Acid `2e11` at `(300, 375, 5)`
- Attractant AC0/AC2 stay at `u=0`
- Clamp `K=2000`, `INITIAL_POP=1800`
- `K_MAX=0.002`, `GROWTH_RATE` compile-time `4π/1800`
- Official ridge: 408 (`Receiver_R_*` 20×10, `Lum_Mean_*` 20×10,
  `Input_Driven_Death_*` 4×2)
- One-way `addQuantity`. No vesicles, no glucose, no Danino

CSV: 48 windows / 768 samples / 768 voxel rows (15 samples/window).
Last sample prefix frozen after smoke: `47;15;299.95`.

## Module 0 — inputs and CRP-only ceilings (no BSim)

Write one 48-line AHL file per scout patient from `Xnorm[:,:,0]`.
Write matching horizon-1 target files (47 values: windows 0→1 .. 46→47).
SHA-256 each file. Verify scout index hashes.

Compute, patient-pooled, delays not crossing patients, **before Java**:

- persistence `t=0..46`
- persistence aligned `t=9..46`
- CRP-only delay-10 (10-D of channel 0, `t=9..46`)
- 5-channel delay-10 on the same 12+8 patients (analysis ceiling)

Lambda for delay-10: leave-one-patient-out on the 12 development
patients, grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`, intercept
unregularized, population standardize on the 11, then refit 12, score
8 once. Confirmation 8 last.

Report in `results/e51_input_ceilings.md`.

## Module 1 — BSim

Smoke (not evidence): one development patient. Print `u[0]`, u map,
CENTER, FLOW=0, 48 windows, finite non-negative AHL.

New BSim:

| Run | n | Notes |
|---|---|---|
| development driven | 12 | seed 111, isolated, patient-specific `u` and warmup `u[0]` |
| confirmation driven | 8 | seed 111, isolated; skip only if development occupancy is DEAD on ≥ 3 of 12. Do not skip because persistence wins |
| silent | 1 | sources off, warmup `u=0`, 48 zero windows, seed 111. Reuse for every patient |
| brownian | 1 | particles do not sense; AHL off; acid held 0.5; 48 windows; seed 111. Reuse `Den` for every patient |

Do not sequential-batch patients.

Occupancy DEAD per patient if `mean_R < 0.05` or `|r(mean_R, u)| < 0.5`.
Keep DEAD rows. Do not raise `Jmax`.

## Module 2 — analysis

Official pairs: `t = 9 .. 46` (38 per patient), aligned to delay-10.
Also report `t = 0 .. 46` as a diagnostic.

Ridge: 408 features, same closed ridge as Narma10b (population
standardize on training patients only; intercept unregularized;
lambda grid as E5.0). Lambda by leave-one-patient-out patient-pooled
NRMSE on the 12 development patients; ties take larger lambda.
Refit all 12. Score confirmation 8 once. Do not refit on 8.

Metric: NRMSE / that patient’s target std, then mean over patients
in the scored set. Concat NRMSE labelled diagnostic only.

Readouts:

- F408 driven
- field-only from driven voxels (`AHL_uM` 20×10)
- occupancy-masked RL surrogate
- reused Brownian `Den`
- reused silent F408 (same silent vector vs each patient’s target
  is legal; it should sit near train-intercept)
- persistence, CRP-only delay-10, 5-channel delay-10
  (from Module 0; 5-ch is ceiling, not a living-layer bar)

## Claims (report, do not retune)

- **System:** confirmation F408 vs Brownian and silent
- **Living-layer:** confirmation F408 vs field
- **Ceiling:** not required to beat 5-ch delay-10
- **Fair input:** report vs persistence and CRP-only delay-10
- **Clinical:** FORBIDDEN. `SYNTHETIC_COHORT` stays

`NO_STORY_MOVE` if `|F408 − field| < 0.03` and the same sign pattern
as “carrier already has the task”. Do not promote a new claim dish.

## Forbidden

- Task A / `Y` classification
- glucose, acid-as-lactate, Track B layout, second AC
- sequential 100-patient dish
- concat NRMSE as the gate
- fitting `K`, `n`, `tau`, `Jmax` to NRMSE
- calling this clinical
- overwriting E5.0 split or scout hashes
- rewriting Track B FAIL

## Deliverables

- `PROTOCOL.md` (this file)
- hashed input/target files
- `results/e51_input_ceilings.md`
- `results/E5_1_SCOUT.md`
- `results/e51_scout.csv` / `.json`
- `check_e51.py`

Checker: `python examples/BSimReservoirPlanE51CRP/check_e51.py`
from the repository root.
