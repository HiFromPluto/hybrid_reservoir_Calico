# E5.2 — occupancy-first one-AHL horizon-1 CRP scout (Task B only)

Frozen: 2026-08-18, encoding selected on occupancy only, before F408.

Track E5 living scout on the claim dish. It does not rewrite Track B,
E4, E5.0, E5.1, Narma10b, Waveform, or Lorenz. `GATE_EVIDENCE.md`
files are not edited. `K`, `n`, `tau_R`, `tau_L`, clamp, mortality,
flow, layout, and `Jmax` are not retuned. Claim dish stays CENTER /
`FLOW=0`.

Task A / `Y` classification is not run. Glucose, Danino, vesicles,
acid-as-lactate, Track B 5-AC sites, sequential batching, C1,
washout twins, two-way, IPC, Lorenz 202/333, Stage99, and E4
retunes are not started.

`COHORT_STATUS` remains `SYNTHETIC_COHORT`. Not MIMIC. Not clinical.
E5.1 `NOT_SCORED` is not rewritten into a PASS.

## Scientific question

Can a predeclared affine map of CRP occupy the Hill on the claim
dish at frozen `Jmax = 1.28e8`, and then can living 408-D forecast
next-window CRP better than Brownian, silent, and field?

Biology is not required to beat 5-channel delay-10. That baseline
sees four channels the dish does not get. Fair input ceilings:
persistence and CRP-only delay-10, same patient-pooled metric,
delays never crossing patient bounds. Target is still raw
`Xnorm` channel 0, not `u`. Concat NRMSE is a labelled diagnostic
only.

## Frozen scout (do not redraw)

From `examples/HybridDish/biomarker_protocol_audit/results/RESET_ISOLATED_DESIGN.md`.
Same 12+8 as E5.0 / E5.1.

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

## Claim dish (copied from E5.1 / Narma10b — do not change)

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
- Warmup `18000 s` at **that patient’s** `u[0]` under the frozen map
- 48 windows, 300 s, pulse 75 s, sample every 20 s, `dt=0.05`

CSV: 48 windows / 768 samples / 768 voxel rows (16 samples/window).
Last sample prefix: `47;15;299.95`.

Silent and Brownian seed 111 are reused from E5.1 if last sample and
48/768/768 match. E5.1 **driven** CSVs are not reused (the command
changed).

## Module 0 — encoding family (occupancy only, before F408)

`Jmax` stays `1.28e8`. One AHL wire. Clip `u` to `[0, 0.5]`.

```
u[n] = clip(u_floor + g * Xnorm[p, n, 0], 0, 0.5)
```

Predeclared maps. A fifth map is not added after seeing occupancy.

| MapID | u_floor | g | Notes |
|---|---|---|---|
| M_E51 | 0.00 | 0.50 | E5.1 replay (mostly DEAD). Reference. |
| M_F15 | 0.15 | 0.35 | Lift floor, keep contrast |
| M_F20 | 0.20 | 0.30 | Mid |
| M_F25 | 0.25 | 0.25 | Highest floor, least CRP range |

Warmup `u[0]` uses the same map.

### Occupancy screen (no living ridge)

For each map × 12 development patients, estimate dish-mean \(R\)
and \(r(R,u)\) over 48 windows.

**Preferred path (used if importable):** E0.2 transport model
(`examples/HybridDish/design_space_preflight/`, status
`VALIDATED_FOR_SCREENING`) plus Stage 3B Hill
(\(K=1.6\), \(n=2\), \(\tau=15\)). Claim condition
`PRI_CENTER_F0p0`: CENTER, `FLOW=0`, `NO_FLUX`, `D=159`,
`k=0.0033`. Warmup is 60 windows at that patient’s `u[0]`, not
the NARMA `0.5` default.

**Density mask (documented before labels):**
`FROZEN_E51_DRIVEN_DEN_WINDOWMEAN`. For development patient \(p\),
window \(n\), `Den_*` is the 20×10 density from E5.1 driven seed 111
(same clamp, acid, layout; AHL command differs). Average `Den` over
the 16 samples in the window. Map E0.2 \(R\) onto the same 20×10
readout. Then

```
R_win[n] = sum_v R_readout[n,v] * Den[n,v] / max(sum_v Den[n,v], 1)
mean_R   = mean_n R_win
r(R,u)   = Pearson(R_win, u)
```

This is the same density weighting used in
`TRANSPORT_MODEL_VALIDATION.md`. It is a screen, not a living score.
Uniform (unweighted) readout-mean \(R\) is recorded as a diagnostic
and is not used for selection.

**Fallback if that path cannot be run:** per-window linear rescale of
E5.1 driven dish-mean AHL for the same patient,

```
AHL_new[n] = AHL_E51[n] * u_new[n] / max(u_E51[n], 1e-6)
R_proxy[n] = AHL_new[n]^2 / (K^2 + AHL_new[n]^2)    with K=1.6
```

then `mean_R = mean(R_proxy)` and `r(R,u) = Pearson(R_proxy, u)`.
This algebraic Hill is an occupancy proxy, not living \(R\).

ALIVE if `mean_R ≥ 0.05` and `|r(mean_R, u)| ≥ 0.5`.

SATURATED (map-level) if median over the 12 patients of `mean_R > 0.40`
or the fraction of the 12×48 windows with `R_win > 0.5` exceeds 0.50.

### Selection (hash before Java)

Among maps that are not SATURATED:

1. Maximize `n_ALIVE` on the 12 development patients.
2. Require `n_ALIVE ≥ 10`. If none qualify, **STOP**. Write
   `ENCODING_FAIL`. Do not raise `Jmax`. Do not invent Map 5.
   Do not run confirmation. E5.1 stays `NOT_SCORED`.
3. Ties: pick the map with **larger** `(max u − min u)` on the
   concatenated 12-patient `u` (more CRP contrast). Remaining tie:
   smaller `u_floor`.

`results/ENCODING_FREEZE.md` records the occupancy table, the winner
or `ENCODING_FAIL`, and SHA-256 of every development `u` file
**before** any F408. `M_E51` appears in the table even if it loses.

CRP-only ceilings (persistence, CRP delay-10, 5-channel delay-10)
are computed on `Xnorm[:,:,0]`, not on `u`, and do not select the map.

## Module 1 — living BSim (winner only)

Smoke: one development patient. Print map id, `u[0]`, CENTER,
`FLOW=0`, 48 windows, finite AHL.

Then:

| Run | n | Notes |
|---|---|---|
| development driven | 12 | seed 111, isolated, frozen map |
| confirmation driven | 8 | **only if** living development `n_ALIVE ≥ 10` |
| silent / Brownian | reuse E5.1 seed 111 | last sample `47;15;299.95`, 48/768/768 |

If living development `n_ALIVE < 10`: `NOT_SCORED`. Do not switch
maps. Do not start confirmation. That is a valid E5.2 outcome.

Keep DEAD rows in the all-patient table. Do not raise `Jmax`.

Occupancy DEAD per patient if `mean_R < 0.05` or `|r(mean_R, u)| < 0.5`.

## Module 2 — Task B ridge (only if occupancy gate passed)

Official pairs: `t = 9 .. 46` (38 per patient), aligned to delay-10.
Also print `t = 0 .. 46` as a diagnostic.

Ridge: 408 features, LOPO lambda on the 12 development patients,
grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`, intercept unregularized,
train-patient standardization. Ties → larger lambda. Refit 12.
Score confirmation 8 once. Do not refit on 8.

Metric: per-patient NRMSE / that patient’s target std, then mean
over patients in the scored set. Concat NRMSE is a labelled
diagnostic only.

Readouts: F408, field-only AHL, occupancy-masked RL surrogate,
reused Brownian, reused silent, persistence, CRP-only delay-10,
5-channel delay-10 (ceiling, not a living-layer bar).

**Primary:** all confirmation patients (DEAD included if any).

**Secondary, predeclared:** in-cohort = `mean_R ≥ 0.05` on that
patient’s living run, applied to development and confirmation
independently, occupancy computed **before** looking at F408.
Label it occupancy-conditioned. Do not make it the only number.

### Standing after scores

| Label | Rule |
|---|---|
| Occupancy | ALIVE count on 12 / 8 |
| System | confirmation F408 < Brownian and < silent (patient-pooled) |
| Living-layer | F408 vs field, reported |
| vs persistence / CRP delay-10 | reported; not required to beat 5-ch delay-10 |
| Clinical | FORBIDDEN |
| Task A / Track B / E5.1 | unchanged |
| Claim dish | not replaced |

## Forbidden

- Task A / `Y` classification
- glucose, acid-as-lactate, Track B layout, second AC
- sequential 100-patient dish
- concat NRMSE as the gate
- fitting `K`, `n`, `tau`, `Jmax` to NRMSE or occupancy
- adding a map after the freeze
- calling this clinical
- overwriting the E5.0 split or scout hashes
- rewriting Track B FAIL or E5.1 `NOT_SCORED`
- C1, washout twins, two-way, IPC

## Deliverables

- `PROTOCOL.md` (this file)
- `README.md`
- hashed input/target files and `results/ENCODING_FREEZE.md`
- `results/E5_2_SCOUT.md`
- `results/e52_occupancy.csv`
- `results/e52_taskb.csv`
- `check_e52.py`

Checker: `python examples/BSimReservoirPlanE52CRP/check_e52.py`
from the repository root.
