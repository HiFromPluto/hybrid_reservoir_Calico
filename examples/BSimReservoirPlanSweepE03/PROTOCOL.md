# Sweep E0.3 — layout-by-flow living scout (development only)

Frozen before scores. This is Track E0 development, not confirmatory
NARMA. It does not edit `GATE_EVIDENCE.md`, historical PASS/FAIL, claim
kinetics, C1, Waveform2, Lorenz, IPC, AC implementation, biomarker BSim,
or Stage99.

Do not fit D, decay, K, n, tau_R, tau_L, source rate, clamp, or acid
mortality after seeing occupancy or scores.

## Scientific question

At fixed commanded AHL payload, do source position and low +x flow
change occupancy, viability, and the living/field/surrogate/direct-input
gaps relative to the frozen claim dish?

A condition may be physically interesting and still lose NARMA. Report
both. Do not select a winner on NRMSE alone. `0.93` stays a weak
predictor in the write-up.

## Frozen living table (exactly these)

| ConditionID | AHL position | Flow um/s | Chemical boundary |
|---|---|---|---|
| PRI_CENTER_F0p0 | (500,250,5) | 0 | NO_FLUX |
| PRI_CENTER_F0p25 | (500,250,5) | 0.25 | OUTFLOW |
| PRI_CENTER_F0p5 | (500,250,5) | 0.5 | OUTFLOW |
| PRI_CENTER_F0p666667 | (500,250,5) | 0.666667 | OUTFLOW |
| PRI_CENTER_F1p0 | (500,250,5) | 1.0 | OUTFLOW |
| PRI_UPSTREAM_CENTER_F0p0 | (200,250,5) | 0 | NO_FLUX |
| PRI_UPSTREAM_CENTER_F0p666667 | (200,250,5) | 0.666667 | OUTFLOW |
| PRI_DOWNSTREAM_CENTER_F0p0 | (800,250,5) | 0 | NO_FLUX |
| PRI_DOWNSTREAM_CENTER_F0p666667 | (800,250,5) | 0.666667 | OUTFLOW |

`HISTORICAL_CENTER_F8_ANCHOR`: do not rerun. Cite SweepS5 seed 111
(`mean_R=0.0029`, occupancy DEAD). Rerun 8 um/s only if the new advect
implementation cannot be shown comparable to SweepS5.

Do not add 2.0 or 3.333 um/s. The preflight min-C exclusion stands.

## What stays frozen

Copied from HybridDish / Narma10b:

- `1000 x 500 x 10` um, `dt=0.05` s, field `50 x 25 x 1`, readout
  `20 x 10` + `4 x 2`
- `D=159`, `k=0.0033`, AHL rate `1.28e8`, 75 s pulse / 300 s window
- `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500` (`alpha=delta=1/1500`)
- acid source `(300,375,5)` held 0.5; attractants silent
- warmup 18000 s, 200 windows, `INITIAL_POP=1800`, clamp `K=2000`
- official 408: `Receiver_R`, `Lum_Mean`, `Input_Driven_Death`
- closed ridge: washout 40 / train 110 / test 50 (`ridge_closed.m` rule)

Primary payload: the same commanded AHL sequence and rate on every
condition. Do not scale rate with flow or position. No matched-occupancy
arm in this scout.

## Advection / OUTFLOW

Reuse SweepS5 Stage-6 upwind advect on AHL, acid, attractant, and
repellent: left inlet 0, no recycling. That is the OUTFLOW chemical
boundary. Stokes +x on `ReservoirBacterium` and `BrownianParticle`.

At `FLOW_SPEED=0`: no chemical advect; NO_FLUX walls as in the claim
dish.

Print every run: `FLOW_SPEED`, `transit_s`,
`chemical_Courant=FLOW*dt/dx`, `CFL_dt_max_flow`, advect fields,
Stokes on/off, AHL `(x,y,z)`, boundary name.

STOP if Courant >= 1. Do not lower dt. Expected Courant at 1 um/s is
0.0025. Source must remain inside bounds.

## Inputs

Frozen Narma10b AHL, acid, and target. Do not regenerate `u`.

`u` SHA-256:
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`

Screen random-drive hash is not a BSim input:
`5919f8fe318308163b1d798967fd8f89e99ca5eaa268af0b8a62e30e0b67307f`

Scout seed is 111 only. Seeds 222/333 wait for a later confirmation
prompt.

## BSim jobs

Reuse, do not rerun:

- `PRI_CENTER_F0p0` driven = Narma10b driven seed 111
- FLOW=0 silent = Narma10b silent seed 111
- FLOW=0 Brownian = Narma10b Brownian seed 111

New production BSim, seed 111 only (16 runs):

Driven (8): CENTER 0.25, 0.5, 0.666667, 1.0; UPSTREAM 0 and 0.666667;
DOWNSTREAM 0 and 0.666667.

Silent (4) and Brownian (4): one each of unique flow > 0
(`0.25, 0.5, 0.666667, 1.0`). Shared by flow, not by layout.

Field-only, unmasked surrogate, occupancy-masked surrogate, and
direct-input baselines are analysis on driven voxels / frozen `u`.
Do not rerun BSim for them.

## Smoke before production

One short HybridDish smoke: CENTER, `FLOW_SPEED=0.666667`, OUTFLOW.
Print geometry, source, flow, Courant, boundary. Finite non-negative
AHL. Source inside bounds. Expected smoke row counts. A smoke pass is
not evidence.

## Claim sanity before new scores

On reused Narma10b seed 111, closed ridge, within `1e-3`:

- driven F408 `0.9312`
- field `1.0289`
- Brownian `1.1625`
- intercept `1.1622`
- persistence `1.0068`

If sanity fails, stop. Do not interpret new conditions.

## Occupancy and viability (before NRMSE claims)

For every driven condition, using window-mean `Receiver_R` and AHL:

- DEAD if `|r(mean_R, u)| < 0.5` or mean of mean_R `< 0.05`
- SATURATED if mean of mean_R `> 0.8` or mean frac(`R>0.5`) `> 0.8`
- else ALIVE

Also report mean_AHL, mean_pop, clamp/acid/OOB deaths, births,
wall/interior AHL ratio, plume coverage.

If DEAD or SATURATED: keep the row, do not promote, do not raise rate,
flow, or K, do not retune clamp or mortality.

Additional scout exclusions (report, do not retune):

- mean R `< 0.02` or `>90%` weak samples (`H<0.05` if H is available;
  else `R<0.05` as the weak proxy on living readout)
- mean H `> 0.90` or `>90%` saturated samples
- population collapse comparable to SweepS4 small-box or SweepS5 flow-8
- Courant/NaN/negative concentration

## Analysis (every living condition and every control that exists)

- official F408 NRMSE / R2
- field-only AHL
- FRL, R-only, L-only (descriptive)
- unmasked and occupancy-masked kinetic surrogates (frozen R/L ODEs)
- `LINEAR_U_DELAY_10` and `NARMA_INFORMED_INPUT`
- train-intercept and persistence
- reused silent/Brownian scored on this same target

Physical table, not just scores: occupancy class, mean_R, mean_AHL,
mean_L, mean_pop, wall/interior ratio, transit time, Courant.

CSV completeness: 200 / 3200 / 3200, last sample `199;15;299.95`.

## Promotion (scout)

`NO_STORY_MOVE` unless a condition is ALIVE and changes the purpose
story by beating Brownian and silent AND beating field by a margin that
is not a 0.03 tick around 0.93, while remaining worse than or only
comparable to direct-input `0.6828` if that baseline still wins.

If a condition beats the claim dish on occupancy/stability but not on
NARMA, say so. That can still be a transport finding.

Do not freeze a new claim dish in this scout. At most name 0–2
candidates for a later confirmation prompt on untouched inputs and
seeds 222/333. Retain every completed, dead, failed, and excluded row.

C1 stays DEFER.
