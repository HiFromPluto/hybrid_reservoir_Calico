# E5.2 — occupancy-first one-AHL CRP scout

Track E5. New named freeze. Select the encoding on occupancy only,
then score Task B (horizon-1 CRP) on isolated HybridDish runs.
E5.1 `NOT_SCORED` is not rewritten. Track B Overall FAIL is not
rewritten. `GATE_EVIDENCE.md` is not edited.

`COHORT_STATUS = SYNTHETIC_COHORT`. Not MIMIC. Not clinical.

## Status

**Module 0 encoding freeze: `M_F15`** (`u_floor=0.15`, `g=0.35`).
E0.2 transport + Stage 3B Hill, density mask
`FROZEN_E51_DRIVEN_DEN_WINDOWMEAN`. Development screen `n_ALIVE`:
M_E51 6/12 (reference), M_F15/M_F20/M_F25 12/12, none SATURATED.
Tie broken by CRP contrast. `Jmax` not raised. No fifth map.

**Module 1 living occupancy (seed 111):** ALIVE **12/12** development
and **8/8** confirmation. Silent / Brownian reused from E5.1 seed 111
(`47;15;299.95`, 48/768/768).

**Module 2 Task B:** confirmation F408 patient-pooled **0.6274**.
System **PASS** vs Brownian **1.4895** and silent **1.7170**.
Living-layer: F408 **0.6274** vs field **0.6412** (reported).
Fair input: persistence aligned **0.8193**, CRP delay-10 **0.6143**.
5-channel delay-10 **0.6137** is a ceiling, not a living-layer bar.
E5.1 stays `NOT_SCORED`. Claim dish is not replaced.
`COHORT_STATUS = SYNTHETIC_COHORT`.

See `results/ENCODING_FREEZE.md` and `results/E5_2_SCOUT.md`.

`Jmax` stays `1.28e8`. `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`.
Claim dish stays CENTER / `FLOW=0`. No fifth map after the freeze.

## Frozen scout (do not redraw)

Development `1,4,6,7,9,13,16,18,20,22,24,25`
(`83acc7102e532112c961bb45a1327504290ffddea1c06e659bfb5e68be032915`).

Confirmation `8,12,14,15,17,19,21,23`
(`c90cca8f5627c737b38024b4863415e46a384d8715e6d06a5823e8cb93ba3787`).

Target: `Xnorm[p, t+1, 0]`. Do not use `Y`.

## Generate and check

From the repository root:

```
python examples/BSimReservoirPlanE52CRP/generate_e52_encoding.py
python examples/BSimReservoirPlanE52CRP/check_e52.py
```

Encoding freeze hashes `u` files before any F408. If
`ENCODING_FAIL`, stop. Do not raise `Jmax`. Do not invent Map 5.

## Compile and smoke (authorized only after a frozen winner)

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanE52CRP.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanE52CRP.BSimReservoirPlanE52CRP sim_config_e52_driven_p0001_seed111.properties
```

Or: `python examples/BSimReservoirPlanE52CRP/run_e52_jobs.py --smoke`

Development driven (12, seed 111) after smoke. Confirmation driven
(8) only if living development `n_ALIVE ≥ 10`. Silent / Brownian
reuse E5.1 seed 111.

Do not retune `K`, `n`, `tau_R`, `tau_L`, source rate, clamp,
mortality, flow, or layout. Do not start C1, washout twins,
two-way, or IPC.
