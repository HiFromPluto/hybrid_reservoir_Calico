# WashoutReset — chemical / hydrodynamic residual scout

Development-only Track E0. Not a task scout. No NARMA, no Mackey–Glass
NRMSE, no waveform AUC, no ridge. Does not edit any `GATE_EVIDENCE.md`.
Does not promote a new claim dish. Does not retune `K`, `n`, `tau`,
rate, `D`, decay, clamp, or mortality.

C1 remains DEFER. E0.3 NARMA remains NO_STORY_MOVE. SweepMGE03 is not
a new claim dish. Waveform1 FAIL and WaveformS stay. Two-way is future
work.

## Scientific question

On the HybridDish clocks, five silent windows already kill extracellular
AHL by decay (`tau_AHL = 1/0.0033 s = 303.0303 s`,
`e^{-1500/303.0303} = 0.00705`). The chip problem is the rest:

1. Does `L` stick (`tau_L = 1500 s`, theory `e^{-1} = 0.367879`)?
2. Does an 8 µm/s flush clear AHL faster than decay, or mainly dump
   cells (and their `L`) out of bounds?
3. Is 0.25 µm/s trickle during reset any different from decay-only?

Hypothesis (predeclared): **chemical washout ≠ reporter washout.**
AHL and `R` die in five silent windows at `FLOW=0`. `L` does not.
Flush can lower dish-mean `L` only if population falls. Trickle should
look like decay-only on AHL/`L`, with a small OOB bump.

## Dish (copy — only epochs change)

- `1000 × 500 × 10` µm, `dt=0.05`, CENTER AHL `(500, 250, 5)`, rate
  `1.28e8` when `u>0`
- Acid `(300, 375, 5)`, rate `2e11` when acid `u>0`
- Attractants silent. `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, `warmup.ahl.input=0.5`
- Window `300 s`, pulse `75 s` when `u>0`, sample every `20 s`
- Seed **101** only unless the W=15 follow-on rule fires
- One-way `addQuantity`. No vesicles, no glucose, no Danino, no extra
  ACs, no ridge

When `FLOW=0`: chemical advect off, `NO_FLUX`.
When `FLOW>0`: SweepS5 upwind on AHL, acid, attractant, repellent;
left inlet 0; no recycling (`OUTFLOW`); Stokes +x on bacterium.
STOP if Courant ≥ 1. Do not lower `dt`. Expected Courant at 8 µm/s is
0.02 (`dx=20` µm). At 0.25 µm/s is 0.000625.

No AHL sequence file. Held `u` from config. Mid-run epoch switch
prints `FLOW_SPEED`, boundary, AHL `(x,y,z)`, Courant, and epoch name
at window 0, wash start, wash end, and last window. Prints `u` and
acid `u` every window.

## Frozen epoch schedule (primary, W=5)

`num.windows=40`. Last sample `39;15;299.95`. CSV 40 / 640 / 640.

| Epoch | Windows | AHL `u` | Acid `u` | Flow | Boundary |
|---|---|---|---|---|---|
| Load | 0..19 | 0.5 held | 0.5 held | 0 | NO_FLUX |
| Wash | 20..24 | 0 | 0 | **arm** | NO_FLUX if 0 else OUTFLOW |
| Post | 25..39 | 0 | 0 | 0 | NO_FLUX |

### Primary arms (seed 101)

| ConditionID | Wash flow |
|---|---|
| WASH_DECAY_W5 | 0 µm/s |
| WASH_TRICKLE_W5 | 0.25 µm/s |
| WASH_FLUSH_W5 | 8 µm/s |

Do not add 0.5 / 1.0 / 3 µm/s, upstream, or downstream.

## Theory (print before BSim, not a retune)

`tau_AHL = 303.0303 s`, `tau_R = 15 s`, `tau_L = 1500 s`.
Wash duration `W=5`: `dt = 1500 s`.

| Quantity | Ideal relative leftover after 5 silent windows |
|---|---|
| AHL | `e^{-1500/303.0303} = 0.00705` |
| `R` | `e^{-1500/15} ≈ 0` |
| `L` | `e^{-1} = 0.367879` |

Ideal leftover after 15 silent windows (`W=15`, `dt = 4500 s`):
AHL ≈ 0, `L = e^{-3} = 0.049787`.

These are spatially uniform exponential predictions. Living residuals
will differ because of walls, clamp, OOB, and occupancy. Report the
difference. Do not change `tau_L` to match.

## Metrics (no ridge)

Last sample of `t_load` = window 19, `t_wash` = window 24 (or 34 if
W=15), and each post window: dish-mean AHL, `R`, `L`, population,
cumulative OOB deaths, cumulative acid deaths, occupancy (`mean_R`,
`r(R,u)` on load only).

Relative leftovers at `t_wash` vs `t_load`:
`AHL_rel`, `R_rel`, `L_rel`, `N_rel`. Also `L_rel` at last post
window vs `t_load`.

Silent floor: Stage 6 / BenchA silent seed 101 analysis-start
`mean_L` if those voxels exist. Do not rerun silent unless missing.

## Predeclared findings (not NRMSE gates)

Evaluate on seed 101, W=5, after wash:

1. **AHL_CLEAR** if `AHL_rel < 0.05` on WASH_DECAY_W5
2. **L_STICKS** if `L_rel > 0.20` on WASH_DECAY_W5
3. **FLUSH_DUMPS_CELLS** if WASH_FLUSH_W5 has `N_rel ≤ 0.80`
4. **FLUSH_DUMPS_L** if WASH_FLUSH_W5 `L_rel` is ≥ 0.10 below
   WASH_DECAY_W5 `L_rel` **and** FLUSH_DUMPS_CELLS
5. **TRICKLE_LIKE_DECAY** if WASH_TRICKLE_W5
   `|L_rel − L_rel(decay)| < 0.05` and
   `|AHL_rel − AHL_rel(decay)| < 0.02`

All five may fire. None promote a claim dish. None retune kinetics.

**STOP / implementation fail (no follow-on, no retune):**

- Courant ≥ 1, or FLUSH printed flow ≠ 8, or source not CENTER
- WASH_DECAY_W5 `AHL_rel ≥ 0.05` — still write the report, do not
  start W=15
- Occupancy DEAD at `t_load` (`mean_R < 0.05`) — do not interpret wash

## Follow-on (only if L_STICKS and AHL_CLEAR)

`num.windows=50`. Last sample `49;15;299.95`. CSV 50 / 800 / 800.
Load 0..19, wash 20..34, post 35..49.
Arms: `WASH_DECAY_W15` (0 µm/s) and `WASH_FLUSH_W15` (8 µm/s).
Do not run trickle at W=15. Do not start seeds 202/303.

## Smoke (not evidence)

`num.windows=6` (2 load, 2 wash at 8 µm/s, 2 post). Print geometry,
CENTER, `FLOW_SPEED=8` during wash, Courant 0.02, finite non-negative
AHL, population finite.

## Checker

```
python examples/BSimReservoirPlanWashoutReset/check_washout_reset.py --theory
python examples/BSimReservoirPlanWashoutReset/check_washout_reset.py
```

`--theory` prints the exponential table before any living score.

## Compile

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanWashoutReset.java VoxelAnalyzer.java
```

## File hashes (recorded before scores)

SHA-256 of each file as stored on disk. Filled after the package is
written; result hashes after scoring.

| File | SHA-256 |
|---|---|
| BSimReservoirPlanWashoutReset.java | 0a576eddd6d34b28b730ab97b1eed6d3c4a642dc9420c64948c7b7ce3c84a9c6 |
| VoxelAnalyzer.java | 4d0abba4397b1ea666c5579ff43a8695670de6b61137a0f172d622b5707be2bc |
| check_washout_reset.py | a3578d50f376d8bba5e5d9b045b83bd09adba06bd86a6c7ec2f0c09e0f278b9a |
| run_washout_jobs.py | 4c478cbb3a317b8a84fc4cd263a486fea9cb7e8340b118b47543c80ae8b5e8fc |
| sim_config_washout.properties | 57716d412d2c4b41242cb499cc260445d7d78c531341006b24d3c59c03881376 |
| sim_config_washout_smoke.properties | 9e51fbe1dd8efcc16200d89c132068304b834001a2ae4965807b9e8089266a8b |
| sim_config_wash_decay_w5_seed101.properties | 7b0a79a30ffb7f42d59e5afbc240098336eec5c8becdbfb48c7a12b6a586e4ae |
| sim_config_wash_trickle_w5_seed101.properties | 946680b7d3d3a46f241df38983f445a6350b87e54ff85884334fb790557100e6 |
| sim_config_wash_flush_w5_seed101.properties | 1231ab1ae61e18a072481871986279da6e59dc36e066bd7485566b601ef76556 |
| sim_config_wash_decay_w15_seed101.properties | 06003260c734a34864760b33717dc4957593af96a359623516c3c8c79282c1cc |
| sim_config_wash_flush_w15_seed101.properties | 24d600a114ad41c90cd6e8c98b982f60678958b1980fbafeb5e8569efaa4f0f7 |

## Result hashes (after scoring)

| File | SHA-256 |
|---|---|
| results/WASHOUT_RESET_SCOUT.md | 7e23cad728ef720f74767cc88546628731345fb37822226c7b5efbfdd1cea7da |
| results/washout_reset_scout.csv | e430c233ee8e237dd35bd83a7927d46b6c7a1b87e8d2ed69ad5f003a2445aa6a |
| results/wash_decay_w5_seed101/run_status.txt | 31ccf60a70d80827324b1f8dc504a61818776ce03e822564381e2920440db46e |
| results/wash_trickle_w5_seed101/run_status.txt | 314aaf84a25a8d9689509a1438749a8a88eadf7d619d29d75669397cc4981414 |
| results/wash_flush_w5_seed101/run_status.txt | ac7c5b753726451a0748d89c1d0a7c1717e2e5e270e0128a18a07c0c0b9c5ce2 |
| results/wash_decay_w15_seed101/run_status.txt | cfd9117910b667d0f3be929c3d42169bcfefa59870567829e3f69de862a294a1 |
| results/wash_flush_w15_seed101/run_status.txt | 92f74b83e05b728b2db3fac77a439be3a0f531c5a6986f569283fa88c121d294 |
| results/washout_smoke/run_status.txt | eb9c03633ad4843bb50d2349a97468dcd5eaa8dd4dd888cf2603e8b9a9ae38c5 |

## Stop list

- Do not run MG, NARMA, waveform, Lorenz, biomarkers, two-way, IPC
- Do not add geometries or extra speeds except the W=15 follow-on
- Do not retune `tau_L` to hit 0.368
- Do not raise `J_max` or lower 8 µm/s if the dish empties
- Do not rewrite BenchA / E0.3 / Waveform gates
- Do not call this a sequential-patient validation
