# SweepMGE03 — Mackey–Glass three-geometry scout (development only)

Frozen before scores. Track E0/E3 development, not a new claim dish. It
does not edit `GATE_EVIDENCE.md`, BenchA Mackey–Glass Overall PASS,
historical PASS/FAIL, claim kinetics, C1, WaveformS, E0.3 NARMA,
Lorenz, IPC, two-way, or Stage99.

Do not fit D, decay, K, n, tau_R, tau_L, source rate, clamp, or acid
mortality after seeing occupancy or scores. Do not scale rate with
flow or position. No matched-occupancy arm. No extra ACs.

## Scientific question

At fixed commanded MG payload, do source position and mild +x flow
change occupancy and the **driven–field gap** relative to the BenchA
claim dish?

Hypothesis (predeclared): MG is carrier-rich (field NRMSE 0.0229 vs
driven 0.3306 on seed 101). Mild flow may smear the spatial delay
line more than it smears L. If field gets worse while driven holds,
the gap shrinks. If both arms move together, layout is transport,
not living computation. Either outcome is a result.

Field-only is a diagnostic, not a kill switch. Persistence of x
(~0.180) is a mandatory trivial baseline. Do not report BenchA
driven MC 15.331 as CHARC memory capacity (smooth correlated drive).

## Frozen geometry table (exactly these three)

| ConditionID | AHL position | Flow | Chemical boundary |
|---|---|---|---|
| PRI_CENTER_F0p0 | (500, 250, 5) µm | 0 µm/s | NO_FLUX |
| PRI_CENTER_F0p25 | (500, 250, 5) µm | 0.25 µm/s | OUTFLOW |
| PRI_UPSTREAM_CENTER_F0p0 | (200, 250, 5) µm | 0 µm/s | NO_FLUX |

Do not add downstream, 0.5 / 0.67 / 1.0 / 8 µm/s, wall, or off-centre
y. 8 µm/s remains the historical occupancy-dead NARMA anchor; do not
rerun it here. Acid stays at `(300, 375, 5)` even when AHL moves.
Attractant AC0/AC2 stay at `u=0`.

## Dish (copy — do not change except the three knobs)

- `1000 × 500 × 10` µm, `dt=0.05`
- Field `50×25×1`, readout `20×10` + `4×2`
- `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL rate `1.28e8`. Acid `2e11` at `(300, 375, 5)`, held 0.5 in
  analysis windows, off in warmup
- Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, window `300 s`, pulse `75 s`, sample every `20 s`
- **200 windows**. Last sample `199;15;299.95`. CSV 200 / 3200 / 3200
- Official 408: `Receiver_R_*`, `Lum_Mean_*`, `Input_Driven_Death_*`
- Brownian: `Den_*` 20×10 only
- Field-only: driven `AHL_uM_*` 20×10
- Scout seed **101** only. Seeds 202/303 only if a story-move rule
  fires
- One-way `addQuantity`. No vesicles, no glucose, no Danino

Advection: SweepS5 / E0.3 Stage-6 upwind on AHL, acid, attractant,
and repellent; left inlet 0; no recycling (`OUTFLOW`) when
`FLOW>0`. Stokes +x on bacterium and `BrownianParticle` when
`FLOW>0`. At `FLOW=0`: chemical advect off, `NO_FLUX`.

Print every run: `FLOW_SPEED`, `transit_s`,
`chemical_Courant=FLOW*dt/dx`, `CFL_dt_max_flow`, advect on/off,
Stokes on/off, AHL `(x,y,z)`, boundary name. STOP if Courant ≥ 1.
Do not lower `dt`. Expected Courant at 0.25 µm/s is 0.000625
(Δx=20 µm).

## Frozen MG task (do not regenerate)

Copied, not rewritten, from BenchA:

- `input_ahl_mg200.txt`
- `input_acid_held05_200.txt`
- `mg_target.csv`

SHA-256 of the 12-decimal `u` sequence:
`e06810ea70f61bbe4baf61939f93b1e1cae9d92d9eee8d757865c4935dbac780`

Target for window `n`: raw Mackey–Glass `x[n+1]` (not affine `u`).
Teacher-forced one-step. Not Jaeger free-run.

Ridge: washout 40 / train 110 / test 50. Lambda on windows 128–149
only. Grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`. Unregularized
intercept. Train-only standardization. NRMSE = RMSE / pop-std of the
**test** target (same as BenchA). Independent lambda per
condition×arm.

## Kinetic surrogates

Unmasked and occupancy-masked surrogates use the same Stage 3B /
Stage 5 ODEs as the zero-simulation audit (`tau_R=15`, `K=1.6`,
`n=2`, `tau_L=1500`) on each driven AHL voxel file. No extra BSim.

Piecewise-constant AHL approximation: voxel samples are 20 s apart.
Between sample `t_{k-1}` and `t_k` the field is held at `AHL[t_{k-1}]`.
The R/L update on a constant-H interval is exact. Occupancy mask
zeros empty density voxels. Occupancy DEAD if `mean_R < 0.05` or
`|r(mean_R, u)| < 0.5`. If DEAD: still report NRMSE, do not promote,
do not raise rate or `K`.

## Claim sanity (before any new score)

Reuse BenchA MG seed 101 voxels for `PRI_CENTER_F0p0`. Recompute
closed ridge. Must match within `1e-3`:

| Arm | Expected seed 101 |
|---|---|
| driven F408 | 0.3306 |
| field-only | 0.0229 |
| Brownian | 1.0415 |
| silent | 1.0416 |

Also print once (target-only, identical for every geometry):

- persistence `x[n] → x[n+1]` on this split (BenchA pack ~0.1799)
- intercept / train-mean predictor (~1.0416)

If sanity misses, stop. Do not interpret new geometries.

## BSim jobs (seed 101)

Reuse, do not rerun:

- `PRI_CENTER_F0p0` driven = BenchA `results/mg_driven_seed101`
- `PRI_CENTER_F0p0` and `PRI_UPSTREAM_CENTER_F0p0` silent = BenchA
  MG silent seed 101 (`FLOW=0`, sources off; layout of an off source
  does not matter)
- `PRI_CENTER_F0p0` and `PRI_UPSTREAM_CENTER_F0p0` Brownian = BenchA
  `results/mg_brownian_seed101` (`FLOW=0`; particles do not sense AHL)

New production BSim, seed 101 only (4 runs):

1. `PRI_CENTER_F0p25` driven
2. `PRI_CENTER_F0p25` silent
3. `PRI_CENTER_F0p25` Brownian
4. `PRI_UPSTREAM_CENTER_F0p0` driven

## Story-move rules (predeclared, seed 101, ALIVE only)

Let `D` = driven F408 NRMSE, `F` = field NRMSE,
`G = D − F` (positive ⇒ field better).

Claim seed 101: `D0=0.3306`, `F0=0.0229`, `G0=0.3077`.

A new condition is a **story move** if occupancy is ALIVE and any of:

1. **Carrier moved:** `|F − F0| ≥ 0.020`
2. **Living moved:** `|D − D0| ≥ 0.050`
3. **Gap moved:** `|G − G0| ≥ 0.050`
4. **Living-layer flip:** `D < F` (biology beats field)
5. **System lost:** `D ≥` Brownian on that condition

If none fire: **NO_STORY_MOVE**. Do not start seeds 202/303. Do not
promote a new claim dish. `0.26` stays a weak one-step predictor vs
field 0.023 and persistence ~0.18.

If a story-move fires: run 202 and 303 **only** for the conditions
that fired, plus the claim `CENTER/0` ridge already on disk. No
best-seed. Do not freeze a new claim dish in this job.

Confirmation runner (authorized only after seed-101 STORY_MOVE):

```
python examples/BSimReservoirPlanSweepMGE03/run_mge03_confirm_jobs.py
```

Ticks smaller than those thresholds are “same story, different
0.26.” Write them down.

## Smoke (not evidence)

One short run: CENTER, `FLOW_SPEED=0.25`, OUTFLOW. Print geometry,
AHL `(500,250,5)`, Courant, boundary. Source inside bounds. Finite
non-negative AHL. A smoke pass is not evidence.

## Checker

```
python examples/BSimReservoirPlanSweepMGE03/check_mg_geom.py --sanity
python examples/BSimReservoirPlanSweepMGE03/check_mg_geom.py
```

`--sanity` must pass on reused BenchA seed 101 before new scores
are interpreted.

## Compile

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanSweepMGE03.java VoxelAnalyzer.java
```

## File hashes (recorded before scores)

SHA-256 of each file as stored on disk.

| File | SHA-256 |
|---|---|
| BSimReservoirPlanSweepMGE03.java | c7ff3b1d22000f39debe38740386f5f03db1f0a39d823940b09ab1bb1dface53 |
| VoxelAnalyzer.java | b699ef6f3ecb7e28dc557b1d723636ca99fff64c435c84739be38307fc8d6382 |
| check_mg_geom.py | cb5af3376449d4cffca4115263e99b4845d2405ef1e7fd5497a55690df3fa048 |
| run_mge03_jobs.py | 9f7dbca72a04b5537b8ef5364c45cd66b68c77f6c0b3ccc883cb0c3cf857b473 |
| sim_config_mg.properties | 04dcf534e3c881d73a59249c590a9417b54617651385ba7c3f515871374c0ca4 |
| sim_config_mge03_smoke.properties | 1c6f9ff89b0d89cb45d2bf5d6136cb2545ddf71863665f23b301ea5e1e9d88bb |
| sim_config_mge03_center_f0p25_driven_seed101.properties | 75fccaf139ad8f0d07f68501867429f344d6dfb7c1360e8249471ba95863d6c9 |
| sim_config_mge03_center_f0p25_silent_seed101.properties | 53c4034133186f8e8cb30ce62dbffd6dfd47a2b3145850f5f06a08f6abfcd082 |
| sim_config_mge03_center_f0p25_brownian_seed101.properties | de7ebfb85656c3d1ee3285b95b0d370c2154ca35114be23a9930bd9f2300b4b6 |
| sim_config_mge03_upstream_f0p0_driven_seed101.properties | 0488caf28d31d42e08e3c99214f56871340aa97652af45af89e23469c1581f9a |
| input_ahl_mg200.txt | df6ddff58a02725e54ab3b9027b30d20168eb6813d1439cd34990938560c63c0 |
| input_acid_held05_200.txt | c00b66f1ee5f2f2f098bd657a96a8746b126a0d5dd7bf8e27e525c21369ca4f9 |
| mg_target.csv | d7475474dc6f7fb02fd6021e7b2c904b38210312d193955caddad9ea2faac393 |

## Result hashes (after scoring)

| File | SHA-256 |
|---|---|
| results/MG_GEOM_SCOUT.md | 993691d81018c8798a5a10bc3515759623c4b33fc0d674c815a273de77a798d2 |
| results/mg_geom_scout.csv | 31e9e6d8b0a565ba91eda5497a5e06b287cdcd250f462a3ac7b690e91924c807 |
| results/mge03_center_f0p25_driven_seed101/run_status.txt | b4247e47f7bedcc0f11d80ec233d3a3d29d6722517c0d529dc46ed51c56d7445 |
| results/mge03_center_f0p25_silent_seed101/run_status.txt | 8b788c31c43c72c7512adeb5f11abdb98978e925c673eb81c23cbba3d4396a90 |
| results/mge03_center_f0p25_brownian_seed101/run_status.txt | 90553679c8bb7002cb231d9214445c710f8225bccc5fbcb1a4cc3d819a536ae9 |
| results/mge03_upstream_f0p0_driven_seed101/run_status.txt | 2e1d55b3928008c8c2f45d8624797be1b2d79b06c62352ba722dd362d9ff1d0f |
| results/mge03_center_f0p25_driven_seed202/run_status.txt | 1308114496f5002fe064c955c497211689fa15a029f563368d3c336ce2b0bf48 |
| results/mge03_center_f0p25_silent_seed202/run_status.txt | 9a3b50202fb0d8ec9904971f68bd501afa2c7180c0ac5950334ec54ab34c24b1 |
| results/mge03_center_f0p25_brownian_seed202/run_status.txt | 55a1aeb6d7e0823aaa4b3a91fccb17dc2cce15abb377c26af2bced116371066e |
| results/mge03_upstream_f0p0_driven_seed202/run_status.txt | b9ca2835d8733082a5f61cee76e77b67019167c3e711ff97eabb640708a0677f |
| results/mge03_center_f0p25_driven_seed303/run_status.txt | 80fa7353a58a4b96a84542fdbd7be8e06e668776038e87ef531966d3c1d9025c |
| results/mge03_center_f0p25_silent_seed303/run_status.txt | 25327c64949beae7dcadb91cda9cbb0dd82684783ba73fba314777cafd2a6fe9 |
| results/mge03_center_f0p25_brownian_seed303/run_status.txt | 51870094fc89f3c26b94ceabb21ad3a1e7dd6054eb6f579ce01f01137cedac1d |
| results/mge03_upstream_f0p0_driven_seed303/run_status.txt | 73c6c93505b9d5d3e561065d33a0666b2af9d0f1e558597192c74b3d6731139d |
| results/mge03_smoke/run_status.txt | 2eeb294660fa52a21577d7a1672081d009973bd1d06eb13dc0afcc94b7832893 |

Post-score checker SHA-256 (confirmation writer and own-CENTER/0 rules): 549bff2e143837515f0141f49081cb10d1af75e9d6e2ae1170554e67c2e813ec.
Pre-score checker hash remains in the table above.

Confirmation runner SHA-256: c971a1edfbccecd85be02cb631b590efebdd152bafcddc83b8e1676d41953bc8.



## Stop list

- Do not run Lorenz, WaveformS, NARMA, biomarkers, two-way, IPC
- Do not add geometries
- Do not retune kinetics after NRMSE
- Do not raise flow, rate, or `K` if occupancy drops
- Do not rewrite BenchA GATE_EVIDENCE
- Do not quote MG MC 15.331 as a capacity result
- Do not promote `CENTER/0.25` or UPSTREAM to the claim dish
