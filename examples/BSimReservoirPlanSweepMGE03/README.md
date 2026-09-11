# SweepMGE03 — Mackey–Glass three-geometry scout

Development-only Track E0/E3 scout. Same MG `u` and target as BenchA.
Same claim kinetics. Layout and low flow are the only knobs.

This does not rewrite BenchA Mackey–Glass Overall PASS. It does not
replace the claim dish. It is not E0.3 NARMA and not the nine-condition
factorial. C1 remains DEFER.

## Frozen geometries

| ConditionID | AHL | Flow | Boundary |
|---|---|---|---|
| PRI_CENTER_F0p0 | (500, 250, 5) | 0 | NO_FLUX |
| PRI_CENTER_F0p25 | (500, 250, 5) | 0.25 | OUTFLOW |
| PRI_UPSTREAM_CENTER_F0p0 | (200, 250, 5) | 0 | NO_FLUX |

Seed 101 fired **STORY_MOVE** on CENTER/0.25 and UPSTREAM versus claim
D0=0.3306. Confirmation seeds 202/303 were run. Versus each seed's
own CENTER/0 those geometries are **NO_STORY_MOVE**. Not a new claim
dish. See `results/MG_GEOM_SCOUT.md`.

## Check (before BSim)

From the repo root:

```
python examples/BSimReservoirPlanSweepMGE03/check_mg_geom.py --sanity
```

Must pass on reused BenchA MG seed 101 voxels.

## Compile and smoke (not evidence)

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanSweepMGE03.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepMGE03.BSimReservoirPlanSweepMGE03 sim_config_mge03_smoke.properties
```

## Production (seed 101, four new BSim jobs)

```
python run_mge03_jobs.py
python ../BSimReservoirPlanSweepMGE03/check_mg_geom.py
```

From the repo root the checker is:

```
python examples/BSimReservoirPlanSweepMGE03/check_mg_geom.py
```

See `PROTOCOL.md` and `results/MG_GEOM_SCOUT.md`. Do not retune K, n,
tau_R, tau_L, source rate, clamp, mortality, D, or decay. Do not
promote a new claim dish.
