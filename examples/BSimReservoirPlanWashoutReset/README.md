# WashoutReset

Declared chemical / hydrodynamic reset after an occupied quiescent
load. Measures what is left. Not a task scout. No ridge.

See `PROTOCOL.md`.

## Theory (before BSim)

```
python examples/BSimReservoirPlanWashoutReset/check_washout_reset.py --theory
```

## Compile and smoke (not evidence)

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanWashoutReset.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWashoutReset.BSimReservoirPlanWashoutReset sim_config_washout_smoke.properties
```

Smoke is 6 windows (2 load / 2 wash at 8 µm/s / 2 post). Do not score
it as residual evidence.

## Production W=5 (seed 101)

```
python run_washout_jobs.py
python check_washout_reset.py
```

From the repo root:

```
python examples/BSimReservoirPlanWashoutReset/check_washout_reset.py --theory
python examples/BSimReservoirPlanWashoutReset/check_washout_reset.py
```

W=15 follow-on (`WASH_DECAY_W15`, `WASH_FLUSH_W15`) runs only if
`L_STICKS` and `AHL_CLEAR` fire. Do not retune kinetics. Do not
promote a claim dish.
