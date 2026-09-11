# Waveform2 — equal-histogram temporal-order classification

Track E1. Strict temporal-order classification on the frozen HybridDish
/ Narma10b / Waveform claim dish. Classes are `ORDER_A`, `ORDER_B`,
`ORDER_C`. Do not call them sine, square, or triangle. This package
does not rewrite Waveform Overall FAIL.

Protocol, hashes, templates, splits, and the u-only stop are frozen in
`PROTOCOL.md` before any biology ridge.

## Status

**u-only audit FAIL.** MEAN/POWER/VARIANCE/MOMENTS/START_U sit at
chance. RAW_U_8 test macro AUC is `0.6058` (need `≥ 0.90`). ORDER_A
and ORDER_C are circular shifts of one another, so the three orders
are not distinct. The task is void. **BSim was not run.** Seeds
202/303 were not started. Kinetics were not retuned.

See `results/WAVEFORM2_SCOUT.md` and `results/u_only_baselines.md`.

## Generate and check (no BSim)

From the repo root:

```
python examples/BSimReservoirPlanWaveform2/generate_waveform2_inputs.py
python examples/BSimReservoirPlanWaveform2/check_waveform2.py --u-only
```

Do not regenerate after seeing reservoir AUC. Do not start a new
class/phase draw until this development set is retained.

## Compile (not authorized until u-only PASS)

From this directory after a future go-ahead:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanWaveform2.java VoxelAnalyzer.java
```

Smoke (HybridDish timing, not evidence):

```
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveform2.BSimReservoirPlanWaveform2 sim_config_waveform2_smoke.properties
```

Scout would be driven + Brownian + silent at seed 101 only. Do not run
it on this development set.
