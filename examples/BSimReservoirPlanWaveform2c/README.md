# Waveform2c — aligned equal-histogram temporal order

Track E1 aligned follow-on. Waveform2 remains TASK_VOID (ORDER_A ≡
ORDER_C orbit). Waveform2b remains the disjoint-orbit / random-phase
negative (RAW_U_8 `0.5298`). This package keeps the same claim dish
and shared amplitude multiset, freezes phase at 0, and tests whether
three aligned 8-vectors are linearly distinct.

Classes are `ORDER_UNI`, `ORDER_DOWN`, `ORDER_UP`. Do not call them
sine, square, or triangle. This is not a cyclic-invariant task.
Waveform `GATE_EVIDENCE.md` is not edited. C1 remains DEFER. E0.3
remains NO_STORY_MOVE.

## Status

**u-only PASS. Occupancy ALIVE.** Driven block macro AUC
`0.8352 ± 0.0138` beats Brownian `0.5369 ± 0.0463`, silent `0.5000`,
and MOMENTS `0.5000`. Field-only `0.8864` is ≥ driven: the aligned
order is already in the plume. Do not add ACs. Do not rewrite
Waveform1 Overall FAIL. This is not shift-invariant classification
evidence. Kinetics were not retuned.

See `PROTOCOL.md` and `results/WAVEFORM2C_SCOUT.md`.

## Generate and check

From the repo root:

```
python examples/BSimReservoirPlanWaveform2c/generate_waveform2c_inputs.py
python examples/BSimReservoirPlanWaveform2c/check_waveform2c.py --u-only
```

Hashes in `PROTOCOL.md` were recorded before ridge. Do not regenerate
after seeing reservoir AUC. Do not rewrite templates. Do not modify
Waveform2 or Waveform2b.

## Compile and smoke (authorized after u-only PASS)

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanWaveform2c.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveform2c.BSimReservoirPlanWaveform2c sim_config_waveform2c_smoke.properties
```

Scout seed 101:

```
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveform2c.BSimReservoirPlanWaveform2c sim_config_waveform2c_driven_seed101.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveform2c.BSimReservoirPlanWaveform2c sim_config_waveform2c_brownian_seed101.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveform2c.BSimReservoirPlanWaveform2c sim_config_waveform2c_silent_seed101.properties
```

Do not retune K, n, tau_R, tau_L, source rate, clamp, mortality, flow,
or layout. Claim dish stays CENTER / FLOW=0.
