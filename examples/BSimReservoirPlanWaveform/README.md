# Waveform classification — sine / square / triangle on Stage 6

Waveform classification copies frozen BenchA Java onto the Stage 6
dish. Three AHL templates. Brownian, silent, and field-only AHL are
required. If the plume ranks the waveform, Overall is FAIL. Do not
then add ACs or switch to accuracy.

Stage 6 NARMA-10 remains PASS. BenchA Mackey–Glass remains PASS.
BenchA Lorenz remains FAIL. Track B remains FAIL. Kinetics are not
retuned.

Protocol, class vector, SHA-256s, splits, templates, and ridge rule
are frozen in `PROTOCOL.md` before AUC is seen.

## Compile and run

From this directory after compiling:

```
javac -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanWaveform.java VoxelAnalyzer.java
```

Six production runs, seeds `101/202/303`. Silent is copied from
Stage 6 and is not rerun.

```
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveform.BSimReservoirPlanWaveform sim_config_waveform_driven_seed101.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveform.BSimReservoirPlanWaveform sim_config_waveform_brownian_seed101.properties
```

Repeat for seeds 202 and 303 with the matching config files.

Evaluate from the repo root:

```
python examples/BSimReservoirPlanWaveform/check_waveform.py
```

If gate 2 fails, write FAIL and stop. Do not start Track B or a wider
dish.
