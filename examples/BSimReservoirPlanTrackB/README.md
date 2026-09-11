# Track B — 5 one-way ACs, synthetic 5-channel patients

Track B copies frozen Stage 6 and tests H5: five mixed chemical inputs
versus a scalar AHL mix on a predeclared patient-classification task.
A Stage 9-style field-only gate is required. If the plumes classify the
patient, Overall is FAIL. Do not then widen the dish, add ACs, or switch
to accuracy.

Stage 6 NARMA-10 remains PASS. Track A Mackey–Glass remains PASS. Track A
Lorenz remains FAIL. Stage 9 product-bit FAIL is not rewritten. Kinetics
are not retuned.

Protocol, layout, class vector, SHA-256s, splits, and ridge rule are
frozen in `PROTOCOL.md` before AUC is seen.

## Compile and run

From this directory after compiling:

```
javac -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanTrackB.java VoxelAnalyzer.java
```

Nine production runs, seeds `101/202/303`. Silent is copied from Stage 6
and is not rerun.

```
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanTrackB.BSimReservoirPlanTrackB sim_config_trackb_driven_seed101.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanTrackB.BSimReservoirPlanTrackB sim_config_trackb_singlesite_seed101.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanTrackB.BSimReservoirPlanTrackB sim_config_trackb_brownian_seed101.properties
```

Repeat for seeds 202 and 303 with the matching config files.

Evaluate from the repo root:

```
python examples/BSimReservoirPlanTrackB/check_trackb.py
```

If gate 2 or 3 fails, write FAIL and stop. Do not start a 2000×1000 copy.
