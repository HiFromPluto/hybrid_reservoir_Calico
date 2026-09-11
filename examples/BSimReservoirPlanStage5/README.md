# Plan Stage 5 — luminescence and feature set

Stage 5 copies frozen Stage 4 and adds a real per-cell luminescence state
downstream of the frozen Stage 3B receiver. Stage 4 remains FAIL. Stage 6
is not started.

Run from this directory after compiling:

```
javac -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanStage5.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanStage5.BSimReservoirPlanStage5 sim_config_stage5_holdout1.properties
```

Evaluate:

```
python analyze_stage5.py --replicates results/stage5_holdout1_seed101 results/stage5_holdout2_seed202 results/stage5_holdout3_seed303 --silent results/stage5_silent_seed404 --evidence results/stage5_gate_evidence.json
```

See `MANIFEST.md` for the feature contract and dropped columns.
