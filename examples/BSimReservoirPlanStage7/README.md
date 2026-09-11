# Plan Stage 7 — input layout rank check

Stage 7 copies frozen Stage 6 and turns on four grounded sources at
predeclared positions. Stage 4 remains FAIL. Stage 5 remains PASS.
Stage 6 remains PASS and is not modified. NARMA-10 is not a Stage 7 gate.
Stage 8 is not started.

This is not the plan's 5-AC / flow-8 / glucose table. See `PROTOCOL.md`.

Run from this directory after compiling:

```
javac -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanStage7.java VoxelAnalyzer.java
```

```
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanStage7.BSimReservoirPlanStage7 sim_config_stage7_driven_seed101.properties
```

Repeat for seeds 202 and 303, then silent seed 404.

Evaluate:

```
python analyze_stage7.py --driven results/stage7_driven_seed101 results/stage7_driven_seed202 results/stage7_driven_seed303 --silent results/stage7_silent_seed404 --pairwise input_pairwise_r.csv --evidence results/stage7_gate_evidence.json --markdown results/GATE_EVIDENCE.md
```

If the gate fails, stop. Do not retune rates or positions after SVD.
Do not start Stage 8.
