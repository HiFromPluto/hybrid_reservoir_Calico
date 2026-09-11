# Plan Stage 4 — input-driven death

Stage 4 copies frozen Stage 3B and adds a dedicated mixed-acid field with a
literature-anchored Hill kill. AC1 remains AHL. The Stage 3B receiver is not
retuned. Stage 5 is not started.

Run from this directory after compiling:

```
javac -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanStage4.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanStage4.BSimReservoirPlanStage4 sim_config_stage4_holdout1.properties
```

Evaluate:

```
python analyze_stage4.py --replicates results/stage4_holdout1_seed101 results/stage4_holdout2_seed202 results/stage4_holdout3_seed303 --silent results/stage4_silent_seed404 --evidence results/stage4_gate_evidence.json
```

See `MANIFEST.md` for mechanism provenance and gates.
See `results/GATE_EVIDENCE.md` for the documented FAIL: keep the 2x
`Total_Deaths` fold; do not treat truncated BSim stdout as incomplete CSVs.
