# Plan Stage 6 — null-model gate (NARMA-10)

Stage 6 copies frozen Stage 5 and tests whether that biology beats a
passive Brownian density null on NARMA-10. Stage 4 remains FAIL. Stage 5
remains PASS and is not modified. Stage 7 is not started.

Protocol, split, ridge rule, and input encoding are frozen in
`PROTOCOL.md` before NRMSE is seen.

Run from this directory after compiling:

```
javac -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanStage6.java VoxelAnalyzer.java
```

Nine production runs, seeds `101/202/303` shared across arms:

```
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanStage6.BSimReservoirPlanStage6 sim_config_stage6_brownian_seed101.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanStage6.BSimReservoirPlanStage6 sim_config_stage6_silent_seed101.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanStage6.BSimReservoirPlanStage6 sim_config_stage6_driven_seed101.properties
```

Repeat for seeds 202 and 303 with the matching config files.

Evaluate:

```
python analyze_stage6.py --brownian results/stage6_brownian_seed101 results/stage6_brownian_seed202 results/stage6_brownian_seed303 --silent results/stage6_silent_seed101 results/stage6_silent_seed202 results/stage6_silent_seed303 --driven results/stage6_driven_seed101 results/stage6_driven_seed202 results/stage6_driven_seed303 --target narma10_target.csv --evidence results/stage6_gate_evidence.json --markdown results/GATE_EVIDENCE.md
```

If the gate fails, stop. Do not add features, retune kinetics, or start Stage 7.
