# Plan Stage 9 — hybrid classification on the Stage 6 dish

Stage 9 tests whether the one-way hybrid (AHL AC → E. coli) classifies a
predeclared nonlinear product bit better than Brownian density and a
linear readout of the AC carrier field. Silent must sit at chance.

This is not a 4-AC or 5-AC claim. Stage 7 remains FAIL and is not used.
Stage 8 is skipped (Danino already archived as Stage 3). Stage 6 remains
PASS and is not modified. NARMA-10 is already done; this stage is
classification on those same CSVs.

Protocol, split, ridge rule, and the product-bit task are frozen in
`PROTOCOL.md` before AUC is seen.

## Evaluate on existing Stage 6 CSVs

Do not rerun BSim unless a CSV is missing or ragged.

From this directory:

```
python analyze_stage9.py --brownian ../BSimReservoirPlanStage6/results/stage6_brownian_seed101 ../BSimReservoirPlanStage6/results/stage6_brownian_seed202 ../BSimReservoirPlanStage6/results/stage6_brownian_seed303 --silent ../BSimReservoirPlanStage6/results/stage6_silent_seed101 ../BSimReservoirPlanStage6/results/stage6_silent_seed202 ../BSimReservoirPlanStage6/results/stage6_silent_seed303 --driven ../BSimReservoirPlanStage6/results/stage6_driven_seed101 ../BSimReservoirPlanStage6/results/stage6_driven_seed202 ../BSimReservoirPlanStage6/results/stage6_driven_seed303 --ahl ../BSimReservoirPlanStage6/input_ahl_narma200.txt --evidence results/stage9_gate_evidence.json --markdown results/GATE_EVIDENCE.md
```

Java is not copied here. Preview, if needed later, must keep Stage 6
physics, restore the Stage 7 camera sequence, and draw marker spheres at
the AHL `(500, 250, 5)` and acid `(300, 375, 5)` sites only. Do not draw
the Stage 7 four-corner layout.

If the gate fails, document why. Do not add features, retune kinetics,
move ACs, import vesicles/Danino, or start two-way coupling.
