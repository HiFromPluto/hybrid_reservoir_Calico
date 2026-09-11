# PocketFill-T1

Reduced occupancy screen for the PocketFill recipe (packed volumetric
LuxI on Danino clocks in the T0 garage). Not a task scout.

```
python examples/BSimReservoirPlanPocketFillT1/run_occupancy_screen.py
python examples/BSimReservoirPlanPocketFillT1/check_pocketfill_t1.py
```

Frozen prompt: `examples/PocketDish/POCKETFILL_T1_FROZEN_BUILDER_PROMPT.md`.
Protocol: `PROTOCOL.md`. Report: `results/OCCUPANCY_SCREEN.md`.

Do not retune HybridDish \(K\) or T0 \(J_{\max}\). Do not start NARMA
if `VOLUMETRIC_LUXI` is DEAD.
