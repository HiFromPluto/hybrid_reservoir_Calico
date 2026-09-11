# PocketFill-T1b

Seed occupancy screen on the T1 volumetric recipe. Not a task scout.
T1 `VOLUMETRIC_LUXI` is DEAD; this package does not overwrite T1.

```
python examples/BSimReservoirPlanPocketFillT1b/run_occupancy_screen.py
python examples/BSimReservoirPlanPocketFillT1b/check_pocketfill_t1b.py
```

Frozen prompt: `examples/PocketDish/POCKETFILL_T1B_FROZEN_BUILDER_PROMPT.md`.
Protocol: `PROTOCOL.md`. Report: `results/OCCUPANCY_SCREEN.md`.

Do not retune `QS_KMLA` or \(N_{\mathrm{pack}}\). Do not start NARMA.
`LIVING_DECISION=STOP_AFTER_T1B_SEED`.
