# Waveform2b — distinct cyclic-orbit temporal order

Track E1 follow-on. Waveform2 remains TASK_VOID (ORDER_A ≡ ORDER_C).
This package uses `ORDER_UNI`, `ORDER_DOWN`, `ORDER_UP` on the same
claim dish and the same shared amplitude multiset. Do not call them
sine, square, or triangle.

## Status

**u-only audit FAIL.** Orbit gate PASS (pairwise intersections empty).
MEAN/POWER/VARIANCE/MOMENTS/START_U at chance. RAW_U_8 test macro AUC
`0.5298` (need `≥ 0.90`). **BSim was not run. Java was not copied.**
Seeds 202/303 were not started. Kinetics were not retuned.

See `results/WAVEFORM2B_SCOUT.md` and `results/u_only_baselines.md`.

## Generate and check (no BSim)

From the repo root:

```
python examples/BSimReservoirPlanWaveform2b/generate_waveform2b_inputs.py
python examples/BSimReservoirPlanWaveform2b/check_waveform2b.py --u-only
```

Do not regenerate after seeing reservoir AUC. Do not rewrite templates
after the RAW_U_8 score. Waveform2 stays the retained TASK_VOID draw.
