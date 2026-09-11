# C1 — packed spatial Object B (NOT_FIG4B)

Frozen packed patch of `DaninoSI_OccupiedDF` cells coupled through a
conservative spatial \(H_e\) field on the P0 pocket+bus mask.

This is **not** Fig. 4b, not Object A, not ACs, not NARMA, not W0,
not a C0 well-mixed identity, and not permission to shrink
\(\mu=0.32\)–\(0.40\).

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).

## Run

```
ant c1-packed
python examples/BSimReservoirPlanDaninoPocketC1/check_c1.py
```
