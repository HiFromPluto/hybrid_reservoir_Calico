# C0 — well-mixed Object B coupling (NOT_FIG4B)

Explicit cells with `DaninoSI_OccupiedDF` intracellular state exchanging
AHL with one well-mixed \(H_e\) compartment. Identity at \(d=0.5\) must
reduce to the frozen D1/D0b 4-DDE.

This is **not** Fig. 4b, not Object A, not C1, not P0 packing, not ACs,
not NARMA.

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).

## Run

```
ant c0-wellmixed
python examples/BSimReservoirPlanDaninoPocketC0/check_c0.py
```
