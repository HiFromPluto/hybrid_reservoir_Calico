# C1b — island C0-equivalent occupancy (NOT_FIG4B)

Chemically isolated \(V_e\) (C0 volumes, \(D_1=0\)). The open P0 chip
(`C1_OPEN_DILUTE`) remains FAIL and is not identity.

C1 standing is not rewritten as PASS. A0 is not automatic.

Standing: [`C1B_ISLAND_STANDING.md`](../PocketDish/C1B_ISLAND_STANDING.md).

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).

```
ant c1b-island
python examples/BSimReservoirPlanDaninoPocketC1b/check_c1b.py
```
