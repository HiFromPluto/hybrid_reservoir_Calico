# C1c — filled-pocket spatial occupancy (NOT_FIG4B)

Cells throughout the 100×100×1.65 µm pocket so extracellular AHL volume
is the colony’s media. Bus voxels are chemical solids on identity.
`D1_spatial=800`. \(\mu=0.40\) on pocket \(V_e\) only.

C1 standing remains FAIL (`C1_OPEN_DILUTE`). C1b remains island PASS.
This is **not** Fig. 4b, not A0, not W0, and not a retune of Object B
onto the empty chip+bus.

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).

```
ant c1c-filled
python examples/BSimReservoirPlanDaninoPocketC1c/check_c1c.py
```
