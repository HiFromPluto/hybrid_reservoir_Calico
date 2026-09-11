# D1 — Java circuit parity of DaninoSI_OccupiedDF (NOT_FIG4B)

Java constant-delay bulk DDE matched to frozen D0b Python fixtures.

This ports **Object B** (`DaninoSI_OccupiedDF`), not Object A
(`Danino2010_Fig4b_bulk_twin`). D0 remains FAIL. D0b remains
FAIL_NO_IDENTITY. This is not a Fig. 4b reproduction.

Claim freeze:
[`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).

## Run

```
ant d1-parity
python examples/BSimReservoirPlanPocketOscSI_D1/check_d1.py
```

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).
