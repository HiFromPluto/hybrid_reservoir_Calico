# A0 — ideal AC source on C1c_FILLED_POCKET (NOT_FIG4B)

Named, immobilized, one-way module. Output is a prescribed molecular
current \(J(t)\) deposited conservatively into the C1c filled-pocket
AHL field. Status label: **A0_IDEAL_SOURCE**. Not a digital twin. Not
Lentini TX–TL. Not Fig. 4b.

Device is `C1c_FILLED_POCKET` only. Do not run on `C1_OPEN_DILUTE`.
C1 packed-spatial standing remains FAIL. C1b island PASS is not this
device.

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).

```
ant a0-ideal
python examples/BSimReservoirPlanDaninoPocketA0/check_a0.py
```
