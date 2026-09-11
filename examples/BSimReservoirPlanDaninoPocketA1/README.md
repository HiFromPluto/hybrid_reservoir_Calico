# A1 — bounded transducer on C1c_FILLED_POCKET (NOT_FIG4B)

One-way stateful releasing envelope. Delayed, leaky, saturating,
payload-limited. Payload = extracellular AHL. Status:
**A1_BOUNDED_TRANSDUCER**, **HYPOTHETICAL_DESIGN_ENVELOPE**. Not
calibrated. Not Lentini TX–TL. Not Fig. 4b.

Device is `C1c_FILLED_POCKET` only. Do not run on `C1_OPEN_DILUTE`.

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).

```
ant a1-bounded
python examples/BSimReservoirPlanDaninoPocketA1/check_a1.py
```
