# PocketFill-T1c — door vs circuit (job protocol)

Frozen 2026-08-21. Same QS_* as T1/T1b. Checkpoint:
[`examples/PocketDish/CHECKPOINT.md`](../PocketDish/CHECKPOINT.md).

No NARMA, no `QS_KMLA` hunt, no extra \(W\) after occupancy, no
membrane BSim LuxI.

## Question

Does the cited volumetric pack occupy a **closed** pocket? If not,
the street is not the QS bug.

| Arm | Geometry | Expect |
|---|---|---|
| `VOLUMETRIC_FLUSH` | T1 `OPEN_BUS_3` | DEAD (replay) |
| `VOLUMETRIC_CLOSED` | pocket only, no-flux | unknown (primary) |
| `VOLUMETRIC_W20` | T2 `W20_L20` | unknown; predeclared |

Kill: closed DEAD → stop Track C. Flush ALIVE → port broke.

## Frozen numbers

T1/T1b QS_*, \(N_{\mathrm{pack}}=5000\), Danino \(k=2.76\times10^{-3}/60\),
14400 s, last-7200 s \(H(\mathrm{LA})\ge 0.05\).
