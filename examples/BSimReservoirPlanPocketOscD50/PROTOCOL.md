# PocketOsc-D50 — cited monolayer occupancy (job protocol)

Frozen 2026-08-21. Architecture:
[`examples/PocketDish/POCKETOSC_D50_FROZEN_BUILDER_PROMPT.md`](../PocketDish/POCKETOSC_D50_FROZEN_BUILDER_PROMPT.md),
[`examples/PocketDish/QS_PACK_PROVENANCE.md`](../PocketDish/QS_PACK_PROVENANCE.md).
QS_*: [`examples/BSimDaninoD1g/PROTOCOL.md`](../BSimDaninoD1g/PROTOCOL.md)
verbatim. ODE clone:
[`examples/BSimReservoirPlanPocketFillT1/fill_model.py`](../BSimReservoirPlanPocketFillT1/fill_model.py)
(`danino_rhs` only). Garage and coupling are **not** T1’s 10 µm write.

No NARMA, ridge, living Java, bus, period-vs-flow, GFP maps, AC,
bath, 0.95 µm, \(d\) sweep, `QS_KMLA` hunt. HybridDish /
`GATE_EVIDENCE.md` / T1–T1c results untouched.

## Question

At Danino Fig. 4b **\(d=0.5\)**, closed SI bulk trap **1.65 µm**,
same D1g QS_*, zero ICs, no bath, no AC: does optical LA occupy?

T1c DEAD was \(d=0.05\) in 10 µm. That is the replay, not the claim.

| Arm | \(d\) | Height | Expect |
|---|---|---|---|
| `D05_H165_CLOSED` | 0.5 | 1.65 µm | unknown (primary) |
| `D005_H10_REPLAY` | 0.05 | 10 µm | T1c-class **DEAD** |

If replay ALIVE: port broke — stop. If claim DEAD:
`LIVING_DECISION=STOP_AFTER_D50`. Do not lower `QS_KMLA`. Do not
raise \(d\).

## Frozen numbers

- Footprint \(100\times100\) µm. Claim height **1.65 µm** (SI bulk).
  Replay height **10 µm** (T1c garage). No 0.95 µm.
- \(d=0.5\) so \(d/(1-d)=1\) on the claim arm. Occupancy state is
  **\(d\)**, not an \(N\) hunt. `CELL_VOL=1` µm³ is ENGINEERING
  converter only (implied \(N=8250\) at 1.65 µm).
- Coupling (Danino SI), extracellular volume \((1-d)V_{\mathrm{trap}}\):
  \(\partial_t C \supset \frac{d}{1-d}\,D_{\mathrm{wall}}(H_i-C)\).
  Do **not** write \(N=50\,000\) into PocketDish-A.
- QS_* , \(\mu=\ln 2/1800\), ICs \(\{0,0,0,0\}\), CONV 602.2: D1g/T1.
- Field \(k=2.76\times10^{-3}/60\) s⁻¹. Closed: leak \(=0\).
- Reduced **0-D** well-mixed (closed; \(L^2/D\sim63\) s \(\ll\) hours).
  No living Java.
- 14400 s, last-7200 s means. Sample 60 s.

## Gates

- Optical: last-7200 s mean \(H(\mathrm{LA})\),
  \(H(x)=x^2/(K_{\mathrm{MLA}}^2+x^2)\), `QS_KMLA=0.01`.
- ALIVE if \(H(\mathrm{LA})\ge0.05\); SATURATED if \(>0.95\).
- Print mean \(C\), LA, LuxI, \(d\), \(d/(1-d)\), \(V_{\mathrm{trap}}\),
  mass residual. Non-negative \(C\). Leak \(\approx0\).
- Replay must stay DEAD (\(H(\mathrm{LA})\) tiny).

`TRANSPORT_MODEL_STATUS`: mass residual \(\le1\%\) on both arms.
