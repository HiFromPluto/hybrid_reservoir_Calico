# PocketOsc-Flow — period vs channel flow (job protocol)

Frozen 2026-08-21 **before traces**. Architecture:
[`examples/PocketDish/POCKETOSC_FLOW_FROZEN_BUILDER_PROMPT.md`](../PocketDish/POCKETOSC_FLOW_FROZEN_BUILDER_PROMPT.md).
Clone ODE:
[`examples/BSimReservoirPlanPocketOscD50/osc_model.py`](../BSimReservoirPlanPocketOscD50/osc_model.py)
(`danino_rhs` verbatim). Add a bus leak. Do not retune QS_*.

No NARMA, GFP maps, 0.95 µm, \(N\) hunt, `W20`, AC, bath,
height/\(d\) sweep, HybridDish \(k=0.0033\). D50 / T1 / HybridDish
/ `GATE_EVIDENCE.md` untouched.

## Question

At \(d=0.5\), SI bulk **1.65 µm**, same D1g QS_*, zero ICs: does
optical LA stay occupied with an open nutrient channel, and does
period **lengthen** as channel flow rises (Danino 180–296 µm/min →
52–90 min; high flow 90 ± 6 min, low flow 55 ± 6 min)?

Matching 55/90 min is **not** a retune ticket.

| Arm | \(d\) | \(h\) | \(v\) | Role |
|---|---|---|---|---|
| `CLOSED_D50` | 0.5 | 1.65 µm | 0 | occupancy control; expect ALIVE |
| `FLOW_180` | 0.5 | 1.65 µm | 180 µm/min = 3.00 µm/s | low-flow identity |
| `FLOW_296` | 0.5 | 1.65 µm | 296 µm/min = 4.93 µm/s | high-flow identity |

## Bus leak (frozen; SI has no \(\mu(v)\))

Danino SI: \(-\mu H_e\) “models dilution of external AHL by external
fluid flow.” They **vary** \(\mu\) as a free rate. They do **not**
state \(\mu(v)\). Experimental \(v\) is channel speed
(180–296 µm/min). The SI “~100 µm/s” is the wave-trap movie, not
this bulk identity.

**ENGINEERING geometric exchange** (two resistances: diffusion across
trap depth, advection at the full-width mouth). Not fitted to 55/90.

\[
\mu_{\mathrm{bus}}(v)=\frac{A_{\mathrm{open}}}{V_{\mathrm{ext}}}
\frac{v D}{v L_y+D},\qquad \mu_{\mathrm{bus}}(0)=0.
\]

| Symbol | Value | Class |
|---|---|---|
| \(L_x\times L_y\times H\) | \(100\times100\times1.65\) µm | TAKEN height/footprint |
| \(d\) | 0.5 | TAKEN (Fig. 4b) |
| \(V_{\mathrm{ext}}=(1-d)V_{\mathrm{trap}}\) | 8250 µm³ | derived |
| \(A_{\mathrm{open}}=L_x H\) | 165 µm² | full-width bulk opening |
| \(D\) | 159 µm²/s | same AHL field as D50 |
| \(v\) | 0 / 3.00 / 4.93 µm/s | TAKEN channel speeds |
| Hydrolysis \(k\) | \(2.76\times10^{-3}/60\) s⁻¹ | D50; **not** the flow knob |

Field (0-D, well-mixed trap):

\[
\partial_t C=\frac{d}{1-d}D_{\mathrm{wall}}(H_i-C)-k C-\mu_{\mathrm{bus}}(v)\,C.
\]

Print \(\mu_{\mathrm{bus}}\) on every arm. Do not change this formula
after occupancy or periods. Reduced **0-D** (closed D50 is the \(v=0\)
limit; port check). No living Java.

## Frozen numbers (QS_*)

D50/D1g verbatim: QS_*, \(\mu_{\mathrm{growth}}=\ln 2/1800\), ICs
zero, CONV 602.2, `QS_KMLA=0.01`. Implied \(N=8250\) ENGINEERING.

## Drive and period rule (frozen before traces)

- Integrate **43200 s** (12 h). Sample 60 s.
- Discard **10800 s** (3 h) before period.
- Occupancy: last-7200 s mean \(H(\mathrm{LA})\ge 0.05\) ALIVE;
  \(>0.95\) SATURATED. Same \(H\) as D50.
- `CLOSED_D50` port: also print D50 window (last-7200 s of the first
  14400 s). That window must stay ALIVE. If DEAD, the port broke —
  stop.
- Period from **LA peaks** (same rule both flow arms): local maxima
  after discard, minimum spacing **600 s**, prominence \(\ge 10^{-4}\)
  µM. Period = mean inter-peak interval. Fewer than **3** peaks →
  `NO_PERIOD`.
- Identity (only if both flow arms ALIVE **and** both have a period):
  \(T(296)>T(180)\). Absolute 55 ± 6 / 90 ± 6 min is **not** a gate.

Open occupancy DEAD → `STOP_AFTER_FLOW_DEAD`.
Occupied but no period or wrong order → `STOP_AFTER_FLOW_NO_IDENTITY`.
Do not ease QS_*. Do not raise \(d\). Do not close the door and call
it identity.
