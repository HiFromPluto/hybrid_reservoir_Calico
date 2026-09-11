# A1 protocol — bounded transducer on C1c_FILLED_POCKET (NOT_FIG4B)

**Gate:** A1 one-way stateful releasing envelope.  
**Status label:** `A1_BOUNDED_TRANSDUCER`.  
**Provenance:** `HYPOTHETICAL_DESIGN_ENVELOPE`. **Not calibrated. Not A1C.**  
**Device:** `C1c_FILLED_POCKET` only. **NOT_FIG4B.**  
**Organism (bacteria-ON):** `DaninoSI_OccupiedDF` (Object B). **NOT_FIG4B.**  
**Not this gate:** Fig. 4b; A0 identity; A2 NARMA/ridge/CHARC; A3 two-way \(h(P)\);
W0; Lentini TX–TL ODE; IPTG/lac receiver; A0 on `C1_OPEN_DILUTE`.

A0 standing remains **PASS** on the memoryless source. This gate replaces
that source with an explicit AC state. It does not rewrite A0 as A1.

C1 packed-spatial standing remains **FAIL**. Do not place this AC on
`C1_OPEN_DILUTE`. C1b island PASS is not the A1 device. C1c filled-pocket
PASS is the device.

Claim freeze:
[`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).  
A0 standing:
[`examples/PocketDish/A0_IDEAL_SOURCE_STANDING.md`](../PocketDish/A0_IDEAL_SOURCE_STANDING.md).  
Design §5.1–5.4 / §14 A1 / §15–16:
[`examples/PocketDish/AC_DANINOPOCKET_PREBUILD_DESIGN.md`](../PocketDish/AC_DANINOPOCKET_PREBUILD_DESIGN.md).

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after A0 PASS, before any A1 integration.

This protocol answers, on an envelope frozen **before traces**:

> Can a named, delayed, leaky, saturating, payload-limited AC on
> `C1c_FILLED_POCKET` release extracellular AHL with a closed mass
> ledger, a \(J(t)\) that is not A0’s \(J_{\max}u(t)\), basal leak at
> \(u=0\), finite payload exhaust, and Object B still occupying the
> C1c class when the AC is off?

Not a digital twin. Not Lentini TX–TL. Not Fig. 4b.

Do not shrink \(\mu=0.32\)–\(0.40\). Do not `TIME_ADJ=60`. Do not retune
Object B. Do not fit \(\tau_{\mathrm{AC}},K,n,J_{\max},J_{\mathrm{leak}},M_0,M_{1/2}\)
after seeing \(I\) or any task score. Do not start A2/W0 from a FAIL.

---

## Named objects (frozen now)

| Object | Role | Status |
|---|---|---|
| `C1c_FILLED_POCKET` | occupancy identity; **A1 device** | PASS (unchanged). NOT_FIG4B |
| `A0_IDEAL_SOURCE` | memoryless \(J=J_{\max}u(t)\) contrast | PASS (unchanged). NOT_FIG4B |
| `A1_BOUNDED_TRANSDUCER` | stateful envelope | this gate. HYPOTHETICAL. NOT_FIG4B |
| `C1_OPEN_DILUTE` | empty chip+bus | FAIL. **Forbidden A1 device.** NOT_FIG4B |
| `C1_PACKED_SPATIAL` | C1 standing | FAIL (unchanged). NOT_FIG4B |
| `C1b_ISLAND` | C0-equivalent voxel | PASS. Not the A1 device. NOT_FIG4B |

Payload species: **extracellular AHL** (same \(H_e\) as Object B).  
Not IPTG. Not an unnamed “signal.”

Architecture citation only (not kinetics, not this ODE, not this
payload):

- Lentini 2014, *Nat. Commun.* **5**, 4012: “AC translates a command the
  bacterium does not need to hear as \(u\).” Not TX–TL. Not hours-scale
  cell-free expression. Not theophylline → αHL → IPTG.
- Lentini 2017, *ACS Cent. Sci.* **3**, 133–140: two-way QS class. **Not
  A1.** No \(h(P)\). That is A3.

Until a lab \(u\to J\) curve exists, this gate is
`HYPOTHETICAL_DESIGN_ENVELOPE`. A1C is later.

---

## Device (unchanged from A0 / C1c)

| Symbol | Frozen value |
|---|---|
| Device | `C1c_FILLED_POCKET` |
| Pocket | \(100\times100\times1.65\) µm, \(N=11000\), \(d=0.5\), \(V_e=16500\) µm³ |
| `D1_spatial` | **800**, \(\mu=0.40\) pocket fluid only |
| Bus | chemical solids |
| Open edge | no-flux |
| AC world | \((50,\ 50,\ 0.825)\) µm, voxel \((25,\ 25,\ 0)\) |
| Far voxel | \((0,\ 0,\ 0)\) |
| Mechanics | OFF |
| Clock | SI-scaled, `rk_dt=0.001`. No `TIME_ADJ` |
| Kernel | N0 `BSimTransportField`; AC `addQuantity` in `depositFluxes` |
| Object B | unchanged |

---

## Envelope (frozen now; HYPOTHETICAL / ENGINEERING)

One immobilized AC. AHL path. No internal TX–TL. No vesicle ODE. No
\(h(P)\). No synthesis (\(J_{\mathrm{syn}}=0\)). No silent refill.
Deterministic; \(N_{\mathrm{AC}}=1\). Injected `BSimRandom` unused. No
`Math.random()`.

State:

| Symbol | Meaning | IC | Units |
|---|---|---|---|
| \(x_1\) | processing / delay | \(0\) | dimensionless, \([0,1]\) target |
| \(M\) | remaining payload | \(M_0\) | field quantity (AHL “molecules”) |

Input \(u(t)\in\{0,1\}\) is the same rectangular command class as A0.

\[
\frac{\mathrm{d}x_1}{\mathrm{d}t}=\frac{u-x_1}{\tau_{\mathrm{AC}}}
\]

\[
J_S=J_{\mathrm{leak}}+J_{\max}\frac{x_1^{n}}{K^{n}+x_1^{n}}\frac{M}{M+M_{1/2}}
\]

\[
\frac{\mathrm{d}M}{\mathrm{d}t}=-J_S,\qquad J_{\mathrm{syn}}=0.
\]

\(M\) is clipped at 0: if \(J_S\Delta t>M\), deposit the remainder and
set \(M=0\), \(J_S=M_{\mathrm{before}}/\Delta t\). No negative payload.
No refill.

On each scheduler interval \([t,t+\Delta t)\), \(u\) is the A0 overlap
(1 if the window intersects, else 0). \(x_1\) is advanced with the
exact linear solution at constant \(u\):

\[
x_1(t+\Delta t)=u+\bigl(x_1(t)-u\bigr)e^{-\Delta t/\tau_{\mathrm{AC}}}.
\]

Then \(J_S\) is evaluated from the updated \(x_1\) and the pre-step
\(M\), and \(J_S\Delta t\) is deposited in **`depositFluxes`** into
voxel \((25,25,0)\). The AC ODE is **not** integrated inside
`bacterium.action`.

### ENGINEERING numbers (frozen now; not from \(I\); not Lentini hours)

| Parameter | Value | Provenance |
|---|---|---|
| \(\tau_{\mathrm{AC}}\) | **2** | SI-scaled ENGINEERING. Visible vs A0 pulse \(\tau=1\) and step window 8. **Not** Lentini hours-scale TX–TL |
| \(n\) | **2** | saturating Hill; ENGINEERING |
| \(K\) | **0.5** | half-activation of \(x_1\); ENGINEERING |
| \(J_{\max}\) | **825** | same A0 mass-budget number (\(\Delta M_{\mathrm{kick}}/\tau_{\mathrm{pulse}}\)). **Not** fitted to colony \(I\) |
| \(J_{\mathrm{leak}}\) | **0.1** | \(\ll J_{\max}\) (\(1.2\times10^{-4}\,J_{\max}\)). ENGINEERING basal |
| \(M_0\) | **5000** | finite payload \(>\) A0 unit pulse 825, \(\lt\) a long `u=1` exhaust |
| \(M_{1/2}\) | **500** | payload saturation; ENGINEERING |
| \(J_{\mathrm{syn}}\) | **0** | no synthesis, no refill |

\(\tau_{\mathrm{AC}}=2\) lives on the Object B SI clock (\(\tau=10\),
\(t_{\mathrm{end}}=1000\)). It is **not** a claim that Lentini 2014
TX–TL runs in two SI-time units.

### Rectangular commands (same windows as A0, plus exhaust)

| Command | \(u(t)\) | \(t_{\mathrm{on}}\) | \(t_{\mathrm{off}}\) | A0 \(\int J\) on that window |
|---|---|---|---|---|
| `U_PULSE` | 1 on the half-open window | 0 | 1 | 825 |
| `U_STEP` | 1 on the half-open window | 0 | 8 | 6600 |
| `U_ZERO` | 0 | — | — | 0 |
| `U_EXHAUST` | 1 on the half-open window | 0 | 50 | n/a (A1 only) |

A1 \(\int J_S\) on `U_PULSE` is **not** required to equal 825. Report
it. Delay, Hill, and payload make it smaller and smeared. That is the
A0 contrast.

---

## Ledger (frozen now)

\[
\Delta M_{\mathrm{release}}=\int_0^{t_{\mathrm{end}}} J_S\,\mathrm{d}t
\]

\[
\Delta M_{\mathrm{field}}=\texttt{sourceAdded}
\]

\[
\Delta M_{\mathrm{payload}}=M_0-M(t_{\mathrm{end}})
\]

Gates \(\le 10^{-3}\) relative:

\[
\frac{|\Delta M_{\mathrm{release}}-\Delta M_{\mathrm{field}}|}{M_\ast}\le 10^{-3},
\quad
\frac{|\Delta M_{\mathrm{payload}}-\Delta M_{\mathrm{release}}|}{M_\ast}\le 10^{-3}.
\]

\(M_\ast=\max(|\Delta M_{\mathrm{release}}|,|\Delta M_{\mathrm{field}}|,|\Delta M_{\mathrm{payload}}|,825,10^{-15})\).

Bacteria-ON AHL residual includes AC sources as in A0:

\[
R=M(0)+\mathrm{synth}+\Delta M_{\mathrm{field}}
-\mu_{\mathrm{loss}}-\gamma_H-\mathrm{outlet}-\mathrm{boundary}-M(t).
\]

---

## A0 contrast (frozen now)

On the **same** \(u(t)\) as A0 `U_PULSE` and `U_STEP`, A1 \(J_S(t)\)
must not be identical to \(J_{\mathrm{A0}}(t)=J_{\max}u(t)\).

**PASS** if, on `A1_FIELD_PULSE_U` and `A1_FIELD_STEP_U`,

\[
\max_t\frac{|J_S(t)-J_{\mathrm{A0}}(t)|}{J_{\max}}\ge 0.10
\]

and \(J_S\) at the first sample with \(t>0\) is **less** than
\(J_{\mathrm{A0}}\) (delay). After `U_PULSE` off (\(t>1\)), \(J_S\)
remains \(>10\,J_{\mathrm{leak}}\) at some recorded sample (memory /
delay tail). Do not retune to enlarge or shrink this gap.

---

## Predeclared arms

Do **not** add A1 on `C1_OPEN_DILUTE`.

| id | bacteria | \(u\) | \(t_{\mathrm{end}}\) | `sample_dt` | required |
|---|---|---|---|---|---|
| `A1_FIELD_STEP_U` | OFF | `U_STEP` | 10 | 0.1 | mass+payload ledger; A0 contrast; He_ac \(>\) He_far while \(0<t\le 8\) |
| `A1_FIELD_PULSE_U` | OFF | `U_PULSE` | 10 | 0.1 | ledger; A0 contrast; report \(\int J_S\) vs A0 825; delay tail |
| `A1_U0_LEAK` | OFF | `U_ZERO` | 10 | 0.1 | \(J_S\approx J_{\mathrm{leak}}\) (\(|\mathrm{mean}J-J_{\mathrm{leak}}|\le 0.01\)); \(\Delta M_{\mathrm{payload}}\approx J_{\mathrm{leak}}t_{\mathrm{end}}\) |
| `A1_PAYLOAD_EXHAUST` | OFF | `U_EXHAUST` | 50 | 0.1 | \(M(t_{\mathrm{end}})/M_0\le 0.10\); \(J_S(t_{\mathrm{end}})\le 0.15\max_t J_S\); ledger; no refill |
| `A1_J0_BACTERIA` | Object B ON, **AC off** (no deposit) | — | 1000 | 0.5 | OSC C1c class \(\|T-63.5833\|/63.5833\le 20\%\). Report vs 61.25 |
| `A1_PULSE_BACTERIA` | Object B ON, envelope on `U_PULSE` | `U_PULSE` | 1000 | 0.5 | ledger; report \(I\) and He vs `A0_J_PULSE_BACTERIA`. **Not** a NARMA/task score. Do not retune |

`rk_dt=0.001`. Mechanics OFF. `AHL_KICK_005`, \(\mu=0.40\) on bacteria-ON.

Optional extra **not identity**, skipped unless invoked:

| id | note |
|---|---|
| `A1_MATCHED_MASS` | HYPOTHETICAL payload-match so \(\int J_S=825\) on the pulse. Not calibrated. Not this identity. **Skip.** |

---

## PASS / FAIL

A1 **PASS** only if all six required arms hold with **no** Object B /
envelope-parameter change after traces, every log/plot/CSV/standing
line says **NOT_FIG4B**, **A1_BOUNDED_TRANSDUCER**, and
**HYPOTHETICAL_DESIGN_ENVELOPE**, and the AC was not run on
`C1_OPEN_DILUTE`.

FAIL: envelope fitted to \(I\); mass mismatch; A1 collapses to A0
without saying so; claiming calibrated or Lentini TX–TL.

If A1 PASSes: **A1C may not start** (no lab \(u\to J\) curve). **A2 may
not start** unless the user later asks for a frozen capacity protocol.
**W0 may not start.** Default A2 = no.

If A1 FAILs: stop. Do not start A1C/A2/W0. Do not retune. Do not move
the AC onto `C1_OPEN_DILUTE`.

---

## Commands

```
ant a1-bounded
python examples/BSimReservoirPlanDaninoPocketA1/check_a1.py
```

Standing: `examples/PocketDish/A1_BOUNDED_TRANSDUCER_STANDING.md`.
