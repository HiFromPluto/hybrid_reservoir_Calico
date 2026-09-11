# Lane A protocol — occupied millimetre dish on N0 (`LaneA_N0_OCCUPANCY`)

**Gate:** `LaneA_N0_OCCUPANCY`  
**Status label:** `LANE_A_OCCUPIED_MILLIMETRE`  
**Object:** HybridDish-class millimetre monolayer, rebuilt on the N0 kernel.  
**Not this gate:** paper 1 / `reservoir_new`; HybridDish Overall / `GATE_EVIDENCE`; NARMA identity; CHARC; Lane B chemotaxis; C1c / A0 / A1; Object B / Fig. 4b; P0 Hertzian; calibrated or Lentini AC.

**frozen_before_traces:** true  
**Frozen:** 2026-08-28, after reading the paper-objects / HybridDish / N0 notes, **before** any Lane A integration.

This protocol answers, on an identity frozen **before traces**:

> Can an occupied one-way AHL → Hill \(R\) → slow \(L\) monolayer exist on
> the HybridDish millimetre geometry with conservative N0 transport and an
> exact mass ledger, on biological windows, with motility off?

Beating a linear map of voxel AHL is a **carrier diagnostic**, not this
job. Occupancy ALIVE is required before any later task number exists.
This job does **not** score NARMA, rewrite Narma10b F408 ~0.928, or
paste scores into paper 1.

Not a digital twin. Not Danino. Not Fig. 4b. Not C1c.

---

## Named objects (do not glue)

| Object | Role here |
|---|---|
| `reservoir_new` (paper 1) | **untouched**. 6 s CHARC dish. Do not rewrite tables. |
| HybridDish / Narma10b | **frozen predecessor**. Claim geometry copied. Overall lines not edited. |
| C1c + A0 + A1 | **not this chat**. 100×100×1.65 µm pocket. NOT_FIG4B. |
| P0 | **not this chat**. Packed rods, chemistry off. |
| Lane B | **not started**. Motility-as-computer. |
| `LaneA_OccupiedMillimetre` | **this job**. N0 successor of HybridDish’s occupancy question. |

Every log, CSV, plot, and standing line says `LANE_A_OCCUPIED_MILLIMETRE`.  
Not Fig. 4b. Not C1c. Not paper 1.

---

## Scientific question (this gate only)

On the frozen millimetre dish, N0 field, A0-class ideal current
\(J=J_{\max}u\), immobilized Hill receivers:

1. Does a closed conservative arm close in the N0 mass band?
2. Does a predeclared short pulse train occupy the Hill (mean \(R \ge 0.05\))?
3. Does commanded / deposited / remaining / decay close on the driven arm?
4. Do Brownian and silent arms run, and is silent **not** occupied like driven?

**Not asked:** NARMA NRMSE; CHARC; living-layer vs field; 408-D ridge;
motility-on; raising \(J_{\max}\) after seeing \(R\).

---

## Frozen identity (HybridDish claim geometry — do not “improve”)

| Item | Frozen value | Provenance |
|---|---|---|
| Dish | \(1000\times 500\times 10\) µm | TAKEN HybridDish |
| Source | CENTER \((500,\ 250,\ 5)\) µm | TAKEN HybridDish |
| Containing voxel | \((i,j,k)=(25,\ 12,\ 0)\) on the \(50\times 25\times 1\) grid | ENGINEERING, from box size |
| FLOW | **0** | TAKEN HybridDish claim. Through-chamber 8 µm/s is not this gate. |
| BC | no-flux (solid walls), not leaky, not periodic | TAKEN HybridDish `NO_FLUX` |
| PDE grid | \(50\times 25\times 1\) (20×20×10 µm boxes) | TAKEN HybridDish |
| Readout grid | \(20\times 10\) **state class**, **not** the PDE grid | TAKEN HybridDish class |
| 408-D / 4×2 death bins | **not copied**. Occupancy is cell-mean \(R\), not a ridge. | this job |
| Hill \(K,n\) | \(1.6\) µM, \(2\) | TAKEN HybridDish |
| \(\tau_R,\tau_L\) | \(15\) s, \(1500\) s | TAKEN HybridDish |
| Windows | \(300\) s; pulse \(75\) s at the start of each window | TAKEN HybridDish duty |
| \(u\) during pulse | **0.5** | TAKEN HybridDish warmup command; NARMA range upper end. Not fitted to \(R\). |
| \(J_{\max}\) | \(1.28\times 10^{8}\) molecules/s | TAKEN HybridDish `field.ahl.source.rate`. ENGINEERING copy. **Not** fitted after \(R\). |
| AHL \(D\) | \(159\) µm²/s | TAKEN HybridDish. Not Dilanji 550. Sensitivity is later. |
| AHL \(k\) | \(0.0033\) s⁻¹ | TAKEN HybridDish |
| Conversion | \(1\) µM \(= 602.2\) molecules/µm³ | TAKEN HybridDish |
| Kernel | N0 `BSimTransportField` + `BSimStepScheduler` | N0 standing |
| Motility | **OFF** on living arms | HybridDish score was \(R,L\), not chemotaxis |
| Growth / division / death | **OFF** | occupancy of \(R\), not clamp/acid ridge |
| AC | A0-class \(J=J_{\max}u\). No internal ODE. | HybridDish tap / A0 law on this dish, not C1c \(J_{\max}=825\) |
| Stateful envelope | **not this job**. If added later: HYPOTHETICAL, not fitted to NRMSE. | |
| `dt` | \(0.05\) s | TAKEN HybridDish |
| `TIME_ADJ` | **unused** | |
| Danino SI \(t_{\mathrm{end}}=1000\), \(d=0.5\) | **unused** | |
| NARMA / CHARC / IPC / KR | **off** | |

Do not retune \(J_{\max}\), \(K\), \(n\), \(\tau\), \(D\), \(k\), layout, or
FLOW after traces. A new named extra may ask a different question; it
cannot rewrite this protocol or a closed HybridDish Overall line.

N0 spreads at **nominal** \(D\), not the legacy ~2× sweep-dependent kernel.
Mean \(R\) may differ from HybridDish ~0.18–0.21. That is allowed. Do
**not** raise \(J_{\max}\) to recover a HybridDish number.

---

## Why warmup is not 18000 s

HybridDish task dishes warm \(L\) (\(\tau_L=1500\) s) for 18000 s.
This gate is occupancy of \(R\) (\(\tau_R=15\) s). AHL decay
\(\tau=1/k\approx 303\) s.

**Frozen:** empty dish at \(t=0\); no 18000 s \(L\) warmup; pulse train
starts immediately. Occupancy is averaged on the **second half** of the
short train, after several AHL decay times of drive.

---

## Predeclared short drive (not NARMA)

Command name: `U_PULSE_TRAIN`.

| Item | Value |
|---|---|
| \(N_{\mathrm{win}}\) | **8** |
| Window | 300 s |
| Pulse | 75 s, \(u=0.5\), then \(u=0\) for 225 s |
| \(t_{\mathrm{end}}\) | **2400** s |
| Samples | HybridDish-class **16** per window: \(t_{\mathrm{in}}=0,20,\ldots,280\) s and the last step (\(299.95\) s). Window index from the scheduler step, not a \(t=300\) fencepost. |
| Occupancy windows | **4–7** (0-based), 16 samples each, **64** samples |
| Expected commanded mass | \(J_{\max}\times 0.5\times 8\times 75=\mathbf{3.84\times 10^{10}}\) |

This is not NARMA-10 `u`. A later extra may **replay** a frozen HybridDish
`u` as a labelled comparison. It is not “paper 2 improved 0.928.”

Silent command: `U_ZERO` (same duration, same cells, \(J=0\)).  
Closed command: `U_CLOSED_IMPULSE` — \(u=1\) on \([0,1)\), \(\Delta M=J_{\max}\),
decay **off**, \(t_{\mathrm{end}}=10\) s, no cells.

---

## Occupancy threshold (frozen now, before traces)

HybridDish NOT_SCORED rule, copied:

\[
\text{ALIVE if } \overline{R} \ge 0.05,\qquad
\text{DEAD if } \overline{R} < 0.05.
\]

\(\overline{R}\) is the **cell-mean** receiver, then time-averaged over
the 64 occupancy-window samples (windows 4–7, HybridDish-class 16
samples each). Not a 20×10 voxel mean. Not \(L\). Not field AHL.

If driven is DEAD: this gate **FAIL**. No living-layer task number exists.
Do **not** raise \(J_{\max}\). Do not drop windows after seeing \(R\).

Silent must be reported and must **not** occupy like driven:
\(\overline{R}_{\mathrm{silent}} < 0.05\) and
\(\overline{R}_{\mathrm{driven}} > \overline{R}_{\mathrm{silent}}\).

Brownian has **no Hill**. Report \(\overline{R}=0\). Density may be
logged; it is not an occupancy competitor.

---

## Population (frozen)

| Item | Value |
|---|---|
| \(N\) | 1800 immobilized receivers (living / silent) |
| Seed box | \(x\in[300,700]\), \(y\in[150,350]\), \(z=5\) | TAKEN HybridDish |
| RNG seed | **101** | TAKEN HybridDish driven seed |
| Brownian | 1800 `BSimParticle` spheres, thermal drift only, same seed box |

Receivers are **not** `BSimBacterium` run–tumble. Positions frozen after
placement. No chemotaxis. That is Lane B, not this gate.

---

## AC law (A0-class on this dish)

\[
J(t)=J_{\max}\,u(t)
\]

Exact overlap of each scheduler interval \([t,t+\Delta t)\) with the
frozen on-window. Deposit \(\Delta M = J_{\max} u \Delta t_{\mathrm{on}}\)
with `addQuantity` on voxel \((25,12,0)\) in **`depositFluxes`**. Not
inside a cell action. Not C1c \(J_{\max}=825\). Not a vesicle/Hill-gate
envelope.

Payload: extracellular AHL. Not IPTG. Not an unnamed “signal.”

---

## Receiver ODEs (frozen)

\[
\frac{\mathrm{d}R}{\mathrm{d}t}=\frac{1}{\tau_R}\left(\frac{C^n}{K^n+C^n}-R\right),\qquad
\frac{\mathrm{d}L}{\mathrm{d}t}=\frac{R-L}{\tau_L}
\]

\(C\) in µM from N0 concentration / 602.2. Exact exponential step for
\(R\) (HybridDish). Exact exponential step for \(L\) with \(R\) held
over \(\Delta t\) (same \(\tau_L\); not an Euler retune). \(L\) is
**reported**, not a gate: 2400 s is not \(L\) steady state.

Cells do **not** consume AHL. One-way receiver.

---

## Ledger (frozen now)

N0 residual:

\[
R_{\mathrm{mass}} = M(0) + \texttt{sourceAdded} - \texttt{decayLoss}
- \texttt{outletLoss} - \texttt{boundaryLoss} - M(t)
\]

Closed dish, FLOW=0, no-flux: outlet = boundary = 0. \(M(0)=0\) after
`ledger.reset()`.

### Closed conservative arm

Decay **off**. Impulse \(\Delta M = J_{\max}\) over \([0,1)\). \(t=10\) s.

**PASS** if

\[
\frac{\lvert \Delta M_{\mathrm{command}} - \texttt{sourceAdded}\rvert}{M_\ast}\le 10^{-3}
\]

and

\[
\frac{\lvert R_{\mathrm{mass}}\rvert}{M_\ast} \le 10^{-11}
\]

(N0 standing band; the 0.1% gate was the N0 *requirement*, the suite
used \(10^{-11}\)).  
\(M_\ast=\max(\lvert\Delta M_{\mathrm{command}}\rvert,\lvert\texttt{sourceAdded}\rvert,\lvert M(t)\rvert,10^{-15})\).

### Driven pulse train

Decay **on** (\(k=0.0033\)). Commanded mass \(3.84\times 10^{10}\).

**PASS** if commanded vs deposited \(\le 10^{-3}\) relative, and
\(R_{\mathrm{mass}}\) in the N0 band \(\le 10^{-11}\). Report remaining
and decay. Do not retune \(k\) to shrink decay.

Silent: commanded = deposited = 0.

---

## Predeclared arms

| id | cells | \(u\) | decay | \(t_{\mathrm{end}}\) | required |
|---|---|---|---|---|---|
| `LANE_A_CLOSED_CONSERVATIVE` | none | `U_CLOSED_IMPULSE` | **0** | 10 | mass replay; N0 residual band |
| `LANE_A_DRIVEN` | 1800 immobilized Hill | `U_PULSE_TRAIN` | 0.0033 | 2400 | occupancy ALIVE; ledger closes; report mean \(L\) |
| `LANE_A_SILENT` | 1800 immobilized Hill | `U_ZERO` | 0.0033 | 2400 | mean \(R < 0.05\); below driven |
| `LANE_A_BROWNIAN` | 1800 particles, motility thermal | `U_PULSE_TRAIN` field (unread) | 0.0033 | 2400 | mean \(R = 0\) (no receiver); run completes |

Scheduler: `BSimStepScheduler`. Not legacy `BSim.export()` inclusive loop.

---

## PASS / FAIL

`LANE_A_OCCUPIED_MILLIMETRE` **PASS** only if all four arms hold with
**no** \(J_{\max}/K/n/\tau/D/k\) change after traces, PROTOCOL was
frozen first, every artifact is labelled `LANE_A_OCCUPIED_MILLIMETRE`,
and none of the forbidden edits occurred.

FAIL: occupancy chased by raising \(J\); HybridDish evidence edited;
motility-on to move a future NRMSE; C1c geometry; CHARC names without a
correct metric protocol; NARMA as identity.

If PASS: a later named extra may replay frozen HybridDish \(u\) or ask
a task, still occupancy-first, still not a rewrite of Narma10b Overall.
Lane B is **not** started from this standing.

If FAIL: stop. Do not raise \(J_{\max}\). Do not edit HybridDish
`GATE_EVIDENCE`. Do not start NARMA.

---

## Explicit non-claims

- HybridDish Overall unchanged (Narma10b F408 ~0.928 stays).
- Paper 1 unchanged.
- Lane B not started.
- Not Fig. 4b. Not C1c. Not Object B. Not P0.
- AC is not calibrated and not Lentini TX–TL.
- Occupancy ALIVE is not a task PASS.

---

## Commands

```
ant lane-a-occupancy
python examples/LaneA_OccupiedMillimetre/check_lane_a.py
```

Standing (after traces):
`examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md`.
