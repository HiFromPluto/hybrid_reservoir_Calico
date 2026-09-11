# LC0 protocol — replication visible on a Warren sizer clock

**Gate:** `LC0_REPLICATION_VISIBLE`  
**Status label:** `LC0_REPLICATION_VISIBLE`  
**Object:** `LANE_C_LIVING_CLOCKS`  
**Device:** `LC0_OPEN_BATH` (named well-mixed open box; **not** the
100 µm P0 pocket; **not** Lane A’s occupied AHL dish)  
**Parents:** Lane A closed (cite only). Lane B occupation FAIL kept
(cite only). P0 packed-device PASS kept (different machine). Job 6
kept (death is **off** here).  
**frozen_before_traces:** true  
**Frozen:** 2026-08-29, before any LC0 Java.

This protocol answers, on an identity frozen **before traces**:

> On a named open bath, Warren sizer clock frozen before traces,
> do births become countable on 300 s windows and match a predeclared
> \(N_{\mathrm{ever}}(t)\) over one generation — and does the
> growth-off arm stay flat?

**One process:** replication only.

Name collision, kept visible:
[`../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md`](../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md)
§8 “Lane C” is the packed-pocket stack (P0 / C1c / A1). This object is
`LANE_C_LIVING_CLOCKS`. Do not merge them.

Investigation (corrected knobs below):
[`../PocketDish/LANE_C_INVESTIGATION.md`](../PocketDish/LANE_C_INVESTIGATION.md).

---

## 0. What this extra is not

- Not Lane A with motility or death on. Not a Lane A NRMSE / AXIS_HOLDS.
- Not B5, not Lane B occupation, not millimetre \(\kappa\).
- Not Paper 1 `reservoir_new` on N0. Not CHARC / IPC / NARMA.
- Not P0 Hertzian rods, not extrusion, not a chemostat.
- Not Schink death. Death is **OFF**. Do not name Schink as identity.
- Not 1800 s dressed as literature.
- Not a glucose PDE, not AHL, not Hill, not att/rep, not FLOW, not \(J\).

If this extra starts to look like a new millimetre chemical dish, it
has left the identity. Stop.

---

## 1. Frozen identity

| Knob | Freeze | Class |
|---|---|---|
| Mechanism | Scalar length + binary fission at \(\ell_{\mathrm{div}}=3\) µm; both daughters reset to \(\ell_0=1\) µm; clone-at-parent **position**; **no Hertzian** | TAKEN clock (Warren / P0 / Scratch Job 2); scalar body is ENGINEERING |
| \(\lambda_S\) | \(1.0\,\mathrm{h}^{-1}=1/3600\,\mathrm{s}^{-1}\) | TAKEN Warren |
| Bath \(C_s\) | \(0.5\) mM | TAKEN Warren |
| \(K_S\) | \(0.02\) mM | TAKEN Warren |
| Monod | \(C_s/(C_s+K_S)=0.5/0.52\) | DERIVED |
| \(T_{\mathrm{div}}\) | \(\ln 2/(\lambda_S\cdot\mathrm{Monod})=2595.1430440164354\) s (spoken \(2595.1\) s) | DERIVED; **not** 1800 s |
| Elongation | \(\ell(a)=\ell_0\exp(\alpha a)\), \(\alpha=\ln(\ell_{\mathrm{div}}/\ell_0)/T_{\mathrm{div}}\) | TAKEN Valdez / Scratch `analyticLengthUm` |
| Birth-rate law | \(\ln 2/T_{\mathrm{div}}\), **not** \(1/T_{\mathrm{div}}\) | IDENTITY (second reader) |
| Predeclared ODE | \(N_{\mathrm{pred}}(t)=N_0\,2^{t/T_{\mathrm{div}}}\) | mean; death off \(\Rightarrow N_{\mathrm{ever,pred}}=N_{\mathrm{pred}}\) |
| Death | **OFF** | identity |
| Motility / chemotaxis | **OFF** | identity |
| Hill / AHL / att / glucose PDE | **OFF** | identity. No field \(\Rightarrow\) no silent / field arms. \(J\) is not a knob. |
| Clamp | **No `bacteria.max`.** Compute cap \(N_{\mathrm{ever}}=512\). If the cap fires: `SCOPE_NOTE`, not a death, not PASS | ENGINEERING stop |
| Geometry | `LC0_OPEN_BATH`: \(200\times 200\times 10\) µm open box. Founders on an \(8\times 8\) grid, \(z=5\). Overlap allowed (no Hertzian). | ENGINEERING. Not P0 \(100\times 100\times 1.65\). Not Lane A \(1000\times 500\times 10\). |
| \(N_0\) | **64** founders | IDENTITY (second reader). Not 16. Not 32. |
| Seed | **101**, injected `BSimRandom`. No `Math.random()`. | identity |
| Horizon | \(T=T_{\mathrm{div}}\) (one generation) | identity |
| Windows | 300 s. Eight **full** windows \([0,300),\ldots,[2100,2400)\). Stub \([2400,T]\) is logged, **not** a visibility window. | identity |
| \(dt\) | \(1.0\) s | ENGINEERING (P0 / Scratch) |
| Initial ages | i.i.d. **stable age** on \([0,T_{\mathrm{div}})\): \(p(a)=(2\ln 2/T)\,2^{-a/T}\) (normalised; \(\int_0^T p=1\)). Inverse CDF \(a=-T\log_2(1-u/2)\), \(u\sim U[0,1)\), one draw per founder in index order \(0\ldots 63\). Then \(\ell=\ell_0\exp(\alpha a)\). | ENGINEERING desynchronisation, so the ODE is the mean. Not a post-hoc fit. Not all-at-\(\ell_0\) (that is a spike at \(T\)). An earlier draft used \(p(a)=(\ln 2/T)\,2^{-a/T}\) / \(a=-T\log_2(1-u)\), whose CDF only reaches \(1/2\) on \([0,T)\) and is **not** this freeze. |

RNG is used **only** for those 64 initial ages, and only after the
protocol is frozen. Positions do **not** consume the stream.

Do **not** raise \(\lambda_S\), \(C_s\), \(N_0\), or \(T\) after seeing
\(N(t)\).

---

## 2. Arms

| Arm | Role |
|---|---|
| `LC0_GROW` | Sizer on. Same 64 stable-age draws (seed 101). |
| `LC0_OFF` | Same draws, then length and age **frozen**, fission off. \(N(t)=N_0\), births \(=0\). |

No silent source. No field-only. \(J\) is not a knob.

---

## 3. Occupancy (frozen before traces)

Not \(\overline{R}\). Not \(\kappa\).

| Observable | Role |
|---|---|
| Births per full 300 s window | Visibility |
| \(N(t)\) at \(t=0,300,\ldots,2400,T\) | Ledger vs \(N_{\mathrm{pred}}(t)\) |
| \(N_{\mathrm{ever}}(t)\) at the same times | Ledger vs \(N_{\mathrm{ever,pred}}(t)\). Death off \(\Rightarrow\) \(N_{\mathrm{ever}}=N\). |

Desynchronised sizer \(\neq\) exact \(2N_0\) at every sample. **Do not
gate on** \(N_{\mathrm{ever}}(T)=128\).

### Birth-rate law (do not use \(1/T\))

Instantaneous: \(\mathbb{E}[\mathrm{births}]\approx N\ln 2\,\Delta t/T_{\mathrm{div}}\).  
Finite window: \(\mathbb{E}[\mathrm{births}\text{ on }[t,t+\Delta t)]=N_{\mathrm{pred}}(t)\,(2^{\Delta t/T_{\mathrm{div}}}-1)\).

At \(N=64\), \(\Delta t=300\) s, \(T_{\mathrm{div}}=2595.1430440164354\) s:

- Instantaneous at \(N=64\): \(64\ln 2\cdot 300/T_{\mathrm{div}}=5.128\).
- Finite first window: \(64\,(2^{300/T_{\mathrm{div}}}-1)=5.339\).

The investigation table’s 23-at-\(N=200\) used \(N\Delta t/T\) and
overstates by \(1/\ln 2\). Visibility still holds at \(N_0=64\). Clock
gates use \(\ln 2/T\).

### Predeclared clock band

\[
N_{\mathrm{pred}}(t)=64\cdot 2^{t/T_{\mathrm{div}}},\qquad
\delta=0.20.
\]

At every observation \(t\in\{0,300,600,900,1200,1500,1800,2100,2400,T\}\):

\[
\left\lvert\frac{N(t)}{N_{\mathrm{pred}}(t)}-1\right\rvert\le\delta
\quad\text{and}\quad
\left\lvert\frac{N_{\mathrm{ever}}(t)}{N_{\mathrm{pred}}(t)}-1\right\rvert\le\delta.
\]

| \(t\) (s) | \(N_{\mathrm{pred}}\) |
|---|---|
| 0 | 64.000 |
| 300 | 69.339 |
| 600 | 75.124 |
| 900 | 81.391 |
| 1200 | 88.181 |
| 1500 | 95.538 |
| 1800 | 103.508 |
| 2100 | 112.144 |
| 2400 | 121.499 |
| \(T_{\mathrm{div}}\) | 128.000 |

\(\delta=0.20\) is ENGINEERING, frozen before traces. It is wider than
\(1/\sqrt{N}\) sampling noise and **narrower than an all-\(\ell_0\)
spike** (that run stays \(N=64\) until \(T\) and fails the band by
\(t=900\)). It is **not** fitted after \(N(t)\).

A fully synchronized newborn cohort is **not** this identity.

---

## 4. Gates (do not relax after traces)

1. **Honesty.** Logs print `growth=ON` (`LC0_GROW`) or `growth=OFF`
   (`LC0_OFF`), and `death=OFF motility=OFF chemistry=OFF`. Seed 101.
   No `Math.random()` on the path. Protocol JSON
   `frozen_before_traces: true`.
2. **Process-off.** On `LC0_OFF`: total births \(=0\),
   \(N(t)=64\) at every observation, \(N_{\mathrm{end}}=64\),
   \(N_{\mathrm{ever}}(T)=64\).
3. **Visibility.** On `LC0_GROW`, at least one **full** 300 s window
   has births \(\ge 5\). If none do: standing is a **bound** that
   frozen \(N_0=64\) was still too small. Stop. Do not raise
   \(\lambda_S\) or \(N_0\).
4. **Clock.** On `LC0_GROW`, \(N(t)\) and \(N_{\mathrm{ever}}(t)\) lie
   in the §3 band (\(\delta=0.20\)) at every observation time. Not a
   gate on exact \(128\) at \(T\).
5. **No silent kill.** On both arms, \(N_{\mathrm{end}}=N_{\mathrm{ever}}\).
   Compute-cap fire (\(N_{\mathrm{ever}}\) would exceed 512)
   \(\Rightarrow\) `SCOPE_NOTE`, not PASS.

Overall PASS only if 1–5 all PASS. FAIL is allowed. A visibility bound
is a bound.

---

## 5. What a PASS may claim

Births are countable on 300 s windows at a Warren generation clock,
with a growth-off control, without a death clamp.

## 6. What a PASS may not claim

A chemostat. P0 packing. Schink death. Chemotactic memory. Lane A
NRMSE / AXIS_HOLDS. Paper 1 CHARC. That 1800 s is literature. That
packed-pocket Lane C and living-clocks Lane C are the same machine.

---

## 7. Forbidden

- Raise \(\lambda_S\), \(C_s\), \(N_0\), \(T\), or \(\delta\) after \(N(t)\)
- Copy Paper 1 `TOXIC_DEATH_RATE`, starvation `0.05 /s`,
  `BASE_DEATH_RATE`, or `bacteria.max`
- Death, motility, Hill, fields, glucose PDE, FLOW, multi-AC
- Glue P0 rods onto a millimetre dish
- Edit Lane A or Lane B standings / protocols / IEEE draft
- CHARC / IPC / NARMA
- A second extra in this identity
- Commit unless Ceylin asks
