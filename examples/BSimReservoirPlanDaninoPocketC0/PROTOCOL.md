# C0 protocol — well-mixed Object B coupling (NOT_FIG4B)

**Gate:** C0 closed-loop well-mixed cells.  
**Organism:** `DaninoSI_OccupiedDF` (Object B). **NOT_FIG4B.**  
**Not ported:** Object A `Danino2010_Fig4b_bulk_twin` (FAIL / FAIL_NO_IDENTITY).  
**Not this gate:** C1 spatial \(H_e\), P0 Hertzian packing, ACs, NARMA, Fig. 4b.

Claim freeze:
[`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).

D0 remains FAIL. D0b remains FAIL_NO_IDENTITY. D1 remains PASS on
circuit parity only. P0 remains PASS on mechanics only. This protocol
does not reopen Fig. 4b identity, does not shrink \(\mu=0.32\)–\(0.40\),
does not retune \(\alpha,\tau,\gamma_A,\gamma_I,\gamma_H,k_1,d,D\), and
does not hunt ICs.

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after D1/P0 standings, before any C0 integration.

No space, rods, growth, division, Hertzian, motility, extrusion, ACs,
NARMA, `TIME_ADJ=60`, \(\mu_{\mathrm{bus}}\), D1g, `QS_KMLA`, or
HybridDish. Mechanics **OFF**. Channel flow is \(\mu\) on one
compartment, not a PDE.

## Scientific question

Does an \(N\)-cell device, each cell with Object B intracellular state
\((A,I,H_i)\) and a delay tape on \(H_i\), exchanging AHL with **one**
well-mixed extracellular compartment \(H_e\), reduce to the frozen
Object B 4-DDE (including the SI factor \(d/(1-d)\)) when all cells are
identical and the volume fraction is the TAKEN \(d=0.5\)?

This is **not** a second copy of `DaninoSIOccupiedDF` with \(N\) painted
on. Membrane flux is per cell. \(H_e\) is updated from the **sum** of
those fluxes plus \(-\mu H_e\). Identity with D1/D0b is a reduction
check, not a Fig. 4b test.

## Volume-fraction algebra (frozen)

SI bulk (Object B, \(D_1=0\)):

\[
\frac{dH_i}{dt}
= \frac{b I}{1+k I}
  -\frac{\gamma_H A H_i}{1+g A}
  + D(H_e-H_i),
\qquad
\frac{dH_e}{dt}
= -\frac{d}{1-d}\,D(H_e-H_i)-\mu H_e.
\]

C0 explicit cells. Cell \(i=1,\ldots,N\) has intracellular volume
\(v_{\mathrm{cell}}\) and state \((A_i,I_i,H_{i,i})\). One compartment
\(H_e\) has volume \(V_e\). Amounts (SI-scaled concentration \(\times\)
µm³):

\[
M_{\mathrm{in}} = \sum_{i=1}^{N} H_{i,i}\,v_{\mathrm{cell}},
\qquad
M_e = H_e\,V_e.
\]

Membrane moves equal-and-opposite AHL **mass**. The intracellular
membrane term \(D(H_e-H_{i,i})\) is a concentration rate, so the mass
rate into cell \(i\) is \(v_{\mathrm{cell}} D(H_e-H_{i,i})\). The
extracellular mass rate is the negative sum:

\[
\left.\frac{dM_e}{dt}\right|_{\mathrm{membrane}}
= -\sum_{i=1}^{N} v_{\mathrm{cell}} D(H_e-H_{i,i}).
\]

Dividing by \(V_e\) and adding SI flow/decay \(-\mu H_e\):

\[
\frac{dH_e}{dt}
= -\frac{v_{\mathrm{cell}}}{V_e}
  \sum_{i=1}^{N} D(H_e-H_{i,i})
  -\mu H_e.
\]

No extra fitted leak. Enzymatic \(\gamma_H\) acts on intracellular
\(H_i\) only. LuxI synthesis \(b I/(1+k I)\) is an intracellular AHL
source.

Volume fraction:

\[
d = \frac{N v_{\mathrm{cell}}}{N v_{\mathrm{cell}}+V_e}
\qquad\Rightarrow\qquad
\frac{N v_{\mathrm{cell}}}{V_e}=\frac{d}{1-d}.
\]

If every cell is identical (\(H_{i,i}=H_i\)):

\[
\frac{dH_e}{dt}
= -\frac{N v_{\mathrm{cell}}}{V_e}\,D(H_e-H_i)-\mu H_e
= -\frac{d}{1-d}\,D(H_e-H_i)-\mu H_e,
\]

which is Object B. Protein equations are already per-cell and use the
same TAKEN table, including \(C_{A,I}[1-(d/d_0)^4]P\) at this \(d\).

## Frozen volumes (before traces)

P0 pocket \(100\times 100\times 1.65\) µm is **not** packed here. It
supplies only a length unit and a typical rod volume. Growth, Hertzian,
and extrusion stay OFF.

| Symbol | Value | Provenance |
|---|---|---|
| \(v_{\mathrm{cell}}\) | \(1.5\) µm³ | Scratch / Trueba & Koppes \(V_0\); not a fit |
| \(N\) | \(8\) | chosen so \(d=0.5\) is exact and \(N>1\) (sum of fluxes) |
| \(V_e\) | \(N v_{\mathrm{cell}}=12\) µm³ | identity \(d=0.5\) |
| \(d\) | \(8\times 1.5/(8\times 1.5+12)=0.5\) | TAKEN / Fig. 4b caption density; **NOT_FIG4B** |
| \(V_{\mathrm{P0}}\) | \(100\times 100\times 1.65=16500\) µm³ | P0 box; **unused** as C0 \(V_e\) (that would force \(N\sim 5500\) packing, which is C1) |

Check: \(N v_{\mathrm{cell}}/V_e=1=d/(1-d)\) at \(d=0.5\).

## TAKEN table (copied; do not retune)

Same as D1/D0b: \(C_A=1\), \(C_I=4\), \(\delta=10^{-3}\), \(\alpha=2500\),
\(\tau=10\), \(k=1\), \(k_1=0.1\), \(b=0.06\), \(\gamma_A=15\),
\(\gamma_I=24\), \(\gamma_H=0.01\), \(f=0.3\), \(g=0.01\), \(d_0=0.88\),
\(D=2.5\), \(d=0.5\), \(D_1=0\). SI-scaled time. No `TIME_ADJ=60`.

\(P=(\delta+\alpha H_\tau^2)/(1+k_1 H_\tau^2)\), \(H_\tau=H_i(t-\tau)\),
delay only in \(H_i\), per cell.

## Time, scheduler, RNG

- Horizon and sample grid copied from D1/D0b: \(t_{\mathrm{end}}=1000\),
  `sample_dt=0.5`, discard 180.
- RK4, `rk_dt=0.001`, linear \(H_i\) delay tape per cell (same numerical
  method as D1, not a biological retune).
- `BSimStepScheduler` with `dt=rk_dt`. Mechanics events: **none**.
  Transport phases: **no-op** (\(H_e\) is a compartment, not a masked
  PDE; that is C1). Coupled RK4 of all \((A_i,I_i,H_{i,i})\) and \(H_e\)
  lives in `integrateModels` so identical cells match D1’s simultaneous
  4-DDE RK4. `depositFluxes` does not apply a second Euler membrane
  step (that would destroy identity).
- Injected `BSimRandom` seed `0` is **unused**. Circuit is deterministic.
- Deprecated `bsim.dde.BSimDdeSolver` is not used.

## AHL mass ledger

Every RK4 step, using the same stage weights as the state:

- \(M_{\mathrm{in}}=\sum_i H_{i,i} v_{\mathrm{cell}}\)
- \(M_e=H_e V_e\)
- synthesis source \(\sum_i v_{\mathrm{cell}} b I_i/(1+k I_i)\)
- \(\gamma_H\) loss \(\sum_i v_{\mathrm{cell}} \gamma_H A_i H_{i,i}/(1+g A_i)\)
- \(\mu\) loss \(V_e \mu H_e\)
- membrane into cells \(\sum_i v_{\mathrm{cell}} D(H_e-H_{i,i})\);
  membrane into extra is the exact negative (bookkeeping cancel)

Residual (BSimTransportLedger-style):

\[
R = M(0) + \text{sources} - \mu_{\mathrm{loss}} - \gamma_{H,\mathrm{loss}} - M(t).
\]

Membrane is omitted from \(R\) because it cancels. No clipping of \(H\).

**Characteristic mass** \(M_\ast=\max(|M(0)|,|M(t)|,|\text{sources}|,10^{-15})\).  
Gate: \(|R|/M_\ast\le 10^{-3}\) (0.1%) on \(\mu=0\) arms. Identity
\(\mu>0\) arms must still **report** \(R\) and \(\mu_{\mathrm{loss}}\);
they are not allowed a second fitted leak.

Closed membrane-only subtest (documented, not an Object B retune):
\(b=0\), \(\gamma_H=0\), \(\mu=0\), so \(H_i+H_e\) mass is conservative.
IVP: all cells \(A=I=0\), \(H_i=0.05\), \(H_e=0\), history \(0.05\),
\(t_{\mathrm{end}}=50\). Gate the same 0.1% (expect near roundoff).

## Predeclared arms (frozen now)

Period rule: D0b/D1 copy (discard 180, spacing 25, rel-prom 0.20,
min_peaks 4, persist last 30%, amp persist 0.40, rel amp \(\ge 0.25\)).
Readout: mean intracellular \(I\) (identical cells: every cell’s \(I\)).
NARMA-blind.

Trajectory band vs frozen D1/D0b fixtures (identity only):
\(|x_{\mathrm{C0}}-x_{\mathrm{P}}|\le 10^{-4}+0.02|x_{\mathrm{P}}|\) on
the sample grid; tighter \(10^{-6}+10^{-3}|P|\) on \([0,\tau]\).
Mean \((A,I,H_i)\) and shared \(H_e\) versus the 4-DDE fixture.

| id | role | IVP | \(\mu\) | \(\gamma_H,b\) | required |
|---|---|---|---:|---|---|
| `C0_IDENTICAL_D05_MU040` | **identity** | all cells `AHL_KICK_005` | 0.40 | TAKEN | OSC, \(\|T-63.5833\|/63.5833<2\%\), traj band vs D1 `primary_mu_0p40` |
| `C0_IDENTICAL_D05_MU032` | basin flag | all cells `AHL_KICK_005` | 0.32 | TAKEN | OSC, \(\|T-56.3077\|/56.3077<2\%\) |
| `C0_IDENTICAL_D05_MU150` | OFF well | all cells `AHL_KICK_005` | 1.50 | TAKEN | `NO_PERIOD` |
| `C0_WRONG_KICK_BASAL` | wrong kick | all cells `SI_BASAL_PERTURB` | 0.40 | TAKEN | `NO_PERIOD` |
| `C0_MU0_LEDGER` | ledger | all cells `AHL_KICK_005` | 0 | TAKEN | \(\|R\|/M_\ast\le 10^{-3}\); report \(\mu_{\mathrm{loss}}=0\) and \(\gamma_H\) split |
| `C0_MU0_CLOSED_MEMBRANE` | conservative | \(H_i=0.05\), \(H_e=0\), \(A=I=0\) | 0 | \(b=\gamma_H=0\) | \(\|R\|/M_\ast\le 10^{-3}\); membrane cancel |
| `C0_HETERO_ONE_KICK` | smoke, **not identity** | cell 0 `AHL_KICK_005`; cells \(1..N-1\): \(A=I=H_i=0\), history 0; shared \(H_e(0)=0.05\) | 0.40 | TAKEN | shared \(H_e\) must move (\(\max|H_e-H_e(0)|>10^{-6}\)); report only; do not retune |

`AHL_KICK_005`: \(A=I=0\), \(H_i=H_e=0.05\), \(H_i(t<0)=0.05\).  
`SI_BASAL_PERTURB`: \(A=0\), \(I=1\), \(H_i=H_e=0\), history 0.

Do not add \(\mu\) points. Do not claim Fig. 4b. Do not claim Fig. 4e
from the heterogeneous smoke.

## Predeclared PASS / FAIL

C0 **PASS** only if all hold, with **no** TAKEN-table change:

1. Claim freeze names Object B as NOT_FIG4B and Object A as FAIL.
2. Algebra check: \(d\) from volumes equals \(0.5\); identical-cell
   \(H_e\) RHS matches Object B \(d/(1-d)\) factor at a predeclared
   test state.
3. Identity arm `C0_IDENTICAL_D05_MU040` matches D1/D0b fixture within
   the band; Java period OSC with relative error \(<2\%\).
4. Flag arms: \(\mu=0.32\) OSC with \(T\) relative error \(<2\%\);
   \(\mu=1.50\) `NO_PERIOD`; wrong kick `NO_PERIOD`.
5. Ledger arms close to \(\le 0.1\%\).
6. Heterogeneous smoke: \(H_e\) moves. Not scored as Fig. 4b/4e.
7. Every log/plot/standing sentence about this organism says NOT_FIG4B.
8. Identical-cell max \(|H_{i,i}-H_{i,0}|\) stays \(\le 10^{-10}\)
   (no painted-\(N\) drift).

If C0 cannot match without retuning Object B, or the ledger does not
close: **C0=FAIL**. Stop. Do not start C1.

## Commands

```
ant c0-wellmixed
python examples/BSimReservoirPlanDaninoPocketC0/check_c0.py
```

Standing: `examples/PocketDish/C0_WELL_MIXED_STANDING.md`.
Must repeat NOT_FIG4B and that D0/D0b/D1/P0 standings are unchanged.
C1 may start only if this gate PASSes.
