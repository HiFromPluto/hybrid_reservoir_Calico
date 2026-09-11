# C1c protocol — filled-pocket spatial occupancy (NOT_FIG4B)

**Gate:** C1c filled-pocket occupancy of Object B.  
**Organism:** `DaninoSI_OccupiedDF` (Object B). **NOT_FIG4B.**  
**Not ported:** Object A `Danino2010_Fig4b_bulk_twin`.  
**Not this gate:** Fig. 4b; A0; W0; retuning Object B to occupy `C1_OPEN_DILUTE`.

C1 packed-spatial standing remains **FAIL**. Do not rewrite it as PASS.  
C1b island standing remains **PASS** on the C0-equivalent voxel only. A0 is
not automatic from C1b. The open chip is still OFF.

Diagnosis: [`examples/PocketDish/C1_FAIL_DIAGNOSIS.md`](../PocketDish/C1_FAIL_DIAGNOSIS.md).  
Claim freeze: [`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after C1b island PASS, before any C1c integration.

This protocol answers, on an identity frozen **before traces**:

> Can Object B occupy a cycle when cells fill the pocket so extracellular
> AHL volume is the colony’s media, not the empty chip+bus?

That is the SI filled-array class (SI Modeling: \(N=200\) interval, periodic
ends, \(D_1\) on \(H_e\)), realized here as a **filled 100×100×1.65 µm
pocket**, not a 10 µm patch in a 400 µm bus.

Do not shrink \(\mu=0.32\)–\(0.40\). Do not `TIME_ADJ=60`. Do not hunt
ICs. Do not retune \(\alpha,\tau,\gamma^\ast,k_1,d,D\). Do not shrink
`D1_spatial` after traces. Do not start A0 from a FAIL.

---

## Named objects (frozen now)

Keep:

| Object | Role | Status |
|---|---|---|
| `C1b_ISLAND` | C0-equivalent voxel occupancy | PASS (unchanged). NOT_FIG4B |
| `C1_OPEN_DILUTE` | empty chip+bus, \(\mu\) on all fluid | FAIL / `NO_PERIOD`. Report-only. Do **not** re-run. Do **not** retune Object B to occupy it. NOT_FIG4B |
| `C1b_ISLAND_D1_800` | extra; occupied closed patch | extra, not the open chip, not this identity. NOT_FIG4B |

New identity:

### `C1c_FILLED_POCKET`

Cells throughout the **100×100×1.65 µm** pocket (not a 10×10 patch).
Bus voxels are **chemical solids on identity** (not in \(V_e\); \(\mu\)
and \(D_1\) never act on them). Pick recorded: solids, not “bus exists
as fluid with pocket-only \(\mu\)”.

| Symbol | Frozen value |
|---|---|
| Pocket | \(100\times100\times1.65\) µm |
| Voxel | \(dx=dy=2\) µm, \(dz=1.65\) µm, \(V_{\mathrm{vox}}=6.6\) µm³ (P0) |
| Pocket grid | \(50\times50\times1\) (all fluid) |
| \(N_{\mathrm{vox}}\) | **2500** |
| \(V_{e,\mathrm{pocket}}\) | \(2500\times 6.6=\mathbf{16500}\) µm³ |
| \(v_{\mathrm{cell}}\) | 1.5 µm³ (point centres; do not subtract from \(V_e\)) |
| \(N\) | **11000** (\(N v_{\mathrm{cell}}=V_{e,\mathrm{pocket}}\)) |
| \(d\) | \(N v_{\mathrm{cell}}/(N v_{\mathrm{cell}}+V_{e,\mathrm{pocket}})=\mathbf{0.5}\) exactly |
| Occupancy map | first 1000 voxels (row-major \(j\) outer, \(i\) inner): 5 centres; remaining 1500: 4 centres. \(1000\times5+1500\times4=11000\) |
| Crowding \(d\) in \([1-(d/d_0)^4]\) | **pocket-global \(d=0.5\)** for every cell. Not per-voxel \(d\) from 4 vs 5 centres |
| Membrane | local voxel \(H_e\); equal-and-opposite `transferQuantity`. Not global \(d/(1-d)\) on \(H_e\) |
| `D1_spatial` | **800**, `SI_TAKEN_SPATIAL` (SI examples 0, 200, 800, 4000). **Not 0.** Not fitted after traces |
| Open-edge BC | **no-flux** on \(x=0\), \(x=100\), \(y=0\), **and** on the \(+y\) mouth. Not Dirichlet. Not island-at-mouth. Bus is solid, so the mouth is a no-flux wall |
| \(\mu\) identity | **0.40** on pocket fluid only |
| IVP | `AHL_KICK_005` on the 11000 packed cells **and pocket \(H_e\) only**. Not a 69300 µm³ whole-domain bath |
| Mechanics | OFF |
| Clock | SI-scaled, \(\tau=10\), \(t_{\mathrm{end}}=1000\), `rk_dt=0.001`. No `TIME_ADJ` |
| Kernel | N0 `BSimTransportField` (signed `transferQuantity`). Not ChipField. Not legacy `diffuse()` |

Integrator reduction (frozen, not a retune): all centres in a voxel share
that voxel’s \(H_e\) and identical ICs, so they remain identical. The
job integrates **one delay-DDE per voxel** with occupancy weight
\(n_{\mathrm{here}}\in\{4,5\}\). Algebraically identical to \(N=11000\)
individual cells. Colony-mean \(I\) is occupancy-weighted.

This is **not** `C1_OPEN_DILUTE`. This is **not** `C1b_ISLAND`. This is
**not** Fig. 4b.

---

## Why this identity (picked now)

C1 failed because \(\mu=0.40\) acted on every fluid voxel of the open
chip+bus (\(V_{\mathrm{fluid}}=69300\) µm³) while 110 cells fed a 10 µm
patch. The bath washed out before \(\tau\). That operator is named
`C1_OPEN_DILUTE` and stays FAIL.

C1b recovered occupancy on C0 volumes (\(N=8\), \(V_e=12\), `D1=0`).
That does not occupy a filled pocket.

SI spatial \(H_e\) is a **filled interval** (\(N=200\), periodic,
\(D_1\) on \(H_e\)), not a dilute patch in a bus. C1c places cells
throughout the pocket so \(V_e\) is the colony’s media at documented
\(d=0.5\).

`D1_spatial=800` is the SI mid-class value already used on C1 / C1b
extra. Identity uses space. `D1=0` would be a well-mixed filled pocket
(C0 with painted \(N\)); that is **not** this identity.

---

## Density (C0 definition)

\[
d=\frac{N v_{\mathrm{cell}}}{N v_{\mathrm{cell}}+V_{e,\mathrm{pocket}}}.
\]

\(V_{e,\mathrm{pocket}}\) is the interstitial pocket fluid: all 2500
pocket voxels, cells as point centres (same convention as C0/C1). The
P0 box is **not** used as \(V_e\) of an empty chip; it **is** the
pocket fluid that the colony fills.

Protein crowding uses this **same** \(d=0.5\) in
\([1-(d/d_0)^4]\). Local 4-vs-5 occupancy is only a membrane-source
checkerboard so spatial \(H_e\) is used; it does not replace documented
\(d\).

---

## Spatial \(H_e\) and open edge

Field lives on the pocket-only \(50\times50\times1\) grid. Domain faces
`setSolid(true)`: no-flux on all six faces, including the geometric
\(+y\) open edge. Bus is **not constructed**.

FTCS / Courant (identity half-step \(\Delta t=5\times10^{-4}\)):

\[
2\,D_1\,\Delta t\bigl(1/dx^2+1/dy^2+1/dz^2\bigr)
\approx 0.694 < 1.
\]

N0 subcycle count is 1. If a later arm exceeds the gate, subcycle; do
not rescale \(D_1\).

Bus advection OFF. Pocket flow 0. \(v=0\).

---

## AHL kick (not the C1 bath)

- packed cells: \(A=I=0\), \(H_i=0.05\), \(H_i(t<0)=0.05\)
- pocket fluid only: \(H_e=0.05\)
- bus: not present (solids)

A whole-domain bath on 69300 µm³ is a **different extra**, not identity.
`SI_BASAL_PERTURB`: packed cells \(A=0\), \(I=1\), \(H_i=0\); pocket
\(H_e=0\).

---

## AHL ledger

Same residual class as C0/N0/C1:

\[
R = M(0)+\mathrm{sources}-\mu_{\mathrm{loss}}-\gamma_{H,\mathrm{loss}}
-\mathrm{outlet}-\mathrm{boundary}-M(t).
\]

Gate \(\le 0.1\%\) on `C1c_MU0_LEDGER` (\(\mu=0\), TAKEN \(\gamma_H\)).
Identity \(\mu>0\) arms still report \(R\).

---

## Period checker

D0b/C0/C1 rules on **colony-mean** \(I\) (occupancy-weighted).

Space may shift \(T\). Occupancy band vs frozen Object B \(T=63.5833\):

\[
|T-63.5833|/63.5833 \le 0.20.
\]

Report the relative error. Do **not** retune to hit 63.58. Do **not**
convert to 55/90 min. Do **not** score increasing-\(T\) with \(\mu\).

NARMA-blind. No TIME_ADJ.

---

## Space used (frozen now)

C0 has a single \(H_e\). C1c must use space.

At recorded samples, track \(\max H_e-\min H_e\) over pocket fluid
(`he_pocket_range`).

**PASS** if \(\max_t(\texttt{he_pocket_range})\ge 10^{-6}\).  
If the range stays at solver noise, C1c has not used space: FAIL. Do
not retune \(D_1\).

Colony-mean \(I\) OSC is the occupancy readout. Per-cell C1 sync of
11000 traces is **not** required for C1c PASS.

---

## Predeclared arms

| id | role | required |
|---|---|---|
| `C1c_FILLED_POCKET_D05_MU040` | **identity** | OSC, \(\|T-63.5833\|/63.5833\le 20\%\), He range \(\ge 10^{-6}\) |
| `C1c_WRONG_KICK` | control | `NO_PERIOD` |
| `C1c_MU150` | OFF well | `NO_PERIOD` |
| `C1c_MU0_LEDGER` | ledger | \(\|R\|/M_\ast\le 10^{-3}\) |
| `C1_OPEN_DILUTE` | already run | report-only `NO_PERIOD` from C1 standing; do not re-run |
| `C1c_BUS_CONNECTED` | optional extra, not identity | pocket+bus one fluid, \(\mu\) on all, `AHL_KICK_005` as a domain bath. Predeclared **`NO_PERIOD`**. Not a retune target. Not required for PASS. Skip unless identity PASSes **and** the job is invoked with this arm |

Do **not** add a \(\mu\) sweep. Do **not** occupy `C1_OPEN_DILUTE` by
shrinking \(D_1\) or \(\mu\) after traces.

---

## PASS / FAIL

C1c **PASS** only if identity + controls + ledger + space-used hold
with **no** Object B / `D1_spatial` / \(d\) / \(\mu\) change after
traces.

C1 FAIL is unchanged. C1b island PASS is unchanged.

If C1c PASSes: A0 may start only as an ideal AC on
`C1c_FILLED_POCKET`, labeled **NOT_FIG4B**, **not** on
`C1_OPEN_DILUTE`. W0 is later. This is **not** a Fig. 4b pass.

If identity stays basal/OFF: **C1c=FAIL**. Stop. Do not start A0/W0.
Do not start A0. Do not occupy the empty-bus chip by shrinking \(D_1\)
or \(\mu\).

---

## Commands

```
ant c1c-filled
python examples/BSimReservoirPlanDaninoPocketC1c/check_c1c.py
```

Standing: `examples/PocketDish/C1C_FILLED_POCKET_STANDING.md`.  
Must say NOT_FIG4B on every result sentence, that
`C1_PACKED_SPATIAL_STANDING.md` remains FAIL, and that
`C1B_ISLAND_STANDING.md` remains island PASS.
