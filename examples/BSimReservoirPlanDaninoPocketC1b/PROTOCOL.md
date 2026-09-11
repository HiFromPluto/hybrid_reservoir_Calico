# C1b protocol — island C0-equivalent occupancy (NOT_FIG4B)

**Gate:** C1b island reduction of Object B.  
**Organism:** `DaninoSI_OccupiedDF` (Object B). **NOT_FIG4B.**  
**Not ported:** Object A `Danino2010_Fig4b_bulk_twin`.  
**Not this gate:** Fig. 4b; A0; W0; retuning Object B to occupy the open chip.

C1 packed-spatial standing remains **FAIL**. Do not rewrite it as PASS.  
Diagnosis: [`examples/PocketDish/C1_FAIL_DIAGNOSIS.md`](../PocketDish/C1_FAIL_DIAGNOSIS.md).  
Claim freeze: [`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after C1 FAIL diagnosis, before any C1b integration.

This protocol freezes a **reduction ladder**. Occupancy identity is a
C0-equivalent extracellular operator, not “open the D1=800 chip and
retune.”

Do not shrink \(\mu=0.32\)–\(0.40\). Do not `TIME_ADJ=60`. Do not hunt
ICs. Do not retune \(\alpha,\tau,\gamma^\ast,k_1,d,D\).

---

## Named objects

### `C1_OPEN_DILUTE` (already run; not identity)

Full P0 mask, \(\mu\) on every fluid voxel, \(H_e(0)=0.05\) everywhere,
`D1_spatial=800`, 110-cell 10 µm patch.

**Status:** FAIL occupancy (`C1_PACKED_SPATIAL_STANDING.md`).  
**Predeclared here:** stays `NO_PERIOD` / report-only. Do **not**
re-integrate. Do **not** retune Object B to occupy it.

### `C1b_ISLAND` (new occupancy identity)

\(H_e\), \(\mu\), and (if \(D_1>0\)) diffusion live **only** on the
island extracellular volume. Other P0 mask voxels are chemical solids /
not in \(V_e\).

Identity arm uses **C0 volumes** and **`D1_spatial=0`** (well-mixed
island = C0):

| Symbol | Frozen value |
|---|---|
| \(N\) | **8** (C0; not the C1 110-cell patch) |
| \(v_{\mathrm{cell}}\) | 1.5 µm³ |
| \(V_e\) | **12** µm³ (one chemical voxel; box \(2\times2\times3\) µm) |
| \(d\) | **0.5** exactly |
| `D1_spatial` | **0** |
| \(\mu\) identity | **0.40** on that \(V_e\) only |
| IVP | `AHL_KICK_005` on those 8 cells **and that \(V_e\) only** |
| Mechanics | OFF |
| Clock | SI-scaled, \(\tau=10\), \(t_{\mathrm{end}}=1000\), `rk_dt=0.001`. No `TIME_ADJ`. |

Integrator: C1 membrane / field path restricted to that one fluid
voxel — not a re-run of C0’s scalar \(H_e\) ODE as a rubber stamp, and
not the open P0 chip. \(H_e\) is N0 `BSimTransportField` quantity in
the single fluid voxel (`setConc` / `decay` hit that voxel only).
Diffuse is a no-op (`D1_spatial=0`). \(\mu\) is field decay on that
\(V_e\) only, not domain-scale decay of 69300 µm³. Intracellular RK4
uses frozen local \(H_e\) as C1; membrane is signed `transferQuantity`
(not source-ledgered). Algebra check still uses C0’s identical-cell
RHS vs Object B.

**Must recover** C0/D1 OSC at \(\mu=0.40\),
\(|T-63.5833|/63.5833 < 2\%\), ledger closed, identical-cell \(H_i\)
spread \(\le 10^{-10}\). Trajectory band vs D1 `primary_mu_0p40` as C0.

If `C1b_ISLAND` fails: membrane / C0-reduction coupling is wrong.
**FAIL. Stop. Do not open `D1_spatial`.**

### `C1b_ISLAND_D1_800` (optional extra; not identity)

Same **chemical island** as the C1 packed patch (25 voxels, \(N=110\),
\(V_{e,\mathrm{local}}=165\) µm³, \(d=0.5\)), constructed as a **5×5×1**
grid (voxel 2×2×1.65 µm). Other P0 mask voxels are **not constructed**
(chemical solids / not in \(V_e\)). `D1_spatial=800` inside the patch
only (no-flux at island faces). \(\mu=0.40\) on island fluid only.
`AHL_KICK_005` on packed cells and island \(H_e\) only.

Report \(T\) (or `NO_PERIOD`). Do **not** retune. Not scored as C1b
occupancy identity. Not Fig. 4b.

---

## Why C0 volumes for identity (picked now)

C1’s 25 voxels at `D1=0` are **25 disconnected compartments** with 4 or
5 cells each (local \(d\neq 0.5\)). That is not C0. C0 is one \(V_e=12\)
shared by \(N=8\) identical cells at \(d=0.5\).

Identity therefore uses **C0’s 8-cell volumes**, not the 110-cell C1
patch. The 110-cell island appears only on the optional `D1=800` extra:
25 fluid voxels on a **5×5×1** island grid (voxel 2×2×1.65 µm,
\(V_{e,\mathrm{local}}=165\)), no-flux faces, other P0 voxels **not
constructed**.

---

## AHL ledger

Same residual class as C0/N0. Gate \(\le 0.1\%\) on
`C1b_ISLAND_MU0_LEDGER` (\(\mu=0\), TAKEN \(\gamma_H\)).

---

## Period checker

D0b/C0 rules on mean \(I\). Identity band **2%** (not C1’s 20%).
NARMA-blind. No TIME_ADJ. No 55/90 min conversion.

---

## Predeclared arms

| id | role | required |
|---|---|---|
| `C1b_ISLAND_D05_MU040` | **identity** | OSC, \(\|T-63.5833\|/63.5833<2\%\), traj band vs D1, \(H_i\) spread \(\le 10^{-10}\) |
| `C1b_ISLAND_WRONG_KICK` | control | `NO_PERIOD` |
| `C1b_ISLAND_MU150` | OFF well | `NO_PERIOD` |
| `C1b_ISLAND_MU0_LEDGER` | ledger | \(\|R\|/M_\ast\le 10^{-3}\) |
| `C1b_ISLAND_D1_800` | extra, not identity | report \(T\) or `NO_PERIOD`; do not retune |
| `C1_OPEN_DILUTE` | already run | report-only `NO_PERIOD` from C1 standing; do not re-run |

---

## PASS / FAIL

C1b **PASS** only if identity + controls + ledger hold with **no**
Object B retune.

C1 FAIL is unchanged. Even if C1b PASSes, the **open chip is still OFF**.
A0 is **not** automatic. W0 is later. Standing must say so.

If identity fails: C1b=FAIL. Stop. Do not start A0/W0. Do not open
`D1_spatial` as a repair.

---

## Commands

```
ant c1b-island
python examples/BSimReservoirPlanDaninoPocketC1b/check_c1b.py
```

Standing: `examples/PocketDish/C1B_ISLAND_STANDING.md`.  
Must say NOT_FIG4B and that `C1_PACKED_SPATIAL_STANDING.md` remains FAIL.
