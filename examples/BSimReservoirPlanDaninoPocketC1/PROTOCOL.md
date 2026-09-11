# C1 protocol — packed spatial Object B (NOT_FIG4B)

**Gate:** C1 packed spatial `DaninoSI_OccupiedDF`.  
**Organism:** `DaninoSI_OccupiedDF` (Object B). **NOT_FIG4B.**  
**Not ported:** Object A `Danino2010_Fig4b_bulk_twin` (FAIL / FAIL_NO_IDENTITY).  
**Not this gate:** Fig. 4b period-versus-flow; Object A; ACs; NARMA; W0
traveling waves; P0 Hertzian growth; C0 well-mixed identity.

Claim freeze:
[`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).

Design §7–9, 11, 12, 14/C1 (except the forbidden “flow-dependent period
in the experimental 52–90 min class” bullet — that is Object A / Fig. 4b),
and §15:
[`examples/PocketDish/AC_DANINOPOCKET_PREBUILD_DESIGN.md`](../PocketDish/AC_DANINOPOCKET_PREBUILD_DESIGN.md).

D0 remains FAIL. D0b remains FAIL_NO_IDENTITY. D1 remains PASS on
circuit parity only. P0 remains PASS on mechanics only. C0 remains PASS
on well-mixed Object B only. This protocol does not reopen Fig. 4b
identity, does not shrink \(\mu=0.32\)–\(0.40\), does not retune
\(\alpha,\tau,\gamma_A,\gamma_I,\gamma_H,k_1,d,D\), does not retune
`D1_spatial` after traces, and does not hunt ICs.

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after C0 PASS, before any C1 integration.

A packed job that oscillates is evidence about Object B in space. It is
**not** evidence that Object A was repaired.

---

## Scientific question

Do Object B intracellular delay-DDEs on a **frozen packed patch** of
cells, coupled through a conservative spatial \(H_e\) field on the P0
Danino-class pocket+bus mask, occupy a cycle at frozen \(\mu=0.40\),
synchronize, and close an AHL ledger — all **NOT_FIG4B**?

This is **not** Fig. 4b period vs flow. This is **not** a C0 painted-\(N\)
copy. Membrane flux uses the **local** \(H_e\) voxel. C0 algebra is
local, not global.

---

## Two clocks (mixing them silently is FAIL)

| Clock | Used by | \(dt\) | \(T\) | Notes |
|---|---|---|---|---|
| **SI-scaled** | Object B / C0 / D1 / **C1 identity** | `rk_dt=0.001` | \(t_{\mathrm{end}}=1000\) | \(\tau=10\). No `TIME_ADJ=60`. |
| **Physical seconds** | P0 mechanics | \(1\) s | \(25950\) s | doubling \(2595\) s. **Not** the C1 identity clock. |

C1 identity arms use **SI-scaled time**, same as D1/C0.  
Mechanics (growth, Hertzian, division, extrusion) are **OFF** on identity
arms. Do **not** apply `TIME_ADJ=60`. Do **not** run C1 identity on the
P0 \(25950\) s clock.

Named extra `C1_P0_SNAPSHOT_SMOKE` may turn chemistry on a P0 end state
and **report** occupancy. It is not identity. P0 did not write end
centres; this smoke is **not required** for C1 PASS and is skipped
unless `p0_end_centres.csv` exists. If local \(d\) is not \(0.5\), say
so. Do not retune to occupy it.

Optional extra `C1_BUS_ADV_SMOKE` may enable N0 upwind in the bus;
report outlet loss; not identity; not Fig. 4b. **Not required** for PASS.
Not run in the identity suite.

---

## Density freeze (before traces)

C0 forced \(d=0.5\) with \(N=8\), \(v_{\mathrm{cell}}=1.5\) µm³,
\(V_e=12\) µm³. The P0 box \(100\times100\times1.65\) µm is
\(16500\) µm³. \(N=503\) rods at \(v_{\mathrm{cell}}=1.5\) would be
\(d\approx 0.046\), **not** TAKEN \(d=0.5\). This protocol does **not**
use the P0 box as \(V_e\).

Identity arm `C1_HELD_PACK_D05_MU040` places a **frozen packed patch**
of equivalent centres so the **documented local** volume fraction in
that patch is \(d=0.5\):

\[
d=\frac{N v_{\mathrm{cell}}}{N v_{\mathrm{cell}}+V_{e,\mathrm{local}}}.
\]

Empty pocket voxels remain media. Membrane flux uses the local \(H_e\)
voxel, equal-and-opposite mass. Mechanics OFF.

| Symbol | Frozen value | Provenance |
|---|---|---|
| \(v_{\mathrm{cell}}\) | \(1.5\) µm³ | same as C0; Scratch / Trueba & Koppes \(V_0\); not a fit |
| Patch bbox (world µm) | \(x\in[44,54]\), \(y\in[44,54]\), \(z\in[0,1.65]\) | pocket interior, not the \(+y\) door |
| Voxel | \(dx=dy=2\) µm, \(dz=1.65\) µm, \(V_{\mathrm{vox}}=6.6\) µm³ | P0 / N0 mask |
| Patch voxels | \(i=97..101\), \(j=22..26\), \(k=0\) (**25** voxels) | mask origin \((-150,0,0)\) |
| \(V_{e,\mathrm{local}}\) | \(25\times 6.6=165\) µm³ | fluid volume of those voxels; cells are point centres |
| \(N\) | **110** | \(N v_{\mathrm{cell}}=V_{e,\mathrm{local}}\) |
| \(d\) | \(165/(165+165)=\mathbf{0.5}\) exactly | TAKEN density; **NOT_FIG4B** |
| Occupancy map | first 10 patch voxels (row-major \(i\) then \(j\)): 5 centres; remaining 15: 4 centres | \(10\times5+15\times4=110\) |
| \(z\) | \(0.825\) µm midplane | monolayer slab |
| P0 box as \(V_e\) | **false** | would change identity \(d\) |
| Mechanics | **OFF** | centres, not Hertzian rods |

BSim coordinates: origin at the mask origin. \(x_{\mathrm{BSim}}=x_{\mathrm{world}}+150\),
\(y_{\mathrm{BSim}}=y_{\mathrm{world}}\), \(z_{\mathrm{BSim}}=z_{\mathrm{world}}\).

Protein crowding uses this frozen patch \(d=0.5\) for every packed cell
(\(C_{A,I}[1-(d/d_0)^4]P\)). Do not substitute whole-pocket \(d\).

Do not retune \(\alpha,\tau,\gamma^\ast,k_1\), membrane \(D\), or \(d\)
after seeing traces.

---

## Spatial \(H_e\) (before traces)

SI spatial \(H_e\) has \(D_1\nabla^2 H_e\); bulk C0/D1 dropped it.

C1 uses N0 `BSimTransportField` (conservative masked FTCS).  
**Do not** use legacy in-place `BSimChemicalField.diffuse()`.  
**Do not** use PocketNeck `ChipField`.

| Item | Frozen value |
|---|---|
| `D1_spatial` | **800** (SI units, µm² / SI-time) |
| Label | `SI_TAKEN_SPATIAL` (mid SI Fig. 6 class: examples 0, 200, 800, 4000) |
| Fitted to waves? | **no** |
| \(\mu\) identity | **0.40** (Object B extracellular decay). Not channel µm/min. Not \(\mu_{\mathrm{bus}}\). Not Fig. 4b identity. |
| Bus advection | **OFF** on identity (\(v=0\)) |
| Pocket flow | **0** |
| Mask | P0 pocket + bus. SHA-256 `2b53bb0b872477ab28afab8b20da48e1c63daa834614caf55aab63043081f99b` |
| Chemical BC | no-flux on the six domain faces; fluid–fluid faces only inside the mask |
| Kernel | `BSimTransportField.diffuse` + `decay`. Advection skipped while \(v=0\) (zero field; not a different operator). |

FTCS / Courant (identity half-step \(\Delta t=5\times10^{-4}\)):

\[
2\,D_1\,\Delta t\bigl(1/dx^2+1/dy^2+1/dz^2\bigr)
= 2\cdot800\cdot5\times10^{-4}\cdot(0.25+0.25+1/1.65^2)
\approx 0.694 < 1
\]

N0 subcycle count is 1 at this load. If a later arm exceeds the gate,
subcycle; do not rescale \(D_1\).

Local membrane (not global \(d/(1-d)\)):

\[
\left.\frac{dM_{\mathrm{cell}}}{dt}\right|_{\mathrm{mem}}
= v_{\mathrm{cell}} D (H_{e,\mathrm{voxel}}-H_i),
\qquad
\left.\frac{dM_{\mathrm{voxel}}}{dt}\right|_{\mathrm{mem}}
= -v_{\mathrm{cell}} D (H_{e,\mathrm{voxel}}-H_i).
\]

AHL_KICK_005 in space (frozen IVP, **not** a local-only kick):

- packed cells: \(A=I=0\), \(H_i=0.05\), \(H_i(t<0)=0.05\)
- **all fluid voxels** (patch, empty pocket, bus): \(H_e=0.05\)

A kick confined to 25 patch voxels would dilute into
\(V_{\mathrm{fluid}}=10500\times6.6=69300\) µm³ on the \(D_1=800\)
timescale \(\ll\tau\), which is a different IVP. Whole-field
\(H_e(0)=0.05\) is the spatial analog of the frozen bulk kick.
Empty voxels remain media. Spatial non-uniformity is required to
develop after \(t=0\) from localized membrane coupling.

`SI_BASAL_PERTURB`: packed cells \(A=0\), \(I=1\), \(H_i=0\), history 0;
all fluid \(H_e=0\).

---

## Scheduler (N0 / design §9)

`BSimStepScheduler` with `dt=rk_dt`. Canonical order:

1. apply boundary (no-op on identity: no-flux, \(v=0\))
2. half transport (`diffuse` + `decay` for \(dt/2\))
3. sample local \(H_e\) at each packed centre
4. integrate delayed intracellular Object B states (RK4, \(H_e\) frozen
   at the sample; delay tape on \(H_i\) only)
5. compute membrane mass (RK4-weighted)
6. deposit equal-and-opposite mass onto the local voxel
   (`transferQuantity`; **not** `addQuantity` source-ledger)
7. complete transport
8. mechanics skipped
9. observe at the new timestamp

RNG: injected `BSimRandom` seed `0`, unused. No `Math.random()`.  
Deprecated `bsim.dde.BSimDdeSolver` is not used.

Horizon and sample grid copied from D1/C0: \(t_{\mathrm{end}}=1000\),
`sample_dt=0.5`, discard 180. `rk_dt=0.001` identity;
`rk_dt=0.002` on the predeclared LORES arm.

---

## TAKEN table (copied; do not retune)

Same as D1/C0/D0b: \(C_A=1\), \(C_I=4\), \(\delta=10^{-3}\),
\(\alpha=2500\), \(\tau=10\), \(k=1\), \(k_1=0.1\), \(b=0.06\),
\(\gamma_A=15\), \(\gamma_I=24\), \(\gamma_H=0.01\), \(f=0.3\),
\(g=0.01\), \(d_0=0.88\), \(D=2.5\), \(d=0.5\).  
`D1_spatial=800`. SI-scaled time. No `TIME_ADJ=60`.

\(P=(\delta+\alpha H_\tau^2)/(1+k_1 H_\tau^2)\), \(H_\tau=H_i(t-\tau)\),
delay only in \(H_i\), per cell. Primary readout: intracellular \(I\)
(LuxI). No HybridDish \(R/L\). No GFP maturation.

---

## AHL mass ledger

Residual class as C0/N0:

\[
R = M(0)+\mathrm{sources}-\mu_{\mathrm{loss}}-\gamma_{H,\mathrm{loss}}
-\mathrm{outlet}-\mathrm{boundary}-M(t).
\]

- \(M=M_{\mathrm{in}}+M_e\), \(M_{\mathrm{in}}=\sum H_i v_{\mathrm{cell}}\),
  \(M_e=\sum_{\mathrm{fluid}} q_{\mathrm{voxel}}\)
- sources: intracellular synthesis \(\sum v_{\mathrm{cell}} b I/(1+k I)\)
- \(\mu_{\mathrm{loss}}\): N0 field exponential decay (Object B \(\mu\))
- \(\gamma_H\) loss: intracellular enzymatic only
- outlet / boundary: N0 transport ledger (identity: expect \(\approx 0\))
- membrane into cells vs into voxels: equal-and-opposite; omitted from
  \(R\); report cancel \(\approx 0\)

No clipping of \(H\). Characteristic mass
\(M_\ast=\max(|M(0)|,|M(t)|,|\mathrm{sources}|,10^{-15})\).  
Gate: \(|R|/M_\ast\le 10^{-3}\) (0.1%) on `C1_HELD_PACK_MU0_LEDGER`
(\(\mu=0\), TAKEN \(\gamma_H\), documented split). Identity \(\mu>0\)
arms still **report** \(R\).

---

## Period checker

Reuse D0b/C0 rules on **colony-mean** \(I\) of the packed patch:

discard 180, spacing 25, rel-prom 0.20, min_peaks 4, persist last 30%,
amp persist 0.40, rel amp \(\ge 0.25\). NARMA-blind.

Space may shift \(T\). Predeclared occupancy band vs frozen Object B
\(T=63.5833\):

\[
|T-63.5833|/63.5833 \le 0.20.
\]

Report the relative error. Do **not** retune to hit 63.58. Do **not**
convert to 55/90 min. Do **not** score increasing-\(T\) with \(\mu\) as
a pass.

---

## Synchronization gate (frozen now)

On `C1_HELD_PACK_D05_MU040` after discard, using the same peak rule on
**each** packed cell's \(I\):

1. Colony-mean \(I\) is `OSC` with the 20% band above.
2. **SYNC_FRAC_OSC:** fraction of the \(N=110\) cells with flag `OSC`
   \(\ge 0.90\).
3. **SYNC_PEAK_WINDOW:** of those OSC cells, fraction whose last peak
   time lies within \(\tau/2=5\) of the median last-peak time of OSC
   cells \(\ge 0.90\).

Metric name: `SYNC_PEAK_WINDOW`. Not circular-mean; not `std(I)/mean(I)`.

## Spatial \(H_e\) non-uniformity (frozen now)

C0 has a single \(H_e\). C1 must use space.

At recorded samples, track
\(\max H_e-\min H_e\) over all fluid voxels (`he_fluid_range`) and over
patch voxels (`he_patch_range`).

**PASS** if \(\max_t(\texttt{he_fluid_range}) \ge 10^{-6}\).  
Report `he_patch_range` too. If both ranges stay at solver noise, C1
has not used space: FAIL. Do not retune \(D_1\).

---

## Predeclared arms (frozen now)

Identity:

| id | \(d\) | \(\mu\) | IVP | `D1_spatial` | \(v\) | mechanics | required |
|---|---|---:|---|---:|---|---|---|
| `C1_HELD_PACK_D05_MU040` | 0.5 local | 0.40 | AHL_KICK_005 all packed cells; He=0.05 all fluid | 800 | 0 | OFF | OSC, 20% \(T\) band, sync gates, He non-uniform |

Controls:

| id | notes | required |
|---|---|---|
| `C1_HELD_PACK_WRONG_KICK` | `SI_BASAL_PERTURB`, \(\mu=0.40\) | `NO_PERIOD` |
| `C1_HELD_PACK_MU150` | AHL_KICK_005, \(\mu=1.50\) | `NO_PERIOD` / OFF |
| `C1_HELD_PACK_MU0_LEDGER` | AHL_KICK_005, \(\mu=0\), TAKEN \(\gamma_H\) | \(\|R\|/M_\ast\le 10^{-3}\); report \(\mu_{\mathrm{loss}}=0\) and \(\gamma_H\) split |

Robustness (second resolution predeclared: **2× coarser time**):

| id | change | required |
|---|---|---|
| `C1_HELD_PACK_MU040_LORES` | `rk_dt=0.002`, same grid, same \(N,d,D_1,\mu\), AHL_KICK_005 | **OSC call unchanged** (still `OSC`; \(T\) may move inside the 20% band) |

Smoke only (not identity, not required for PASS):

| id | notes |
|---|---|
| `C1_P0_SNAPSHOT_SMOKE` | P0 end centres if present; chemistry ON; SI clock; report occupancy and local \(d\); do not retune |
| `C1_BUS_ADV_SMOKE` | optional N0 upwind in the bus; report outlet; not Fig. 4b |

Do **not** add a \(\mu\) sweep to “fix” Fig. 4b. If \(T\) at
\(\mu=0.32\) vs \(0.40\) vs \(0.60\) is ever printed, label NOT_FIG4B
and do not score increasing-\(T\) as a pass.

---

## Predeclared PASS / FAIL

C1 **PASS** only if all hold, with **no** Object B / `D1_spatial` /
\(d\) change after traces:

1. Claim freeze names Object B as NOT_FIG4B and Object A as FAIL.
2. Identity `C1_HELD_PACK_D05_MU040`: colony-mean \(I\) is OSC in the
   20% band around \(63.5833\); sync gates pass; He is spatially
   non-uniform as frozen.
3. Spatial He + membrane ledger: `C1_HELD_PACK_MU0_LEDGER` \(\le 0.1\%\).
4. `C1_HELD_PACK_MU040_LORES` OSC call unchanged.
5. Wrong-kick and \(\mu=1.50\) remain `NO_PERIOD`.
6. Every log, plot, and standing sentence says NOT_FIG4B.
7. Mechanics OFF; `TIME_ADJ` unused; P0 box not used as \(V_e\);
   mask SHA matches P0; N0 kernel; no ChipField; no NARMA.

If occupancy, sync, or ledger holds only after retuning Object B / D1 /
\(\mu\) / \(d\): **C1=FAIL**. Stop. Do not start A0/W0. Do not claim
Fig. 4b.

A0 may start only if C1 PASSes. A0 is an ideal AC source on this device,
still NOT_FIG4B, still not Fig. 4b. W0 is later.

---

## Commands

```
ant c1-packed
python examples/BSimReservoirPlanDaninoPocketC1/check_c1.py
```

Standing: `examples/PocketDish/C1_PACKED_SPATIAL_STANDING.md`.  
Must repeat NOT_FIG4B on every result sentence, and that D0/D0b/D1/P0/C0
standings are unchanged.
