# B0 protocol — motility-as-memory on a millimetre dish

**Gate:** B0 motility-as-memory feasibility. No task score.  
**Status label:** `B0_MOTILITY_MEMORY`  
**Object:** `LANE_B_CHEMOTACTIC_SPATIAL`  
**Device:** `LANE_B_MM_DISH_FLOW0` (1000×500×10 µm, FLOW=0)  
**Not this gate:** NARMA, CHARC, IPC, delay-line of \(u\), Danino / Object B /
Fig. 4b, C1c, A1, HybridDish Hill \(R\to L\), paper-1 CHARC rewrite,
calibrated Lentini AC, growth/division/death as readout.

**frozen_before_traces:** true  
**Frozen:** 2026-08-28, before any B0 integration.

This protocol answers, on an identity frozen **before traces**:

> On a frozen millimetre dish with N0 transport, does motile *E. coli*
> hold spatial density information about a painted attractant/repellent
> scene for a window that is short compared with mixing at the
> **readout wavelength**, and does that information vanish in
> motility-off and remain distinct from the chemical field itself?

Label: `LANE_B_CHEMOTACTIC_SPATIAL`, `B0_MOTILITY_MEMORY`.  
Not Fig. 4b. Not Lane A. Not C1c.

Related files:

- Lane comparison: [`../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md`](../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md)
- Costing review (**process only**; §4 table is **not** this freeze):
  [`../PocketDish/LANES_REVIEW_AND_LANE_B_DESIGN.md`](../PocketDish/LANES_REVIEW_AND_LANE_B_DESIGN.md)
- N0 kernel: [`../PocketDish/N0_NUMERICAL_KERNEL_STANDING.md`](../PocketDish/N0_NUMERICAL_KERNEL_STANDING.md)
- Claim freeze (do not port Object A/B):
  [`../PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md)
- HybridDish field-vs-cells prior:
  [`../HybridDish/WHAT_IS_ESTABLISHED.md`](../HybridDish/WHAT_IS_ESTABLISHED.md)

---

## 1. What this job is not

- Not Danino, not `DaninoSI_OccupiedDF`, not Fig. 4b.
- Not C1c occupancy voxels, not A1, not P0 packed rods.
- Not HybridDish Hill \(R\to L\) and not a HybridDish NRMSE rewrite.
- Not a rewrite of paper-1 CHARC / voxel-feature tables.
- Not NARMA, CHARC, IPC, or “bacteria beat a delay-line of \(u\).”
- Not a calibrated AC (A0-class current only; `HYPOTHETICAL_DESIGN_ENVELOPE`
  if a later stateful envelope is added).

If B0 **FAIL**s versus field-only, Lane B does **not** start capacity work.

---

## 2. What was discarded from the costing review

[`LANES_REVIEW_AND_LANE_B_DESIGN.md`](../PocketDish/LANES_REVIEW_AND_LANE_B_DESIGN.md)
§4 is **not identity**. Do not copy these into a gate:

| Discarded | Why |
|---|---|
| \(\tau=L_{\mathrm{dish}}^2/(4D)\) with \(L=500\) or \(1000\) µm | Pattern forgetting on a voxel map uses the **resolved wavelength** \(\lambda\), not the dish width. Eigenmode decay is \(\sim\lambda^2/(\pi^2 D)\), not \(L^2/(4D)\). |
| 6 s “frozen” / 300 s “usable” / 521 s mix | 300 s was an audit window for **counting births**, the opposite observable. “Both want 300 s” is not evidence. |
| Stokes \(D=0.66\) µm²/s and **182×** contrast | Unit error (\(\eta\) as Pa·s not mPa·s). Water 37 °C, \(r=1\) µm, is \(\sim 330\) µm²/s. Control = **motility-off**, not Stokes 0.66. |
| Gate mixing to \(L^2/(4\mu_0)\) | Circular. \(\mu_0\) is unstimulated random motility. Directed chemotactic drift can be faster. |
| Default paper-1 `flow.x=8` as the no-flow costing | Fields advect at 8 µm/s; cells used `FLOW_COUPLING_FACTOR=0.05`. HybridDish claim is FLOW=0. Freeze flow **on or off as named arms**. |

Keep from the review: pre-register feasibility, **no task score**, N0 ledger,
chemotaxis citations, motility-off and field-only as the kill comparison.

---

## 3. Named objects (frozen now)

| Object | Role |
|---|---|
| `LANE_B_CHEMOTACTIC_SPATIAL` | The dish class: ACs paint att/rep; motile *E. coli* read by Barkai–Leibler chemotaxis; colony **density** is the reservoir candidate. |
| `LANE_B_MM_DISH_FLOW0` | Identity device for B0: 1000×500×10 µm, **FLOW=0**, N0 kernel, solid no-flux walls. |
| `LANE_B_MM_DISH_FLOW8` | **Named extra, not identity.** Field advection 8 µm/s. Not required for B0 PASS. Do not mix into FLOW=0 numbers. |
| `B0_MOTILITY_MEMORY` | This gate. Feasibility only. |
| `A0_LANE_B_PAINT` | Immobilized A0-class current \(J=J_{\max}u(t)\) into **attractant** or **repellent**. No internal AC ODE. |

Paper-1 `reservoir_new` remains an external copy. This job **ports the
Barkai–Leibler att−rep living layer and A0-class paint**, on N0 clocks.
It does not import CHARC MATLAB, QS, Monod, or death.

---

## 4. Device (frozen now)

| Symbol | Frozen value | Provenance |
|---|---|---|
| Dish | \(1000\times 500\times 10\) µm | paper-1 / HybridDish millimetre class, **named** |
| FLOW | **0** | identity. Matches Lane A / HybridDish so field vs cells is comparable |
| Walls | solid, no-flux (chemicals and cells) | N0 BSim solid = no-flux |
| Kernel | N0 `BSimTransportField` | [`N0_NUMERICAL_KERNEL_STANDING.md`](../PocketDish/N0_NUMERICAL_KERNEL_STANDING.md) |
| PDE grid | \(100\times 50\times 1\) | \(dx=dy=dz=10\) µm. Finer than readout. |
| Readout grid | \(20\times 10\) | \(\lambda_x=\lambda_y=\mathbf{50}\) µm. Coarser than the PDE. |
| \(\lambda\) | **50 µm** | resolved wavelength of the density/field maps used in gates |
| AHL field | **absent** | B0: no QS confound |
| Glucose / growth / division / death | **OFF** | B0 is not a generation-scale memory test |
| Scheduler | `BSimStepScheduler` | not legacy `BSim.export()` inclusive loop |
| RNG | one injected `BSimRandom`, seed **0** | N0 contract; no `Math.random()` on this path |
| \(T\) | 305 K | BSim / Berg motility medium |
| \(\eta\) | \(2.7\times 10^{-3}\) Pa·s | BSim TAKEN (Berg *E. coli* motility conditions) |
| \(dt\) | 0.02 s | ENGINEERING; \(\ll\) run/tumble means |
| Sample | 0.10 s | maps and ledgers |

`FLOW=8` is listed in §10 as a named extra. It is **not** run for B0 PASS.

---

## 5. Chemicals (frozen now)

Payload species: **attractant** and **repellent** (Tar-class ligands).  
Not AHL. Not IPTG. Not an unnamed “signal.”

| Quantity | Frozen value | Provenance |
|---|---|---|
| \(D_{\mathrm{att}}=D_{\mathrm{rep}}\) | **800 µm²/s** | ENGINEERING, small-molecule-in-water **order**. Not claimed as Middlebrooks TAKEN (PDF \(D_{\mathrm{Asp}}\) not re-read here). Paper-1 used 50/100 — **not** copied. |
| \(k_{\mathrm{att}}=k_{\mathrm{rep}}\) | **0** | ENGINEERING. B0 forgetting is **transport**, not hydrolysis. |
| Units | molecules / voxel; conc = molecules/µm³ | N0 ledger units |
| Conversion note | 1 molecule/µm³ \(\approx\) 1.7 nM | BSim comment; **not** used as a gate |

Why 800, not paper-1 50: aqueous amino-acid diffusivities are \(O(10^{-5})\) cm²/s
\(= O(10^3)\) µm²/s. A slow ENGINEERING field would make B0.4 a field-only
test (HybridDish already showed that trap). The number is frozen before
traces and is **not** fitted after a map score.

Diagnostic only (not a gate): eigenmode
\(\tau_{\mathrm{field}}(\lambda)\sim\lambda^2/(\pi^2 D)\approx 0.32\) s at
\(\lambda=50\) µm, \(D=800\). Actual \(\tau\) is **measured**.

---

## 6. Artificial cells (frozen now)

A0-class, immobilized, no internal ODE.

\[
J(t)=J_{\max}\,u(t)
\quad\text{(molecules/s into the containing PDE voxel)}
\]

| Item | Frozen value |
|---|---|
| Count | 4 |
| Mobility | immobilized |
| \(J_{\max}\) | **\(2\times 10^4\)** molecules/s (ENGINEERING occupancy of Barkai–Leibler \(L\sim O(1)\) after spreading; **not** fitted) |
| Command | rectangular pulse \(u=1\) on \([0,t_{\mathrm{off}})\), else 0 |
| \(t_{\mathrm{off}}\) | **10 s** | ENGINEERING; directed swim at \(\sim 20\) µm/s can move \(\sim 200\) µm, not \(L_{\mathrm{dish}}^2/4\mu_0\) |
| Deposit | scheduler `depositFluxes`, `addQuantity` on containing voxel |

Layout (world µm, \(z=5\)):

| AC | \((x,y,z)\) | payload |
|---|---|---|
| 1 | (200, 150, 5) | attractant |
| 2 | (200, 350, 5) | attractant |
| 3 | (800, 150, 5) | repellent |
| 4 | (800, 350, 5) | repellent |

Painted scene: **attractant left, repellent right.**  
Silent control: same AC objects, \(J=0\).

---

## 7. Living layer (frozen now)

Port of paper-1 `ChemotaxisDynamics` / run–tumble, **stripped**:

**On:** Barkai–Leibler 3-state ODE; ligand \(L=c_{\mathrm{att}}-c_{\mathrm{rep}}\);
run–tumble; CheY-P modulates \(p_{\mathrm{end\,run}}\).

**Off:** AHL, QS, gene-expression ODE, Monod, glucose uptake, toxic/starvation
death, crowding, 35 µm lookahead `dirFactor`, fast/slow adapters,
`FLOW_COUPLING_FACTOR`, growth/division.

ODE (paper-1, per-second units; **not re-fit**):

\[
\begin{aligned}
a_\infty(L,m) &= \bigl(1+\exp[\alpha(L-m)]\bigr)^{-1}, \\
\dot a &= (a_\infty-a)/\tau_a, \\
\dot m &= k_R(1-a)-k_B a, \\
\dot y &= k_A a - \gamma_y y.
\end{aligned}
\]

| Symbol | Frozen value | Provenance |
|---|---|---|
| \(\alpha,k_R,k_B,k_A,\gamma_y,\tau_a\) | 2, 0.05, 0.05, 1, 0.5, 0.5 | paper-1 `sim_config` / `ChemotaxisDynamics` |
| IC | \(m=0\), \(a=0.5\), \(y=1\) plus seeded \(U(-0.05,0.05)\) on \(m,a\) | L=0 Barkai steady state; no `Math.random()` |
| \(p_{\mathrm{end\,run}}\) (unstimulated) | \(1/0.86\) s⁻¹ | Berg, BSim TAKEN |
| \(p_{\mathrm{end\,tumble}}\) | \(1/0.14\) s⁻¹ | Berg, BSim TAKEN |
| Mapping | \(p_{\mathrm{end\,run}}=p_{\mathrm{else}}\cdot\mathrm{clip}(y/y_0,0.2,3)\) with \(y_0=k_A\cdot 0.5/\gamma_y=1\) | ENGINEERING. Unstimulated recovers Berg. No lookahead. |
| Flagellar force | 1 pN | BSim TAKEN \(\Rightarrow v\approx 20\) µm/s at frozen \(\eta,r\) |
| Radius | 1 µm | BSim default |
| Brownian on motile arm | **OFF** | MSD is swimming \(\mu\), not Stokes |
| Point particles | yes | not Hertzian rods; not C1c voxels |

AHL **does not enter** this ODE (same as paper-1 `ChemotaxisDynamics` comment).

---

## 8. Readout (frozen now)

Density and field maps on the **20×10** grid (\(\lambda=50\) µm), **not** the
PDE grid.

Left–right **contrast** (amplitude-aware; Pearson vs layout is forbidden as
the B0.4 distance because a faint leftover plume still correlates):

\[
\kappa(t)=\frac{S_L(t)-S_R(t)}{S_L(t)+S_R(t)+\varepsilon},\quad\varepsilon=10^{-12}.
\]

- Cells: \(S_L\) = count with \(x<500\) µm, \(S_R\) with \(x\ge 500\) µm.
- Field: \(S_L\) = attractant mass in \(x<500\) plus repellent mass in
  \(x\ge 500\); \(S_R\) = the swapped masses. Scene-aligned field mass.

Density-map autocorrelation \(C(\tau)\): Pearson correlation of the flattened
20×10 density at time \(t_0\) vs \(t_0+\tau\).

\(\tau_{\mathrm{mix}}(\lambda)\): first lag where blob-arm \(C(\tau)\le 1/e\),
measured, **not** \(L_{\mathrm{dish}}^2/(4\mu_0)\).

Memory window (B0.3):

\[
W=0.3\,\tau_{\mathrm{mix}}(\lambda)
\quad\text{declared band: } W/\tau_{\mathrm{mix}}\in[0.1,1].
\]

Score B0.4 at \(t^\star=t_{\mathrm{off}}+W\). Record \(W/\tau_{\mathrm{mix}}\)
on **every** paint run. If \(\tau_{\mathrm{mix}}\) is unmeasurable (never
crosses \(1/e\)), FAIL B0.2 and stop.

---

## 9. Arms (frozen now)

| Arm | Cells | Motility | ACs | Chemistry |
|---|---|---|---|---|
| `B0_MSD_MOTILE` | 400, interior spawn | run–tumble, Brownian off | none | off |
| `B0_MIX_BLOB` | 1200, Gaussian blob \(\sigma=40\) µm at centre | motile, no gradient | none | off |
| `B0_PAINT_MOTILE` | 1200, uniform | motile | pulse | N0 att/rep |
| `B0_PAINT_OFF` | 1200, uniform | **motility-off** (frozen positions) | pulse | N0 att/rep |
| `B0_PAINT_FIELD` | **none** | — | pulse | N0 att/rep |
| `B0_PAINT_SILENT` | 1200, uniform | motile | \(J=0\) | N0 fields exist, no source |
| `B0_THERMAL` | 400, interior | **Brownian only** (1 µm, frozen \(\eta,T\)) | none | off |

`B0_THERMAL` is **reported separately** from motility-off. Do not call it
Stokes \(D=0.66\). Do not use it as the B0.4 null.

Motility-off = same particles, run/tumble off, flagellar force off,
Brownian off, positions frozen.

---

## 10. Flow arms (named, not identity)

| Arm | Status |
|---|---|
| `LANE_B_MM_DISH_FLOW0` | **B0 identity** |
| `LANE_B_MM_DISH_FLOW8` | Named extra. Field \(v_x=8\) µm/s (\(\sim 125\) s to clear 1000 µm). Cells in paper-1 used coupling 0.05; B0 does **not** silently copy that. Not run for B0 PASS. |

---

## 11. Gates (pre-registered; do not use 521 s / 0.58)

### B0.1 Unstimulated MSD → \(\mu\)

From `B0_MSD_MOTILE`, 2-D xy MSD. Fit \(\langle r_{xy}^2\rangle=4\mu t\) on
\(t\in[2,10]\) s (after ballistic / one run). Report 2-D \(\mu\) and, if
computed, a 3-D diagnostic \(\langle r^2\rangle/6t\).

Compare to literature **after unit conversion** (not a fit):

- Middlebrooks et al. 2021: \(\mu_0=1.2\times 10^{-6}\) cm²/s \(= \mathbf{120}\) µm²/s
  (already read in the Scratch dossier).
- Zhao & Ford 2022 / bioRxiv 2021.12.30.474376: \(\mu_0=1.3\pm 0.21\times 10^{-10}\) m²/s
  \(= \mathbf{130\pm 21}\) µm²/s (gradient-assay random motility, not this MSD).

**FAIL only if** \(\mu<1\) or \(\mu>5000\) µm²/s (dead or unphysical), or MSD
is nonfinite. Do **not** require Middlebrooks identity. Quasi-2-D confinement
(\(H=10\) µm) is declared.

### B0.2 Measured \(\tau_{\mathrm{mix}}(\lambda)\)

From `B0_MIX_BLOB`. \(\tau_{\mathrm{mix}}\) = time to \(C=1/e\) on the 20×10
map. **Must not** require equality to \(L_{\mathrm{dish}}^2/(4\mu_0)\). That
quantity may be printed as `NOT_A_GATE`.

**FAIL if** \(C(\tau)\) never drops to \(1/e\) before \(T=16\) s, or is
nonfinite.

### B0.3 Window band

\(W=0.3\,\tau_{\mathrm{mix}}\). **PASS if** \(W/\tau_{\mathrm{mix}}\in[0.1,1]\)
(true by construction unless \(\tau_{\mathrm{mix}}\) is unusable). Record the
ratio. Do not copy 300 s.

### B0.4 Ablation (kill gate)

At \(t^\star=t_{\mathrm{off}}+W\):

| Inequality | Frozen threshold |
|---|---|
| \(\lvert\kappa_{\mathrm{motile}}\rvert\) | \(\ge 0.12\) |
| \(\lvert\kappa_{\mathrm{motile}}\rvert-\lvert\kappa_{\mathrm{field}}\rvert\) | \(\ge 0.08\) |
| \(\lvert\kappa_{\mathrm{motile}}\rvert-\lvert\kappa_{\mathrm{off}}\rvert\) | \(\ge 0.08\) |
| \(\lvert\kappa_{\mathrm{silent}}\rvert\) | \(< 0.08\) |

If motile does **not** separate from field-only, **STOP**. Lane B is dead for
capacity work. That is the useful result. **Do not add NARMA to rescue it.**

Pearson correlation of z-scored maps with the AC layout is **not** the gate
(scale-invariant leftover plumes). Contrast is.

### B0.5 N0 chemical mass ledger

Closed residual on each painted field:

\[
R = M_0 + M_{\mathrm{src}} - M_{\mathrm{decay}} - M_{\mathrm{outlet}} - M_{\mathrm{boundary}} - M_{\mathrm{now}}.
\]

**PASS if** \(\lvert R\rvert / \max(M_{\mathrm{src}}, 1) < 10^{-3}\) on att and
rep for every chemistry arm, and remaining mass is nonnegative/finite.
Decay and outlet are 0 on this closed dish.

---

## 12. Thermal Stokes (report only)

In-dish Stokes–Einstein at frozen Berg \(\eta,T,r=1\) µm:

\[
D_{\mathrm{SE}}=\frac{k_B T}{6\pi\eta r}.
\]

**Arithmetic (formula frozen; digit corrected 2026-08-28 after the thermal arm):**
\(k_B T/(6\pi\eta r)=8.27\times 10^{-14}\) m²/s \(=\mathbf{0.083}\) µm²/s, not 83
(the \(10^{12}\) m²→µm² factor was applied as \(10^{15}\) in the first draft).
Water 37 °C yardstick (\(r=1\) µm, \(\eta\approx 0.69\) mPa·s): \(\sim\mathbf{0.33}\) µm²/s, not 330.

`B0_THERMAL` MSD is compared to this \(D_{\mathrm{SE}}\) as a sanity print.
It is **not** the B0.4 control. Motility-off is.

**Forbidden as a gate:** Stokes \(D=0.66\) µm²/s and 182× “Brownian is frozen.”
Those belong to the discarded costing table, not this identity.

---

## 13. Citations (verified numbers only)

| Source | Use in B0 |
|---|---|
| Barkai & Leibler 1997 *Nature* | BUILD — ODE structure already in paper-1 `ChemotaxisDynamics` |
| Berg *Chemotaxis in Escherichia coli* (BSim TAKEN run/tumble / Stokes swim) | BUILD — motility kinematics |
| Middlebrooks et al. 2021 10.1002/bit.27930 | BUILD **candidate** for \(\mu_0=120\) µm²/s comparison; **not** a fit target; \(D_{\mathrm{Asp}}\) not TAKEN |
| Zhao & Ford 2022 10.1002/bit.28161 | BUILD **candidate** for independent \(\mu_0=130\pm 21\) µm²/s |
| Dilanji et al. 2012 10.1021/ja211593q | **CITE-ONLY** for Lane B. Immobilized cells; AHL \(D\) spread. Not motility identity |

See [`../PocketDish/CITATION_DOSSIER.md`](../PocketDish/CITATION_DOSSIER.md) §7.

---

## 14. Durations (frozen now)

| Arm | \(T\) (s) |
|---|---|
| MSD / thermal | 12 |
| Mix blob | 16 |
| Paint arms | \(t_{\mathrm{off}}+12=22\) (must cover \(t^\star\); if \(W>12\), extend to \(t_{\mathrm{off}}+W+2\) in code without retuning biology) |

---

## 15. Forbidden moves

1. Freeze LANES_REVIEW §4 as this PROTOCOL.
2. Use Stokes \(D=0.66\) or 182×.
3. Gate \(\tau_{\mathrm{mix}}\) to \(L_{\mathrm{dish}}^2/(4\mu_0)\).
4. Copy 300 s / 521 s / 0.58.
5. Port Object A/B, C1c, A1, Fig. 4b, HybridDish \(R,L\).
6. Turn growth/death on as a B0 readout.
7. Fit \(\chi,\mu,J,D\) after seeing a task NRMSE (there is no task NRMSE).
8. Add NARMA/CHARC after a B0.4 FAIL.
9. Run FLOW=8 as identity.
10. Commit unless asked.

---

## 16. Claim boundary

**May claim after traces, if gates pass:** motile density on this frozen dish
holds a painted att/rep scene for a measured window \(W\) with
\(W/\tau_{\mathrm{mix}}(\lambda)\) in the declared band; that contrast
vanishes in motility-off and is not already in the coarsened field at
\(t^\star\); N0 ledgers close.

**May not claim:** Lane B computes a temporal task; bacteria beat a delay
line of \(u\); CHARC/IPC; Danino; Lentini twin; paper-1 score improvement;
that 300 s is the operating window; that Brownian is 182× slower than
swimming.
