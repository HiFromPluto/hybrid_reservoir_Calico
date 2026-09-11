# PocketNeck-T2 — occupancy / leak vs door width (job protocol)

Frozen 2026-08-20. Job-level only. Architecture freeze:
[`examples/PocketDish/PROTOCOL.md`](../PocketDish/PROTOCOL.md).
Neck geometry (this freeze):
[`examples/PocketDish/GEOMETRY_NECK.md`](../PocketDish/GEOMETRY_NECK.md).
T0 job protocol (flush clone, stencil, plateau language):
[`examples/BSimReservoirPlanPocketDishT0/PROTOCOL.md`](../BSimReservoirPlanPocketDishT0/PROTOCOL.md).
T0 standing (cite only; do not rewrite):
[`examples/PocketDish/T0_STANDING.md`](../PocketDish/T0_STANDING.md).
Builder authorization:
[`examples/PocketDish/POCKETNECK_T2_FROZEN_BUILDER_PROMPT.md`](../PocketDish/POCKETNECK_T2_FROZEN_BUILDER_PROMPT.md).

This file does **not** authorize NARMA, Mackey–Glass, waveform AUC,
ridge, a λ grid, week-2 living Java, PocketFill volumetric LuxI, a D1g
seed bath, packed-plug \(D_{\mathrm{eff}}\), Danino 4-ODE, vesicles,
glucose, Grober pucks, Stokes, two-way \(h(P)\), a circle pocket,
hexagonal traps, or mother-machine trenches.

Do not edit any `GATE_EVIDENCE.md`. Do not retune HybridDish
\(K,n,\tau_R,\tau_L\) or HybridDish / T0 claim \(J_{\max}\). Do not
edit T0 or T1 results. Do not pick \(W\) after seeing \(R\).

## Scientific question

On HybridDish clocks (\(D=159\) µm²/s, \(k=0.0033\) s⁻¹,
\(K=1.6\), \(n=2\), \(\tau_R=15\), \(\tau_L=1500\)) and T0
\(J_{\max}=6.36\times10^5\) molecules/s at pocket centre, does a
**narrow neck** leave PocketHill occupied when the flush 100 µm door
does not?

Hypothesis (predeclared, **before** occupancy):

- `W100_FLUSH` is **DEAD** and must match T0 `OPEN_BUS_3` class
  (`mean_R ≈ 0.0318`). If it is ALIVE, the clone is wrong — **stop
  and fix geometry**, do not celebrate a neck.
- Narrower \(W\) **may** ALIVE because opening conductance (and leak)
  scales with door width. Order-of-magnitude only: T0 \(C\approx 0.29\) µM;
  if leak \(\propto W\), \(W=50\) µm is near the occupancy gate. This
  is a **hypothesis**, not a fit target.
- If **all** neck arms are DEAD at frozen \(W\): the door is still not
  enough on these clocks. **Stop.** Do not retune \(K\), \(k\),
  \(J_{\max}\), or \(L_n\). Do not add \(W=5\) µm after looking.

This is **not** “more Danino.” Danino’s bulk trap is a full-edge feed.
The neck is a Kim 2016 crevice / crypt-class ENGINEERING door
(flow-shielded mouth). See [`GEOMETRY_NECK.md`](../PocketDish/GEOMETRY_NECK.md).

## Path

Preferred: **reduced 2-D field + Hill/\(L\) surrogate**, clone T0
`transport_model.py` (SciPy `solve_ivp` BDF + sparse stencil). No
living Java this job.

Occupancy means = **pocket voxels only** (neck and bus are diagnostics,
not in `mean_R`).

## Frozen dish (PocketDish-A clocks, T0 source, T2 door)

Copied clocks and source **verbatim** from T0. Door from
`GEOMETRY_NECK.md`. **Do not steal neck volume from the 100×100
pocket** (closed-box \(C_{\mathrm{ss}}\) must stay T0).

| Quantity | Value | Class |
|---|---|---|
| Pocket | \(100\times 100\times 10\) µm | PocketDish-A |
| Pocket flow | 0 | frozen |
| Bus | \(400\times 80\) µm, \(v=3.0\) µm/s, inlet \(C=0\), conservative outlet | PocketDish-A |
| Grid | \(dx=dy=5\) µm, one \(z\) box (depth 10 µm) | T0 numerics |
| \(D\) | 159 µm²/s | HybridDish |
| \(k\) | 0.0033 s⁻¹ | HybridDish ENGINEERING |
| \(K,n,\tau_R\) | 1.6 µM, 2, 15 s | HybridDish; do not retune |
| \(\tau_L\) | 1500 s | HybridDish; do not retune |
| Conversion | 1 µM = 602.2 molecules/µm³ | DERIVED |
| \(J_{\max}\) | \(6.36\times 10^5\) molecules/s at pocket centre | T0 verbatim |
| Command | \(J=J_{\max} u\) | one-way; no vesicles |
| AC | (50, 50) µm in pocket coords | T0 verbatim |

**Flush attachment (`W100_FLUSH`):** T0 `OPEN_BUS_3` replay. Pocket
\(x\in[150,250]\), \(y\in[0,100]\); bus \(y\in[100,180]\), \(x\in[0,400]\);
**no** neck voxels.

**Neck attachment (necked arms only):** pocket \(x\in[150,250]\),
\(y\in[0,100]\); neck centred on pocket `+y`, width \(W\),
\(y\in[100,120]\); bus \(y\in[120,200]\), \(x\in[0,400]\). Bus width
stays 80 µm. Neck `FLOW=0`. Bus flow does not enter the neck as a
through-pocket velocity. Diffusion couples pocket ↔ neck ↔ bus.

Neck length \(L_n=20\) µm (`ENGINEERING`). \(W=10\) µm is **2 voxels**
at \(dx=5\). That is acceptable. Do **not** change \(W\) to 15 µm for
prettier cells.

Closed-pocket well-mixed algebra (unchanged from T0; pocket \(V\) is
still \(10^5\) µm³):

\[
C_{\mathrm{ss}}(u=1)=2K=3.2~\mu\mathrm{M},\qquad
J_{\max}=k\,C_{\mathrm{ss}}\,V,\qquad
V=10^5~\mu\mathrm{m}^3.
\]

At \(u=0.5\), leak-free \(C_{\mathrm{ss}}=K=1.6\) µM. Open-pocket
occupancy is a measurement of **interface leak**. Do not raise
\(J_{\max}\) if necks are `DEAD`. Do not lower \(K\) if they are
`SATURATED`.

Reduced model: cell-centred finite volume on the 5 µm grid; centred
diffusive fluxes; first-order upwind advection on the **bus only**;
SciPy `solve_ivp` BDF with `rtol=1e-7`, `atol=1e-10`. Mass identity
(molecules):

`residual = initial + injected − remaining − decay_loss − boundary_loss`

Each drive starts from \(C=R=L=0\). Unresolved
\(|\mathrm{residual}| > 1\%\) of the dominant flux (production or leak),
or \(C < -10^{-9}\) µM, is a numerical failure. This surrogate has
**no cells**. \(N(t)\) is `NA`.

## Predeclared arms (do not add after occupancy)

| Arm | Door | Role |
|---|---|---|
| `W100_FLUSH` | full 100 µm `+y` face, **no** neck, T0 bus attachment | T0 `OPEN_BUS_3` replay. Must **DEAD**, `mean_R` within ~10% of **0.0318** |
| `W50_L20` | neck \(W=50\) µm, \(L_n=20\) µm | modest constriction |
| `W20_L20` | neck \(W=20\) µm, \(L_n=20\) µm | literature constriction class |
| `W10_L20` | neck \(W=10\) µm, \(L_n=20\) µm | aggressive (not 3 µm hex) |

Do not add `W=80`, `W=5`, \(L_n=40\), Robin \(h\), or a circle after
the first occupancy number.

Voxel check (must print): \(W/dx\) and \(L_n/dx\) are integers for
every necked arm (10, 4, 2 widths; \(L_n/dx=4\)).

## Predeclared drive (no task labels)

Command range \([0,1]\). No NARMA file. No sequence labels.

**`STEP_ON` only:** \(u=0.5\) held for **\(T=9000\) s** (T0 plateau
rule). Sample \(\le 60\) s. Plateau: terminal pocket `mean_R` within
5% of the mean of `mean_R` over the last 10% of time (8100–9000 s).
Occupancy flag is evaluated on the terminal pocket field. \(r(R,u)\)
is **NA** (\(u\) does not vary).

No `SINGLE_PULSE` (T0 already showed chemical ≠ reporter).
No `THROUGH_8` (already DEAD, not this door).

## Occupancy flags (not NRMSE)

| Flag | Rule |
|---|---|
| `ALIVE` | pocket `mean_R ≥ 0.05` and not `SATURATED` |
| `DEAD` | pocket `mean_R < 0.05` |
| `SATURATED` | pocket `mean_R > 0.95` (unexpected on an open neck; if it happens, units/source are wrong — fix conversion, not \(K\)) |

Also print at plateau, per arm:

- pocket-mean AHL, \(R\), \(L\)
- neck-mean AHL (NA on `W100_FLUSH`)
- bus-mean AHL
- corner/interior AHL ratio in the **pocket** (T0: closed-end ~0.98)
- mass budget: production, decay, pocket→neck (or pocket→bus) interface
  flux, bus outlet, residual. T0 lesson: `W100_FLUSH` interface should
  dominate decay (~2.3×10⁹ vs ~1.1×10⁹ class)
- \(L^2/D\), Damköhler, bus Pe
- non-negative \(C\)

`TRANSPORT_MODEL_STATUS`: mass residual ≤1% of the dominant flux
(production or leak). Clone T0 validation discipline.

If `W100_FLUSH` is ALIVE or `mean_R` is not T0-class: **stop**. Fix
the flush clone before interpreting necks.

If any neck arm is ALIVE: that is the **leak story**. Write it. **Do
not** start NARMA, LuxI, or a living pack. Do not add more widths.

If all three neck arms are DEAD: **stop**. Neck is not enough on
HybridDish clocks. Do not retune. A later named prompt may ask hours
+ seed (PocketFill), not a quieter \(K\).

## Optional grid note (not occupancy)

After the four-arm table is written, one extra `W10_L20` at
\(dx=2.5\) µm may be run, labelled `CONVERGENCE` not occupancy.
Must keep \(W=10\), \(L_n=20\). Do not use a finer-grid ALIVE to
replace a coarse DEAD in the official table. If skipped, say so.

## Gates that are not this job

No NRMSE, AUC, NARMA target, ridge, or λ. The checker **must refuse**
NARMA paths and must not import `BSimReservoirPlanNarma10b`.
Do not write Overall PASS/FAIL on a task. Lines:

- `OCCUPANCY_W100_FLUSH=…`
- `OCCUPANCY_W50_L20=…`
- `OCCUPANCY_W20_L20=…`
- `OCCUPANCY_W10_L20=…`
- `TRANSPORT_MODEL_STATUS=…`
- `LIVING_DECISION=STOP_AFTER_NECK_SCREEN`

## Outputs

- `results/OCCUPANCY_SCREEN.md`
- `results/occupancy_screen.csv`
- maps of pocket (+ neck + bus) \(C\) and \(R\) for every arm
- SHA-256 of configs and checker

## Commands

```
python examples/BSimReservoirPlanPocketNeckT2/check_pocketneck_t2.py --theory
python examples/BSimReservoirPlanPocketNeckT2/test_mass_budget.py
python examples/BSimReservoirPlanPocketNeckT2/run_occupancy_screen.py
python examples/BSimReservoirPlanPocketNeckT2/check_pocketneck_t2.py
```

`--smoke` on the runner integrates 50 s only and must not be scored as
occupancy evidence.

## Stop list

- Do not add arms, \(h\), \(J_{\max}\), \(K\), or extra \(W\) after occupancy
- Do not retune HybridDish clocks to occupy the pocket
- Do not start NARMA, LuxI, or living Java from this job
- Do not call this an exact Danino/Prindle blueprint
- Do not rewrite T0 / T1 / WashoutReset / Narma10b results
- Do not start week-2 living Java, circle twin, or QS seed unless a
  **new** frozen prompt says so
