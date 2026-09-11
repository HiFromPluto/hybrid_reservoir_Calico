# PocketDish protocol (architecture freeze)

Frozen 2026-08-19. This file does **not** authorize NRMSE hunting,
ridge λ grids, or a HybridDish rewrite.

Copied clocks: HybridDish \(D,k,K,n,\tau_R,\tau_L\).
**New** geometry, \(J_{\max}\), population story, and primary readout.

## Integrity (do not violate)

Do **not** edit:

- any `GATE_EVIDENCE.md`
- `examples/BSimReservoirPlanNarma10b/` results
- HybridDish claim-dish Java / configs as a “quiet fix”
- `reservoir_new` Results / CHARC tables
- E4.1 / E4.2 scout numbers (read-only)

Do **not**:

- retune HybridDish \(K,n,\tau_R,\tau_L,J_{\max}\) so PocketDish occupies
- resurrect Stage 3 Danino 4-ODE as the 300 s receiver
- call this an “exact Danino/Prindle blueprint”
- call vesicles TX–TL unless the ODE and payload match Lentini 2014/2017
- implement Grober Stokes/rotlets in the first job
- start pucks in Java before `examples/PocketDish/` has a reduced-model
  memo (expected **KILL**: PocketHill has no QS→motility chain)
- mix glucose / Monod into this dish
- use Stage 8 ~4.2 h AHL half-life as HybridDish or PocketHill

Paper 1 = `reservoir_new`. Paper 2 = HybridDish claim dish (already
scored). PocketDish = **paper-3 device / outlook object**, not a
retcon of paper 2.

## Architecture names

| Name | Geometry | Receiver | Status |
|---|---|---|---|
| **PocketHill** | PocketDish-A (flush) | HybridDish Hill + \(L\) | T0: **OPEN_BUS_3 occupancy DEAD**. Backbone cannot score. See [`T0_STANDING.md`](T0_STANDING.md) |
| PocketOsc | PocketDish-A, hours, filled | Danino-class oscillator | later; **new** occupancy protocol (not a T0 \(k\) retune) |
| PocketRect | image 1, 1000×500 + vertical bus | PocketHill | not this protocol |
| **PocketNeck** | PocketDish-A + ENGINEERING neck \(W,L_n\) | PocketHill (T0 clocks) | T3 Hill **ALIVE** on `W20`. T4 population **EMPTY** ([`T4_STANDING.md`](T4_STANDING.md)). **No NARMA**. No T5 colony movie |

## Dish (PocketDish-A)

See [`GEOMETRY_FREEZE.md`](GEOMETRY_FREEZE.md).

- Pocket `100 × 100 × 10` µm. Bus along +y open edge, width 80 µm,
  \(v_{\mathrm{bus}}=3.0\) µm/s (Danino 180 µm/min class).
- Pocket `FLOW_SPEED = 0`. Through-pocket `8` µm/s is **control only**.
- `dt = 0.05` s. Chemical grid: **dx = dy = 5 µm** in the pocket
  (20 × 20 × 1). HybridDish 20 µm voxels are too coarse here.
  Numerics, not a computing axis.
- AHL \(D=159\) µm²/s, \(k=0.0033\) s⁻¹.
- Hill: \(K=1.6\) µM, \(n=2\), \(\tau_R=15\) s.
- \(L\): \(\tau_L=1500\) s.
- Conversion: 1 µM = 602.2 molecules/µm³.

## Source freeze (before occupancy, not after)

HybridDish \(J_{\max}=1.28\times10^8\) is **not** copied. SweepS4
showed a box change at frozen rate occupancy-kills.

Closed-pocket well-mixed target at command 1:

\[
C_{\mathrm{ss}}(u=1)=2K=3.2~\mu\mathrm{M},\qquad
J_{\max}=k\,C_{\mathrm{ss}}\,V
\]

\(V=100\times100\times10=10^5\) µm³ →
\(J_{\max}=6.36\times10^5\) molecules/s.

Class: `ENGINEERING_DESIGN`. Not a LuxI rate. Not HybridDish \(J_{\max}\).

Warmup / step command uses \(u=0.5\) → closed-pocket \(C_{\mathrm{ss}}=K\)
if leak is zero. **Open-pocket occupancy is then a measurement of
interface leak.** Do not raise \(J_{\max}\) if it is DEAD. Do not
lower \(K\) if it is SATURATED.

AC: one-way `addQuantity` at pocket centre. No vesicles. Inventory
12/12 `MISSING` → not `CALIBRATED_A1`.

## Population

Image 2 is a **filled** pocket. Do not tell the story with 1800 cells
in a 1000×500 box.

- Seed **small** (Danino “once seeded”). Suggested ENGINEERING:
  \(N_0=50\) in the pocket.
- Growth: keep HybridDish surface-area convention \(4\pi/1800\)
  (30 min doubling algebra) unless a later named freeze cites Senn.
- Clamp \(K_{\mathrm{clamp}}=4000\) as a **safety cap** (order of a
  jammed 100×100 µm monolayer), not as the claim mechanism.
- Spillover = particle removal at the open edge.
- **Do not claim the clamp is removed** until \(N(t)\) equilibrates
  on a clamp-off (or clamp-never-hit) run.

Acid / attractant ACs: **off**. PocketDish has no HybridDish acid-bias
story unless a later named protocol adds it.

## Readout (no ridge in the first jobs)

Primary: spatial \(L\) field + occupancy map \(R\) + density.
See [`READOUT_STACK.md`](READOUT_STACK.md).

Occupancy:

| Flag | Rule |
|---|---|
| ALIVE | `mean_R ≥ 0.05` |
| DEAD | `mean_R < 0.05` → **NOT_SCORED** |
| SATURATED | `mean_R > 0.95` |

Also report: domain-mean AHL, corner/interior AHL ratio,
\(N(t)\), Pe, \(L^2/D\), Damköhler \(k L^2/D\).

**No NRMSE, AUC, NARMA target, or λ grid** until occupancy ALIVE
and a later named task protocol exists. Week-6 may add **one** task;
it is not this file.

## Chemical boundaries

| Arm (names only) | Pocket walls | Open edge | Bus | Use |
|---|---|---|---|---|
| `CLOSED_NOFLUX` | NO_FLUX ×4 | none (closed) | none | occupancy sanity; should ALIVE at \(u=0.5\) |
| `OPEN_ABSORBING` | NO_FLUX ×3 | Dirichlet \(C=0\) | none | harsh leak diagnostic |
| `OPEN_BUS_3` | NO_FLUX ×3 | coupled to bus | \(v=3\) µm/s, inlet \(C=0\) | **backbone geometry** |
| `THROUGH_8` | — | through-pocket advection 8 µm/s | — | **negative control** (expect DEAD) |

Robin \(h\) only if labelled `HYPOTHETICAL_DESIGN_ENVELOPE`. Do not
import millimetre-dish \(h=0.159\) µm/s as calibration.

## Standardized drives (field / occupancy — no labels)

1. `STEP_ON`: \(u=0.5\) for 18 000 s (HybridDish warmup length) or
   until AHL/\(R\) plateau, whichever is declared in the job PROTOCOL.
2. `SINGLE_PULSE`: command 1 for 75 s, then zeros (HybridDish pulse
   length) — leak vs decay.
3. Through-flow control: `THROUGH_8` on the same step.

WashoutReset lesson: chemical ≠ reporter. After a pulse, AHL can die
in ~5×300 s; \(L\) does not. Do not call AHL drain a sequential-patient
reset.

## First authorized job — **done**

[`POCKETDISH_TRANSPORT_FROZEN_BUILDER_PROMPT.md`](POCKETDISH_TRANSPORT_FROZEN_BUILDER_PROMPT.md)
ran as `examples/BSimReservoirPlanPocketDishT0/`.
Standing: [`T0_STANDING.md`](T0_STANDING.md).

`OPEN_BUS_3` STEP_ON `mean_R = 0.0318` **DEAD**. Stop list in that
prompt is now in force: no NARMA, no \(K/k/J_{\max}\) retune, no
week-3 population Java.

Not authorized here: PocketOsc, PocketRect, two-way, pucks, NARMA
ridge, IPC, glucose, second strain, Grober BEM, any extra T0 arm.

## PocketNeck-T2 — **done** (named door, not a T0 retune)

[`POCKETNECK_T2_FROZEN_BUILDER_PROMPT.md`](POCKETNECK_T2_FROZEN_BUILDER_PROMPT.md)
ran as `examples/BSimReservoirPlanPocketNeckT2/`.
Standing: [`T2_STANDING.md`](T2_STANDING.md).

`W100_FLUSH` `mean_R=0.03179` **DEAD** (T0 class). `W50_L20`,
`W20_L20`, `W10_L20` **ALIVE**. Stop list in that prompt is in
force: no NARMA, no extra widths, no LuxI, no living Java from T2.
