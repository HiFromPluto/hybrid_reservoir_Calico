# PocketNeck-T3 — living BSim occupancy (job protocol)

Frozen 2026-08-20. Job-level only. Architecture freeze:
[`examples/PocketDish/PROTOCOL.md`](../PocketDish/PROTOCOL.md).
Chip direction:
[`examples/PocketDish/CHIP_PLAN.md`](../PocketDish/CHIP_PLAN.md).
Neck geometry (T2 freeze, do not invent):
[`examples/PocketDish/GEOMETRY_NECK.md`](../PocketDish/GEOMETRY_NECK.md).
T2 class numbers (cite only; do not rewrite):
[`examples/PocketDish/T2_STANDING.md`](../PocketDish/T2_STANDING.md),
[`examples/BSimReservoirPlanPocketNeckT2/PROTOCOL.md`](../BSimReservoirPlanPocketNeckT2/PROTOCOL.md).
T0 standing (cite only; do not rewrite):
[`examples/PocketDish/T0_STANDING.md`](../PocketDish/T0_STANDING.md).
Clone source (patterns only; **do not edit**):
[`examples/HybridDish/bsim/BSimHybridDish.java`](../HybridDish/bsim/BSimHybridDish.java).

This file does **not** authorize NARMA, Mackey–Glass, waveform AUC,
ridge, a λ grid, LuxI, vesicles, acid death, attractant chemotaxis,
Grober, two-way, T4 growth-on, extra widths \(W=50\) or \(W=10\), a
circle/hex twin, or a HybridDish rewrite.

Do not edit any `GATE_EVIDENCE.md`. Do not retune HybridDish
\(K,n,\tau_R,\tau_L\) or T0 \(J_{\max}\). Do not edit T0 / T1 / T2
results. Do not put LuxI in this package.

## Scientific question

On HybridDish clocks and T0 \(J_{\max}\), in **BSim**:

1. Does a Java wall-mask + bus advection reproduce the T2 leak
   story **without cells** (field port)?
2. With \(N_0=50\) receivers in the pocket (growth **off**), does
   pocket occupancy stay in the same class?

Hypothesis (predeclared, **before** occupancy):

- `FIELD_W100_FLUSH` **DEAD**, T0/T2 class (`mean_R` within ~15% of
  **0.0318**). If ALIVE, the mask/bus is wrong — **stop and fix**.
- `FIELD_W20_L20` **ALIVE**, T2 class (`mean_R` near **0.125**, not
  SATURATED). BSim explicit stencil need not match scipy BDF digit
  for digit. If DEAD with a clean mass budget: write DEAD, do **not**
  retune \(K\) or \(J_{\max}\).
- Living arms: cells **read** AHL (Hill); they do **not** write AHL.
  Growth **off** so T3 does not silently become T4 packing. Expect
  occupancy class to follow the field arms. Report \(N(t)\) (should
  stay 50 except spillover).

This is **not** HybridDish `preview` on a millimetre dish. HybridDish
`FLOW_SPEED` is a particle Stokes force, not a chemical bus.

## Path

BSim Java in this package (`PocketNeckT3.BSimPocketNeckT3`). Clone
ticker / drawer / `ReservoirBacterium` Hill \(R\)/\(L\) patterns from
HybridDish. Wrap **interior walls** and **bus AHL advection** in the
**example ticker**. Do not edit `src/bsim/BSimChemicalField.java`.

Occupancy means = **pocket boxes only** (neck and bus are diagnostics).

## Frozen dish (port T2, do not invent)

| Quantity | Value | Class |
|---|---|---|
| Pocket | \(100\times100\times10\) µm, FLOW=0 | PocketDish-A |
| Bus | \(400\times80\) µm, \(v=3.0\) µm/s, inlet \(C=0\), conservative outlet | PocketDish-A |
| Flush | T2 `W100_FLUSH` attachment (no neck) | T0 `OPEN_BUS_3` replay |
| Chip | T2 `W20_L20`: \(W=20\), \(L_n=20\), neck between pocket and bus | literature constriction |
| Grid | \(dx=dy=5\) µm, \(z=10\) µm (one box) | T0/T2 numerics |
| \(D,k,K,n,\tau_R,\tau_L\) | 159, 0.0033, 1.6, 2, 15, 1500 | HybridDish; do not retune |
| Conversion | 602.2 molecules/µm³/µM | DERIVED |
| \(J_{\max}\) | \(6.36\times10^5\) molecules/s at pocket centre | T0 verbatim |
| Occupancy voxels | **pocket only** | same rule as T0/T2 |

**`dt`:** HybridDish `0.05` s is **unstable** for explicit
`BSimChemicalField.diffuse` at \(dx=5\)
(\(k=D\,dt/dx^2=0.318\); 2-D/3-D FTCS). Freeze **`dt=0.02` s**
(`ENGINEERING` numerics so \(6D\,dt/dx^2\lesssim 1\)). This is **not**
a clock retune. Do not change \(k\) or \(K\) to stabilize.

**Flush attachment:** pocket \(x\in[150,250]\), \(y\in[0,100]\); bus
\(y\in[100,180]\), \(x\in[0,400]\); **no** neck voxels. Sim bound is
the 400×200×10 µm box that also holds the necked chip; \(y\in[180,200]\)
on the flush arm is interior plastic (inactive).

**Neck attachment:** pocket \(x\in[150,250]\), \(y\in[0,100]\); neck
centred on pocket `+y`, \(W=20\), \(y\in[100,120]\); bus \(y\in[120,200]\),
\(x\in[0,400]\). Neck `FLOW=0`. Do not steal neck volume from the pocket.

Walls: outer rectangle `setSolid`. Interior plastic (not pocket, not
neck, not bus) = inactive boxes: no diffusion in/out; conc held at 0
each step after the engine-style update. Implemented in the ticker
(masked conservative FTCS on active boxes only, then decay, then
zero walls, then first-order upwind advection \(v_x=3\) µm/s **in bus
voxels only**). Label `ENGINEERING` port of T2.

Particles (living arms only): seed **inside the pocket**. Mirror on
pocket solid walls (and neck walls). **Remove** if a cell centre
enters the bus (spillover). Do not mirror on the outer 400 µm box
the way HybridDish mirrors \(y\).

Acid field: **off**. Attractant / repellent: unused (do not `setGoal`).
Clamp death: **off**. Growth: **off** (`setSurfaceAreaGrowthRate(0)`).

## Predeclared arms (do not add after occupancy)

| Arm | Cells | Door | Role |
|---|---|---|---|
| `FIELD_W100_FLUSH` | none | flush 100 µm | T2 field port. Must **DEAD**, T0 class |
| `FIELD_W20_L20` | none | \(W=20\), \(L_n=20\) | T2 field port. Expect **ALIVE** |
| `LIVE_W100_FLUSH` | \(N_0=50\), seed **101**, growth **off** | flush | living negative |
| `LIVE_W20_L20` | \(N_0=50\), seed **101**, growth **off** | \(W=20\) | living chip + **`preview`** |

Do not add `W=50`, `W=10`, growth-on, or a circle after the first
occupancy number. Growth-on is T4.

## Drive

`STEP_ON` \(u=0.5\) for **9000 s** (T0/T2 plateau). Sample ≤60 s
(this job: 20 s). Occupancy on terminal pocket field: mean Hill \(R\)
of **pocket boxes** (same gate as T2). Also print mean cell \(R\) and
mean cell \(L\) on living arms (diagnostic; gate is the field).

Plateau: terminal pocket `mean_R` within 5% of the last 10% of time.

No `dt` sensitivity is declared. Occupancy uses frozen `dt=0.02` only.

## Occupancy flags (not NRMSE)

| Flag | Rule |
|---|---|
| `ALIVE` | pocket-field `mean_R ≥ 0.05` and not `SATURATED` |
| `DEAD` | pocket-field `mean_R < 0.05` |
| `SATURATED` | pocket-field `mean_R > 0.95` (unexpected; fix units, not \(K\)) |

Also print: pocket-mean AHL; neck-mean AHL (NA on flush); bus-mean AHL;
\(N(t)\) start/end; spillover count; mass budget (production, decay,
bus outlet, residual); \(C\ge 0\); wall boxes stay ~0.

`TRANSPORT_MODEL_STATUS`: residual ≤1% of dominant flux on **field**
arms. Living arms: same budget plus “cells did not write AHL.”

If `FIELD_W100_FLUSH` is ALIVE or not T0-class: **stop**. Fix mask.
If `FIELD_W20_L20` is DEAD: **stop**. Do not retune. Compare stencil
to T2 (explicit vs BDF, `dt`).
If field arms match class and a living arm disagrees wildly: cells
are leaking through walls or sitting in the bus — fix BC, do not
retune \(K\).

## Preview (IntelliJ)

```
compile_and_run.cmd config\live_w20.properties preview
```

`preview` overrides `headless=true`. Camera: this box is ~400×200×10
µm, **not** HybridDish 1000×500. Pull the eye back by **bound Y**
(HybridDish lesson), not by thickness 10. Draw AHL + green cells.
Pocket/neck/bus outline in the drawer so the garage is obvious.

Headless occupancy is the score. Preview is for seeing it.
`BSim.preview()` loops until the window is closed (real-time sleep
`dt`); it is **not** occupancy evidence. Close the window.

## Gates that are not this job

No NRMSE, AUC, NARMA target, ridge, or λ. The checker **must refuse**
NARMA paths and must not import `BSimReservoirPlanNarma10b`.
Do not write Overall PASS/FAIL on a task. Lines:

- `OCCUPANCY_FIELD_W100_FLUSH=…`
- `OCCUPANCY_FIELD_W20_L20=…`
- `OCCUPANCY_LIVE_W100_FLUSH=…`
- `OCCUPANCY_LIVE_W20_L20=…`
- `TRANSPORT_MODEL_STATUS=…`
- `LIVING_DECISION=STOP_AFTER_T3_OCCUPANCY`

## Commands

```
compile_and_run.cmd
compile_and_run.cmd config\field_w100.properties
compile_and_run.cmd config\field_w20.properties
compile_and_run.cmd config\live_w100.properties
compile_and_run.cmd config\live_w20.properties
compile_and_run.cmd config\live_w20.properties preview
python check_pocketneck_t3.py --theory
python run_occupancy_screen.py
python check_pocketneck_t3.py
```

`--smoke` on the runner integrates 50 s only and must not be scored as
occupancy evidence.

## Stop list

- Do not add arms, \(W\), growth-on, or a circle after occupancy
- Do not retune HybridDish clocks or T0 \(J_{\max}\) to occupy
- Do not start NARMA, LuxI, or T4 from this job
- Do not call this an exact Danino/Prindle blueprint
- Do not rewrite T0 / T1 / T2 / Narma10b results
- Do not edit HybridDish Java or any `GATE_EVIDENCE.md`
- Do not edit `src/bsim/BSimChemicalField.java`
