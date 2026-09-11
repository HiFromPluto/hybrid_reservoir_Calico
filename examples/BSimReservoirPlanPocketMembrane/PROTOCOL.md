# PocketMembrane — particle-closed pocket, chemically open neck

Frozen 2026-08-20. Job-level only. Architecture freeze:
`examples/PocketDish/PROTOCOL.md`. Chip direction:
`examples/PocketDish/CHIP_PLAN.md`. Fill order:
`examples/PocketDish/FILL_DESIGN.md`. T5 standing (cite only; do not
rewrite population numbers): `examples/PocketDish/T5_STANDING.md`.
T4 standing (cite only): `examples/PocketDish/T4_STANDING.md`.
Clone source (do not edit T4 results):
`examples/BSimReservoirPlanPocketNeckT4/`.
Clamp-off / `replicate()` hard cap (cite only; do not copy weir or
\(k_{\mathrm{ov}}\); do not edit T5 results):
`examples/BSimReservoirPlanPocketMonolayerT5/`.

This file does not authorize NARMA, Mackey-Glass, waveform AUC,
ridge, a lambda grid, LuxI, vesicles, acid death, attractant
chemotaxis, Grober, two-way, extra widths, a circle/hex twin,
PocketCascade, shrinking the AHL box, a \(W=10\) hunt, T5 spring
retune, extra AHL membrane resistance (pore tortuosity), or a
HybridDish rewrite. Do not edit any GATE_EVIDENCE.md. Do not retune
K, J_max, W, L_n, growth, or the T5 spring after seeing N.

## Scientific question

With T3 AHL occupancy frozen, if cells **cannot** leave through the
neck, does HybridDish surface-area growth keep a colony on
`W20_L20`?

Hypothesis (predeclared, before looking at N):

- Occupancy stays ALIVE (T3 class `mean_R ~ 0.125`). Same chemical
  leak as T3–T5. If DEAD, the field port broke — stop, do not raise
  \(J_{\max}\).
- Spillover = 0. Any spill is a bounce bug. Stop and fix the bounce.
  Do not retune \(K\) or growth.
- Population unknown among CLAMPED (hit 4000 hard cap) |
  EQUILIBRIUM (plateau off the cap) | TRANSIENT (still rising at
  18000 s) | EMPTY (membrane failed). With clamp death off and a
  sealed garage, CLAMPED is the expected fill flag. Do not call
  CLAMPED a chemostat. EQUILIBRIUM is the only flag that allows
  that sentence.

## Path

Package `PocketMembrane.BSimPocketMembrane`. Clone T4 Java (T3 wall
mask, bus AHL advection v=3 um/s, dt=0.02 s, dx=5 um, cyan AC
visual). Same W20_L20, seed 101, N_0=50, pocket only. Growth on:
GROWTH_RATE = 4*pi/1800 /s on surface area. Acid off. Clamp *death*
off (`clamp.off=true`). Safety hard cap N=4000 in `replicate()`
only.

Particle membrane (ENGINEERING, Groisman / HEMA–EDMA / CNF class):
after `updatePosition`, if the centre is outside the pocket
rectangle [150,250] x [0,100] um, **mirror back into the pocket**.
The neck is chemically open and **cell-closed**. Do not allow the
T5 `inNeckX` corridor. Do not HybridDish-mirror the outer 400 um
box. Do not delete at the bus — cells should never get there. If a
centre is still in neck or bus after bounce, throw. No extra AHL
pore tortuosity; AHL stencil unchanged.

Not in this file: z weir, overlap spring, clamp death, LuxI.

Cells read AHL; they do not write AHL. Occupancy voxels = pocket
only.

## Frozen dish (T3 chip; clamp death off)

Chip W20_L20 only. No flush arm (a cell-closed flush is a different
closed box). No second open-neck control (T4/T5 already). Clocks /
J_max / Hill: T3 / T0 verbatim. N_0=50, seed 101, positions in the
pocket. Growth on. Acid off. AC immobilized point source, visual
cyan 10 um; no new dynamics. dt=0.02 s (T3 FTCS freeze). Chemical
box 400 x 200 x 10 um. Pocket particle height stays 10 um (not T5
weir).

Clamp death off. The 4000 cap is a hard safety, not HybridDish
p_removal. Do not call this a chemostat unless the population flag
is EQUILIBRIUM.

Preview AC: cyan sphere at (200, 50, 5) um, radius 5 um (10 um
diameter, GUV-class VISUAL_ONLY). Green = bacteria (drawn 8 um;
true BSim radius 1 um). Do not make the AC a BSimParticle. Do not
change J.

## Predeclared arms (do not add after seeing N)

LIVE_W20_MEM: T3 W20_L20 AHL; particle-closed pocket; growth on;
clamp death off. Fill try.

Do not add W=10, W=50, flush, or LuxI.

## Drive

STEP_ON u=0.5 for 18000 s. Sample N, pocket R, spillover at 20 s
(<=60 s). If 18000 s is not a plateau, write that. Do not extend
after looking. Smoke is 50 s and is not evidence (smoke <=200 s
must not be scored).

## Gates (not NRMSE)

Occupancy (pocket-field mean_R, same as T3): ALIVE if mean_R >= 0.05;
DEAD if mean_R < 0.05; SATURATED if mean_R > 0.95.
Expect ALIVE. If DEAD: stop, do not raise J_max.

Population flags, in this order:
EMPTY: N_end=0 (membrane implementation failed — fix bounce, do not
raise growth).
CLAMPED: N hit 4000 or 0.95*4000=3800 at any time (hard cap only;
death is off).
EQUILIBRIUM: last 20% of time, every sampled N within 10% of the
window mean, and never CLAMPED.
TRANSIENT: none of the above.

Also print N(0), N_max, N_end, cumulative spillover (must be 0),
births, clamp-removals (must be 0), pocket mean_R, mass residual
class.

Chemostat sentence allowed only if EQUILIBRIUM. If CLAMPED or
EQUILIBRIUM: fill recipe exists; still do not start NARMA or an
image-2 movie from this job.

LIVING_DECISION=STOP_AFTER_MEMBRANE_SCREEN

TRANSPORT_MODEL_STATUS: residual <=1% of dominant flux; C>=0;
walls ~0; cells did not write AHL.

## Preview

compile_and_run.cmd config\live_w20_mem.properties preview

Cyan AC, green cells, yellow AHL, garage outline. Cells must stay
in the garage outline. Close the window. Headless 18000 s is the
score. BSim.preview() does not stop at 18000 s.

## Output lines (no task Overall PASS)

OCCUPANCY_LIVE_W20_MEM=...
POPULATION_LIVE_W20_MEM=EMPTY|CLAMPED|EQUILIBRIUM|TRANSIENT
SPILLOVER_LIVE_W20_MEM=...
TRANSPORT_MODEL_STATUS=...
LIVING_DECISION=STOP_AFTER_MEMBRANE_SCREEN

## Commands

compile_and_run.cmd
compile_and_run.cmd config\live_w20_mem.properties
compile_and_run.cmd config\live_w20_mem.properties preview
python check_pocketmembrane.py --theory
python run_population_screen.py
python check_pocketmembrane.py

## Stop list

Do not add arms or widths after seeing N. Do not retune growth,
clocks, or J_max. Do not start NARMA, LuxI, cascade, or extra W.
Do not rewrite T0/T1/T2/T3/T4/T5 results. Do not edit HybridDish
Java, GATE_EVIDENCE.md, src/bsim/BSimChemicalField.java, or T4/T5
results.
