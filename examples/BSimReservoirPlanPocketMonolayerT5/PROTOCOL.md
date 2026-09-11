# PocketMonolayer-T5 — Danino-class height weir (job protocol)

Frozen 2026-08-20. Job-level only. Architecture freeze:
`examples/PocketDish/PROTOCOL.md`. Chip direction:
`examples/PocketDish/CHIP_PLAN.md`. Fill order:
`examples/PocketDish/FILL_DESIGN.md`. T4 standing (cite only; do not
rewrite population numbers): `examples/PocketDish/T4_STANDING.md`.
Clone source (do not edit T4 results):
`examples/BSimReservoirPlanPocketNeckT4/`.

This file does not authorize NARMA, Mackey-Glass, waveform AUC,
ridge, a lambda grid, LuxI, vesicles, acid death, attractant
chemotaxis, Grober, two-way, extra widths, a circle/hex twin,
PocketMembrane, PocketCascade, shrinking the AHL box to 2 um, or a
HybridDish rewrite. Do not edit any GATE_EVIDENCE.md. Do not retune
K, J_max, W, L_n, or growth after seeing N. Do not retune k_ov
after seeing N.

## Scientific question

With T3 AHL occupancy frozen (pocket 100 x 100 x 10 um, J_max
unchanged), does a 2 um particle weir plus labelled excluded volume
keep a growing colony on W20_L20, or do swimmers still empty the
garage?

Hypothesis (predeclared, before looking at N):

- Occupancy on both arms stays ALIVE (T3 class mean_R ~ 0.125).
  If DEAD, the field port broke — stop, do not raise J_max.
- Z2_NOREP EMPTY (height-only). If it EQUILIBRIUM, say so — then
  BSim motility in a thin slab was enough and repulsion is extra.
- Z2_REP unknown. Do not retune the spring after seeing N.
  Do not add "motility off when overlapping" after occupancy.

## Path

Package `PocketMonolayerT5.BSimPocketMonolayerT5`. Clone T4 Java
(T3 wall mask, bus AHL advection v=3 um/s, dt=0.02 s, dx=5 um,
cyan AC visual). Same W20_L20, seed 101, N_0=50, pocket only.
Growth on: GROWTH_RATE = 4*pi/1800 /s on surface area. Acid off.
Clamp *death* off (clamp.off=true). Safety hard cap N=4000 still
recorded. Spillover still = remove at the bus.

Particle z weir (both arms): after updatePosition, clamp z to
[radius, h_weir - radius] with h_weir=2.0 um (ENGINEERING,
Fluigi/Danino ~1.65 um class, rounded). If 2R > h_weir, clamp z
to weir-box centre. Do not change true radius 1 um. Do not shrink
the AHL box to 2 um.

Excluded volume (Z2_REP only), predeclared: each tick, all pairs
with d < 2R get equal-and-opposite spring
F = k_ov (2R - d) r_hat in xy (k_ov=50 pN/um ENGINEERING, not
fitted). Cap force 100 pN. O(N^2) is fine at N <= 4000. Do not
change k_ov. Tutorial interaction()/outerDistance is an overlap
flag only; this job adds the force.

Cells read AHL; they do not write AHL. Occupancy voxels = pocket
only.

## Frozen dish (T3 chip; clamp death off)

Chip W20_L20 only. No flush arm (T4 already). No membrane.
Clocks / J_max / Hill: T3 / T0 verbatim. N_0=50, seed 101,
positions in the pocket. Growth on. Acid off. AC immobilized
point source, visual cyan 10 um; no new dynamics. dt=0.02 s
(T3 FTCS freeze). Chemical box 400 x 200 x 10 um.

Clamp death off. The 4000 cap is a hard safety (jammed 100x100
monolayer order), not HybridDish p_removal. Do not call this a
chemostat unless the population flag is EQUILIBRIUM.

Preview AC: cyan sphere at (200, 50, 5) um, radius 5 um (10 um
diameter, GUV-class VISUAL_ONLY). Green = bacteria (drawn 8 um;
true BSim radius 1 um). Do not make the AC a BSimParticle. Do not
change J.

## Predeclared arms (do not add after seeing N)

Z2_NOREP: z weir 2 um, no repulsion, growth on, clamp death off.
  Height-only control; expect EMPTY.
Z2_REP: z weir 2 um, k_ov=50 pN/um, growth on, clamp death off.
  Danino-class try.

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
EMPTY: N_end=0.
CLAMPED: N hit 4000 or 0.95*4000=3800 at any time (hard cap only;
death is off).
EQUILIBRIUM: last 20% of time, every sampled N within 10% of the
window mean, and never CLAMPED.
TRANSIENT: none of the above.

Also print N(0), N_max, N_end, cumulative spillover, births,
clamp-removals (must be 0), pocket mean_R, mass residual class.

Chemostat sentence allowed only if EQUILIBRIUM. If both arms EMPTY:
stop. Danino weir failed in this engine. Next named job is
PocketMembrane, not a spring hunt.

LIVING_DECISION=STOP_AFTER_MONOLAYER_SCREEN

TRANSPORT_MODEL_STATUS: residual <=1% of dominant flux; C>=0;
walls ~0; cells did not write AHL.

## Preview

compile_and_run.cmd config\z2_rep.properties preview

Cyan AC, green cells, yellow AHL, garage outline. Chemical box
still 400 x 200 x 10. Close the window. Headless 18000 s is the
score. BSim.preview() does not stop at 18000 s.

## Output lines (no task Overall PASS)

OCCUPANCY_Z2_NOREP=...
POPULATION_Z2_NOREP=...
OCCUPANCY_Z2_REP=...
POPULATION_Z2_REP=...
TRANSPORT_MODEL_STATUS=...
LIVING_DECISION=STOP_AFTER_MONOLAYER_SCREEN

## Commands

compile_and_run.cmd
compile_and_run.cmd config\z2_norep.properties
compile_and_run.cmd config\z2_rep.properties
compile_and_run.cmd config\z2_rep.properties preview
python check_pocketmonolayer_t5.py --theory
python run_population_screen.py
python check_pocketmonolayer_t5.py

## Stop list

Do not add arms or widths after seeing N. Do not retune growth,
clocks, J_max, or k_ov. Do not start NARMA, LuxI, membrane, or
cascade. Do not rewrite T0/T1/T2/T3/T4 results. Do not edit
HybridDish Java, GATE_EVIDENCE.md, src/bsim/BSimChemicalField.java,
or T4 results.
