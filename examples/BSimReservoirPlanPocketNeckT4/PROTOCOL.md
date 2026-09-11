# PocketNeck-T4 — growth-on N(t) / spillover (job protocol)

Frozen 2026-08-20. Job-level only. Architecture freeze:
`examples/PocketDish/PROTOCOL.md`. Chip direction:
`examples/PocketDish/CHIP_PLAN.md`. T3 standing (cite only; do not
rewrite occupancy numbers): `examples/PocketDish/T3_STANDING.md`.
Clone source (do not edit T3 results):
`examples/BSimReservoirPlanPocketNeckT3/`. Growth / clamp patterns
only (do not edit HybridDish Java): `examples/HybridDish/bsim/`.

This file does not authorize NARMA, Mackey-Glass, waveform AUC,
ridge, a lambda grid, LuxI, vesicles, acid death, attractant
chemotaxis, Grober, two-way, extra widths, a circle/hex twin, T5
optical pictures, or a HybridDish rewrite. Do not edit any
GATE_EVIDENCE.md. Do not retune K or J_max. Do not retune growth
after seeing N.

## Scientific question

On the living W20_L20 chip (T3 occupancy already ALIVE), HybridDish
surface-area growth 4*pi/1800 (~30 min doubling algebra), seed 101,
N_0=50 in the pocket: does N(t) plateau off the clamp, hit the
safety cap, or still go empty through the neck?

Hypothesis (predeclared, before looking at N):

- Pocket-field occupancy on W20 stays ALIVE (AC still writes;
  growth does not change J_max). If DEAD, stop — do not raise J_max.
- Population is unknown. T3 empty-garage is the growth-off baseline.
  Write the flag. Do not retune growth rate after seeing N.
- Flush growth-on (diagnostic) should lose cells faster than the neck.

## Path

Package `PocketNeckT4.BSimPocketNeckT4`. Clone T3 wall mask, bus AHL
advection v=3 um/s in the ticker, dt=0.02 s, dx=5 um, pocket
occupancy voxels, spillover = remove when a cell centre enters the
bus. Copy the T3 cyan AC visual. Enable HybridDish
setSurfaceAreaGrowthRate and replicate(). Cells read AHL; they do
not write AHL. Occupancy voxels = pocket only.

## Frozen dish (T3 chip; growth on)

Chip W20_L20 only as the living geometry. Flush is a diagnostic arm,
not a width hunt. Clocks / J_max / Hill: T3 / T0 verbatim. N_0=50,
seed 101, positions in the pocket. Growth on:
GROWTH_RATE = 4*pi/1800 /s on surface area (ENGINEERING, same as
PocketDish PROTOCOL). Acid off. AC immobilized point source, visual
cyan 10 um; no new dynamics. dt=0.02 s (T3 FTCS freeze).

Clamp (declared before occupancy): HybridDish stochastic p_removal
with CARRYING_CAPACITY=4000. Not a hard cap. Acid still off.

T_removal = T_gen = 1800 s
p_removal = (dt / T_removal) * 2^(-(1 - N/K_clamp))

Each living cell draws this Bernoulli trial in action(). First
removal cause wins (clamp in action, then spillover in
updatePosition). The cap is a safety (jammed 100x100 monolayer
order), not the claim mechanism. Do not call this a chemostat unless
the population flag is EQUILIBRIUM.

Preview AC: cyan sphere at (200, 50, 5) um, radius 5 um (10 um
diameter, GUV-class VISUAL_ONLY). Green = bacteria (drawn 8 um;
true BSim radius 1 um). Do not make the AC a BSimParticle. Do not
change J.

## Predeclared arms (do not add after seeing N)

LIVE_W20_GROWTH: growth on, W=20 L_n=20, primary N(t) + occupancy.
LIVE_W100_GROWTH: growth on, flush 100 um, diagnostic (door faster?).

Do not add W=10, W=50, field-only, or LuxI.

## Drive

STEP_ON u=0.5 for 18000 s. Sample N, pocket R, spillover at 20 s
(<=60 s). If 18000 s is not a plateau, write that. Do not extend
after looking. Smoke is 50 s and is not evidence (smoke <=200 s
must not be scored).

## Gates (not NRMSE)

Occupancy (pocket-field mean_R, same as T3): ALIVE if mean_R >= 0.05
on LIVE_W20_GROWTH; DEAD if mean_R < 0.05; SATURATED if mean_R > 0.95.
Expect ALIVE. If DEAD: stop, do not raise J_max.

Population flags on LIVE_W20_GROWTH, in this order:
EMPTY: N_end=0 (door still wins with growth).
CLAMPED: N reached 4000 or 0.95*4000=3800 at any time.
EQUILIBRIUM: last 20% of time, every sampled N within 10% of the
window mean, and never CLAMPED.
TRANSIENT: none of the above.

Also print N(0), N_max, N_end, cumulative spillover, births,
clamp-removals, pocket mean_R, mass residual class.

Chemostat sentence allowed only if EQUILIBRIUM. If CLAMPED, the cap
is still the story. If EMPTY, growth lost to the neck — same physics
as T3 plus birth, not a K retune.

LIVING_DECISION=STOP_AFTER_T4_POPULATION

TRANSPORT_MODEL_STATUS: residual <=1% of dominant flux; C>=0;
walls ~0; cells did not write AHL.

## Preview

compile_and_run.cmd config\live_w20_growth.properties preview

Cyan AC, green cells, yellow AHL, garage outline. Close the window.
Headless 18000 s is the score. BSim.preview() does not stop at 18000 s.

## Output lines (no task Overall PASS)

OCCUPANCY_LIVE_W20_GROWTH=...
POPULATION_LIVE_W20_GROWTH=EMPTY|CLAMPED|EQUILIBRIUM|TRANSIENT
OCCUPANCY_LIVE_W100_GROWTH=... (if run)
POPULATION_LIVE_W100_GROWTH=... (if run)
TRANSPORT_MODEL_STATUS=...
LIVING_DECISION=STOP_AFTER_T4_POPULATION

## Commands

compile_and_run.cmd
compile_and_run.cmd config\live_w20_growth.properties
compile_and_run.cmd config\live_w20_growth.properties preview
compile_and_run.cmd config\live_w100_growth.properties
python check_pocketneck_t4.py --theory
python run_population_screen.py
python check_pocketneck_t4.py

## Stop list

Do not add arms or widths after seeing N. Do not retune growth,
clocks, or J_max. Do not start NARMA, LuxI, T5, or AC payload
dynamics. Do not rewrite T0/T1/T2/T3 results. Do not edit HybridDish
Java, GATE_EVIDENCE.md, src/bsim/BSimChemicalField.java, or T3 results.

