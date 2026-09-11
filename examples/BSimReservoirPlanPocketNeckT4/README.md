# PocketNeck-T4

Growth-on N(t) / spillover on the frozen T3 `W20_L20` chip. Not a task
scout. No ridge. No NARMA. No LuxI.

Question: with HybridDish surface-area growth on, does a population
stay in the garage, or does the door still win?

Architecture: `examples/PocketDish/PROTOCOL.md`.
Chip plan: `examples/PocketDish/CHIP_PLAN.md`.
T3 standing (cite only): `examples/PocketDish/T3_STANDING.md`.
This job: `PROTOCOL.md`.

Clone source: T3 Java (wall-mask, bus advection, cyan AC visual).
Growth / clamp algebra: HybridDish `GROWTH_RATE=4*pi/1800` and
`p_removal` at `CARRYING_CAPACITY=4000`. Acid off.
Do **not** edit `examples/HybridDish/` or T3 results.

## Rerun

From this folder:

```
compile_and_run.cmd
python check_pocketneck_t4.py --theory
python run_population_screen.py
python check_pocketneck_t4.py
```

Preview (IntelliJ / seeing the garage; **not** the score):

```
compile_and_run.cmd config\live_w20_growth.properties preview
```

The immobilized AC is a **cyan** sphere at pocket centre, **10 um
diameter** (GUV-class visual only). Green = bacteria. Yellow-orange =
AHL. Close the preview window. Headless 18000 s is the score.
`BSim.preview()` sleeps `dt` per step and does not stop at 18000 s.

The checker cannot see a NARMA target and will refuse NARMA paths.

Smoke (not population evidence):

```
python run_population_screen.py --smoke
```

## Package map

| File | Role |
|---|---|
| `PROTOCOL.md` | Job freeze; clamp algebra declared before occupancy |
| `config/sim_config.properties` | Frozen dish / clocks / dt=0.02 / growth on |
| `config/live_w20_growth.properties` | `LIVE_W20_GROWTH` primary + preview |
| `config/live_w100_growth.properties` | `LIVE_W100_GROWTH` flush diagnostic |
| `bsim/BSimPocketNeckT4.java` | BSim job (`package PocketNeckT4`) |
| `compile_and_run.cmd` | HybridDish-style compile/run |
| `run_population_screen.py` | Two-arm headless driver |
| `check_pocketneck_t4.py` | Occupancy + population flags; NARMA-blind |
| `results/POPULATION_SCREEN.md` | Population + occupancy status |

Do not retune K, J_max, or growth after seeing N.
Do not start T5, NARMA, or LuxI from this package.
