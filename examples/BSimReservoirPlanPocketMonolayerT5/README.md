# PocketMonolayer-T5

Danino-class **height weir** on the frozen T3 `W20_L20` chip. Particle
\(z\le 2\) µm; chemical box stays 10 µm. Not a task scout. No ridge.
No NARMA. No LuxI. No membrane.

Question: with T3 AHL occupancy frozen, does a 2 µm particle weir
plus labelled excluded volume keep a growing colony, or do swimmers
still empty the garage?

Architecture: `examples/PocketDish/PROTOCOL.md`.
Chip plan: `examples/PocketDish/CHIP_PLAN.md`.
Fill order: `examples/PocketDish/FILL_DESIGN.md`.
T4 standing (cite only): `examples/PocketDish/T4_STANDING.md`.
This job: `PROTOCOL.md`.

Clone source: T4 Java (wall-mask, bus advection, cyan AC visual).
Growth: HybridDish `GROWTH_RATE=4*pi/1800`. Clamp *death* off.
Safety hard cap \(N=4000\) still recorded.
Do **not** edit `examples/HybridDish/` or T3/T4 results.

## Rerun

From this folder:

```
compile_and_run.cmd
python check_pocketmonolayer_t5.py --theory
python run_population_screen.py
python check_pocketmonolayer_t5.py
```

Preview (IntelliJ / seeing the garage; **not** the score):

```
compile_and_run.cmd config\z2_rep.properties preview
```

The immobilized AC is a **cyan** sphere at pocket centre, **10 um
diameter** (GUV-class visual only). Green = bacteria. Yellow-orange =
AHL. Chemical box still 400×200×**10**. Close the preview window.
Headless 18000 s is the score. `BSim.preview()` sleeps `dt` per step
and does not stop at 18000 s.

The checker cannot see a NARMA target and will refuse NARMA paths.

Smoke (not population evidence):

```
python run_population_screen.py --smoke
```

## Package map

| File | Role |
|---|---|
| `PROTOCOL.md` | Job freeze; weir + spring declared before occupancy |
| `config/sim_config.properties` | Frozen dish / clocks / dt=0.02 / clamp.off |
| `config/z2_norep.properties` | `Z2_NOREP` height-only control |
| `config/z2_rep.properties` | `Z2_REP` Danino-class try + preview |
| `bsim/BSimPocketMonolayerT5.java` | BSim job (`package PocketMonolayerT5`) |
| `compile_and_run.cmd` | HybridDish-style compile/run |
| `run_population_screen.py` | Two-arm headless driver |
| `check_pocketmonolayer_t5.py` | Occupancy + population flags; NARMA-blind |
| `results/POPULATION_SCREEN.md` | Population + occupancy status |

Do not retune K, J_max, growth, or \(k_{\mathrm{ov}}\) after seeing N.
Do not start membrane, cascade, NARMA, or LuxI from this package.
