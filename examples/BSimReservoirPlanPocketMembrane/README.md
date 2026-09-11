# PocketMembrane

Particle-closed pocket, chemically open neck on the frozen T3
`W20_L20` chip. Groisman / nanoporous-wall class: cells stay, AHL
still leaks. Not a task scout. No ridge. No NARMA. No LuxI. No
cascade. Not more Danino (Danino *does* spill).

Question: with T3 AHL occupancy frozen, if cells cannot leave
through the neck, does HybridDish surface-area growth keep a colony?

Architecture: `examples/PocketDish/PROTOCOL.md`.
Chip plan: `examples/PocketDish/CHIP_PLAN.md`.
Fill order: `examples/PocketDish/FILL_DESIGN.md`.
T5 standing (cite only): `examples/PocketDish/T5_STANDING.md`.
T4 standing (cite only): `examples/PocketDish/T4_STANDING.md`.
This job: `PROTOCOL.md`.

Clone source: T4 Java (wall-mask, bus advection, cyan AC visual).
Clamp-off / `replicate()` hard cap: T5 (cite only). Do **not** copy
the 2 µm \(z\) weir or \(k_{\mathrm{ov}}\).
Growth: HybridDish `GROWTH_RATE=4*pi/1800`. Clamp *death* off.
Safety hard cap \(N=4000\) in `replicate()` only.
Do **not** edit `examples/HybridDish/` or T3/T4/T5 results.

## Rerun

From this folder:

```
compile_and_run.cmd
python check_pocketmembrane.py --theory
python run_population_screen.py
python check_pocketmembrane.py
```

Preview (IntelliJ / seeing the garage; **not** the score):

```
compile_and_run.cmd config\live_w20_mem.properties preview
```

The immobilized AC is a **cyan** sphere at pocket centre, **10 um
diameter** (GUV-class visual only). Green = bacteria. Yellow-orange =
AHL. Cells must stay in the garage outline. Close the preview
window. Headless 18000 s is the score. `BSim.preview()` sleeps `dt`
per step and does not stop at 18000 s.

The checker cannot see a NARMA target and will refuse NARMA paths.

Smoke (not population evidence):

```
python run_population_screen.py --smoke
```

## Package map

| File | Role |
|---|---|
| `PROTOCOL.md` | Job freeze; membrane declared before occupancy |
| `config/sim_config.properties` | Frozen dish / clocks / dt=0.02 / clamp.off |
| `config/live_w20_mem.properties` | `LIVE_W20_MEM` fill try + preview |
| `bsim/BSimPocketMembrane.java` | BSim job (`package PocketMembrane`) |
| `compile_and_run.cmd` | HybridDish-style compile/run |
| `run_population_screen.py` | One-arm headless driver |
| `check_pocketmembrane.py` | Occupancy + population flags; NARMA-blind |
| `results/POPULATION_SCREEN.md` | Population + occupancy status |

Do not retune K, J_max, or growth after seeing N.
Do not start cascade, NARMA, or LuxI from this package.
