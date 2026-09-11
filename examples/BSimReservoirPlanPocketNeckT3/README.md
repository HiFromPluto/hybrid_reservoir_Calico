# PocketNeck-T3

Living BSim occupancy + `preview` for the PocketNeck chip. Not a task
scout. No ridge. No NARMA. No LuxI.

Question: does the T2 leak story (flush 100 µm door **DEAD**, Kim-class
\(W=20\) neck **ALIVE**) survive the real BSim field and bacteria?

Architecture: `examples/PocketDish/PROTOCOL.md`.
Chip plan: `examples/PocketDish/CHIP_PLAN.md`.
Neck freeze: `examples/PocketDish/GEOMETRY_NECK.md`.
This job: `PROTOCOL.md`.

Clone source: HybridDish ticker / drawer / Hill receiver patterns.
Do **not** edit `examples/HybridDish/`. Interior plastic and bus AHL
advection are wrapped in this example's ticker.

## Rerun

From this folder:

```
compile_and_run.cmd
python check_pocketneck_t3.py --theory
python run_occupancy_screen.py
python check_pocketneck_t3.py
```

Preview (IntelliJ / seeing the garage; **not** the score):

```
compile_and_run.cmd config\live_w20.properties preview
```

The immobilized AC is a **cyan** sphere at pocket centre, **10 µm
diameter** (GUV-class visual only). Green = bacteria. Yellow-orange =
AHL. Recompile after pulling if the AC is missing.

Close the preview window. Headless occupancy is the score.

Close the preview window. Headless occupancy is the score.
`BSim.preview()` sleeps `dt` per step and does not stop at 9000 s.

The checker cannot see a NARMA target and will refuse NARMA paths.

Smoke (not occupancy evidence):

```
python run_occupancy_screen.py --smoke
```

## Package map

| File | Role |
|---|---|
| `PROTOCOL.md` | Job freeze; cites `CHIP_PLAN.md`, `GEOMETRY_NECK.md`, HybridDish |
| `config/sim_config.properties` | Frozen dish / clocks / `dt=0.02` |
| `config/field_w100.properties` | `FIELD_W100_FLUSH` |
| `config/field_w20.properties` | `FIELD_W20_L20` |
| `config/live_w100.properties` | `LIVE_W100_FLUSH` |
| `config/live_w20.properties` | `LIVE_W20_L20` + preview |
| `bsim/BSimPocketNeckT3.java` | BSim job (`package PocketNeckT3`) |
| `compile_and_run.cmd` | HybridDish-style compile/run |
| `run_occupancy_screen.py` | Four-arm headless driver |
| `check_pocketneck_t3.py` | Occupancy flags, mass status; NARMA-blind |
| `results/OCCUPANCY_SCREEN.md` | Occupancy + leak status |

Do not retune \(K\) or \(J_{\max}\) if the neck is DEAD.
Do not start T4 growth-on from this package.
