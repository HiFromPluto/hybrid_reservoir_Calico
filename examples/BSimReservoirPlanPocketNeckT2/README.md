# PocketNeck-T2

Occupancy / leak screen for a **narrow neck** on frozen HybridDish
clocks and the T0 point source. Not a task scout. No ridge. No living
Java.

Question: does a Kim-class door leave PocketHill occupied when the
flush 100 µm `OPEN_BUS_3` door does not?

Architecture: `examples/PocketDish/PROTOCOL.md`.  
Neck freeze: `examples/PocketDish/GEOMETRY_NECK.md`.  
T0 job (flush clone): `examples/BSimReservoirPlanPocketDishT0/PROTOCOL.md`.  
This job: `PROTOCOL.md`.

## Rerun

From the repository root:

```
python examples/BSimReservoirPlanPocketNeckT2/check_pocketneck_t2.py --theory
python examples/BSimReservoirPlanPocketNeckT2/test_mass_budget.py
python examples/BSimReservoirPlanPocketNeckT2/run_occupancy_screen.py
python examples/BSimReservoirPlanPocketNeckT2/check_pocketneck_t2.py
```

The checker cannot see a NARMA target and will refuse NARMA paths.

Smoke (not occupancy evidence):

```
python examples/BSimReservoirPlanPocketNeckT2/run_occupancy_screen.py --smoke
```

Official occupancy uses \(dx=5\) µm. An optional `W10_L20` run at
\(dx=2.5\) µm is labelled `CONVERGENCE` and does not replace the
table.

## Package map

| File | Role |
|---|---|
| `PROTOCOL.md` | Job freeze; cites `GEOMETRY_NECK.md` and T0 PROTOCOL |
| `configs/dish.json` | Geometry and T0 clocks / \(J_{\max}\) |
| `configs/arms.json` | Four predeclared arms |
| `configs/drives.json` | `STEP_ON` only |
| `transport_model.py` | Reduced 2-D FV + Hill/\(L\); T0 stencil clone |
| `run_occupancy_screen.py` | Screen driver; writes CSV and maps |
| `check_pocketneck_t2.py` | Occupancy flags, mass status, markdown; NARMA-blind |
| `test_mass_budget.py` | Mass / voxel / flush-clone tests |
| `results/OCCUPANCY_SCREEN.md` | Occupancy + leak status |

Reduced model only. Week-2 living Java is not this package.
Do not retune \(K\) or \(J_{\max}\) if necks are DEAD.
