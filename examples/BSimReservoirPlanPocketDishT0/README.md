# PocketDish-T0

Transport / occupancy screen for frozen **PocketDish-A** (Danino-class
square pocket + horizontal bus) with HybridDish \(D,k,K,n,\tau\) and the
**predeclared** PocketDish \(J_{\max}\). Not a task scout. No ridge.

Architecture: `examples/PocketDish/PROTOCOL.md`.  
This job: `PROTOCOL.md`.

## Rerun

From the repository root:

```
python examples/BSimReservoirPlanPocketDishT0/check_pocketdish_t0.py --theory
python examples/BSimReservoirPlanPocketDishT0/test_mass_budget.py
python examples/BSimReservoirPlanPocketDishT0/run_transport_screen.py
python examples/BSimReservoirPlanPocketDishT0/check_pocketdish_t0.py
```

The checker cannot see a NARMA target and will refuse NARMA paths.

Smoke (not occupancy evidence):

```
python examples/BSimReservoirPlanPocketDishT0/run_transport_screen.py --smoke
```

## Package map

| File | Role |
|---|---|
| `PROTOCOL.md` | Job freeze; cites PocketDish architecture PROTOCOL |
| `configs/dish.json` | Geometry and clocks |
| `configs/arms.json` | Four predeclared arms |
| `configs/drives.json` | STEP_ON / SINGLE_PULSE |
| `configs/sim_config_*.properties` | Living-ready dish record (not executed by T0 unless the living rule fires) |
| `transport_model.py` | Reduced 2-D FV + Hill/\(L\) |
| `run_transport_screen.py` | Screen driver; writes CSV and maps |
| `check_pocketdish_t0.py` | Occupancy flags, mass status, markdown; NARMA-blind |
| `test_mass_budget.py` | Mass-identity tests A–E |
| `results/TRANSPORT_SCREEN.md` | Occupancy + transport status |

Reduced model only. Living BSim seed 101 is gated in `PROTOCOL.md`.
Week-3 population Java is not this package.
