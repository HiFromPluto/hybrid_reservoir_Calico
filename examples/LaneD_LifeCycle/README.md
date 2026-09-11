# Lane D — life cycle (phase-1 device only)

**Object:** `LANE_D_LIFE_CYCLE`  
**Device:** `LANE_D_CLOSED_BATH`  
**Extras:** D0 PASS; D1 REPORT; D2 ARGUMENT_HOLDS; D3.1 PORT_HOLDS
(Fig. 1g \(\lambda(t)\) overlay remains a `SCOPE_NOTE` **fact**);
D3.2 CONVERGES; D3.3 PORT_HOLDS (lifted that overlay for D4);
`D4_THREE_PHASE` PASS; `D5_GATES` PASS.  
**Status:** three-phase object exists and its composition gates hold.
Not Erickson BUILD. LC1 FAIL stays FAIL.

Read the plan erratum first:
[`../PocketDish/LANE_D_PLAN.md`](../PocketDish/LANE_D_PLAN.md) (F1–F8).
Investigation erratum:
[`../PocketDish/LIFE_CYCLE_COMPOSITION_INVESTIGATION.md`](../PocketDish/LIFE_CYCLE_COMPOSITION_INVESTIGATION.md)
(E1–E3).

This tree is **not** `LC0_OPEN_BATH` with a depleting \(C\), **not** a
`NutrientField` mode, **not** Erickson, **not** Job 6 famine, and **not**
the packed-pocket stack (P0 / C1c / A1) that
[`../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md`](../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md)
§8 called “Lane C”. `LANE_C_LIVING_CLOCKS` is the birth/death clocks.
Do not merge them.

| File | Role |
|---|---|
| [`PROTOCOL_D0.md`](PROTOCOL_D0.md) | Frozen D0 identity and gates |
| [`configs/protocol_d0.json`](configs/protocol_d0.json) | Machine-readable freeze |
| [`check_d0.py`](check_d0.py) | D0 checker |
| [`PROTOCOL_D1.md`](PROTOCOL_D1.md) | Frozen D1 windows and integral |
| [`configs/protocol_d1.json`](configs/protocol_d1.json) | D1 freeze |
| [`check_d1.py`](check_d1.py) | D1 metrics on `d0_grow.csv` |
| [`PROTOCOL_D33.md`](PROTOCOL_D33.md) | Frozen Table 2 / Fig. 4b,c overlay |
| [`check_d33.py`](check_d33.py) | D3.3 checker |
| [`PROTOCOL_D4.md`](PROTOCOL_D4.md) | Frozen D4 three-phase identity |
| [`configs/protocol_d4.json`](configs/protocol_d4.json) | D4 freeze |
| [`check_d4.py`](check_d4.py) | D4 checker |
| [`../PocketDish/LANE_D_PLAN.md`](../PocketDish/LANE_D_PLAN.md) | Why this object exists |
| [`../PocketDish/LANE_D_D0_CLOSED_BATH_STANDING.md`](../PocketDish/LANE_D_D0_CLOSED_BATH_STANDING.md) | D0 standing (PASS; kept) |
| [`../PocketDish/LANE_D_D1_MEASURED_DECLINE_STANDING.md`](../PocketDish/LANE_D_D1_MEASURED_DECLINE_STANDING.md) | D1 standing (REPORT) |
| [`PROTOCOL_D5.md`](PROTOCOL_D5.md) | Frozen D5 composition gates |
| [`check_d5.py`](check_d5.py) | D5 checker |
| [`../PocketDish/LANE_D_D4_THREE_PHASE_STANDING.md`](../PocketDish/LANE_D_D4_THREE_PHASE_STANDING.md) | D4 standing |
| [`../PocketDish/LANE_D_D5_GATES_STANDING.md`](../PocketDish/LANE_D_D5_GATES_STANDING.md) | D5 standing |

D0: `ant laned-d0` then `python examples/LaneD_LifeCycle/check_d0.py`.  
D1: `python examples/LaneD_LifeCycle/check_d1.py` (reads `d0_grow.csv`; no new Java).  
No `Math.random()`. Seed 101. Death off. Sink uses \(\rho_{\mathrm{cell}}\),
not \(\rho_0\). Elongation is incremental; `elongateCited` is not called.

D3.1: `python examples/LaneD_LifeCycle/check_d31.py`  
D3.2: `python examples/LaneD_LifeCycle/check_d32.py`  
D3.3: `python examples/LaneD_LifeCycle/check_d33.py`  
D4: `ant laned-d4` then `python examples/LaneD_LifeCycle/check_d4.py`.  
D5: `ant laned-d5` then `python examples/LaneD_LifeCycle/check_d5.py`.
