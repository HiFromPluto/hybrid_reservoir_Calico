# Lane C — living clocks

**Object:** `LANE_C_LIVING_CLOCKS`  
**Device:** `LC0_OPEN_BATH`  
**Status:** extras closed. Claim note is the object.

Read this first:
[`../PocketDish/LANE_C_LIVING_CLOCKS_CLAIM.md`](../PocketDish/LANE_C_LIVING_CLOCKS_CLAIM.md).

This tree is **not** Lane A with motility on, **not** a fifth Lane B
occupation job, **not** `reservoir_new` on N0, and **not** the
packed-pocket stack (P0 / C1c / A1) that
[`../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md`](../PocketDish/PAPER_OBJECTS_AND_DESIGN_LANES.md)
§8 also called “Lane C”.

| File | Role |
|---|---|
| [`PROTOCOL.md`](PROTOCOL.md) | Frozen identity and gates |
| [`configs/protocol_lc0.json`](configs/protocol_lc0.json) | Machine-readable freeze |
| [`check_lc0.py`](check_lc0.py) | Checker |
| [`../PocketDish/LANE_C_INVESTIGATION.md`](../PocketDish/LANE_C_INVESTIGATION.md) | Why this object exists |
| [`../PocketDish/LANE_C_LC0_REPLICATION_STANDING.md`](../PocketDish/LANE_C_LC0_REPLICATION_STANDING.md) | LC0 standing (PASS; kept) |
| [`PROTOCOL_TWO_GEN.md`](PROTOCOL_TWO_GEN.md) | Two-generation polish identity |
| [`configs/protocol_lc0_two_gen.json`](configs/protocol_lc0_two_gen.json) | Two-gen freeze |
| [`check_lc0_two_gen.py`](check_lc0_two_gen.py) | Two-gen checker |
| [`../PocketDish/LANE_C_LC0_TWO_GEN_STANDING.md`](../PocketDish/LANE_C_LC0_TWO_GEN_STANDING.md) | Two-gen standing |
| [`PROTOCOL_DEATH.md`](PROTOCOL_DEATH.md) | LC1 death identity (Schink; growth off) |
| [`configs/protocol_lc1.json`](configs/protocol_lc1.json) | LC1 freeze |
| [`check_lc1.py`](check_lc1.py) | LC1 checker |
| [`../PocketDish/LANE_C_LC1_DEATH_STANDING.md`](../PocketDish/LANE_C_LC1_DEATH_STANDING.md) | LC1 standing (FAIL; kept) |
| [`PROTOCOL_DEATH_BINOMIAL.md`](PROTOCOL_DEATH_BINOMIAL.md) | LC1b binomial \(2\sigma\) identity |
| [`configs/protocol_lc1b.json`](configs/protocol_lc1b.json) | LC1b freeze |
| [`check_lc1b.py`](check_lc1b.py) | LC1b checker |
| [`../PocketDish/LANE_C_LC1B_DEATH_BINOMIAL_STANDING.md`](../PocketDish/LANE_C_LC1B_DEATH_BINOMIAL_STANDING.md) | LC1b standing |

Command: `ant lanec-lc0` then `python examples/LaneC_LivingClocks/check_lc0.py`.  
Two-gen polish: `ant lanec-lc0-two-gen` then
`python examples/LaneC_LivingClocks/check_lc0_two_gen.py`.  
Death extra: `ant lanec-lc1` then
`python examples/LaneC_LivingClocks/check_lc1.py`.  
Binomial clock: `ant lanec-lc1b` then
`python examples/LaneC_LivingClocks/check_lc1b.py`.  
No `Math.random()`. LC0/two-gen: death off. LC1/LC1b: growth off.
LC1 seed 101. LC1b seed 303.
