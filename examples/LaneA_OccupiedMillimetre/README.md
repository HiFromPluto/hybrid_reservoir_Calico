# Lane A — occupied millimetre dish (N0)

Named successor of HybridDish: one-way AHL → Hill \(R\) → slow \(L\) on
`1000×500×10` µm, `FLOW=0`, N0 transport, occupancy first.

**Not** paper 1. **Not** a rewrite of Narma10b 0.928. **Not** C1c / Fig. 4b.
**Not** Lane B.

Protocol (frozen before traces): [`PROTOCOL.md`](PROTOCOL.md).

```
ant lane-a-occupancy
python examples/LaneA_OccupiedMillimetre/check_lane_a.py
ant lane-a-carrier
python examples/LaneA_OccupiedMillimetre/check_lane_a_carrier.py
python examples/LaneA_OccupiedMillimetre/check_lane_a_carrier_audit.py
ant lane-a-narma10
python examples/LaneA_OccupiedMillimetre/check_lane_a_narma10.py
ant lane-a-mg
python examples/LaneA_OccupiedMillimetre/check_lane_a_mg.py
ant lane-a-lorenz
python examples/LaneA_OccupiedMillimetre/check_lane_a_lorenz.py
python examples/LaneA_OccupiedMillimetre/check_lane_a_waveform.py --u-only
ant lane-a-waveform
python examples/LaneA_OccupiedMillimetre/check_lane_a_waveform.py
```

Occupancy standing: [`../PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md`](../PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md).  
Carrier extra: [`PROTOCOL_CARRIER.md`](PROTOCOL_CARRIER.md), standing **FAIL vs field** [`../PocketDish/LANE_A_CARRIER_STANDING.md`](../PocketDish/LANE_A_CARRIER_STANDING.md) (occupancy still PASS).

Path after carrier FAIL: [`../PocketDish/LANE_A_PATH_AFTER_CARRIER.md`](../PocketDish/LANE_A_PATH_AFTER_CARRIER.md).
Carrier maps audit (**SCOPE_NOTE**; carrier still FAIL): [`PROTOCOL_CARRIER_AUDIT.md`](PROTOCOL_CARRIER_AUDIT.md), standing [`../PocketDish/LANE_A_CARRIER_AUDIT_STANDING.md`](../PocketDish/LANE_A_CARRIER_AUDIT_STANDING.md).
NARMA-10 extra (**system PASS**; not 0.928; carrier still FAIL): [`PROTOCOL_NARMA10.md`](PROTOCOL_NARMA10.md), standing [`../PocketDish/LANE_A_NARMA10_STANDING.md`](../PocketDish/LANE_A_NARMA10_STANDING.md).
MG extra (**system PASS**; field wins living; not 0.260): [`PROTOCOL_MG.md`](PROTOCOL_MG.md), standing [`../PocketDish/LANE_A_MG_STANDING.md`](../PocketDish/LANE_A_MG_STANDING.md).
Lorenz extra (**system PASS**; named living-layer note on this \(u\); not 0.826; not CHARC): [`PROTOCOL_LORENZ.md`](PROTOCOL_LORENZ.md), standing [`../PocketDish/LANE_A_LORENZ_STANDING.md`](../PocketDish/LANE_A_LORENZ_STANDING.md).
Waveform extra (**system PASS**; named living-layer note on this \(u\); not 0.835; path step 2 closed): [`PROTOCOL_WAVEFORM.md`](PROTOCOL_WAVEFORM.md), standing [`../PocketDish/LANE_A_WAVEFORM_STANDING.md`](../PocketDish/LANE_A_WAVEFORM_STANDING.md).

Closed-benchmark report (HybridDish-style, with figures):
[`../PocketDish/LANE_A_BENCHMARK_REPORT.md`](../PocketDish/LANE_A_BENCHMARK_REPORT.md).
Regenerate figures: `python examples/LaneA_OccupiedMillimetre/plot_lane_a_benchmarks.py`.

Path memo §4 readiness (multi-AC / two-way / new Lane B identity):
[`../PocketDish/LANE_A_SECTION4_READINESS.md`](../PocketDish/LANE_A_SECTION4_READINESS.md).

Three-seed replicate (202/303 driven-only; waveform axis is Brier, not AUC 1.000).
One command, resume-safe:

```
python examples/LaneA_OccupiedMillimetre/run_seed_replicate.py
```

Metric lock (`LANE_A_METRIC_LOCK`; Jaeger locked; CHARC / Dambre IPC
out of scope on the frozen NARMA maps; not a paper):
[`PROTOCOL_METRIC.md`](PROTOCOL_METRIC.md),
`python examples/LaneA_OccupiedMillimetre/check_lane_a_metric_lock.py`,
standing [`../PocketDish/LANE_A_METRIC_LOCK_STANDING.md`](../PocketDish/LANE_A_METRIC_LOCK_STANDING.md).

New i.i.d. / constant rank drives (`LANE_A_KR_GR`; **SCORED**; T=840;
not a NARMA-map rewrite; not a parameter sweep):
[`PROTOCOL_KR_GR.md`](PROTOCOL_KR_GR.md),
`python examples/LaneA_OccupiedMillimetre/generate_lane_a_kr_gr_u.py`,
`ant lane-a-kr-gr`,
`python examples/LaneA_OccupiedMillimetre/check_lane_a_kr_gr.py`,
standing [`../PocketDish/LANE_A_KR_GR_STANDING.md`](../PocketDish/LANE_A_KR_GR_STANDING.md).

Protocol: [`PROTOCOL_SEED_REPLICATE.md`](PROTOCOL_SEED_REPLICATE.md).
Clean-chat prompt: [`../PocketDish/CHAT_LANE_A_SEED_REPLICATE_BOOTSTRAP.md`](../PocketDish/CHAT_LANE_A_SEED_REPLICATE_BOOTSTRAP.md).
