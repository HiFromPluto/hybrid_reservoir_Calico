# BSimReservoirPlanStage3B mechanism manifest

## Scope and provenance

This package is the spatial integration of the predeclared Stage 3B receiver
benchmark in `benchmark_receiver.py`. The benchmark passed before integration
with `K=1.6 uM`, Hill exponent `n=2`, `tau=15 s`, and
`t95=-ln(0.05)*tau=44.936 s`. Those values are copied exactly and are not
configurable or retuned here.

`K=1.6 uM` is trace-matched to the archived Stage 3 AHL field traces. It is not
claimed to be a literature EC50. The archived
`examples/BSimReservoirPlanStage3` package and its Danino constants were not
altered; Stage 3 remains preserved negative evidence.

## Exact receiver

Every bacterium carries one scalar state `R`. Once per simulation tick, its
local AHL field concentration is converted from molecules/um3 to uM using
`C = field_concentration / 602.2`. The target and exact constant-target update
over the tick are:

`target = C^2 / (1.6^2 + C^2)`

`R = target + (R - target) * exp(-dt / 15)`

`Mean_q` is a compatibility label for mean `R`, and
`Fraction_q_gt_0_5` is the fraction of cell observations with `R>0.5`.
Receiver-specific outputs are otherwise named `Receiver_R`; they do not denote
Danino LuxR-AHL or LuxI.

There is no bacterial AHL exchange, source, sink, intracellular AHL, Danino
ODE, or Danino kinetic constant in Stage 3B. The direct center-AC1 source and
AHL field decay remain unchanged.

## Preserved architecture and exclusions

Stage 3B retains the Stage 3 spatial simulation, Stage 2 timing, 50x25x1
chemical grid, readout grids, population growth/removal mechanisms, zero flow,
18,000 s warmup, matched 0.5 AHL input for 75 s of every 300 s warmup cycle,
40 production windows of 300 s, 75 s measurement pulses, whole-window
sampling every 20 s plus the final tick, source strength 128,000,000
molecules/s, silent AC0/AC2 inputs, direct AC1 source, and field decay
`0.0033 /s`.

The three production configs use the three balanced permutations, seeds
11/22/33, and independent output directories. Each complete run must produce
40 summary rows and 640 rows in each sample and voxel CSV.

There is no nutrient model, ArtificialCell, death extension, toxicity,
luminescence, or Stage 4 mechanism in this package.

## Evidence gate

`analyze_stage3b.py` requires all of the following:

1. Every replicate has `|r(AHL_Input, Mean_q)| > 0.5`.
2. Mean mid-input `Fraction_q_gt_0_5` is in `[0.2, 0.8]`.
3. Pooled response is ordered low < mid < high.
4. Summary/sample/voxel CSVs are rectangular and have 40/640/640 rows.
5. The explicitly recorded response-timescale gate has `t95 <= 300 s`.

A silent control is optional diagnostic input and is not required by the
Stage 3B full-rerun gate.

## Final outcome

The deterministic benchmark passed first with correlations
`0.8183/0.8372/0.8356`, mid occupancy `0.3038`, ordered response
`0.2392 < 0.4212 < 0.5561`, and `t95=44.936 s`.

The unchanged candidate then passed the spatial three-permutation calibration rerun:
correlations `0.7870/0.8177/0.8171`, mid occupancy `0.3374`, and ordered
response `0.2536 < 0.4327 < 0.5605`. All CSV validation and timescale gates
pass. Because those permutations also supplied the benchmark traces, this is a
conditional/calibration pass rather than independent closure evidence. See
`results/GATE_EVIDENCE.md` and
`results/stage3b_gate_evidence.json`.

`t95=44.936 s` is analytic for the exact first-order equation, not an
empirically measured response time. Biological claims are limited to a
phenomenological, trace-matched receiver.

The frozen receiver then passed three untouched balanced holdouts with fresh
seeds `44/55/66`: correlations `0.8321/0.8223/0.8029`, mid occupancy `0.3450`,
and ordered response `0.2364 < 0.4355 < 0.5707`. The seed-77 source-off
control had exactly zero mean response and occupancy. See
`results/stage3b_holdout_gate_evidence.json`.

Stage 3B is closed as PASS on independent holdouts. Biological claims remain
limited to a phenomenological, trace-matched receiver. Stage 4 remains blocked.
