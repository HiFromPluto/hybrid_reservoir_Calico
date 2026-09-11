# BSimReservoirPlanStage7 mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanStage6`. Stage 6 is
not modified and remains **PASS**. Stage 5 remains **PASS**. Stage 4
remains **FAIL**; the 2x `Total_Deaths` fold is not rewritten. Receiver
`K=1.6`, `n=2`, `tau=15`, luminescence, acid source `2e11`, `K_MAX`,
clamp, `GROWTH_RATE`, field D/decay, and `FLOW_SPEED=0` are not retuned.
Glucose, Monod, nutrient AC, Danino, pH taxis, sign-flip repellent, and
reservoir_new inlet BCs are not imported. Stage 8 is not started.

This is **not** a copy of the plan's 5-AC table. That table is a different
dish (flow 8 µm/s, glucose AC5). This rebuild has four grounded sources.

## Four sources

| Source | Frozen position | Field |
|---|---|---|
| Attractant A | (150, 100, 5) | attractantField |
| AHL | (500, 250, 5) | ahlField |
| Acid | (850, 100, 5) | acidField |
| Attractant B | (850, 400, 5) | attractantField |

Attractant A and B are two sources on one attractant field. AC0/AC2 are
no longer silent architecture. Acid is no longer held at 0.5.

## Rank matrix

`Den_*` (20×10), `Receiver_R_*` (20×10), `Lum_Mean_*` (20×10),
`Input_Driven_Death_*` (4×2) = 608 columns. Not Births, Total_Deaths,
window AHL, occupancy, pH, or raw field voxels.

## Protocol

40 windows, same Stage 5/6 timing. Four independent input sequences with
pairwise |r| < 0.3. See `PROTOCOL.md`.

## Gates

All required. NARMA-10 is not a gate. Production outcome: **FAIL**
(silent seed 404 also has 39 significant components; rank is intrinsic
drift). Positions, rates, and decays were not retuned after SVD.
Stage 8 was not started. See `results/GATE_EVIDENCE.md`.
