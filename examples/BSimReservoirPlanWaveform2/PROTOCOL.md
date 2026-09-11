# Waveform2 protocol (frozen before ridge)

Track E1. Equal-histogram temporal-order classification on the frozen
HybridDish / Narma10b / Waveform claim dish. This package does **not**
convert the frozen Waveform Overall FAIL into PASS. Waveform
`GATE_EVIDENCE.md` is not edited. Kinetics, layout, flow, source rate,
clamp, mortality, K, n, tau_R, and tau_L are not retuned.

Choice (do not reopen): **strict temporal-order classification**.
Classes are `0=ORDER_A`, `1=ORDER_B`, `2=ORDER_C`. Do not name them
sine, square, or triangle. The plan’s 8-point sine/square/triangle
family is forbidden here: those match mean/min/max but not
variance/power. Waveform2 uses one shared amplitude multiset so every
scalar block moment is identical across classes. Class identity is
temporal order plus phase, nothing else.

C1 remains DEFER. E0.3 remains NO_STORY_MOVE. This development set is
retained in full. Do not start confirmation on a new class/phase draw
until this package is complete.

## Dish, timing, seeds (copied — do not change)

- `1000 × 500 × 10` µm, `dt=0.05`, `FLOW_SPEED=0`, `NO_FLUX`
  (`setSolid(true,true,true)`)
- Field grid `50×25×1`, readout `20×10` and `4×2`
- Receiver `K=1.6`, `n=2`, `tau_R=15`
- `ALPHA_LUX=DELTA_LUX=1/1500` (`tau_L=1500`)
- AHL source `1.28e8` at `(500, 250, 5)`. Acid source `2e11` at
  `(300, 375, 5)` held at `0.5` every analysis window, off during warmup
- Attractant AC0/AC2 stay at `u=0`
- `K_MAX=0.002`, `GROWTH_RATE` compile-time `4π/1800`
- Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, windows `300 s`, pulses `75 s`, whole-window
  sampling every `20 s`, **400 windows**
- One-way `addQuantity` ACs. No vesicles, no glucose, no Danino, no
  `setGoal(glucose)`, no extra ACs
- Seeds `101 / 202 / 303`. Seed 202/303 run only if seed 101 scout
  occupancy is ALIVE and driven block AUC beats Brownian and MOMENTS

Official 408: `Receiver_R_*` 20×10, `Lum_Mean_*` 20×10,
`Input_Driven_Death_*` 4×2. Brownian: `Den_*` 20×10 only. Field-only:
driven `AHL_uM_*` 20×10. Not in the biology ridge: voxel AHL, Den, pH,
Att, raw `u`.

Silent cannot be reused from 200-window runs. Brownian cannot be reused
(duration changed).

CSV completeness: 400 / 6400 / 6400 windows/samples. Last-sample
convention copied from Narma10b: `399;15;299.95`. Freeze this triple
from the smoke/production header if the exporter prints it.

## Frozen task geometry

8 windows per block = one period. 50 blocks = 400 windows.
Chance accuracy 1/3. Chance macro OVR AUC 0.5.

| Split | Blocks | Windows |
|---|---|---|
| Washout | 0..7 | 0..63 |
| Train | 8..33 | 64..271 |
| Test | 34..49 | 272..399 |

Inner validation = last 6 train blocks (28..33), windows 224..271.
Lambda and standardization never see test. Inner standardization uses
inner-train blocks only (blocks 8..27, windows 64..223). Final fit uses
all 26 train blocks.

Class counts, frozen:

- washout 8: 3, 3, 2
- train 26: 9, 9, 8
- test 16: 6, 5, 5

That is 16 independent test blocks, not 10. Block boundaries are never
broken.

## Frozen amplitude alphabet

Let `s2 = 0.25 * sqrt(2) / 2`. Exactly, to 12 decimals after rounding:

```
LO   = 0.073223304703
MID  = 0.250000000000
HI   = 0.426776695297
MINV = 0.000000000000
MAXV = 0.500000000000
```

Shared 8-value multiset (order does not matter here):

```
{ MINV, LO, LO, MID, MID, HI, HI, MAXV }
```

Phase-0 templates, 12-decimal `u` clipped to `[0, 0.5]`:

```
ORDER_A: MID, HI, MAXV, HI, MID, LO, MINV, LO
ORDER_B: MAXV, HI, HI, MID, MID, LO, LO, MINV
ORDER_C: MINV, LO, MID, HI, MAXV, HI, MID, LO
```

A block with class `c` and integer phase `p` in `{0..7}` is the circular
shift of that template by `p` (index `(w+p) mod 8`). Full-period
histogram is invariant to `p`. No extra noise. Do not change templates
after seeing AUC.

Generator diagnostic (not a retune): ORDER_A and ORDER_C occupy the
same cyclic orbit (`ORDER_C` is ORDER_A shifted by 6). ORDER_B is a
distinct orbit. This is recorded before ridge. It is not a license to
edit the frozen alphabet.

## Phase / class schedule

Class counts were placed in each split, then shuffled with
`random.Random(20260817)`. Phases were assigned with
`random.Random(20260817 + 1 + subseed)` to meet:

1. Class counts match the table inside each split.
2. Within each split, phases are as balanced as integer counts allow.
   No class in train or test uses a single phase for all of its blocks.
3. Starting amplitude alone must not separate classes in train or test.

Chosen `phase_subseed=0` (first admissible draw; no hand-edit).

Length-50 class vector:

```
1,1,2,0,1,2,0,0,2,1,0,2,0,0,0,1,1,0,1,0,1,2,1,0,1,1,1,2,0,2,2,2,0,2,1,1,0,1,0,1,2,2,2,2,1,0,0,2,0,0
```

Length-50 phase vector:

```
0,3,2,2,7,4,7,6,1,3,5,3,7,0,2,1,6,3,0,4,7,2,5,1,4,1,2,4,6,6,7,5,1,0,0,4,4,5,1,6,5,0,3,7,7,7,2,6,6,0
```

SHA-256 of the 12-decimal `u` sequence (Narma10b / Waveform recipe):
`8f47d9f30d888b0e223e23b8804b8d875187b800f283a166d15ab47a44663d31`.

SHA-256 of the ASCII class vector:
`83aab1be5b5499566c802d3be58f918c572fb77b2a50baf6b9c6de911ca9499d`.

SHA-256 of the ASCII phase vector:
`5cf285dc94883d4b892fa00b22dd751a4f4d5c3ea672787f2e0d02b57f98c73a`.

Files: `input_ahl_waveform2_400.txt`, `input_acid_held05_400.txt`,
`waveform2_labels.csv` (`n;block;class;phase;u`),
`class_phase_schedule.txt`.

## Block moments by class (must match to 1e-12)

| class | mean | variance | pop SD | power | min | max | range | n |
|---|---|---|---|---|---|---|---|---|
| ORDER_A | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |
| ORDER_B | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |
| ORDER_C | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |

## u-only audit (mandatory stop, before BSim)

Block-level one-vs-rest ridge on:

1. MEAN_ONLY
2. POWER_ONLY
3. VARIANCE_ONLY
4. MOMENTS: mean, variance, power, min, max, range
5. RAW_U_8: the eight ordered `u` values in the block
6. START_U: first window `u` only

Same block splits and lambda rules.

Pass criteria, frozen now:

- MEAN, POWER, VARIANCE, MOMENTS, START_U: test macro AUC within 0.15
  of 0.5. If any of these is ≥ 0.70, the templates leaked a scalar cue.
  STOP. Do not run BSim.
- RAW_U_8: test macro AUC must be ≥ 0.90. This is the input-separability
  ceiling, not reservoir skill. If it is near chance, the three orders
  are not distinct and the task is void.

Write `results/u_only_baselines.md` and `results/u_only_baselines.csv`.
Do not proceed to BSim until this audit is recorded below as PASS.

### u-only audit result

**FAIL.** Block-level test macro OVR AUC:

- MEAN_ONLY 0.5000
- POWER_ONLY 0.5000
- VARIANCE_ONLY 0.5000
- MOMENTS 0.5000
- RAW_U_8 0.6058
- START_U 0.4919

See `results/u_only_baselines.md`.

STOP. RAW_U_8 is below 0.90. ORDER_A and ORDER_C are cyclic
shifts of one another, so the three orders are not distinct.
The task is void. Do not run BSim. Do not retune K, n, tau,
source, clamp, mortality, flow, or layout. Do not start a new
class/phase draw until this development set is retained.


## Ridge (biology / field / Brownian)

Preserve Waveform OVR behaviour:

- three binary heads
- lambda grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`
- pick lambda by highest macro OVR validation AUC; ties → larger lambda
- unregularized intercept
- train-only standardization
- no sklearn required

Primary metrics (block-level):

- mean the 8 window OVR scores inside each block, then macro OVR AUC
  and accuracy on the 16 test blocks
- per-class AUC and confusion

Secondary:

- window-level macro AUC / accuracy on the 128 test windows
- block-mean 408-D feature ridge (26 train blocks; descriptive only)

| Arm | Features |
|---|---|
| Driven / silent | 408 official |
| Brownian | 200 Den |
| Field-only | 200 AHL_uM from driven |

Also compute unmasked and occupancy-masked kinetic surrogates with the
frozen R/L ODEs from the zero-simulation audit, scored at block level.

## Occupancy

Report driven `mean_R`, `r(mean_R,u)`, `mean_AHL`, `mean_pop`.
DEAD if `mean_R < 0.05` or `|r(mean_R,u)| < 0.5`.
Do not raise rate or K if DEAD. Occupancy-dead rows are not Waveform2
skill.

## BSim schedule (internal stops)

Smoke first: HybridDish smoke timing, 400-window config not required.
Print dish, source, FLOW=0, finite non-negative AHL. Not evidence.

Scout, seed 101 only: driven + Brownian + silent (3 production runs).

If scout occupancy is DEAD, or driven block AUC does not beat Brownian
and MOMENTS, stop. Do not run 202/303. Do not retune.

If scout is ALIVE and driven beats Brownian and MOMENTS on block AUC,
run seeds 202 and 303 (6 more runs). No best-seed selection.

Total if fully authorized: 9 production BSim runs.

## Gates (Waveform2 only; do not touch old Overall)

Report, do not retune:

1. Driven block macro AUC > Brownian, every completed seed.
2. Driven block macro AUC > MOMENTS and MEAN/POWER/VARIANCE.
3. Silent block macro AUC not clearly above 0.5.
4. Field-only block macro AUC as a carrier diagnostic. If field ≥
   driven, say the plume already contains the order. That is not a
   license to add ACs or to rewrite the old waveform FAIL.
5. RAW_U_8 remains the input ceiling; biology is not required to beat
   it. If biology is far below RAW_U_8 and also below field, the living
   readout is not the temporal-order substrate.

0.73 / 0.82 from Waveform1 stay historical. Do not mix them into
Waveform2 means. The 0.93 NARMA story is unchanged.

## Checker

```
python examples/BSimReservoirPlanWaveform2/check_waveform2.py --u-only
python examples/BSimReservoirPlanWaveform2/check_waveform2.py
```
