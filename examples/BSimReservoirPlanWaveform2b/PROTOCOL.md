# Waveform2b protocol (frozen before ridge)

Track E1 follow-on. Distinct cyclic-orbit temporal-order classification
on the frozen HybridDish / Narma10b / Waveform claim dish.

Waveform2 remains **TASK_VOID** (`ORDER_A` ≡ `ORDER_C` orbit; RAW_U_8
`0.6058`). That package is not rewritten. Waveform `GATE_EVIDENCE.md`
is not edited. Kinetics, layout, flow, source rate, clamp, mortality,
K, n, tau_R, and tau_L are not retuned. Claim dish stays CENTER /
FLOW=0.

Classes are `0=ORDER_UNI`, `1=ORDER_DOWN`, `2=ORDER_UP`. Do not name
them sine, square, or triangle. Shared amplitude multiset; class
identity is temporal order plus phase.

C1 remains DEFER. E0.3 remains NO_STORY_MOVE.

## Why Waveform2 was void

```
ORDER_A = MID, HI, MAXV, HI, MID, LO, MINV, LO
ORDER_C = MINV, LO, MID, HI, MAXV, HI, MID, LO
```

ORDER_C is ORDER_A circularly shifted by 6. ORDER_A is also closed
under reversal (reverse = shift 5). Reverse-A is not used as a third
class here.

## Dish, timing, seeds (copied from Waveform2 — do not change)

- `1000 × 500 × 10` µm, `dt=0.05`, `FLOW_SPEED=0`, `NO_FLUX`
- Field grid `50×25×1`, readout `20×10` and `4×2`
- Receiver `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL source `1.28e8` at `(500, 250, 5)`. Acid `(300, 375, 5)` held at
  `0.5` in analysis windows, off in warmup
- Attractants silent. Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, windows `300 s`, pulses `75 s`, sample every `20 s`,
  **400 windows**
- Official 408: `Receiver_R`, `Lum_Mean`, `Input_Driven_Death`
- Brownian: `Den` 20×10 only. Field-only: driven `AHL_uM` 20×10
- Seeds `101 / 202 / 303`. 202/303 only if seed 101 scout is ALIVE and
  driven block AUC beats Brownian and MOMENTS
- Last-sample convention: `399;15;299.95`. CSV 400 / 6400 / 6400

## Frozen task geometry

8 windows per block. 50 blocks = 400 windows. Chance accuracy 1/3.
Chance macro OVR AUC 0.5.

| Split | Blocks | Windows |
|---|---|---|
| Washout | 0..7 | 0..63 |
| Train | 8..33 | 64..271 |
| Test | 34..49 | 272..399 |

Inner validation = last 6 train blocks (28..33), windows 224..271.
Lambda and standardization never see test. Inner standardization uses
inner-train blocks only. Final fit uses all 26 train blocks.

Class counts: washout 3,3,2 ; train 9,9,8 ; test 6,5,5.
Sixteen independent test blocks. Block boundaries never broken.

## Frozen amplitude alphabet

```
LO   = 0.073223304703
MID  = 0.250000000000
HI   = 0.426776695297
MINV = 0.000000000000
MAXV = 0.500000000000
```

Shared multiset: `{ MINV, LO, LO, MID, MID, HI, HI, MAXV }`

Phase-0 templates:

```
ORDER_UNI  = MID, HI, MAXV, HI, MID, LO, MINV, LO
ORDER_DOWN = MAXV, HI, HI, MID, MID, LO, LO, MINV
ORDER_UP   = MINV, LO, LO, MID, MID, HI, HI, MAXV
```

A block is `template[(w+p) mod 8]`, `p` in `{0..7}`. No extra noise.
Do not change templates after seeing AUC.

## Mandatory pre-hash orbit gate

Computed before writing inputs or hashes. Set of 8 circular shifts of
each phase-0 template. Pairwise intersections empty. No template is a
circular shift of another. ORDER_UNI reverse lies in the UNI orbit and
is not used as another class.

3×3 orbit-intersection matrix (counts of shared circular shifts):

|  | ORDER_UNI | ORDER_DOWN | ORDER_UP |
|---|---|---|---|
| ORDER_UNI | 8 | 0 | 0 |
| ORDER_DOWN | 0 | 8 | 0 |
| ORDER_UP | 0 | 0 | 8 |

**ORBIT_GATE_PASS.**

## Block moments by class (must match to 1e-12)

| class | mean | variance | pop SD | power | min | max | range | n |
|---|---|---|---|---|---|---|---|---|
| ORDER_UNI | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |
| ORDER_DOWN | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |
| ORDER_UP | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |

## Phase / class schedule

Class counts placed in each split, then shuffled with
`random.Random(20260818)`. Phases assigned with
`random.Random(20260818 + 1 + subseed)`. Chosen `phase_subseed=0`.

Length-50 class vector:

```
0,2,2,0,1,1,0,1,0,1,1,0,1,2,2,0,1,0,1,0,2,0,1,0,2,0,1,0,2,1,2,1,2,2,0,1,0,0,1,2,0,1,2,2,2,0,0,1,1,2
```

Length-50 phase vector:

```
2,5,3,0,2,4,6,3,2,4,2,0,7,7,0,5,0,7,6,3,1,1,1,7,4,4,5,6,3,3,6,6,2,5,6,6,7,3,4,3,4,1,7,4,6,5,2,2,3,2
```

SHA-256 of the 12-decimal `u` sequence:
`339e8eb8c6a89510bfe1ef0f38938cc4407c6470e526e5fcf7c2ec734b9d5cdf`.

SHA-256 of the ASCII class vector:
`26b44c1c2ee3932becc7fdebfc72179a2872bb87baca975624355c0f315772e6`.

SHA-256 of the ASCII phase vector:
`01f972847d1dad3c3b6591ff7fab7537b8f4b020da40f6f35376dfa02e97b7bc`.

Files: `input_ahl_waveform2b_400.txt`, `input_acid_held05_400.txt`,
`waveform2b_labels.csv` (`n;block;class;phase;u`),
`class_phase_schedule.txt`.

## u-only audit (mandatory stop, before BSim)

Block-level one-vs-rest ridge:

1. MEAN_ONLY
2. POWER_ONLY
3. VARIANCE_ONLY
4. MOMENTS: mean, variance, power, min, max, range
5. RAW_U_8
6. START_U

Same block splits and lambda rules as Waveform2.

- MEAN, POWER, VARIANCE, MOMENTS, START_U: test macro AUC within 0.15
  of 0.5. If any ≥ 0.70, scalar leak; STOP; do not BSim.
- RAW_U_8: test macro AUC ≥ 0.90. If below, orders are not distinct;
  STOP; do not rewrite templates after that score.

If PASS: smoke, then seed 101 driven/Brownian/silent. Seeds 202/303
only if seed 101 is ALIVE and beats Brownian and MOMENTS on block AUC.
If FAIL: retain this draw, write WAVEFORM2B_SCOUT.md, no Java.

### u-only audit result

**FAIL.** Block-level test macro OVR AUC:

- MEAN_ONLY 0.5000
- POWER_ONLY 0.5000
- VARIANCE_ONLY 0.5000
- MOMENTS 0.5000
- RAW_U_8 0.5298
- START_U 0.5571

See `results/u_only_baselines.md`.

STOP. RAW_U_8 is below 0.90. The three orders are not
distinct. The task is void. Do not run BSim. Do not rewrite
templates after that score. Do not retune K, n, tau, source,
clamp, mortality, flow, or layout.


## Ridge (biology / field / Brownian)

Preserve Waveform OVR behaviour: three binary heads; lambda grid
`{1e-6,1e-4,1e-2,1,1e2,1e4,1e6}`; highest macro OVR validation AUC;
ties → larger lambda; unregularized intercept; train-only
standardization; no sklearn.

Primary: mean the 8 window OVR scores inside each block, then macro
OVR AUC and accuracy on the 16 test blocks. Secondary: window-level
on 128 test windows; block-mean 408-D feature ridge (descriptive).

Occupancy: driven `mean_R`, `r(mean_R,u)`, `mean_AHL`, `mean_pop`.
DEAD if `mean_R < 0.05` or `|r(mean_R,u)| < 0.5`. Do not raise rate
or K if DEAD.

## Gates (Waveform2b only; do not touch old Overall)

1. Driven block macro AUC > Brownian, every completed seed.
2. Driven block macro AUC > MOMENTS and MEAN/POWER/VARIANCE.
3. Silent block macro AUC not clearly above 0.5.
4. Field-only as a carrier diagnostic.
5. RAW_U_8 remains the input ceiling; biology is not required to beat it.

0.73 / 0.82 from Waveform1 stay historical. 0.93 NARMA unchanged.

## Checker

```
python examples/BSimReservoirPlanWaveform2b/check_waveform2b.py --u-only
python examples/BSimReservoirPlanWaveform2b/check_waveform2b.py
```
