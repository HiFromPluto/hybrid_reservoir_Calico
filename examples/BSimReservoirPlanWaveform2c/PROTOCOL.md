# Waveform2c protocol (frozen before ridge)

Track E1 aligned follow-on. Equal-histogram temporal-order
classification at a frozen alignment on the HybridDish / Narma10b /
Waveform claim dish.

Waveform2 remains **TASK_VOID** (ORDER_A ≡ ORDER_C orbit; RAW_U_8
`0.6058`). Waveform2b remains the disjoint-orbit / random-phase
negative (RAW_U_8 `0.5298`). Those packages are not rewritten.
Waveform `GATE_EVIDENCE.md` is not edited. Kinetics, layout, flow,
source rate, clamp, and mortality are not retuned. Claim dish stays
CENTER / FLOW=0.

This is **not** a cyclic-invariant task. Phase is identically 0.
Classes are `0=ORDER_UNI`, `1=ORDER_DOWN`, `2=ORDER_UP`. Do not name
them sine, square, or triangle.

C1 remains DEFER. E0.3 remains NO_STORY_MOVE.

## Why 2 and 2b stopped

Waveform2: ORDER_A and ORDER_C were one cyclic orbit.
Waveform2b: orbits disjoint, but random phase made linear RAW_U_8
near chance (`0.5298`). ORDER_UP is reverse(ORDER_DOWN), so power
spectra cannot separate them either.

Waveform2c tests temporal order at a frozen alignment with matched
histograms.

## Dish, timing, seeds (copied — do not change)

- `1000 × 500 × 10` µm, `dt=0.05`, `FLOW_SPEED=0`, `NO_FLUX`
- Field `50×25×1`, readout `20×10` + `4×2`
- `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL `1.28e8` at `(500, 250, 5)`. Acid `(300, 375, 5)` held 0.5 in
  analysis windows, off in warmup
- Attractants silent. Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, window `300 s`, pulse `75 s`, sample every `20 s`,
  **400 windows**
- Official 408. Brownian Den 20×10. Field-only driven AHL_uM 20×10
- Seeds `101 / 202 / 303`. 202/303 only if seed 101 is ALIVE and
  driven block AUC beats Brownian and MOMENTS
- Last sample `399;15;299.95`. CSV 400 / 6400 / 6400

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

## Frozen alphabet

```
LO   = 0.073223304703
MID  = 0.250000000000
HI   = 0.426776695297
MINV = 0.000000000000
MAXV = 0.500000000000
```

Shared multiset: `{ MINV, LO, LO, MID, MID, HI, HI, MAXV }`

Canonical representatives, all starting at MID, phase `p=0` on every
block (no circular shift):

```
ORDER_UNI  = MID, HI, MAXV, HI, MID, LO, MINV, LO
ORDER_DOWN = MID, MID, LO, LO, MINV, MAXV, HI, HI
ORDER_UP   = MID, MID, HI, HI, MAXV, MINV, LO, LO
```

ORDER_DOWN is Waveform2b DOWN circularly rotated by 3.
ORDER_UP is Waveform2b UP circularly rotated by 3.
No extra noise. Do not change templates after seeing AUC.

## Pre-hash gates (all before SHA-256)

1. Pairwise circular-shift sets have empty intersections.
2. The three canonical 8-vectors are not equal.
3. All three start with MID, so START_U is identical.
4. Block moments match to 1e-12.
5. Phase column is 0 for all 50 blocks.

3×3 orbit-intersection matrix:

|  | ORDER_UNI | ORDER_DOWN | ORDER_UP |
|---|---|---|---|
| ORDER_UNI | 8 | 0 | 0 |
| ORDER_DOWN | 0 | 8 | 0 |
| ORDER_UP | 0 | 0 | 8 |

**PREHASH_GATES_PASS.**

## Block moments (match to 1e-12)

| class | mean | variance | pop SD | power | min | max | range | n |
|---|---|---|---|---|---|---|---|---|
| ORDER_UNI | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |
| ORDER_DOWN | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |
| ORDER_UP | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 | 0.500000000000 | 8 |

## Class schedule

Class counts placed in each split, then shuffled with
`random.Random(20260819)`. No phase subseed. Phase identically 0.
Do not reuse Waveform2/2b class vectors.

Length-50 class vector:

```
0,0,2,0,1,1,1,2,2,0,1,2,0,1,1,0,0,2,1,1,2,2,1,0,1,1,2,0,0,0,1,0,2,2,1,0,0,0,2,0,1,0,2,2,1,2,2,1,0,1
```

Length-50 phase vector (all zeros):

```
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
```

SHA-256 of the 12-decimal `u` sequence:
`0979090c051c513dce14632e83ae9bf025c4b666e497dc879ae1db4bba30d465`.

SHA-256 of the ASCII class vector:
`465b3ad7132bd1e709f346fa4be2c946d28130ce398b47f41b9c2859baf70add`.

SHA-256 of the ASCII all-zero phase vector:
`ac6e5db2e0d7e63815992a3ca671b395e6bc533228054cffa27c17f35fa89cb8`.

Files: `input_ahl_waveform2c_400.txt`, `input_acid_held05_400.txt`,
`waveform2c_labels.csv` (`n;block;class;phase;u`),
`class_phase_schedule.txt`.

## u-only audit (mandatory stop, before BSim)

Block-level one-vs-rest ridge, same lambda rules:

1. MEAN_ONLY
2. POWER_ONLY
3. VARIANCE_ONLY
4. MOMENTS
5. START_U
6. RAW_U_8
7. CYCLIC_MATCH (diagnostic, not a gate): for each block, for each
   class, min L2 over 8 circular shifts of that class’s canonical
   template; three distances as features. With `p=0` this should also
   be high; it is the ceiling for a later cyclic task, not this gate.

PASS, frozen now:

- MEAN, POWER, VARIANCE, MOMENTS, START_U: test macro AUC within 0.15
  of 0.5. If any ≥ 0.70, scalar leak; STOP; no BSim.
- RAW_U_8: test macro AUC ≥ 0.90. If below, the three aligned
  8-vectors are not linearly distinct; STOP; do not rewrite after
  that score.

If PASS: smoke, then seed 101 driven + Brownian + silent. Seeds
202/303 only if seed 101 is ALIVE and driven block AUC beats Brownian
and MOMENTS. If FAIL: retain this draw, no Java.

### u-only audit result

**PASS.** Block-level test macro OVR AUC:

- MEAN_ONLY 0.5000
- POWER_ONLY 0.5000
- VARIANCE_ONLY 0.5000
- MOMENTS 0.5000
- RAW_U_8 1.0000
- START_U 0.5000
- CYCLIC_MATCH 1.0000 (diagnostic)

See `results/u_only_baselines.md`.

BSim is authorized to the smoke, then seed-101 scout.


## Ridge (biology / field / Brownian)

Preserve Waveform OVR behaviour. Primary: mean the 8 window OVR
scores inside each block, then macro OVR AUC and accuracy on the 16
test blocks. Occupancy DEAD if `mean_R < 0.05` or `|r(mean_R,u)| < 0.5`.
Do not raise rate or K if DEAD.

Biology need not beat RAW_U_8. Biology must beat Brownian, silent,
and MOMENTS at block level. Field-only is a carrier diagnostic. If
field ≥ driven, the aligned order is already in the plume. Do not add
ACs. Do not rewrite Waveform1 FAIL. This is not evidence of
shift-invariant classification.

## Checker

```
python examples/BSimReservoirPlanWaveform2c/check_waveform2c.py --u-only
python examples/BSimReservoirPlanWaveform2c/check_waveform2c.py
```
