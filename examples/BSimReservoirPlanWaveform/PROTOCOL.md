# Waveform classification protocol (frozen before AUC)

This is waveform classification on the frozen Plan Stage 6 hybrid dish.
Mechanisms are not retuned. Stage 6 remains the NARMA-10 provenance
checkpoint (`examples/BSimReservoirPlanStage6/`, Overall: PASS). BenchA
Mackey–Glass remains PASS. BenchA Lorenz remains FAIL. Track B remains
FAIL. Stage 9 product-bit FAIL still stands. Stage 6, BenchA, BenchA2,
Track B, and Stage 9 Java are not edited.

Three classes: sine, square, triangle. Chance accuracy = 1/3. Chance
macro one-vs-rest AUC = 0.5. Track B was binary; 0.5 is not the
accuracy floor here. Valid only with Brownian, silent, and a Stage
9-style field-only AHL baseline.

If field-only ranks the waveform, Overall is FAIL. That means the
plume is the waveform. Do not then mix in Att/Rep, raise source rates,
or make the three templates more similar to hide the field. Do not
switch the gate to accuracy after seeing numbers.

## Dish, timing, seeds (copied from Stage 6 / BenchA — do not change)

- `dt=0.05`, warmup `18000 s`, windows `300 s`, pulses `75 s`,
  whole-window sampling every `20 s`, 200 windows
- Field grid `50×25×1`, readout `20×10` and `4×2`
- Receiver `K=1.6`, `n=2`, `tau=15`
- `ALPHA_LUX=DELTA_LUX=1/1500`
- AHL source `1.28e8`, acid source `2e11`, `K_MAX=0.002`,
  `GROWTH_RATE` compile-time `4π/1800`
- Clamp `K=2000`, `INITIAL_POP=1800`
- `FLOW_SPEED=0`
- AHL AC at `(500, 250, 5)`. Acid at `(300, 375, 5)` held at `0.5`
  every analysis window, off during warmup
- Attractant AC0/AC2 stay at `u=0`
- Seeds `101 / 202 / 303` shared across Brownian, silent, and driven
- One-way `addQuantity` ACs. No vesicles, no glucose, no Danino, no
  `setGoal(glucose)`, no extra ACs

`num.windows=200` is unchanged from Stage 6.

## Arms

1. **Driven.** AHL carries the waveform sequence; acid held at 0.5.
2. **Brownian.** Same pulses; ridge `Den` only. Particles do not sense.
   Brownian is rerun because the AHL file is new.
3. **Silent.** Reuse Stage 6 silent CSVs (sources off). Do not rerun
   silent.

## Frozen task (declared before AUC)

40 blocks × 5 windows = 200. Block `p` uses windows `5p .. 5p+4`.
Label `y[n]` = class of that block. Classes: `0=sine`, `1=square`,
`2=triangle`.

Splits on blocks (same window counts as Stage 6):

| Split | Blocks | Windows | Count |
|---|---|---|---|
| Washout | 0 .. 7 | 0 .. 39 | 40 |
| Train | 8 .. 29 | 40 .. 149 | 110 |
| Test | 30 .. 39 | 150 .. 199 | 50 |

Class counts forced then shuffled with `random.Random(20260822)`
inside each split:

- washout 3 sine, 3 square, 2 triangle
- train 8 sine, 7 square, 7 triangle
- test 3 sine, 3 square, 4 triangle

Length-40 class vector (comma-separated 0/1/2):

```
0,1,0,1,2,1,2,0,2,0,0,1,0,1,1,1,2,1,2,2,0,1,0,0,2,2,0,0,2,1,1,2,2,1,2,0,2,0,0,1
```

SHA-256 of that ASCII string:
`c2155819f0fb512b770db53d44ef89e2577de3025f81141166278a8903446c29`.

Deterministic templates, `w = 0 .. 4` within a block, `u` clipped to
`[0, 0.5]`:

```
sine:     u = 0.25 + 0.25 * sin(2 * pi * w / 5)
square:   u = 0.45 if w < 3 else 0.05
triangle: u = [0.05, 0.25, 0.45, 0.25, 0.05][w]
```

No extra noise. Do not change templates after seeing AUC.

Files: `input_ahl_waveform200.txt`, `input_acid_held05_200.txt`
(Stage 6 copy), `waveform_labels.csv` with columns `n;block;y;u`.

SHA-256 of the 12-decimal `u` sequence (Stage 6 recipe):
`c97a7dcb92c663892d16d50bfe342003b021b173211d9246b3c86c1e9734e0f8`.

## Ridge

Washout drop; 160 rows. Lambda on rows 88..109 only (windows
128..149). Test never used for lambda or standardization. Grid
`{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`. Bias; intercept not
regularized. Training mean/std. Independent lambda per arm×seed.

Multiclass: one-vs-rest ridge (three binary heads, same lambda grid
picked by highest macro OVR validation AUC; ties larger lambda).

Primary test metric: macro one-vs-rest AUC (mean of three
Mann–Whitney AUCs). sklearn is not required. Accuracy (argmax of
three scores) is printed and compared to 1/3. NOT a gate.

Diagnostic: block-pooled (mean of 5 window scores per test block, 10
blocks). Print; do not replace window AUC.

### Features

| Arm | Features | Count |
|---|---|---|
| Driven / silent | `Receiver_R_*` 20×10, `Lum_Mean_*` 20×10, `Input_Driven_Death_*` 4×2 | 408 |
| Brownian | `Den_*` 20×10 | 200 |
| Field-only (required) | driven voxel `AHL_uM_*` | 200 |

Not in the biology ridge: voxel AHL, Den (biology), pH, Att, raw `u`.

## Gates (all required for Overall PASS)

1. Driven macro OVR AUC > Brownian, non-overlapping mean±s.e., every
   seed.
2. Driven macro OVR AUC > field-only AHL, non-overlapping mean±s.e.,
   every seed. If this fails: Overall FAIL. Do not add ACs. Do not
   switch to accuracy vs 0.334.
3. Silent macro OVR AUC not clearly above 0.5 (mean±s.e. overlaps
   0.5). Silent not ≈ driven if driven has skill.
4. CSV 200/3200/3200; last sample `199;15;299.95`.
5. Stage 6 NARMA-10 PASS, BenchA MG PASS, BenchA Lorenz FAIL, Track B
   FAIL unchanged. Java of those packages untouched.

If gate 2 fails, write FAIL and stop.

## Runs

6 BSim runs: driven + brownian × 101/202/303. Copy three silent dirs.

    results/waveform_{driven,brownian,silent}_seedXXX

Checker: `python examples/BSimReservoirPlanWaveform/check_waveform.py`
Evidence: `examples/BSimReservoirPlanWaveform/results/GATE_EVIDENCE.md`

Print per-class OVR AUC, confusion counts, accuracy vs 1/3, and
field-only.

Do not start Track B or a wider dish. Do not start 5-AC Track B,
Lorenz, vesicles, glucose, or Danino.
