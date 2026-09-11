# Track B protocol (frozen before AUC)

This is Track B on a copy of the frozen Plan Stage 6 hybrid dish: five
one-way ACs and a synthetic 5-channel patient classification task.
Mechanisms are not retuned. Stage 6 remains the NARMA-10 provenance
checkpoint (`examples/BSimReservoirPlanStage6/`, Overall: PASS). Track A
Mackey–Glass remains PASS. Track A Lorenz remains FAIL. Stage 9
product-bit FAIL (driven AUC ≈ field AHL) still stands. Stage 7
four-source SVD FAIL still stands. Stage 6, BenchA, BenchA2, Stage
3B/4/5/7/9, D1g, and M1–M3 Java are not edited.

Track B tests H5: five mixed chemical inputs vs a scalar AHL mix, with a
field-only gate copied from Stage 9. If the five plumes classify the
patient, Overall is FAIL. Do not then add ACs, stretch the domain, or
switch the gate to accuracy.

## Frozen 5-AC layout (declared before BSim; do not move after AUC)

| Channel | Chemical | Position | Rate × u |
|---|---|---|---|
| AHL | ahlField | (500, 250, 5) | 1.28e8 |
| acid | acidField | (300, 375, 5) | 2e11 |
| attA | attractantField | (150, 100, 5) | 1e6 |
| attB | attractantField | (850, 400, 5) | 1e6 |
| rep | repellentField | (150, 400, 5) | 1e6 |

Do not inject at Stage 6’s silent sites `(250, 250, 5)` or `(750, 250, 5)`.
attA and attB share one attractant field (two plumes, one chemical). That
is intended. Do not split into two attractant PDEs. `PROD_RATE=1e6` is
engineering; do not raise it after seeing AUC.

`setGoal(attractantField)` only. No `setGoal(repellent)`, no pH taxis, no
glucose. One-way `addQuantity`. No vesicles.

## Dish, timing, seeds (copied from Stage 6 — do not change)

- `dt=0.05`, warmup `18000 s`, windows `300 s`, pulses `75 s`,
  whole-window sampling every `20 s`, 200 windows
- Field grid `50×25×1`, readout `20×10` and `4×2`, domain `1000×500×10`
- Receiver `K=1.6`, `n=2`, `tau=15`
- `ALPHA_LUX=DELTA_LUX=1/1500`
- AHL source `1.28e8`, acid source `2e11`, `K_MAX=0.002`,
  `GROWTH_RATE` compile-time `4π/1800`
- Clamp `K=2000`, `INITIAL_POP=1800`
- `FLOW_SPEED=0`
- Seeds `101 / 202 / 303`

Warmup: AHL `0.5` at the AHL site only, acid off, att/rep off (same as
Stage 6).

## I/O change (not kinetics)

Stage 6 voxels have `AHL_uM`, `pH`, `Den`, `R`, `L`, deaths — no Att/Rep
maps. Track B adds attractant and repellent voxel exports on the 20×10
grid:

- `Att_conc_*` = `attractantField.getConc` (molecules/µm³)
- `Rep_conc_*` = `repellentField.getConc` (molecules/µm³)

`AHL_uM` and `pH` stay as in Stage 6. Field-only needs these maps. Do
not add them to the biology ridge.

## Arms

1. **Driven 5-channel.** All five sites carry their biomarker files.
   Acid is a free input (NOT held at 0.5).
2. **Single-site.** Only the AHL site fires. `u_AHL` = mean of the five
   channels that window, clipped to `[0, 0.5]`. attA=attB=rep=0. Acid
   held at 0.5 (Stage 6 option a). Same patient labels. This is the H5
   control.
3. **Brownian.** Same 5-channel pulses as driven (identical dish). Ridge
   is `Den` only. Particles do not sense.
4. **Silent.** All source rates 0, `warmup.ahl.input=0`. Reuse Stage 6
   silent CSVs
   `examples/BSimReservoirPlanStage6/results/stage6_silent_seed{101,202,303}/`.
   Do not rerun silent.

## Frozen patient task (declared before AUC)

40 patients × 5 windows = 200 windows, patient-major time order.
Patient `p` uses windows `5p .. 5p+4`. Label `y[n]` = class of that
patient (constant on the block).

Splits (patient boundaries, not a random window cut):

| Split | Patients | Windows | Count |
|---|---|---|---|
| Washout | 0 .. 7 | 0 .. 39 | 40 |
| Train | 8 .. 29 | 40 .. 149 | 110 |
| Test | 30 .. 39 | 150 .. 199 | 50 |

Class balance is forced inside each split, then shuffled with
`random.Random(20260820)`:

- washout 4 pos + 4 neg
- train 11 pos + 11 neg
- test 5 pos + 5 neg

Length-40 class vector (comma-separated 0/1):

```
1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0
```

SHA-256 of that ASCII string:
`16de0e28bc19c9f94bd468afb8a3418b29b747ba544e8769855a3499eee2ca48`.

Channel traces, given class, `random.Random(20260821)` — a different
seed from the label shuffle:

- HIGH ~ Uniform`[0.30, 0.50]`
- LOW  ~ Uniform`[0.00, 0.20]`
- DIST ~ Uniform`[0.00, 0.50]`
- Class 1: every window attA=HIGH, acid=HIGH
- Class 0: every window, coin flip: (attA=HIGH, acid=LOW) or (attA=LOW, acid=HIGH)
- AHL, attB, rep: DIST every window, independent of class (distractors)

Frozen draw order per window `n = 0 .. 199`: DIST AHL, DIST attB, DIST
rep, then class-conditional attA/acid (class 1: HIGH attA then HIGH
acid; class 0: coin, then the HIGH/LOW pair in attA then acid order).

Do not change HIGH/LOW, the product structure, or seeds after seeing
AUC.

Files:

- `input_ahl_trackb200.txt`
- `input_acid_trackb200.txt`
- `input_atta_trackb200.txt`
- `input_attb_trackb200.txt`
- `input_rep_trackb200.txt`
- `trackb_labels.csv` columns `n;patient;y;u_ahl;u_acid;u_atta;u_attb;u_rep`

SHA-256 of each `u` sequence (12-decimal, Stage 6 recipe), frozen
before BSim:

| Sequence | SHA-256 |
|---|---|
| AHL | `b8684fafe2452bc2f9330bd7e368e559a99bd3bceea4b219226499dbe7d7825e` |
| acid | `6762700e28f9a4e2087054acc25d933080358829eadc652ea165a1a6d915b0a3` |
| attA | `dbb60b2bb4b264cdc9341e33d58422c7f3f66977889cd6e57437c44751fe34c4` |
| attB | `25d27243479ccb2b2fdf324982811ebad72f4ed21cd6a4aa55de7a7545608cc2` |
| rep | `4d50166d9c199ab517d4bb6593a80f83f451006c7f9ef8b735ff589e239bf420` |
| single-site AHL | `1efc565c29aebe9baa91260584cf3ccfff7f474e981f63725bcd95ea31167df3` |

Single-site AHL file is the per-window mean of the five `u`s, written
as `input_ahl_singlesite200.txt`. Acid hold file is the Stage 6
`input_acid_held05_200.txt` copy.

## Ridge

Washout drop then 160 rows. Lambda on rows 88..109 only (windows
128..149). Test never used for lambda or standardization. Grid
`{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`. Bias; intercept not
regularized. Training mean/std. Independent lambda per arm×seed.
Select lambda by **highest validation AUC**; ties take the larger
lambda.

Primary test metric: **AUC** (Mann–Whitney / Wilcoxon on ridge scores
vs `{0,1}`; average ranks for ties). sklearn is not required.
Accuracy at 0.5 is printed. It is NOT a gate (Stage 9: driven acc can
look good while AUC ties the field).

### Features

| Arm | Features | Count |
|---|---|---|
| Driven / silent / single-site biology | `Receiver_R_*` 20×10, `Lum_Mean_*` 20×10, `Input_Driven_Death_*` 4×2, `Den_*` 20×10 | 608 |
| Brownian | `Den_*` 20×10 | 200 |
| Field-only (driven 5-channel CSVs) | `AHL_uM_*`, `Att_conc_*`, `Rep_conc_*`, `pH_*` | 800 |

Den is in the biology ridge because attractant taxis is now a live
input. Do not add Att/Rep/AHL/pH voxels to biology. Field-only is
“read the four chemicals, skip the cells.” Do not put raw `u` into
this baseline. Do not use `Total_Deaths` as a feature.

Diagnostic, not a gate: patient-pooled AUC (mean of 5 window scores
per test patient, 10 scores). Print it. Do not replace window AUC
with it after seeing numbers.

## Gates (all required for Overall PASS)

1. Driven 5-channel test AUC > Brownian, non-overlapping mean±s.e.
   across 3 seeds, and driven > Brownian in every seed.
2. Driven 5-channel test AUC > field-only (800 chemical voxels),
   non-overlapping mean±s.e., and driven > field in every seed.
   If this fails: Overall FAIL. The plumes classify the patient. Do
   not widen the dish. Do not add ACs. Do not switch to accuracy. Do
   not raise `PROD_RATE`.
3. Driven 5-channel test AUC > single-site test AUC, non-overlapping
   mean±s.e., and 5-channel > single-site in every seed.
   If this fails: H5 is not supported. Overall FAIL. Do not retune.
4. Silent not above chance by a clear margin (mean±s.e. overlaps 0.5,
   or silent not ≈ driven). Silent must not match driven.
5. CSV 200 / 3200 / 3200; last sample `199;15;299.95`. Completeness is
   the CSV, not stdout.
6. Stage 6 NARMA-10 `GATE_EVIDENCE` still Overall PASS. BenchA
   Mackey–Glass still PASS. BenchA Lorenz still FAIL. Stage
   6/BenchA/A2 Java untouched. No glucose / Danino / vesicles / extra
   PDE / Track A2 SKIP change.

If gate 2 or 3 fails, write FAIL and stop. Do not start a 2000×1000
copy.

## Runs

9 production BSim runs: driven 5-channel, single-site, brownian
5-channel × seeds 101/202/303. Copy three silent dirs.

`output.dir`:

- `results/trackb_driven_seedXXX`
- `results/trackb_singlesite_seedXXX`
- `results/trackb_brownian_seedXXX`
- `results/trackb_silent_seedXXX` (copies)

Checker: `python examples/BSimReservoirPlanTrackB/check_trackb.py`
Evidence: `examples/BSimReservoirPlanTrackB/results/GATE_EVIDENCE.md`
