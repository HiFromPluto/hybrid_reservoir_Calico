# Stage 9 protocol (frozen before AUC)

This is Plan Stage 9 on the frozen Stage 6 dish: one-way hybrid
classification (AHL AC → E. coli). NARMA-10 is already closed as Stage 6
PASS. This stage does not claim a 4-AC or 5-AC reservoir. Two-way
AC↔bacteria coupling is the next project, not this one.

Stage 8 is skipped. Plan Stage 8 is the Danino LuxI/AHL/AiiA/LuxR-AHL
ODE; that already ran as Plan Stage 3 and is ARCHIVED FAIL. The frozen
receiver remains Stage 3B: `dR/dt = (C^2/(1.6^2+C^2)-R)/15`, `K=1.6`,
`n=2`, `tau=15`. Do not retune Danino constants. Do not port Stage 10's
ODE. The plan text allows this skip: the Hill receiver is defensible.

Stage 4 remains FAIL. Stage 5 remains PASS. Stage 6 remains PASS and is
not modified. Stage 7 remains FAIL: do not classify on the four-corner
layout, do not loosen SVD, do not raise `PROD_RATE`. Kinetics are not
retuned.

There is no `classification_task/` pipeline in this repo. The task below
is the stand-in for the missing biomarker pipeline. It is declared here
before any AUC is computed. Do not change the threshold after seeing AUC.

## Dish (do not change)

Evaluate on the existing Stage 6 CSVs. Do not rerun BSim unless a CSV is
missing or ragged. Completeness is last sample `199;15;299.95` on disk,
not BSim stdout.

Copied from Stage 6 / Stage 5 and not retuned:

- One AHL AC at `(500, 250, 5)`, acid AC at `(300, 375, 5)` held at `0.5`,
  attractant ACs silent, `FLOW_SPEED=0`
- `dt=0.05`, warmup `18000 s`, windows `300 s`, pulses `75 s`,
  whole-window sampling every `20 s`
- Field grid `50x25x1`, readout `20x10` and `4x2`
- Receiver `K=1.6`, `n=2`, `tau=15`
- `ALPHA_LUX=1/1500`, `DELTA_LUX=1/1500`
- Acid source rate `2e11`, `K_MAX=0.002`, `GROWTH_RATE` compile-time constant
- Clamp `K=2000` unchanged
- Seeds `101 / 202 / 303` shared across Brownian, silent, and driven

ACs in this rebuild are one-way injection points: stimulus `u[n]` →
`field.addQuantity`. Pore scatter is not added. Vesicle store /
leaky-integrator `ArtificialCell` is not added. Glucose, Monod,
ODE-gated vesicle ACs, Danino, and pH taxis are not imported.

If a preview Java copy is added later, physics must stay identical.
Restore the Stage 7 camera sequence and draw marker spheres at the AHL
and acid AC sites (preview I/O only). Do not draw the Stage 7
four-corner layout. Do not change rates.

## Frozen AHL sequence

The same Stage 6 file `input_ahl_narma200.txt`. `u[n] ~ Uniform[0, 0.5]`
from `random.Random(20260814)`. SHA-256 of the 12-decimal sequence:

`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`

Do not regenerate it. Acid remains held at `0.5` (Stage 6 option a).

## Frozen classification task

For window `n`:

```
y[n] = 1 if u[n] * u[n-1] > 0.0625 else 0
u[-1] = 0
```

That is a product bit (nonlinear, one lag). Threshold `0.0625` is frozen.
Do not change it after seeing AUC.

Class balance on the frozen sequence (computed from `u` only, before AUC):

| Split | Windows | Count | Positives | Fraction |
|---|---|---|---|---|
| Washout | 0 .. 39 | 40 | 14 | 0.350 |
| Train | 40 .. 149 | 110 | 38 | 0.345 |
| Test | 150 .. 199 | 50 | 23 | 0.460 |
| All | 0 .. 199 | 200 | 75 | 0.375 |

`analyze_stage9.py` reprints these counts.

## Washout / train / test

| Split | Windows | Count |
|---|---|---|
| Washout | 0 .. 39 | 40 |
| Train | 40 .. 149 | 110 |
| Test | 150 .. 199 | 50 |

Lambda on train only: fit windows 40 .. 127, pick `lambda` on windows
128 .. 149, then refit on all 110 train windows. After stripping
washout, the design matrix `X` has 160 rows. Validation MUST be
`X[inner_train_end:TRAIN]` (windows 128 .. 149), **not**
`X[inner_train_end:]`. The Stage 6 slice `X[inner_train_end:]` included
test. Do not copy it.

Test windows are never used to choose `lambda` or to standardize.

## Ridge rule (same as Stage 6, now on {0,1})

- Bias column of ones; intercept is not regularized
- Features standardized with training mean/std (zero-variance columns kept at 0)
- Grid: `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`
- Select the `lambda` with **highest validation AUC**; ties take the larger `lambda`
- Test metric: **AUC** (primary). Accuracy at threshold 0.5 is printed, not the gate
- Independent `lambda` per arm × seed, identical procedure

AUC is the Mann–Whitney / Wilcoxon rank statistic on ridge scores versus
`{0,1}` labels (average ranks for ties). sklearn is not required.

## Readout state (same 408-feature contract as Stage 6)

Per-window features, after whole-window sampling:

- State channels: mean of the 16 intra-window samples
- Count channel: last sample of the window (full-window accumulation)

| Arm | Features | Count |
|---|---|---|
| Brownian | `Den_*` 20×10 | 200 |
| Silent | `Receiver_R_*` 20×10, `Lum_Mean_*` 20×10, `Input_Driven_Death_*` 4×2 | 408 |
| Driven biology | same 408 as silent | 408 |
| Linear / field baseline | driven-arm voxel `AHL_uM_*` 20×10 only | 200 |

The field baseline is a required gate, not optional. It is "read the AC
carrier, skip the cells." Do not put raw `u[n]` into this baseline (that
is knowing the stimulus without a dish).

A diagnostic delay-line of `(u[n], u[n-1])` may be printed. It is not a
substitute for the field baseline.

Not in any ridge: window AHL summary, occupancy, pH, Births,
Total_Deaths, Clamp, OOB, Population, Lum_Sum, voxel pH/Den (biology
arms), `Att_*`.

Brownian particles do not sense. Silent sources are off. Driven biology
uses the frozen Stage 5/6 list only.

## Gates (all required)

1. Driven test AUC > Brownian in every seed, and non-overlapping
   mean ± s.e. across the 3 seeds.
2. Driven test AUC > field-only linear baseline in every seed, and
   non-overlapping mean ± s.e. If the field wins, Stage 9 is FAIL: the
   AC plume is doing the work, not the hybrid medium.
3. Silent test AUC is not above chance by a clear margin (mean ± s.e.
   overlaps 0.5, or silent worse than driven with no overlap). If
   silent ≈ driven, FAIL.
4. Driven uses only the frozen 408 features. Print the list.
5. CSV rectangularity 200/3200/3200 on the Stage 6 files used.
   Completeness is last sample `199;15;299.95`, not stdout.
6. Stage 4 FAIL, Stage 5 PASS, Stage 6 PASS, Stage 7 FAIL unchanged.
   Stage 8 skipped (Danino already archived). Kinetics not retuned.

If the gate fails, document why (brownian wins, field wins, or
silent≈driven). Do not add features, retune kinetics, move ACs, import
vesicles/Danino, or start two-way coupling.
