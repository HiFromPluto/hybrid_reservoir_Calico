# BSimReservoirPlanTrackB mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanStage6` Java onto
the Stage 6 dish. Stage 6 is not modified and remains the NARMA-10
provenance checkpoint (**PASS**). Track A Mackey–Glass remains **PASS**.
Track A Lorenz remains **FAIL** and is not rewritten. Stage 9
product-bit remains **FAIL**. Stage 7 four-source SVD remains **FAIL**.
Kinetics are not retuned. Clamp, acid source rate, AHL source rate,
`K_MAX`, `GROWTH_RATE`, `PROD_RATE=1e6`, and Stage 3B receiver `K=1.6`,
`n=2`, `tau=15` are not retuned. Glucose, Monod, Danino, ODE-gated
vesicle ACs, pH taxis, `setGoal(repellent)`, Stage 7 four-mode rank,
and a wider dish are not imported.

Track B tests H5: five one-way chemical inputs versus a scalar AHL mix
on a synthetic 5-channel patient classification task. Field-only gate
copied from Stage 9. If the plumes classify the patient, Overall is
FAIL. Do not then add ACs, stretch the domain, or switch the gate to
accuracy.

## Frozen 5-AC layout (do not move after AUC)

| Channel | Chemical | Position | Rate × u |
|---|---|---|---|
| AHL | ahlField | (500, 250, 5) | 1.28e8 |
| acid | acidField | (300, 375, 5) | 2e11 |
| attA | attractantField | (150, 100, 5) | 1e6 |
| attB | attractantField | (850, 400, 5) | 1e6 |
| rep | repellentField | (150, 400, 5) | 1e6 |

attA and attB share one attractant PDE. Stage 6 silent sites
`(250, 250, 5)` and `(750, 250, 5)` are not used. Four chemical fields,
same as Stage 6. No extra PDE.

## Four arms

Identical dish, timing, and seeds `101/202/303` except as named.

1. **Driven 5-channel.** All five sites carry biomarker files. Acid is a
   free input.
2. **Single-site H5 control.** Only AHL fires (`u` = mean of five
   channels that window, clipped to `[0, 0.5]`). attA=attB=rep=0. Acid
   held at 0.5. Same patient labels.
3. **Brownian null.** Same 5-channel pulses as driven. Particles do not
   sense. Ridge is `Den_*` 20×10 = 200.
4. **Biology silent.** Stage 6 silent CSVs reused; not rerun.

## Frozen readout

Driven / silent / single-site biology, and only these:

- voxel `Receiver_R_*` / Mean_q (20×10)
- voxel `Lum_Mean_*` / Mean_L (20×10)
- coarse `Input_Driven_Death_*` (4×2)
- voxel `Den_*` (20×10)

= 608. Den is in because attractant taxis is a live input. Att/Rep/AHL/pH
voxels stay out of biology. `Total_Deaths` is not a feature.

Field-only (required gate, driven 5-channel CSVs): `AHL_uM_*`,
`Att_conc_*`, `Rep_conc_*`, `pH_*` = 800. Do not put raw `u` into this
baseline.

## Protocol

See `PROTOCOL.md`. Class vector, balances, and all SHA-256s are frozen
there before BSim. Washout 40 / train 110 / test 50 at patient
boundaries. Lambda on rows 88..109 only. Highest validation AUC; ties
take the larger lambda. Primary metric is window AUC.

## Gates

All required. If gate 2 or 3 fails, write FAIL and stop. Do not start a
2000×1000 copy. Do not raise `PROD_RATE`. Do not switch to accuracy.
