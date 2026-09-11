# BSimReservoirPlanWaveform mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanBenchA` Java onto
the Stage 6 dish. Stage 6 is not modified and remains the NARMA-10
provenance checkpoint (**PASS**). BenchA Mackey–Glass remains **PASS**.
BenchA Lorenz remains **FAIL** and is not rewritten. Track B remains
**FAIL**. Stage 9 product-bit remains **FAIL**. Kinetics are not
retuned. Clamp, acid source rate, AHL source rate, `K_MAX`,
`GROWTH_RATE`, and Stage 3B receiver `K=1.6`, `n=2`, `tau=15` are not
retuned. Glucose, Monod, Danino, ODE-gated vesicle ACs, attractant/
repellent encoding, Stage 7 four-mode rank, and a wider dish are not
imported.

Waveform classification drives the same dish with deterministic sine /
square / triangle AHL templates. New AHL file and checker only.
408-feature ridge. 200-window budget. Field-only AHL is a required
gate. If the plume ranks the waveform, Overall is FAIL. Do not then
mix in Att/Rep, raise source rates, or make the templates more similar.

## Three arms

Identical dish, timing, and seeds `101/202/303` except as named.

1. **Brownian null.** `BSimParticle` drift and diffusion only. Same AHL
   / acid pulses as driven. Ridge is spatial density `Den_*` on the
   20×10 grid.
2. **Biology silent.** Stage 6 silent CSVs reused; not rerun.
3. **Biology driven.** Full Stage 6 cells. AHL carries the waveform
   sequence. Acid held at 0.5.

## Frozen readout

Silent and driven ridge features, and only these:

- voxel `Receiver_R_*` / Mean_q (20×10)
- voxel `Lum_Mean_*` / Mean_L (20×10)
- coarse `Input_Driven_Death_*` (4×2)

= 408. Voxel AHL, Den (biology), pH, Att, and raw `u` stay out.

Field-only (required gate, driven CSVs): voxel `AHL_uM_*` = 200.

## Protocol

See `PROTOCOL.md`. Class vector, templates, and SHA-256s are frozen
there before BSim. Washout 40 / train 110 / test 50 at block
boundaries. Lambda on rows 88..109 only. Highest macro OVR validation
AUC; ties take the larger lambda. Primary metric is window macro OVR
AUC. Accuracy vs 1/3 is printed, not a gate.

## Gates

All required. If gate 2 fails, write FAIL and stop. Do not add ACs.
Do not switch to accuracy vs 0.334.
