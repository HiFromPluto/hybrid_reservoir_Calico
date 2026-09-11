# E4.1 — one-way AHL AC transducer vs ideal current

Track E4, one-way only. Frozen HybridDish / Narma10b claim dish.
`GATE_EVIDENCE.md` is not edited. Kinetics, layout, and flow are not
retuned. Claim dish stays CENTER / `FLOW=0`.

This package asks whether replacing the ideal current
\(J=J_{\max}u\) with one static AHL Hill gate, at matched commanded
payload, still carries NARMA-10 above Brownian and silent, and
whether the living 408-D readout still beats field-only.

**MODEL_STATUS = `HYPOTHETICAL_DESIGN_ENVELOPE`.** There is no
measured AC dose-response in the project. This is not a digital twin.
`reservoir_new` vesicle `ArtificialCell` is not copied.

## Status

**Module 0 complete.** All 12 device items `MISSING`.
**MODEL_STATUS = `HYPOTHETICAL_DESIGN_ENVELOPE`.** Not a digital twin.

**Module 1 complete: PASS.** Unit tests A–E passed. Payload relative
error `1.62e-16`. Living Java was authorized after this PASS.

**Module 2 seed 111 complete.** Occupancy **ALIVE**
(`mean_R=0.1915`). System **PASS** (F408 `0.9007` vs Brownian
`1.1625` and silent `1.1622`). Living-layer **PASS** vs A1 field
`1.0064`. |A1−A0 F408| = `0.0305`. Direct-input `0.6828` still wins.
AC 10-tap `g` scores `0.6535` (attribution, not a living-layer
claim). **Do not promote A1 as a new claim dish** (`0.03` tick around
the weak `0.93` band). Seeds 222/333 were not started. Two-way was
not started.

See `PROTOCOL.md`, `results/AC_SPEC_INVENTORY.md`,
`results/E4_AC_CALIBRATION.md`, and `results/E4_AC_SCOUT.md`.

## Frozen A1 numbers

See `PROTOCOL.md`. Envelope `g0=0`, `n_AC=2`, `K_AC=0.25`.
`Jmax_A1=7.467750975821304e7` molecules/s matches A0 commanded mass
on the frozen Narma10b `u` file
(`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`).

Do not retune `K`, `n`, `tau_R`, `tau_L`, clamp, mortality, flow, or
layout. Do not fit `K_AC` to NARMA NRMSE. Direct-input `0.6828`
remains the task ceiling. `0.93` stays a weak predictor.
