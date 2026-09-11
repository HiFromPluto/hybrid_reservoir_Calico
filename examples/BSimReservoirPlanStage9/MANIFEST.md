# BSimReservoirPlanStage9 mechanism manifest

## Scope

This package evaluates Plan Stage 9 on frozen Stage 6 CSVs. Stage 6 is
not modified and remains **PASS**. Stage 5 remains **PASS**. Stage 4
remains **FAIL**. Stage 7 remains **FAIL** and is not used (no
four-corner layout, no SVD loosening, no `PROD_RATE` raise). Stage 8 is
**skipped**: Danino already archived as Plan Stage 3. Clamp, acid source
rate, `K_MAX`, `GROWTH_RATE`, and Stage 3B receiver `K=1.6`, `n=2`,
`tau=15` are not retuned. Glucose, Monod, Danino, ODE-gated vesicle ACs,
pH taxis, pore scatter, and two-way AC↔bacteria coupling are not
imported.

ACs remain one-way injection points: `u[n]` → `field.addQuantity`.

## Why not a new BSim run

The dish is the Stage 6 dish. Prefer evaluation on
`examples/BSimReservoirPlanStage6/results/`. Re-run only if a CSV is
missing or ragged. Completeness is last sample `199;15;299.95`.

## Four arms

Identical Stage 6 dish, timing, and seeds `101/202/303` except as named.

1. **Brownian.** `BSimParticle` drift only. Particles do not sense.
   Ridge is spatial density `Den_*` on the 20×10 grid (200).
2. **Silent.** Full Stage 5/6 cells, sources off. Ridge is the frozen
   408 biology features. Must sit near chance.
3. **Driven biology.** AHL carries the frozen Stage 6 sequence. Acid
   held at 0.5. Ridge: `Receiver_R_*`, `Lum_Mean_*`,
   `Input_Driven_Death_*` (408).
4. **Linear / field baseline (required gate).** Driven-arm voxel
   `AHL_uM_*` 20×10 only. Read the AC carrier, skip the cells. Raw
   `u[n]` is not a field-baseline feature.

A diagnostic delay-line of `(u[n], u[n-1])` may be printed. It is not a
substitute for arm 4.

## Frozen readout

Silent and driven ridge features, and only these:

- voxel `Receiver_R_*` / Mean_q (20×10)
- voxel `Lum_Mean_*` / Mean_L (20×10)
- coarse `Input_Driven_Death_*` (4×2)

Not in any ridge: window AHL, occupancy, pH, Births, Total_Deaths,
Clamp, OOB, Population, Lum_Sum, voxel pH/Den (biology arms), `Att_*`.

## Frozen task

See `PROTOCOL.md`. Product bit on the Stage 6 AHL sequence:

`y[n] = 1` if `u[n] * u[n-1] > 0.0625` else `0`, `u[-1] = 0`.

SHA-256 `d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.
Do not regenerate `u`. Do not change 0.0625 after seeing AUC.

## Protocol

Same Stage 6 timing and splits. Lambda fit 40..127, pick on 128..149
only (`X[inner_train_end:TRAIN]` after washout strip), refit 40..149.
Do not copy the Stage 6 `X[inner_train_end:]` slice (it included test).

## Gates

All required. Primary metric is test AUC. Accuracy at 0.5 is printed.

1. Driven test AUC > Brownian every seed, non-overlapping mean ± s.e.
2. Driven test AUC > field-only every seed, non-overlapping mean ± s.e.
   Field win = FAIL (plume, not hybrid).
3. Silent not above chance by a clear margin; silent not ≈ driven.
4. Driven uses only the frozen 408 features. Print the list.
5. CSV 200/3200/3200 rectangular; last sample `199;15;299.95`.
6. Prior stage labels unchanged. Stage 8 skipped. Kinetics not retuned.

If this fails, document why. Do not add features, retune, move ACs,
import vesicles/Danino, or start two-way coupling.
