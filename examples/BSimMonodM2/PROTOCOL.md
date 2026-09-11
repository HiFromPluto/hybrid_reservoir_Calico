# Monod M2 protocol (frozen before traces)

Parent is `examples/BSimMonodM1/` (Monod growth on a clamped µM bath, 602
conversion, `K_s = 0.18 µM`, `µ_max = ln 2 / 1800`). M1 remains a closed
PASS and is not rewritten.

This package asks whether uptake actually removes glucose from the field.
There is **no replenisher** and **no nutrient AC**. Those are M3 / Stage 11
and stay out.

Main reservoir (out of scope): Stage 6 PASS with the Stage 3B Hill receiver,
no glucose, no Monod. If M2 passes, still do not copy glucose into
`BSimReservoirPlanStage6`. Do not start M3, D2, or Plan Stage 3–9 work.
Do not modify `BSimDaninoD1`, `BSimDaninoD1g`, `BSimReservoirPlanStage*`,
or `BSimReservoirStage11`.

## Why M1 is not enough

M1 showed `µ = µ_max G/(K_s+G)` with G glued at the commanded bath. That
cannot tell whether cells write the field. M2 drops the bath reset and
adds a per-cell sink. Field decay is 0 so the only sink is uptake
(otherwise decay looks like uptake).

## Growth law (unchanged from M1)

```
µ = µ_max * G / (K_s + G)
K_s    = 0.18 µM
µ_max  = Math.log(2.0) / 1800.0
GROWTH_RATE_SA = 4.0 * Math.PI / 1800.0
```

Do not retune these after seeing traces. Do not use `reservoir_new`'s
unitless `mu_max = 1`, `K_G = 5`. Do not use Senn's 0.92 h⁻¹.

Surface-area growth, replication on, no clamp, start `N = 50`, desynchronised
`setRadius()`, same as M1.

## Glucose units (unchanged)

```
1 µM = 602 molecules/µm³
G_uM = conc / 602
```

## Change from M1 (only these)

1. Do **not** reset the field to a bath every tick. Initialize once at
   `t = 0` to `G0`, then let diffusion / decay / uptake act.
   - Field decay = 0 (only sink is uptake).
   - Diffusivity = 100 µm²/s (Plan attractant-scale; not tuned after traces).
2. Each cell, each tick, after `grow()`:

   ```
   gUm = conc / 602
   uptake = K_UPTAKE * gUm / (K_s + gUm)    // molecules/s
   field.addQuantity(position, -uptake * dt)
   ```

   Then clamp any negative voxel concentration to 0.
3. `K_UPTAKE` is frozen from a mass-balance estimate, **not** from seeing
   `G(t)`:

   ```
   box volume = 200 * 200 * 10 = 4e5 µm³
   G0 = 1.80 µM = 1.80 * 602 molecules/µm³
   total molecules ≈ 4.33e8
   N = 50, want ~50% depletion in ~3600 s at saturating Monod
   => K_UPTAKE = 1.2e3 molecules/s
   ```

   Print it. Label: engineering identity rate, not Senn 1994.
   Do not change it after seeing traces.
   Fix wiring bugs (sign, dt, 602, clamp) if G is flat because of a code
   error. That is not a retune of `K_UPTAKE`.

## Box (frozen, same as M1)

- Domain 200 × 200 × 10 µm
- Glucose field grid 20 × 20 × 1
- `dt = 0.05` s
- `FLOW_SPEED = 0`
- No AHL, QS, Danino, acid, attractant, luminescence, vesicles, clamp,
  voxels, or 300 s RC windows
- Camera: same ortho / camera / perspective sequence as D1, this box size
- Window 800 × 600
- Seed 101

## Production run (frozen)

One run, seed 101, duration 7200 s. Completeness is the last CSV row
`t = 7200`, not stdout.

| Run | G0 (µM) | notes |
|---|---|---|
| M2a | 1.80 | same G as M1c, but G is free |

Reference (do not rerun M1): M1c at G = 1.80 had `N(7200) = 664` and
`µ_hat/µ_max = 0.880` with G glued at 1.80.

## Export

Every 10 s, one CSV (semicolon):

```
t_s;G_uM_mean;G_uM_min;N;run_label
```

Path: `results/m2a_seed101/m2_timeseries.csv`

`G_uM_*` MUST be field/602. Last row `t = 7200`.

## Optional diagnostic (not a gate)

`µ_hat` on `[600, 7200]` vs M1c. If G collapses early, `µ_hat` will be
smaller; that is expected, not a FAIL.

```
µ_hat = (ln N_end - ln N_start) / (t_end - t_start)
```

## Gates (all required)

1. `G_uM_mean` at `t = 0` is `1.80 ± 0.01` µM.
2. Uptake is live: `G_uM_mean(7200) < 0.90 * G_uM_mean(0)` (at least 10%
   drop). If G is flat, FAIL: uptake is not hitting the field. Stop.
   Do not raise `K_UPTAKE` after seeing it.
3. `G_uM_mean` is lower in `[6000, 7200]` than in `[0, 600]`. Print both.
4. `N_final > N_initial` (uptake did not sterilise the box).
5. `K_s`, `µ_max`, `GROWTH_RATE_SA` unchanged from M1. No replenisher, no AC.
6. M1 still PASS; Plan Stage 6 still has no glucose.

If gate 2 fails, document and stop. Do not add Model A replenisher or an
AC to keep G pretty. Do not start M3.

If M2 passes, stop. M3 (one nutrient point source) is a later prompt and
still does not go into Stage 6.

## Out of scope

Do not modify `examples/BSimMonodM1/`. Do not start M3, D2, or Plan
Stage 3–9 work. Do not copy glucose into the Plan Stage 6 dish.
