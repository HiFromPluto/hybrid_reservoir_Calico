# Monod M3 protocol (frozen before traces)

Parent is `examples/BSimMonodM2/` (Monod growth + uptake, no replenisher).
M1 and M2 remain closed PASSes and are not rewritten.

This package asks whether one constant point source makes a spatial G
contrast against that uptake. It is **not** Stage 7 AC5 and **not**
Stage 11. No vesicle ODE. No Model A replenisher. No chemotaxis on glucose.

Main reservoir (out of scope): Stage 6 PASS with the Stage 3B Hill receiver,
no glucose, no Monod. If M3 passes, still do not copy glucose into
`BSimReservoirPlanStage6`. Do not start D2 or Plan Stage 3–9 work.
Do not modify `BSimDaninoD1`, `BSimDaninoD1g`, `BSimReservoirPlanStage*`,
or `BSimReservoirStage11`.

## Why M2 is not enough

M2 showed uptake writes the field: a uniform G0 = 1.80 µM collapsed. That
cannot tell whether a point source can hold a near/far contrast against
uptake and diffusion. M3 adds one constant source and a silent spatial
control. Field decay stays 0 so the only sink is uptake.

## Unchanged from M2

```
µ = µ_max * G / (K_s + G)
K_s            = 0.18 µM
µ_max          = Math.log(2.0) / 1800.0
GROWTH_RATE_SA = 4.0 * Math.PI / 1800.0
K_UPTAKE       = 1.2e3 molecules/s
1 µM           = 602 molecules/µm³
D              = 100 µm²/s
decay          = 0
```

Surface-area growth, replication on, no clamp, start `N = 50`,
desynchronised `setRadius()`, no `setGoal`. Do not retune `K_s`, `µ_max`,
or `K_UPTAKE` after seeing traces.

## Change from M2 (only these)

1. Initialize the field once at `t = 0` to `G0 = 0.18 µM` (`K_s`), not
   1.80. No per-tick bath reset.
2. One source at `(50, 100, 5)`. Each tick, after `field.update()`:

   ```
   field.addQuantity(source, K_SOURCE * dt)
   ```

   `K_SOURCE` is frozen from M2's uptake scale, **not** from seeing
   contrast:

   ```
   N0 * K_UPTAKE = 50 * 1.2e3 = 6e4 molecules/s at saturating G
   K_SOURCE = 1.0e5 molecules/s   (~1.7× that)
   ```

   Print it. Engineering identity rate. Do not change after traces.
   Fix wiring bugs (position, dt, 602) if both arms look the same because
   the source never fired. That is not a retune of `K_SOURCE`.
3. Two runs, seed 101, 7200 s each:

   | Arm | K_SOURCE | notes |
   |---|---|---|
   | M3on | 1.0e5 | source on |
   | M3off | 0 | silent spatial control |

   Java may force `K_SOURCE = 0` when `arm = off`.

## Near / far (frozen, print)

```
near: x < 70 µm    (source at x = 50)
far:  x > 150 µm
```

`G_near` / `G_far` = mean field/602 over voxels whose centre `x` satisfies
that cut. Births use the child particle's `x`, not the voxel centre.

## Births (diagnostic only)

Count replications by child position into near vs far, cumulative.
Over 7200 s, run-tumble mixing length is larger than this box, so births
may equalize. Do **not** fail M3 on births. Do **not** add glucose taxis
to force a birth contrast.

## Box (frozen, same as M2)

- Domain 200 × 200 × 10 µm
- Glucose field grid 20 × 20 × 1
- `dt = 0.05` s
- `FLOW_SPEED = 0`
- No AHL, QS, Danino, acid, attractant, luminescence, vesicles, clamp,
  voxels, or 300 s RC windows
- Camera: same ortho / camera / perspective sequence as D1, this box size
- Window 800 × 600
- Seed 101

## Export

Every 10 s, one CSV (semicolon) per arm:

```
t_s;G_uM_mean;G_uM_near;G_uM_far;N;births_near_cum;births_far_cum;arm
```

Paths:

- `results/m3on_seed101/m3_timeseries.csv`
- `results/m3off_seed101/m3_timeseries.csv`

`G_uM_*` MUST be field/602. Last row `t = 7200`. Completeness is that row,
not stdout.

## Gates (all required)

1. M3on: mean `G_near` in `[6000, 7200]` > mean `G_far` in `[6000, 7200]`.
   Print both. This is the spatial-supply identity.
2. M3off: mean `G_near` / mean `G_far` in `[6000, 7200]` is in `[0.5, 2.0]`
   (no large leftover gradient when the source is off). Print the ratio.
3. M3on late `G_near` > M3off late `G_near`. Print both.
4. `N_final > N_initial` on both arms.
5. `K_UPTAKE`, `K_s`, `µ_max` unchanged from M2. No replenisher, no vesicle
   AC, no `setGoal(glucose)`.
6. M1 PASS, M2 PASS, Plan Stage 6 has no glucose.

If gate 1 fails, the point source does not beat uptake/diffusion on this
box: **FAIL, stop**. Do not raise `K_SOURCE` after seeing it. Do not add
chemotaxis or a replenisher.

If M3 passes, stop. Still do not merge glucose into the Stage 6 dish.

## Out of scope

Do not modify `examples/BSimMonodM1/` or `examples/BSimMonodM2/`.
Do not start D2 or Plan Stage 3–9 work. Do not copy glucose into the
Plan Stage 6 dish.
