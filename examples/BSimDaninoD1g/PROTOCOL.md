# Danino D1g protocol (frozen before traces)

Citation box only. Parent is `examples/BSimDaninoD1/` (602 conversion, open-loop
bath, QS_* from Stage 10). D1 remains a closed FAIL and is not rewritten.

This package asks whether the Danino 4-ODE can raise LA on a 0.05 µM bath once
dilution and zero ICs are restored. It is not a reservoir window experiment
and is not merged into Plan Stage 6 even if it passes.

Stage 6 PASS with the Stage 3B Hill receiver stays the hybrid dish. Plan Stage 3
remains ARCHIVED FAIL. D2, Monod/M1, vesicle ACs, two-way coupling, and Plan
Stage 3–9 work are not started. This ODE is not put into the 300 s hybrid dish.

## Why D1 failed (not rewritten here)

Open loop, 602 conversion, and micromolar AHL_in all passed. LA fell because
(1) Stage 10 ICs `{0.05,…}` start LA above what this bath can hold, and
(2) growth/dilution was off, so AiiA ran away. BSim replication copies `y` to
the daughter and does not dilute. Danino protein equations need an explicit
linear µ.

## Keep frozen from D1

- Domain 200 × 200 × 10 µm, grid 20 × 20 × 1, D = 159 µm²/s, dt = 0.05 s
- Field decay as in D1: `AHL_DECAY_RATE = 2.76e-3/60` s⁻¹
- N = 50, BSim growth rate 0, clamp off, no replication, no chemotaxis,
  no acid/glucose/L, no voxels, no 300 s windows
- `FLOW_SPEED = 0`, seed 101, one production run
- Open-loop well-mixed bath; cells do not write the field
- Pulse: t = 0..600 off; 600..2400 ON at 0.05 µM = 0.05×602 molecules/µm³;
  2400..5400 off
- QS_* rates, `QS_KMLA=1e-2`, `QS_LTOT=15`, `QS_N=2`, `CELL_WALL_DIFF=3/60`
  copied from D1/Stage 10. Do not retune after traces.
- `AHL_ext_uM = field/602`. No 1e15.
- Camera sequence as D1. Window 800×600.

Code convention at 1 s log resolution (half-open on the falling edge):

- bath ON iff `600 <= t < 2400`
- bath OFF otherwise

Each tick, after `field.update()`, set the field concentration uniformly.
Do not calibrate a source rate after seeing traces.

## Change from D1 (only these)

1. ICs `{0, 0, 0, 0}` for LuxI, AHL_in, AiiA, LA. Not `{0.05,…}`.
2. Dilution `µ = Math.log(2.0) / 1800.0` (ln 2 / 30 min doubling).
   Compile-time constant. Print it.
   Add `-µ*y[i]` to `dy[0]` (LuxI), `dy[2]` (AiiA), `dy[3]` (LA) only.
   Do not add µ to `dy[1]` (AHL already has `QS_T_A` and membrane exchange).
   Do not turn on BSim growth or change N. N stays 50.

## Export

Every 1 s, one CSV (semicolon) to `results/d1g_seed101/d1g_timeseries.csv`:

```
t_s;AHL_ext_uM_mean;LuxI_mean;AHL_in_mean;AiiA_mean;LA_mean;N
```

`AHL_ext_uM_mean` MUST be field/602. Last row `t = 5400`. Completeness is that
row, not stdout.

Sidecar `results/d1g_seed101/d1g_params.csv` prints µ and ICs so an audit can
see they changed from D1.

## Gates (same identity gates as D1; all required)

1. During the pulse (`600 <= t < 2400`), `AHL_ext_uM_mean = 0.05 ± 1e-6`.
2. Mean AHL_in in `[1800, 2400]` ∈ `[0.001, 10]` µM (not 1e-12).
3. Mean LA in `[1800, 2400]` > mean LA in `[0, 600]`. Print both.
4. Mean AHL_in in `[4800, 5400]` < mean AHL_in in `[1800, 2400]`. Print both.
5. N = 50 at every log row.
6. D1 remains FAIL. Plan Stage 3 remains ARCHIVED FAIL. QS_* unchanged
   except the added `-µ` terms. Conversion 602. ICs are zero.

If gate 3 fails, D1g is FAIL. Document it. Do not retune `QS_KPLI`, `QS_KP2`,
`QS_KR1ON`, µ, or the 0.05 µM bath. Do not start D2.

If it passes, stop. Do not merge into `BSimReservoirPlanStage6`.

## Out of scope

Do not modify `examples/BSimDaninoD1/`. Do not start D2, Monod/M1, or Plan
Stage 3–9. Do not modify `BSimReservoirPlanStage*` or `BSimReservoirStage10`.
Do not put this ODE into the 300 s hybrid dish.
