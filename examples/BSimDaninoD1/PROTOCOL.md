# Danino D1 protocol (frozen before traces)

This is an isolated BSim example box, in the style of `examples/BSimQuorumOscillator`:
one mechanism, traces, no reservoir windows. It asks whether the Danino 4-ODE
plus a correct µM ↔ molecules conversion responds to an external AHL step.

Plan Stage 3 already failed that ODE as a 300 s reservoir receiver. That FAIL
stays archived in `examples/BSimReservoirPlanStage3/ARCHIVED_NEGATIVE_RESULT.md`.
D1 does not reopen or retune it. D2, Monod/M1, vesicle ACs, two-way coupling,
and Plan Stage 3–9 work are not started.

## Why Stage 10 is not the parent

`BSimReservoirStage10` copies the Danino ODE but couples `y[1]` (µM) to the
field with `1e15` / `1e-15`. The correct factor is 602 molecules/µm³ per µM
(1 µM = 10⁻⁶ mol/L × 6.022×10²³ × 10⁻¹⁵ L/µm³ = 602 molecules/µm³).

D1 copies the ODE from Stage 10 and replaces **only** that conversion. QS_*
rates, `QS_KMLA`, Hill `n`, `QS_LTOT`, and `CELL_WALL_DIFF` are not retuned.

## Anchor files (constants copied, not retuned)

- `examples/BSimReservoirStage10/BSimReservoirStage10.java`
  QS_* block `TIME_ADJ=60`, `QS_DELTA1` .. `QS_N`, `QS_KMLA=1e-2`, `QS_LTOT=15`,
  `CELL_WALL_DIFF=3/60` (Kaplan ~20 s). ODE in `QSGRN.derivativeSystem`.
- `examples/BSimReservoirStage8/BSimReservoirStage8.java` (same constants, comments)
- `examples/BSimReservoirPlanStage3/ARCHIVED_NEGATIVE_RESULT.md`
- `examples/BSimQuorumOscillator/BSimQuorumOscillator.java` (example shape)

## Box (frozen)

- Domain 200 × 200 × 10 µm
- AHL field grid 20 × 20 × 1, D = 159 µm²/s
- Field decay: Stage 10 / Stage 8 `AHL_DECAY_RATE = 2.76e-3/60` s⁻¹
- `dt = 0.05` s
- N = 50. Growth rate 0. Clamp off. No replication, no death, no chemotaxis,
  no attractant/repellent/acid/glucose, no luminescence L, no voxels,
  no 300 s windows, no NARMA.
- `FLOW_SPEED = 0`
- Seed 101, one production run
- Preview camera (I/O), proven 1000×500 sequence scaled to this box:
  `ortho(0, BOUND_X, BOUND_Y, 0, -1000, 10000);`
  `camera(BOUND_X/2, BOUND_Y/2, BOUND_Y, BOUND_X/2, BOUND_Y/2, 0, 0, 1, 0);`
  `perspective(PI/2, BOUND_X/BOUND_Y, 0.1, 10000);`
  Window 800×600

## Open-loop bath protocol (frozen)

The extracellular field is a well-mixed bath, not a point source and not an
AC vesicle.

| Interval | Bath | AHL_ext |
|---|---|---|
| t = 0 .. 600 s | OFF | 0 |
| t = 600 .. 2400 s | ON | 0.05 µM everywhere |
| t = 2400 .. 5400 s | OFF | 0 |

0.05 µM is 5 × `QS_KMLA` (0.01 µM). In field units:

`0.05 µM = 0.05 * 602.0` molecules/µm³.

Code convention at 1 s log resolution (half-open on the falling edge):

- bath ON iff `600 <= t < 2400`
- bath OFF otherwise (`t < 600` and `t >= 2400`)

Each tick, after `field.update()`, set the field concentration uniformly to
that value (or 0). Do not calibrate a source rate after seeing traces.

## Membrane coupling (open loop)

Cells run the full 4-ODE. Membrane exchange **reads** the bath:

```
AHL_ext_uM = field.getConc(position) / 602.0
extraintradiff_uM = y[1] - AHL_ext_uM
dy[1] includes - CELL_WALL_DIFF * extraintradiff_uM
```

as in Stage 10, with 602 instead of `1e15`/`1e-15`.

Cells do **not** write the field: skip `addQuantity` from the membrane term.
Intracellular production terms stay (`QS_KP2 * y[0]`, etc.).
ICs: Stage 10 `getICs()` `{0.05, 0.05, 0.05, 0.05}`.

## Export

Every 1 s, one CSV (semicolon) to `results/d1_seed101/d1_timeseries.csv`:

```
t_s;AHL_ext_uM_mean;LuxI_mean;AHL_in_mean;AiiA_mean;LA_mean;N
```

`AHL_ext_uM_mean` MUST be field/602, never field×1e-15.
Last row `t = 5400`. Completeness is that row, not stdout.

## Gates (identity, not a reservoir gate)

All required:

1. During the pulse (`600 <= t < 2400`), `AHL_ext_uM_mean = 0.05 ± 1e-6`
   (bath actually on).
2. During the last 10 min of the pulse (`t` in `[1800, 2400]`), mean AHL_in
   (`y[1]`) is in `[0.001, 10]` µM — order of `QS_KMLA`, not 1e-12.
3. Mean LA (`y[3]`) in `[1800, 2400]` is higher than mean LA in the
   pre-pulse window `[0, 600]`. Print both numbers.
4. Mean AHL_in in the last 10 min after the pulse (`t` in `[4800, 5400]`)
   is lower than in `[1800, 2400]`. Print both.
5. N is 50 at every log row. No growth/death.
6. Plan Stage 3 remains labelled FAIL. QS_* rates unchanged from Stage 10.

If gate 2 fails, the unit coupling is still wrong: stop. Do not retune
`QS_KPLI`, `QS_KP2`, or the bath concentration.

If gates 3–4 fail with gate 2 passing, document
“ODE does not respond to a 0.05 µM bath on this timescale” and stop.
Do not start D2.

## Out of scope

Do not start D2 (closed-loop oscillation), Monod, vesicle ACs, or two-way
coupling. Do not copy this box into the Plan Stage 6 dish. Do not modify
`BSimReservoirPlanStage*` or `BSimReservoirStage10`.
