# Monod M1 protocol (frozen before traces)

This is an isolated BSim example box, in the same spirit as
`examples/BSimDaninoD1/`: one mechanism, traces, no reservoir windows.
It asks whether Monod growth

```
µ = µ_max * G / (K_s + G)
```

is actually wired, with literature units. Glucose is a uniform bath, like
D1's AHL bath. There is **no uptake**.

Stage 11's Monod + nutrient-AC merge stays a failed reservoir path. M1 does
not reopen it. Do not start M2 (uptake), M3 (AC source), D2, or Plan
Stage 3–9 work. Do not modify `BSimDaninoD1`, `BSimDaninoD1g`,
`BSimReservoirPlanStage*`, or `BSimReservoirStage11`.

## Why this is not Stage 11 and not the hybrid dish

`BSimReservoirStage11` coupled Monod to a glucose field with uptake, a
Model A replenisher, and a nutrient AC, using `reservoir_new`'s unitless
`mu_max = 1`, `K_G = 5`. That merge is closed as a reservoir path.

The main reservoir stays Plan Stage 6 PASS with the Stage 3B Hill receiver:
no glucose, no Monod. If M1 passes, still do not copy glucose into
`BSimReservoirPlanStage6`.

## Growth law (frozen, also declared in the Java header)

```
µ = µ_max * G / (K_s + G)
K_s    = 0.18 µM
µ_max  = Math.log(2.0) / 1800.0    // 30 min doubling, THIS rebuild
```

- `K_s = 0.18 µM` is Senn et al. 1994 (BNID 111049 range 0.18–0.55 µM).
- `µ_max = ln 2 / 1800 s⁻¹` matches Stage 10 / Plan `GROWTH_RATE` (30 min
  doubling). Do **not** use Senn's 0.92 h⁻¹ here — that would be a 45 min
  doubling and would disagree with the Plan dish.
- Do **not** use `reservoir_new`'s unitless `mu_max = 1`, `K_G = 5`.
- Print both `µ_max` and `K_s` at the start of each run.
- Do not retune `K_s` or `µ_max` after seeing `µ_hat`.

## How growth is applied (frozen)

Use BSim surface-area growth so N can increase.

```
GROWTH_RATE_SA = 4 * π / 1800
```

Same compile-time constant as the Plan dish (`4.0 * Math.PI / 1800.0`).
Scale the instantaneous surface-area rate by the Monod factor
`G / (K_s + G)` with `G` in µM from the bath. If `G` is uniform, every
cell sees the same `µ`.

- No clamp.
- Start `N = 50`.
- Replication on.
- Initial radii desynchronised with BSim `setRadius()` (surface area uniform
  in `[S(rep)/2, S(rep)]`) so some cells can divide inside 7200 s even at
  low `G`.
- No death, no vesicles, no chemotaxis goal field.

## Glucose units (frozen)

Store the field in µM-equivalent molecules, same conversion as D1 AHL:

```
1 µM = 602 molecules/µm³
G_uM = conc / 602
```

Set the field uniformly every tick after `field.update()`. Do not calibrate
a source rate. Do not subtract glucose when cells grow. Do not add a
replenisher. Do not add an AC source. Those are M2/M3.

Diffusivity and decay are unused because the bath is overwritten every tick.
They are set to 0 so the field cannot drift between `setConc` calls.

## Box (frozen)

- Domain 200 × 200 × 10 µm
- Glucose field grid 20 × 20 × 1
- `dt = 0.05` s
- `FLOW_SPEED = 0`
- No AHL, QS, Danino, acid, attractant, luminescence, vesicles, clamp,
  voxels, or 300 s RC windows
- Camera: same ortho / camera / perspective sequence as D1, this box size:
  `ortho(0, BOUND_X, BOUND_Y, 0, -1000, 10000);`
  `camera(BOUND_X/2, BOUND_Y/2, BOUND_Y, BOUND_X/2, BOUND_Y/2, 0, 0, 1, 0);`
  `perspective(PI/2, BOUND_X/BOUND_Y, 0.1, 10000);`
- Window 800 × 600
- Seed 101

## Production runs (frozen)

Three runs, seed 101, duration 7200 s each (2 h). Completeness is the last
CSV row `t = 7200`, not stdout.

| Run | G (µM) | vs K_s | expected µ/µ_max |
|---|---|---|---|
| M1a | 0.02 | 0.1 × K_s | 0.02 / 0.20 = 0.10 |
| M1b | 0.18 | 1 × K_s | 0.50 |
| M1c | 1.80 | 10 × K_s | 1.80 / 1.98 ≈ 0.909 |

## Export

Every 10 s, one CSV (semicolon) per run:

```
t_s;G_uM;N;run_label
```

Paths:

- `results/m1a_seed101/m1_timeseries.csv`
- `results/m1b_seed101/m1_timeseries.csv`
- `results/m1c_seed101/m1_timeseries.csv`

`G_uM` MUST be field/602. Last row `t = 7200`.

## Measured specific growth rate (frozen before seeing N)

On `t` in `[600, 7200]` (skip initial division lag):

```
µ_hat = (ln N_end - ln N_start) / (t_end - t_start)
```

Print `µ_hat` and `µ_hat / µ_max` for each run.

## Gates (all required)

1. `G_uM` is the commanded bath ± 1e-6 on every log row of that run.
2. `N` increases in every run. `N_final > N_initial`.
3. Ordered rates: `µ_hat(0.02) < µ_hat(0.18) < µ_hat(1.80)`.
4. `µ_hat(1.80) / µ_max ∈ [0.85, 1.00]`.
5. `µ_hat(0.18) / µ_max ∈ [0.40, 0.60]` (Monod at `K_s` is 0.5).
6. No uptake: `G` does not fall. Stage 11 / Plan Stage 6 not modified.
   No glucose in the hybrid dish.

If `N` climbs the same at all three `G`, Monod is not wired: **FAIL, stop**.
Do not then add uptake or an AC to “help”. Do not retune `K_s` or `µ_max`
after seeing `µ_hat`.

If M1 passes, stop. M2 (uptake) is a later prompt.

## Out of scope

Do not start M2, M3, D2, or Plan Stage 3–9 work. Do not copy glucose into
the Plan Stage 6 dish. Do not modify `BSimDaninoD1`, `BSimDaninoD1g`,
`BSimReservoirPlanStage*`, or `BSimReservoirStage11`.
