# HybridDish protocol (frozen)

Copied from the working NARMA-10 dish. Mechanisms are not knobs.

## Dish

- Bounds `1000 x 500 x 10` µm, `dt = 0.05` s
- Field grid `50 x 25 x 1`, readout `20 x 10` and `4 x 2`
- Warmup `18000` s, windows `300` s, AHL pulse `75` s
- Sampling every `20` s for the whole window (16 samples)
- Initial population `1800`, clamp `K = 2000`
- AHL AC at `(500, 250, 5)`, acid AC at `(300, 375, 5)`
- Attractant ACs silent; `FLOW_SPEED = 0`

## Receiver and light

```
dR/dt = (C^n / (K^n + C^n) - R) / tau
  K = 1.6 µM, n = 2, tau = 15 s

dL/dt = (R - L) / 1500
```

## Inputs

- NARMA-10: `input/ahl_narma10.txt`, `u ~ Uniform[0, 0.5]`
- Acid: `input/acid_held05.txt`, every window `0.5`
- Target: `input/narma10_target.csv`

```
y[0] = 0
y[n+1] = 0.3 y[n] + 0.05 y[n] * sum_{i=0..9} y[n-i]
         + 1.5 u[n-9] u[n] + 0.1
```

Window `n` state predicts `y[n+1]`.

## Readout

| Arm | Features |
|---|---|
| Driven / silent | `Receiver_R` 20×10 mean, `Lum_Mean` 20×10 mean, `Input_Driven_Death` 4×2 last |
| Brownian | `Den` 20×10 mean |

Ridge: bias column of ones, intercept not regularized, features standardized
on train only, lambda grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`, NRMSE
uses population std. Validation slice is train windows 128–149 only.
