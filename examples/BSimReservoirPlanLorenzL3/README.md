# Lorenz L3 — sample-clock screen, then at most one living scout

Track E3. Frozen HybridDish / BenchA claim dish. BenchA Lorenz FAIL
and BenchA2 SKIP=50 FAIL are not rewritten. C1 remains DEFER. E0.3
remains NO_STORY_MOVE. Waveform1 FAIL and Waveform2c field≥driven
stay as written.

This package asks whether an intermediate Lorenz sample clock exists
between the trivial Δt=0.02 map and SKIP=50 destruction, and if so
whether the living 408-D readout beats the legal linear baselines on
that clock.

## Status

**Module 1 complete.** Anchor sanity PASS. k=1 TRIVIAL, k=50
DESTROYED. Selected living scout: **k=10, Δt=0.20, AUTO_X**.

**Module 2 seed 101: living-layer FAIL.** Occupancy ALIVE
(`mean_R=0.230`, `r(R,u)=0.897`). Driven F408 test NRMSE `0.8256`
does not beat the legal AR `0.8233`. It does beat Brownian and silent
(`1.0076`). Field-only `0.7973` is below driven. Seeds 202/303 were
not run. Kinetics were not retuned. BenchA / BenchA2 Overall lines
were not edited.

See `PROTOCOL.md`, `results/LORENZ_L3_CLOCK_SCREEN.md`, and
`results/LORENZ_L3_SCOUT.md`.

## Generate then screen

From the repo root:

```
python examples/BSimReservoirPlanLorenzL3/generate_lorenz_l3.py
python examples/BSimReservoirPlanLorenzL3/screen_lorenz_l3.py
```

Do not retune K, n, tau, rate, clamp, mortality, flow, or layout.
Claim dish stays CENTER / FLOW=0. One AHL channel. Do not inject y or
z. Do not use acid as a Lorenz channel.

See `PROTOCOL.md` and `results/LORENZ_L3_CLOCK_SCREEN.md`.
