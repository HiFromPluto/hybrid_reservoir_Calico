# BSimReservoirPlanStage5 mechanism manifest

## Scope

This package copies frozen `examples/BSimReservoirPlanStage4` and adds a
per-cell luminescence state. Stage 4 is not modified and remains **FAIL**.
The 2x `Total_Deaths` fold is not rewritten. The clamp, acid source rate,
`K_MAX`, and Stage 3B receiver (`K=1.6`, `n=2`, `tau=15`) are not retuned.
Stage 6 is not started.

## Frozen feature contract

USE as input-tracking / analysis channels after collinear deletions:

- `Mean_q` (frozen Stage 3B receiver)
- `Input_Driven_Deaths`
- `Mean_L`

EXPORT but do not claim as input channels:

- `Births`, `Total_Deaths`, `Clamp_Deaths`, `OOB_Deaths`, `Population`, `Lum_Sum`
- Window `Extracellular_AHL_uM_Mean`, `Fraction_q_gt_0_5`, `pH_Mean` (diagnostic / occupancy / field; dropped from analysis for `|r|>0.9`)
- Voxel AHL, pH, density, receiver, occupancy, Lum mean/sum, and split births/deaths

## Luminescence

New per-cell state `L`, initialized at 0. The 18000 s warmup sets it.
Euler update each tick:

`L += (ALPHA_LUX * R - DELTA_LUX * L) * dt`

`DELTA_LUX = 1/1500 s^-1`, so `tau_L = 1500 s ≈ 5 windows`.
Analytic `t95_L = -ln(0.05)/DELTA_LUX = 4493.599 s` (not measured).
`ALPHA_LUX = 1/1500` so steady-state `L → R`, putting mid-input mean L in
the same usable `(0,1)` band as R. `L` is copied to daughters as an
intensive luciferase level, not halved.

This is not a relabel of `R`, not a boolean of `R>0.5`, and not Danino.

## Dropped columns

- Voxel `Att_*` and `Rep_*`: silent AC0/AC2 architecture, dead fields.
- AHL molecules/um3: linear duplicate of `Extracellular_AHL_uM_Mean`.
- `Receiver_R_Mean`: identical to `Mean_q`.
- `Acid_mM_*`: affine duplicate of pH (`pH = 7.1 - acid_mM/2`).
- Unsplit voxel `Death_*`: replaced by clamp / input-driven / OOB.

Glucose, pH taxis, pH-dependent growth, Danino ODE, and Stage 7 layout
changes are not added.

## Protocol

Same Stage 4/3B timing and the same AHL/acid holdout permutations and
seeds `101/202/303` plus acid source-off seed `404`. AHL drive stays on
in the acid-off control.

## Gates

All required. Production holdouts passed. See `results/GATE_EVIDENCE.md`.
Stage 4 remains labelled FAIL. Stage 6 was not started.

- Analysis channels `Mean_q`, `Input_Driven_Deaths`, `Mean_L` have non-zero
  variance and no pair with `|r|>0.9`.
- Receiver `r=0.787/0.819/0.816`, mid occupancy `0.333`, ordered.
- `Input_Driven_Deaths` `r=0.890/0.930/0.882`.
- `|r(Mean_L, Mean_q)|=0.193/0.400/0.383`; mid-input mean L `0.403`.
- CSV 40/640/640 rectangular.
