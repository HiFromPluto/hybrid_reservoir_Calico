# ChassisHybridmm NARMA I3n protocol

**Frozen before I3n Java exists and before any occupancy or NRMSE is observed.**

I3n asks one new-dish question: on a millimetre monolayer of Hertzian
*E. coli* rods, with HybridDish Hill \(R,L\) and the copied Narma10b
input, does an occupied driven 408-D readout beat Brownian density and a
silent dish?

This is not a rerun or rewrite of paper-2 Narma10b. Narma10b spheres in
the `1000×500×10 µm` claim dish remain Overall PASS at driven F408
`0.928`; those files and Overall lines are not edited. ChassisPocket I0c
proved Hertzian fill/jam at `W=20`; that fill result is not a NARMA
score. Packing is the body here. The task remains AHL → \(R\) → \(L\) →
ridge.

## Predeclared outcomes

- Driven occupancy `mean_R < 0.05`: **NOT_SCORED**. Stop. Do not raise
  `J_max` or lower `K`.
- Driven occupancy `mean_R > 0.95`: **SATURATED / NOT_SCORED**. Stop.
  Do not move `J_max`.
- Occupied and driven mean test NRMSE is lower than both Brownian and
  silent means across seeds `111/222/333`: **System PASS** for
  ChassisHybridmm I3n.
- Occupied and driven loses to Brownian or silent: **System FAIL** for
  this dish. No retuning.
- Driven versus field-only AHL is reported separately as the
  living-layer diagnostic. A field win is not permission to add ACs.

The `10`-tap linear map of \(u\), NRMSE `0.683`, is a task ceiling and
not a living-layer gate.

## Physical scenario freeze

```
1000 × 500 × 1 µm
FLOW = 0
NO_FLUX chemicals
CENTER AC at (500, 250, 0.5)
Hertzian rods; motility OFF
Hill R,L on each rod
408-D ridge on the Narma10b readout grids
```

| Item | Frozen value | Class | Notes |
|---|---:|---|---|
| Domain | `1000 × 500 × 1 µm` | ENGINEERING | Claim xy; Job 3c monolayer height, not claim height 10 µm |
| Flow | `0 µm/s` | TAKEN | Claim-dish flow 8 µm/s was occupancy DEAD |
| Chemical BC | NO_FLUX | ENGINEERING | `setSolid(true,true,true)` |
| Field dt | `0.05 s` | ENGINEERING | Narma10b |
| Growth/mechanics cadence | every `1.0 s` = 20 field ticks | ENGINEERING | `ChassisParameters.DT_S`; do not run Hertzian at 0.05 s |
| Nutrient | uniform `C_s=0.5 mM`, no PDE | TAKEN Warren | `elongateCited`; no sphere `4π/1800` growth |
| Warmup | `18000 s` at `u=0.5` | ENGINEERING | Reporter \(L\) warmup |
| Window / pulse / sample | `300 / 75 / 20 s`, 16 samples | ENGINEERING | Last sample `199;15;299.95` |
| Windows | `200` | ENGINEERING | washout 40 / train 110 / test 50 |
| NARMA input | copied `input_ahl_narma200.txt` | TAKEN | SHA-256 of canonical 12-decimal sequence below |
| Input range | `[0,0.5]` | TAKEN | Uniform `[0,1]` forbidden |
| Acid input | copied `input_acid_held05_200.txt` | TAKEN | Held 0.5; not a second NARMA wire |
| Field grid | `50 × 25 × 1` (`dx=20 µm`) | ENGINEERING | Narma10b |
| State readout | `20 × 10 × 1` | ENGINEERING | Per-window mean \(R\) and \(L\) |
| Death readout | `4 × 2 × 1` | ENGINEERING | Last sample, input-driven deaths |
| Ridge width | `200 R + 200 L + 8 deaths = 408` | ENGINEERING | No field or density columns in biology ridge |
| Hill \(R\) | `K=1.6 µM`, `n=2`, `tau_R=15 s` | TAKEN HybridDish | Plume occupancy, not LuxR EC50 |
| Reporter \(L\) | `alpha=delta=1/1500 s^-1` | TAKEN HybridDish | Per rod |
| AHL | `D=159 µm²/s`, `k=0.0033 s^-1` | HybridDish clocks | Not Dilanji and not Weber |
| `J_max` | `1.28e7 molecules/s` | ENGINEERING | Declared 1/10 volume match before occupancy |
| AHL AC | `(500,250,0.5)` | ENGINEERING | One point source, \(J=J_max u\) |
| Acid source | `(300,375,0.5)`, `2e10 molecules/s` at input 0.5 | ENGINEERING | HybridDish rate scaled 1/10 for the 1-µm field |
| Acid field | `D=200 µm²/s`, `k=0.0067 s^-1` | HybridDish clocks | Effective acid field |
| Acid cell production | `1e5 molecules/s` | ENGINEERING | HybridDish volumetric rate scaled 1/10 |
| Acid death | HybridDish pH map, `K_max=0.002 s^-1` | TAKEN HybridDish | No Biselli; do not tune |
| Clamp | `N=2000` | ENGINEERING / MODEL_CONVENIENCE | Do not tune to NRMSE |
| Clamp removal | per-field-tick HybridDish formula with `T_removal=1800 s` | ENGINEERING | `p=(dt/T_removal) 2^{-(1-N/2000)}`; not Biselli |
| Initial population | `1800` placed rods | ENGINEERING | Not grown from one founder |
| Initial cell-cycle phase | seeded uniform over one Warren bath division cycle | ENGINEERING | Avoids an artificial synchronous 1800-cell division; frozen before Java |
| Division | `SYMMETRY_BROKEN`; arm RNG seed | Job 3c | Seeds `111/222/333` |
| Hertzian | Job 3 `k_cc`, `k_ac` | TAKEN | `ValdezHertzian`; no retune |
| Contact broadphase | conservative `4.5 µm` centre cutoff; unchanged Valdez solve per connected component, rebuilt twice | ENGINEERING | Maximum possible contact reach is 4.0 µm at `L_div`; acceleration only |
| Motility / LuxI / glucose PDE / Danino / Weber / Biselli | OFF | — | No species crossing |

### Declared source scaling

The claim dish was 10 µm high. Sending the same molecules into a 1 µm
field would approximately multiply a well-mixed concentration by ten
and could saturate the Hill receiver. I3n therefore freezes
`J_max=1.28e8 × (1/10) = 1.28e7 molecules/s` before occupancy is seen.
HybridDish volumetric acid production rates are also scaled by `1/10`.
If occupancy is DEAD or SATURATED, that is the I3n finding; no source
hunt follows.

### Initial placement

Place all 1800 rods at `z=0.5 µm` on a sparse central `60×30` placement
grid with seeded jitter, seeded in-plane axes, and no initial overlap.
Seed each rod uniformly over one Warren bath cell-cycle phase so 1800
placed cells do not divide synchronously. Every rod imports
`BacteriumFromScratch.EcoliRodCell`; its body is not copied. Division
uses `SYMMETRY_BROKEN` with the arm RNG. Motility is off. The closed
1-µm z walls enforce the Job 3c monolayer.

The mechanics broadphase is numerical only. At `L≤L_div`, two rods
cannot contact when their centres are more than `4.0 µm` apart. Each
1-s mechanics tick conservatively joins rods within `4.5 µm`, solves
each connected component with the unchanged
`ValdezHertzian.relaxContacts`, rebuilds, and solves once more. No force,
stiffness, residual, or wall rule changes.

## Inputs and identity

Copy, do not regenerate:

- `examples/BSimReservoirPlanNarma10b/input_ahl_narma200.txt`
- `examples/BSimReservoirPlanNarma10b/input_acid_held05_200.txt`
- `examples/BSimReservoirPlanNarma10b/narma10_target.csv`

Canonical NARMA-input SHA-256, computed exactly as Narma10b by joining
the parsed values formatted to 12 decimals with commas:

`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`

NARMA target:

```
y[0] = 0
y[n+1] = 0.3 y[n] + 0.05 y[n] sum(i=0..9, y[n-i])
         + 1.5 u[n-9] u[n] + 0.1
```

Negative-index \(u,y\) values are zero. Window \(n\), sampled after
pulse \(u[n]\), predicts \(y[n+1]\).

## Arms

| Arm | Simulation | Ridge |
|---|---|---|
| Driven | Copied AHL command, acid held 0.5, warmup AHL 0.5 | 408 \(R/L/\)input-driven-death features |
| Brownian | Same driven chemical program; Narma10b passive Brownian-particle density null (no biology or growth) | 200 `Den_*` features only |
| Silent | AHL and acid sources off; warmup AHL 0 | 408 \(R/L/\)input-driven-death features |
| Field | No extra simulation; driven AHL voxels | `AHL_uM_*` only, reported diagnostic |

Required seeds are `111 / 222 / 333`, shared across arms. A seed is not
dropped after a bad NRMSE.

## Logging contract

- `window_summary.csv`: one row per completed analysis window.
- `results.csv` and `voxels.csv`: 16 rows per window, 3200 for a full
  run; final row begins `199;15;299.95`.
- `voxels.csv` exposes exactly named grids used by the checker:
  `Receiver_R_*` (200), `Lum_Mean_*` (200),
  `Input_Driven_Death_*` (8), `Den_*` (200), and `AHL_uM_*` (200,
  sampled on the state-readout grid as in Narma10b).
- State features are means of 16 intra-window samples. Death features
  are the final full-window accumulated values.
- Mechanics rows every `10 s` include `N`, `offplane`,
  `delta_cc_max_um`, `d_centers_min_um`, and in-plane nematic order.
- Run metadata labels the AHL and reporter numbers as HybridDish clocks,
  the source scaling as ENGINEERING, and this package as a new dish.

## Ridge freeze

Washout windows `0..39`; train `40..149`; test `150..199`. After
dropping washout, inner fitting uses rows `0..87` (windows `40..127`)
and lambda validation uses rows `88..109` (windows `128..149`) only.
Test rows never select lambda or standardize features.

- Training-only mean/std standardization; zero-variance columns stay 0.
- Bias column of ones; intercept unregularized.
- Lambda grid `{1e-6,1e-4,1e-2,1,1e2,1e4,1e6}`.
- Lowest validation NRMSE wins; exact ties take the larger lambda.
- Refit all 110 train windows, then score 50 test windows.
- NRMSE is test RMSE divided by population standard deviation of test
  labels.

## Pre-registered gates

| Gate | Pass | Stop/fail condition |
|---|---|---|
| I3n.0 occupancy | Driven seed 111 warmup+smoke `mean_R ≥ 0.05` and `≤0.95` | `<0.05` DEAD or `>0.95` SATURATED → NOT_SCORED, stop |
| I3n.1 monolayer | End-warmup `offplane ≤0.05` | Stacked → stop; do not raise height |
| I3n.2 input identity | Canonical SHA-256 matches Narma10b | Regenerated or mismatched input |
| I3n.3 system | Mean driven test NRMSE `<` Brownian and `<` silent across seeds 111/222/333 | Loses to either null → System FAIL |
| I3n.4 living layer | Report driven versus field; PASS if driven `<` field | Field wins is reported; no extra AC |
| I3n.5 ridge hygiene | 408 biology columns; lambda windows 128..149 only; intercept unregularized | Validation leak or feature leak |
| I3n.6 honesty | New ChassisHybridmm dish; Narma10b 0.928 untouched | “Paper 2 improved” or occupancy called a task score |

Packing gate: on logged mechanics rows with `N≥8`,
`delta_cc_max_um ≤0.25 µm`. End-warmup filament morphology also stops
the run: the predeclared distributed-population guard is in-plane
nematic order `≤0.95` (`1` means all rod axes are parallel). No
mechanics parameter is retuned.

## Occupancy-first execution order

1. Freeze this protocol before Java.
2. Copy all three input/target sidecars verbatim and verify the input
   hash.
3. Implement a new Java package using imported chassis classes and a
   new `check_i3n.py`; do not edit claim Java.
4. Run driven seed 111 through warmup only or warmup plus five windows.
   Inspect occupancy, monolayer, overlap, and filament morphology.
5. DEAD, SATURATED, stacked, over-overlapped, or filamentary: stop and
   write `CHASSISHYBRIDMM_NARMA_STANDING.md` as NOT_SCORED.
6. Only if the smoke is ALIVE, run all driven seeds, then Brownian and
   silent seeds. Field is read from driven voxels.
7. Run `python check_i3n.py`, write the standing honestly, and stop.

## Integrity fences

Do not edit `examples/HybridDish/bsim/BSimHybridDish.java`, any
`GATE_EVIDENCE.md`, Narma10b results or Overall lines, PocketDish /
PocketNeck claim Java, ChassisPocket I0/I0b/I0c Java or CSVs, or Job 3b
CSVs. Do not retune `k_cc`, `K`, `n`, `tau_R`, `tau_L`, AHL decay,
`J_max`, clamp, or acid `K_max` after occupancy or NRMSE. Do not add
LuxI, Danino 4-ODE, Weber, glucose PDE, motility, Biselli death, extra
ACs, Waveform, C1, E5, PocketDish-A, PocketOsc, or I1.

Do not claim spatial QS structure: this remains a `dx=20 µm`
millimetre dish, not Dilanji's 32-mm lane.

## Reference only: Narma10b spheres, not this run

| Arm | Test NRMSE |
|---|---:|
| Driven F408 | 0.928 |
| Field AHL | 1.029 |
| 10-tap of \(u\) | 0.683 |
| Brownian / silent | 1.162 |

I3n is not required to beat `0.928`. A worse absolute rod score can
still be System PASS if it beats this new dish's Brownian and silent
controls.

## Post-run clerical clarification

The first text freeze accidentally described the Brownian row as also
having “rod growth.” That conflicts with the copied Narma10b control
definition, which is explicitly passive particles with `Den_*` only.
The Java used the Narma10b passive density null from its first run. This
row is corrected transparently; no simulation parameter, output, ridge
feature, or score changed.
