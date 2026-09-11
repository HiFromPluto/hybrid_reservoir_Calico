# E4.1 — one-way AHL AC transducer vs ideal current

Frozen before living NRMSE. Track E4 architecture, **one-way only**.
This package does not rewrite Stage 7 or Track B. `GATE_EVIDENCE.md`
files are not edited. K, n, tau_R, tau_L, clamp, mortality, flow, and
layout are not retuned. Claim dish stays CENTER / `FLOW=0`.

Two-way feedback, multi-AC patients, vesicle stores, glucose, Danino,
C1, Lorenz 202/303, and Stage99 are not started.

Do not copy `reservoir_new` `ArtificialCell` (store depletion +
saturating gate). Do not rename `addQuantity` and call it a digital
twin.

## Scientific question

If the ideal HybridDish current A0 \(J=J_{\max}u\) is replaced by one
static AHL gate A1 at matched cumulative payload, does the complete
one-way hybrid still carry NARMA-10 above Brownian and silent, and
does the living 408-D readout still beat field-only?

AC gating is legitimate hybrid computation. A field/AC-output win is
attribution, not “biology failed.” Direct-input `0.6828` remains the
task ceiling; biology is not required to beat it.

Two-way communication is not this experiment.

## MODEL_STATUS

**`HYPOTHETICAL_DESIGN_ENVELOPE`**

Module 0 found no measured steady \(u\to\)flux curve with units
(`results/AC_SPEC_INVENTORY.md`). A1 is therefore not
`CALIBRATED_A1` and is **not a digital twin**. Envelope `g0`, `n_AC`,
and `K_AC` were frozen before occupancy or NRMSE. `Jmax_A1` was
chosen only to match A0 commanded mass on the frozen Narma10b `u`
file. Do not fit `K_AC` to NARMA NRMSE. Do not change `K_AC` after
seeing occupancy.

`tau_on`/`tau_off` were not found. A2 is not implemented. Living BSim
is A1 only.

## Frozen dish and task

Copied from Narma10b. Mechanisms are not knobs.

- `1000 x 500 x 10` µm, `dt=0.05` s, field `50 x 25 x 1`, readout
  `20 x 10` + `4 x 2`
- warmup `18000` s, windows `300` s, AHL pulse `75` s, sampling every
  `20` s
- receiver `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL at `(500,250,5)`, acid at `(300,375,5)` held `0.5` in analysis
  windows, attractant AC0/AC2 silent
- clamp `K=2000`, `INITIAL_POP=1800`, `FLOW_SPEED=0`
- official ridge: 408 channels, closed NARMA split (washout 40 /
  train 110 / test 50, val windows 128–149 only)

Frozen `u` and target, **not regenerated**:

- `input_ahl_narma200.txt` from `examples/BSimReservoirPlanNarma10b/`
- `narma10_target.csv` from the same package
- acid `input_acid_held05_200.txt`

`u` SHA-256 (12-decimal payload, must match Narma10b):

`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`

Reuse, do not rerun:

- Narma10b driven seed 111 → A0 living reference
- Narma10b field-only from that run
- Narma10b silent seed 111
- Narma10b Brownian seed 111 (A0 pulses)

A1 changes \(J(u)\), so A1 driven and A1 Brownian are new BSim.

## A0 reference (not retuned)

\[
J(t)=J_{\max}u(t)
\]

during the 75 s pulse, 0 after. \(J_{\max}=1.28\times10^8\)
molecules/s, \(u\in[0,0.5]\), warmup command `0.5`.

## Frozen A1 envelope

\[
g(u)=g_0+(1-g_0)\frac{u^{n}}{K_{AC}^{n}+u^{n}},
\qquad
J(t)=J_{\max,A1}\,g(u(t))
\]

during the 75 s pulse, 0 after. No store variable.

| Parameter | Frozen value | Provenance |
|---|---|---|
| `MODEL_STATUS` | `HYPOTHETICAL_DESIGN_ENVELOPE` | Module 0: no measured curve |
| `g0` | `0` | envelope; not fitted |
| `n_AC` | `2` | envelope; not fitted |
| `K_AC` | `0.25` | envelope; not fitted |
| `g(0)` | `0` | follows from `g0=0` |
| `g(0.5)` | `0.8` | documented; \(0.5^2/(0.25^2+0.5^2)\) |
| `g(∞)` | `1` | documented saturation |
| `Jmax_A1` | `7.467750975821304e7` molecules/s | payload match only, relative error `1.62e-16` |
| vesicle store | none | forbidden |

`Jmax_A1` was not chosen from occupancy or NRMSE. If payload could
not be matched at this `K_AC`, stop; do not retune `K_AC`.

## Module 1 — cell-free (no bacteria)

```
python examples/BSimReservoirPlanE4AC/generate_a1_flux.py
```

Unit tests (all required before living Java):

| Test | Requirement |
|---|---|
| A | `u=0` ⇒ `J=Jmax*g0` |
| B | document `g(0.5)` and `g(∞)` |
| C | remaining + decay + boundary = injected on a no-bacteria field; leftover warmup is NARMA initial (E0.2 accounting) |
| D | A1 cumulative commanded mass vs A0 within 1% |
| E | no store; repeated `u` yields the same `J` |

Exports: gate, flux, cumulative mass, probes. Plots in
`results/figures/`. Report: `results/E4_AC_CALIBRATION.md`.

PASS Module 1 before any living Java.

## Module 2 — living scout (seed 111 only, if Module 1 PASSes)

New BSim:

- A1 driven seed 111
- A1 Brownian seed 111

Smoke first: print A1 parameters, `g(u)` at 0 and 0.5, `Jmax_A1`,
cumulative-mass match, dish CENTER, `FLOW=0`, finite non-negative
AHL. Not evidence.

CSV 200 / 3200 / 3200; last sample `199;15;299.95`.

Occupancy **DEAD** if `mean_R < 0.05` or `|r(mean_R,u)| < 0.5`. If
DEAD: keep the row, do not raise `Jmax` or `K_AC`, do not run
222/333.

Analysis, closed NARMA ridge, same target:

A0 (reused Narma10b seed 111):

- F408, field, Brownian, silent
- direct-u `0.6828`, informed `0.7071`
- masked/unmasked kinetic surrogates as in zero-sim

A1:

- F408, field-only from A1 voxels, A1 Brownian, reused silent
- same direct-u baselines (task unchanged)
- masked/unmasked surrogates on A1 AHL
- AC-output baseline: ridge on scalar `g[n]` / window-mean `J`, and
  on a 10-tap of `g` (legal AC-only delay line)

Claims (report, do not retune):

- **System:** A1-driven F408 beats A1-Brownian and silent.
- **Living-layer:** A1 F408 beats A1 field. If not, the gate/plume
  already has the task.
- **Transducer:** compare A1 F408 vs A0 F408. A 0.03 tick around
  0.93 is `NO_STORY_MOVE`. Do not promote A1 as a new claim dish on
  that tick.
- **Ceiling:** neither A0 nor A1 is required to beat `0.6828`.
  `0.93` stays a weak predictor.

Seeds 222/333 only if A1 is ALIVE AND (system pass) AND
(`|A1_F408 − A0_F408| ≥ 0.03` OR A1 living-layer sign flips vs A0).
Otherwise stop. No best seed.

## Forbidden

- two-way / cell→AC
- acid or attractant as extra analog wires
- 5-AC patient layout
- fitting `n`, `K_AC`, `g0`, `Jmax` to test NRMSE
- calling this envelope a digital twin
- vesicle `CS_in`

## Standing claims this package must not rewrite

- C1 remains **DEFER**
- E0.3 remains **NO_STORY_MOVE**
- Waveform1 FAIL stays
- Waveform2c field ≥ driven stays
- L3 living-layer FAIL stays
- E0–E3 standing claims in
  `examples/HybridDish/HIGH_TIER_SUBMISSION_EVIDENCE_PLAN.md` stay
