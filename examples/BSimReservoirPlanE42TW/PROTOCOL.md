# E4.2 — hypothetical two-way mediator + Turing/replay controls

Frozen before occupancy or NRMSE. Track E4 architecture, **two-way
mediator only**. This package does not rewrite Stage 7, Track B, E4.1,
Narma10b, Waveform, or Lorenz. `GATE_EVIDENCE.md` files are not edited.
K, n, tau_R, tau_L, clamp, mortality, flow, layout, K_AC, and Jmax_A1
are not retuned. Claim dish stays CENTER / `FLOW=0`.

Do not copy `reservoir_new` `ArtificialCell` (store depletion +
saturating gate). Do not implement vesicles, `CS_in`, glucose,
Danino/LuxI AHL from cells, AiiA, A2, multi-AC patients, acid-as-feedback,
C1, Lorenz 202/303, Stage99, or E4.1 seeds 222/333.

Do not call this a digital twin.

## Scientific question

On the frozen NARMA-10 dish, if a payload-matched A1 AC also senses a
cell-produced inert reporter P and modulates AHL release, does the
closed loop add capacity that is not already in:

1. the one-way A1 hybrid (nonfunctional AC: sense P, do not send);
2. a replay of recorded P(t) into that same AC (broken feedback).

System claim: closed-loop F408 beats Brownian and silent.
Living-layer claim: closed-loop F408 beats field-only.
Interaction claim (the Lentini one): closed-loop F408 beats
nonfunctional A1 AND beats replay.

If replay ties the closed loop, cells are not writing anything the AC
needs. If closed ties A1, there is no two-way on this task.

Direct-input `0.6828` remains the task ceiling. Biology is not required
to beat it. A1 10-tap of `g` at `0.6535` is an AC-output result, not a
living-layer win. Do not fit anything to NRMSE.

## MODEL_STATUS

**`HYPOTHETICAL_DESIGN_ENVELOPE`**

No wet-lab P sensor curve, no TX-TL, no vesicle geometry. Lentini et al.
(*ACS Cent. Sci.* 2017, 3, 117–123, DOI 10.1021/acscentsci.6b00330) is
**topology + Turing/replay controls only**. It is not a parameter source.
Numbers were not copied from that paper into `K`, `n`, `tau`, `K_AC`,
`K_P`, `k_P`, or `alpha`. Lab measurement is not available.

## Topology (Lentini mediator, not same-channel QS)

```
u  →  g(u)           frozen A1 Hill gate
P  →  h(P)           new AC sensor, local P at the AHL source voxel
J  =  Jmax_A1 * g(u) * h(P)

AHL field → living Hill receiver R, L   (frozen)
living R  → secrete P                   (new, labelled engineering)
P field   → AC sensor only
P does not activate the AHL Hill receiver.
P is not acid. P is not a second QS HSL. P is not L.
```

Open-loop nested identity, required:

`h(0) = 1` ⇒ `P=0` recovers A1 exactly, including `Jmax_A1`.

Nonfunctional AC (Turing control) = A1 already on disk:

`J = Jmax_A1 * g(u)` (ignore live P)

Do not close the loop on the same AHL the receiver already sees.

## Frozen dish and task

Copied from E4AC / Narma10b. Mechanisms are not knobs.

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

- `input_ahl_narma200.txt` from `examples/BSimReservoirPlanE4AC/`
  (originally Narma10b)
- `narma10_target.csv` from the same package
- acid `input_acid_held05_200.txt`

`u` SHA-256 (12-decimal payload, must match):

`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`

Reuse, do not rerun, do not overwrite:

- Narma10b driven/field/Brownian/silent seed 111 (A0)
- E4AC A1 driven/Brownian seed 111 (nonfunctional AC)
- frozen u / target / acid files from E4AC/Narma10b

A1 driven seed 111 IS the nonfunctional AC. Brownian secretes `P=0`,
so functional Brownian ≡ A1 Brownian.

## Frozen A1 envelope (untouched)

\[
g(u)=\frac{u^{2}}{K_{AC}^{2}+u^{2}},\qquad
J_{A1}(t)=J_{\max,A1}\,g(u(t))
\]

during the 75 s pulse, 0 after. No store variable.

| Parameter | Frozen value | Provenance |
|---|---|---|
| `MODEL_STATUS` | `HYPOTHETICAL_DESIGN_ENVELOPE` | no measured curve |
| `g0` | `0` | E4.1 envelope; not fitted |
| `n_AC` | `2` | E4.1 envelope; not fitted |
| `K_AC` | `0.25` | E4.1 envelope; not fitted |
| `Jmax_A1` | `7.467750975821e7` molecules/s | E4.1 payload match; not retuned |
| vesicle store | none | forbidden |

## Frozen two-way envelope

\[
h(P)=1+\alpha\frac{P^{n_P}}{K_P^{n_P}+P^{n_P}},\qquad
J(t)=J_{\max,A1}\,g(u)\,h(P)
\]

during the 75 s pulse, 0 after. `h(0)=1` recovers A1.

| Parameter | Frozen value | Provenance |
|---|---|---|
| `n_P` | `2` | envelope; not fitted |
| `alpha` | `0.5` initially, so `h ∈ [1, 1.5]` | envelope; one predeclared cut to `0.25` is allowed only if Module 1.6 is SATURATED/UNSTABLE |
| `K_P` | `1.6` µM | same number as receiver `K`, different species; not from Lentini |
| `k_P` | freeze in Module 1 | well-mixed algebra in `results/p_scale_freeze.json` |
| `D_P` | `159` µm²/s | AHL-like inert reporter |
| `decay_P` | `0.0033` s⁻¹ | AHL-like inert reporter |
| `tau_I` | not introduced | if the loop is invisible vs A1, that is a finding |

P secretion (fast engineering envelope, **not** TX-TL hours): each
living cell secretes `k_P * R` molecules/s into the P field. No
intracellular synthase ODE. Brownian particles secrete 0. Dead /
absent cells secrete 0.

P concentration conversion is identical to AHL:
`µM = getConc / 602.2`.

No vesicle store. Repeated `(u,P)` yields the same `J`.

## Module 0 — inventory and fences

Write `PROTOCOL.md`, `README.md`, `results/AC_TW_INVENTORY.md`.
Vendor Lentini PDF to `examples/HybridDish/external_analysis/oc6b00330.pdf`
and record SHA-256 here. Copy u, target, acid from E4AC. Do not
regenerate.

Forbidden Java strings: glucose, Danino, `CS_in`, vesicle store,
AiiA, LuxI as bacterial AHL source.

## Module 1 — cell-free (no bacteria)

```
python examples/BSimReservoirPlanE42TW/generate_tw_flux.py
```

PASS before any living Java. No living BSim.

Unit tests (all required):

| Test | Requirement |
|---|---|
| A | `P=0` ⇒ `J_TW(u) ≡ J_A1(u)` for the frozen Narma10b `u`; payload relative error ≤ `1e-12` vs A1 commanded mass |
| B | `h(0)=1`, `h(∞)=1+alpha`, `h(K_P)=1+alpha/2` |
| C | P clamp: increasing prescribed P at fixed u increases J, never decreases it, saturates at `1+alpha` |
| D | no store: repeated `(u,P)` yields the same J |
| E | no cells ⇒ P field stays 0; AHL mass identity as in E4.1 (remaining + decay + boundary = injected; leftover warmup is NARMA initial) |
| F | P is absent from the AHL Hill receiver |

Exports: `g(u)`, `h(P)`, `J(u,P)`, cumulative mass, probes.
Report: `results/E4_2_CALIBRATION.md`.
Freeze: `results/tw_frozen.json`.

Do not fit `K_P`, `k_P`, or `alpha` to NARMA.

## Module 1.6 — 0-D loop-gain screen (no BSim)

Well-mixed ODE using the frozen envelope, A1 mean occupancy as the
scale already used for `k_P`, and constant `u` in `{0, 0.25, 0.5}`
plus the frozen NARMA `u` as an open drive.

Classify: `DEAD` / `STABLE_RESPONSIVE` / `SATURATED` (mean R > 0.9
pinned) / `UNSTABLE`.

If SATURATED or UNSTABLE: one predeclared cut `alpha 0.5 → 0.25`,
re-freeze in `tw_frozen.json`, re-run Module 1 A–F. Do not iterate
further. If still SATURATED/UNSTABLE/DEAD: STOP. Keep the row.
No living Java.

If STABLE_RESPONSIVE: proceed.

## Module 2 — living scout, seed 111 only

Smoke first (not evidence): print `tw_frozen.json`, P=0 identity,
`h(K_P)`, `alpha`, dish CENTER, `FLOW=0`, finite non-negative AHL
and P.

New BSim, sequential:

1. CLOSED driven seed 111. Live P → AC; log P at the source voxel
   every 1 s to `results/e42_closed_seed111/p_ac_timeseries.txt`.
2. REPLAY driven seed 111. Living cells still secrete P, but AC uses
   the recorded `P_ac(t)` from (1), time-aligned, not live P.
   Occupancy may differ; that is the point.

Do not rerun A1. Do not rerun Narma10b silent/Brownian.

CSV 200 / 3200 / 3200; last sample `199;15;299.95`.
AHL and P finite and non-negative.

Occupancy DEAD if `mean_R < 0.05` or `|r(mean_R,u)| < 0.5`.
If DEAD: keep the row, do not raise `k_P`, `alpha`, or `Jmax`, do not
run replay if closed is DEAD, do not run shuffle / 222 / 333.

Analysis, closed NARMA ridge, same target, same split.

Reuse:

- A0 F408 / field / Brownian / silent
- A1 F408 / field / Brownian / silent / 10-tap g
- `LINEAR_U_DELAY_10 = 0.6828`
- `NARMA_INFORMED_INPUT = 0.7071`

New:

- CLOSED F408, field-only, occupancy-masked RL surrogate
- REPLAY F408, field-only
- AC-output: scalar `g[n]*h[n]`, window-mean J, 10-tap of `(g*h)`
- report `mean_P_ac`, `mean_h`, `mean_R`, `r(R,u)`, `r(P_ac,u)`

Claims (report, do not retune):

- **System:** CLOSED F408 vs Brownian and silent
- **Living-layer:** CLOSED F408 vs CLOSED field
- **Transducer:** CLOSED F408 vs A1 F408. `|Δ| < 0.03` in the 0.93
  band is `NO_STORY_MOVE`
- **Interaction:** CLOSED F408 vs A1 F408 (nonfunctional) AND CLOSED
  F408 vs REPLAY F408. Interaction PASS only if closed beats both.
  If closed beats A1 but replay ties closed: `FEEDFORWARD_SCHEDULE`,
  not two-way. If closed ties A1: `NO_STORY_MOVE`, two-way
  untested-as-useful.
- **Ceiling:** not required to beat `0.6828`. `0.93` stays weak.

Time-shuffle BSim is authorized only if interaction is not
`NO_STORY_MOVE` (i.e. `|CLOSED−A1| ≥ 0.03` or `|CLOSED−REPLAY| ≥ 0.03`
on F408). Shuffle: permute the 200 analysis-window `P_ac` blocks from
the closed recording; keep warmup P; living cells on; AC uses shuffled
P. One seed 111. If authorized, do it. If not, write `SKIP_SHUFFLE`
with the deltas.

Seeds 222/333: not in this prompt.

Do not promote the two-way dish as a new claim dish on a 0.03 tick.

## Forbidden

- fitting `n_P`, `K_P`, `k_P`, `alpha`, `Jmax` to test NRMSE
- retuning receiver `K`, `n`, `tau` after scores
- same-AHL QS loop / bacterial LuxI
- acid or attractant as the feedback species
- vesicle stores, TX-TL ODEs, Lentini 10 µM / 5 h numbers
- calling `MODEL_STATUS` a digital twin
- overwriting E4AC or Narma10b voxels
- rewriting standing claims, C1 DEFER, E0.3 NO_STORY_MOVE,
  Waveform1 FAIL, Waveform2c field ≥ driven, L3 living-layer FAIL

## Standing claims this package must not rewrite

- C1 remains **DEFER**
- E0.3 remains **NO_STORY_MOVE**
- Waveform1 FAIL stays
- Waveform2c field ≥ driven stays
- L3 living-layer FAIL stays
- E4.1 A1 remains the one-way seed-111 scout; not a new claim dish
- E0–E3 standing claims in
  `examples/HybridDish/HIGH_TIER_SUBMISSION_EVIDENCE_PLAN.md` stay
- two-way communication was untested before this package; this scout
  does not rewrite that standing sentence in the evidence plan
