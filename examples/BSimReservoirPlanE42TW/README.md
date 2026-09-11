# E4.2 — hypothetical two-way mediator + Turing/replay controls

Track E4, two-way mediator only. Frozen HybridDish / Narma10b /
E4.1 claim dish. `GATE_EVIDENCE.md` is not edited. Kinetics, layout,
flow, `K_AC`, and `Jmax_A1` are not retuned. Claim dish stays CENTER /
`FLOW=0`.

This package asks whether a payload-matched A1 AC that also senses a
cell-produced inert reporter P, and modulates AHL release, adds NARMA-10
capacity beyond (i) nonfunctional A1 and (ii) a replay of recorded P(t).

**MODEL_STATUS = `HYPOTHETICAL_DESIGN_ENVELOPE`.** There is no measured
P sensor curve or TX-TL device in the project. This is not a digital
twin. Lentini et al. 2017 is topology + Turing/replay controls only;
its µM / fold-change / hour numbers are not parameters.

## Status

**Module 0 complete.** No wet-lab P sensor curve, no TX-TL, no vesicle
geometry. Lentini PDF vendored
(`5c4beb444515ffc58ba01d333598ba5ddaa47e75447eaa0e0691ac2549411ab6`).
**MODEL_STATUS = `HYPOTHETICAL_DESIGN_ENVELOPE`.** Not a digital twin.

**Module 1 complete: PASS.** Unit tests A–F passed. `P=0` recovers A1
payload to relative error 0. `h(0)=1`, `h(K_P)=1.25`, `h(∞)=1.5`.

**Module 1.6 complete: STABLE_RESPONSIVE** at frozen `alpha=0.5`.
Pulsed constant-u and NARMA 0-D loop-gain did not saturate. Living
Java was authorized after this class. The one predeclared alpha cut
was not used.

**Module 2 seed 111 complete.** Occupancy **ALIVE**
(`mean_R=0.2214`, `r(R,u)=0.791`, `mean_h=1.1244`).
System **PASS** (CLOSED F408 `0.8906` vs Brownian `1.1625` and silent
`1.1622`). Living-layer **PASS** vs CLOSED field `0.9817`.
Transducer |CLOSED−A1| = `0.0101`: **NO_STORY_MOVE**.
Interaction: **NO_STORY_MOVE** (closed ties A1; replay ties closed at
`0.0002`). **SKIP_SHUFFLE.** Do not promote this dish. Seeds 222/333
were not started. Direct-input `0.6828` still wins. `0.93` stays weak.
**MODEL_STATUS remains `HYPOTHETICAL_DESIGN_ENVELOPE`.**

## Frozen numbers

A1 (untouched): `g0=0`, `n_AC=2`, `K_AC=0.25`,
`Jmax_A1=7.467750975821e7` molecules/s.

Feedback: `h(P)=1+alpha*P^2/(K_P^2+P^2)` with `n_P=2`, `alpha=0.5`
unless Module 1.6 takes the one predeclared cut to `0.25`.
`K_P` and `k_P` are frozen in Module 1 from well-mixed algebra, not
from NRMSE and not from Lentini.

`u` SHA-256
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.

Direct-input `0.6828` remains the task ceiling. `0.93` stays a weak
predictor. Do not promote this dish on a 0.03 tick.
