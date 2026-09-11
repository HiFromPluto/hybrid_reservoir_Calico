# PocketOsc-SI — bulk delay-DDE, period vs μ (job protocol)

Frozen 2026-08-21 **before traces**. Architecture:
[`examples/PocketDish/POCKETOSC_SI_FROZEN_BUILDER_PROMPT.md`](../PocketDish/POCKETOSC_SI_FROZEN_BUILDER_PROMPT.md).
Equations and TAKEN table copied verbatim from
[`examples/PocketDish/DANINO_SI_MODEL.md`](../PocketDish/DANINO_SI_MODEL.md).
Job 3 CSTR identity is **DEAD**
([`examples/PocketDish/FLOW_STANDING.md`](../PocketDish/FLOW_STANDING.md)).
D50 is a different ODE (port check only; `danino_rhs` is **not** cloned).

No NARMA, GFP maps, consortium, `QS_KMLA` hunt, D1g QS_*, HybridDish,
`GATE_EVIDENCE.md`, Groisman+AC Java, \(\mu_{\mathrm{bus}}\) formula,
1-D \(N=200\) array, 55/90 min fit. Time stays **SI-scaled**. No
`TIME_ADJ=60`. Readout is **\(I\)** (LuxI), not \(H(\mathrm{LA})\).

## Question

SI delay-DDE, TAKEN table, bulk (\(D_1=0\)), \(d=0.5\): on the
**predeclared** \(\mu\) line below, do we get oscillations on a finite
interval, with \(T\) increasing in \(\mu\) (main text Fig. 4b/c)?

Hypothesis: unknown. Do not retune \(\alpha,\tau,\gamma_\ast\).

## TAKEN equations (SI “Modeling”; bulk)

Delayed production \(H_\tau(t)=H_i(t-\tau)\):

\[
P(\alpha,\tau)=\delta+\frac{\alpha H_\tau^{2}}{1+k_1 H_\tau^{2}}.
\]

\[
\begin{aligned}
\partial_t A &= C_A\bigl[1-(d/d_0)^4\bigr]P
  -\frac{\gamma_A A}{1+f(A+I)},\\
\partial_t I &= C_I\bigl[1-(d/d_0)^4\bigr]P
  -\frac{\gamma_I I}{1+f(A+I)},\\
\partial_t H_i &= \frac{b I}{1+k I}
  -\frac{\gamma_H A H_i}{1+g A}+D(H_e-H_i),\\
\partial_t H_e &= -\frac{d}{1-d}D(H_e-H_i)-\mu H_e.
\end{aligned}
\]

\(D_1\partial_x^2 H_e\) is **dropped** (bulk). Ordinary DDE.

## TAKEN scaled table

| Symbol | Value | Class |
|---|---:|---|
| \(C_A\) | 1 | TAKEN |
| \(C_I\) | 4 | TAKEN |
| \(\delta\) | \(10^{-3}\) | TAKEN |
| \(\alpha\) | 2500 | TAKEN |
| \(\tau\) | 10 | TAKEN (scaled time) |
| \(k\) | 1 | TAKEN |
| \(k_1\) | 0.1 | TAKEN |
| \(b\) | 0.06 | TAKEN |
| \(\gamma_A\) | 15 | TAKEN |
| \(\gamma_I\) | 24 | TAKEN |
| \(\gamma_H\) | 0.01 | TAKEN |
| \(f\) | 0.3 | TAKEN |
| \(g\) | 0.01 | TAKEN |
| \(d_0\) | 0.88 | TAKEN |
| \(D\) (membrane) | 2.5 | TAKEN |
| \(d\) | 0.5 | TAKEN (Fig. 4b) |
| \(\mu\) | **varied** (line below) | TAKEN as a free rate; **no** \(\mu(v)\) |
| \(D_1\) | 0 | bulk drop |

Integrator: scipy method of steps, constant delay \(\tau=10\), BDF,
`rtol=1e-7`, `atol=1e-9`, `max_step=0.1`. Equivalent to `jitcdde`
for this constant-delay bulk DDE. No living Java.

## μ line (frozen; do not add after occupancy)

SI **varies** \(\mu\); it does not list Fig. 4b ticks. Do **not**
use \(\mu_{\mathrm{bus}}=0.0208~\mathrm{s}^{-1}\). Covering slice of
Fig. 4c, labelled COVER except \(d\):

| Arm | \(d\) | \(\mu\) | \(D_1\) | Role |
|---|---|---|---|---|
| `MU_0` | 0.5 | 0 | 0 | COVER; closed analog |
| `MU_025` | 0.5 | 0.25 | 0 | COVER; low flow |
| `MU_05` | 0.5 | 0.5 | 0 | COVER; mid |
| `MU_1` | 0.5 | 1.0 | 0 | COVER; higher |
| `MU_2` | 0.5 | 2.0 | 0 | COVER; high |
| `MU_CSTR_MIN` | 0.5 | 1.25 | 0 | ENGINEERING diagnostic only: Job 3 \(0.0208~\mathrm{s}^{-1}\) as min⁻¹. **Not** identity. |

Do not add \(\mu\) points after occupancy. Do not fit 55/90 min.
Do not retune from `MU_CSTR_MIN`.

## ICs (ENGINEERING; SI silent)

\(A=I=0\), \(H_i=H_e=1\), history \(H_i(t<0)=1\). Print them. Do
not change after traces.

## Drive and period rule (frozen before traces)

- Integrate **1000** SI time units. Sample **0.5**.
- Discard first **200** before means and period.
- Period from **\(I\)** peaks: local maxima after discard, minimum
  spacing **5**, prominence \(\ge 0.10\times(\max I-\min I)\) on the
  post-discard window (relative). Period = mean inter-peak interval.
  Fewer than **3** peaks → `NO_PERIOD`.
- `OSC` if a period exists. Else not `OSC`.
- Print mean \(I,A,H_i,H_e\) (post-discard) and period.

Identity **only if** at least **two** COVER arms are `OSC`: \(T\)
strictly increases with \(\mu\) on those OSC COVER arms (Fig. 4b
direction). Absolute 55/90 min is **not** a gate (units are scaled).
`MU_CSTR_MIN` is never scored for identity.

`LIVING_DECISION=STOP_AFTER_SI_BULK_DEAD` if no COVER arm oscillates.
`LIVING_DECISION=STOP_AFTER_SI_NO_IDENTITY` if they oscillate but
\(T\) does not lengthen with \(\mu\) (or fewer than two COVER `OSC`).
`LIVING_DECISION=STOP_AFTER_SI_BULK_ALIVE` if identity holds —
**stop**. No spatial add-on. No Fig. 2c. No NARMA.
