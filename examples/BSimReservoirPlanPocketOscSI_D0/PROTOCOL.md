# D0 protocol — published bulk delay-DDE, Figure 4

Frozen **before traces**. Equations:
[`EQUATIONS.md`](EQUATIONS.md). Prior twin
`examples/BSimReservoirPlanPocketOscSI/` is **not** this oracle.

No NARMA, GFP maps, Java, spatial \(N=200\), D1g, \(\mu_{\mathrm{bus}}\),
`QS_KMLA`, HybridDish, or parameter retune after seeing periods.

## Scientific question

At \(d=0.5\), bulk \(D_1=0\), TAKEN table, SI Hill \(P\): does the
published SI DDE oscillate, and does period increase with \(\mu\) over
the published flow/\(\mu\) interval, consistent with Fig. 4b/c?

Main-text Fig. 4b (PMC author manuscript caption):

> “Period of oscillations as a function of the flow rate \(\mu\) at cell
> density \(d=0.5\) (top panel).”

Fig. 4c:

> “Oscillations occur over a finite range of cell densities, and period
> increases with \(\mu\) after the bifurcation line is crossed.”

Experimental comparison class (main text): channel flow 180–296 µm/min
gave periods 52–90 min; high flow 90±6 min, low flow 55±6 min. Absolute
experimental minutes are **not** a fudge target. Model periods are
SI-scaled; the published **model class** is about 48–88
scaled/min-comparable units as used in Fig. 4b.

## TAKEN table (do not retune)

| Symbol | Value | Provenance |
|---|---:|---|
| \(C_A\) | 1 | SI “most of our simulations” |
| \(C_I\) | 4 | SI |
| \(\delta\) | \(10^{-3}\) | SI |
| \(\alpha\) | 2500 | SI |
| \(\tau\) | 10 | SI (scaled time) |
| \(k\) | 1 | SI |
| \(k_1\) | 0.1 | SI |
| \(b\) | 0.06 | SI |
| \(\gamma_A\) | 15 | SI |
| \(\gamma_I\) | 24 | SI (page-break continuation) |
| \(\gamma_H\) | 0.01 | SI |
| \(f\) | 0.3 | SI |
| \(g\) | 0.01 | SI |
| \(d_0\) | 0.88 | SI |
| \(D\) (membrane) | 2.5 | SI |
| \(d\) | 0.5 | Fig. 4b caption |
| \(D_1\) | 0 | bulk drop |
| \(\mu\) | varied | frozen grid below |

## \(\mu\) grid (frozen; do not add after occupancy)

The SI does not list Fig. 4b tick values. Identity line follows the
predeclared Fig. 4 operating interval around \(\mu=1\)–\(2\):

\(\mu = 1.0, 1.1, \ldots, 2.0\).

Predeclared extras to report the finite operating interval (Fig. 4c
bifurcation / absence outside a band), **frozen now**:

\(\mu \in \{0, 0.25, 0.5, 0.75, 2.25, 2.5, 3.0, 4.0, 5.0\}\).

Do not add points after occupancy. Do not convert from Pocket Job-3
\(\mu_{\mathrm{bus}}\).

## History / IC ensemble (frozen; `PREDECLARED_UNSTATED_ICS`)

The SI does not state bulk \(A,I,H_i,H_e\) or \(H_i(t<0)\). No
Danino-group code release is used. Ensemble frozen before traces.
Primary is **not** chosen after seeing periods.

All arms: \(d=0.5\), constant delay history.

| id | \(A,I,H_i,H_e\) at \(t=0\) | \(H_i(t<0)\) | Role |
|---|---|---|---|
| `SI_BASAL_PERTURB` | \(0, 1, 0, 0\) | \(0\) (equals state) | **PRIMARY.** SI Fig. 6 analogue \(I_{N/2}=1\), \(A=0\); unstated AHL set to basal 0 |
| `SI_BASAL_REST` | \(0, 0, 0, 0\) | \(0\) | low-AHL rest; SI “basal state with \(A=I=0\)” |
| `MODERATE_AHL` | \(0, 1, 0.2, 0.2\) | \(0.2\) (equals state) | moderate-AHL start; constant history = IC |
| `HISTORY_OFFSET` | \(0, 1, 0, 0\) | \(0.2\) | constant history **different** from \(t=0\) |
| `FAILED_PRIOR_IVP` | \(0, 0, 1, 1\) | \(1\) | negative control; PocketOsc-SI ENGINEERING start |

Do not hunt ICs after failure. Do not promote a post-hoc history into
the primary arm.

## Integrator (frozen)

- Constant delay \(\tau=10\).
- Method of steps: scipy `solve_ivp` BDF on successive intervals of
  length \(\tau\), dense output as the delay interpolant.
- Delay lookup **raises** if \(t-\tau\) falls in a gap. History for
  \(t-\tau\le 0\) is the frozen constant. No silent skip.
- Nominal tolerances: `rtol=1e-8`, `atol=1e-10`, `max_step=0.05`.
- Tightness check: rerun primary \(\mu=1.5\) at `rtol=5e-9`,
  `atol=5e-11` (2× tighter). OSC/NO_PERIOD must match; period relative
  change \(\le 0.02\).
- No hidden clipping. Concentrations are not projected during
  integration. Production uses \(H_\tau\) as returned by the
  interpolant. If any sampled coordinate is negative, the run is
  recorded as `NEGATIVE_STATE` and is not scored as OSC.
- Deterministic: same frozen JSON + same SciPy → identical samples.

`jitcdde` / `ddeint` are not installed in this environment. Method of
steps with BDF is the documented solver.

## Horizon, sampling, peaks (frozen)

Justified from the 48–88 class and SI Fig. 5 (spatial plots to 1000),
**not** copied as PocketOsc-SI occupancy rules.

- \(t_{\mathrm{end}}=1000\) SI-scaled units (several cycles after
  transients at 48–88; SI Fig. 5 horizon).
- Discard / warmup \(t<180\) (two periods of the long ~90 class).
- Sample interval \(0.5\).
- Peaked variable: intracellular LuxI \(I\) (Fig. 4a; SI fluorescence
  is the *luxI*-driven reporter).
- Peak detector (`scipy.signal.find_peaks`): prominence
  \(0.20\times(P_{95}-P_{5})\) of \(I\) on the accepted window;
  minimum spacing 25 (below 48, above \(\tau=10\)).
- Period = mean inter-peak interval on the accepted window.
- Minimum peak count: 4.
- Occupancy / anti-transient gate (all required for `OSC`):
  1. \(\ge 4\) accepted \(I\) peaks after discard;
  2. at least one of those peaks lies in the last 30% of
     \([180,1000]\);
  3. \(P_{95}-P_{5}\) of \(I\) on the last 40% of the accepted window
     is at least \(0.40\) times that of the first 40%;
  4. relative amplitude \((P_{95}-P_{5})/\mathrm{median}(I)\ge 0.25\).
- Else `NO_PERIOD`. A decaying first fire is not OSC.

SI Data Analysis: “Peak-to-peak values were taken for all period
measurements.” Amplitude is reported as mean peak-to-previous-trough
of \(I\) on the accepted window, but is not an identity gate.

## Identity gates (predeclared)

D0 **PASS** only if the **primary** arm satisfies all of:

1. `OSC` at every identity \(\mu\in[1.0,2.0]\).
2. Period strictly increases with \(\mu\) on that interval: consecutive
   identity periods satisfy \(T_{i+1}>T_i\), Spearman
   \(\rho(\mu,T)\ge 0.95\), and \(T(2.0)>T(1.0)\).
3. Every identity period lies in the published model class window
   \([40,110]\) (documented margin around 48–88; 55/90 experimental
   minutes are not a target).
4. Tolerance check passes at primary \(\mu=1.5\).

Finite-interval extras are **reported**. Fig. 4c’s finite band is
primarily in cell density \(d\); at fixed \(d=0.5\) a \(\mu\)
bifurcation is expected (“after the bifurcation line is crossed”) but
is not a hard PASS requirement if identity 1–2 already OSC.

If the primary arm does not recover the oscillatory regime: **FAIL**.
Stop. Do not retune \(\alpha,\tau,\gamma_\ast,k_1,d\). Do not start D1.

## Required evidence

1. Equation-diff of corrected vs PocketOsc-SI \(P\), with SI quotation.
2. Failed-prior replay of `FAILED_PRIOR_IVP` on both \(P\) forms.
3. Primary-arm traces at \(d=0.5\) over the frozen \(\mu\) grid.
4. Period table vs \(\mu\) with OSC / NO_PERIOD.
5. Trend test on the identity interval.
6. Finite-interval extras reported.
7. Tolerance check.
8. History sensitivity: which predeclared histories OSC at \(\mu=1.5\);
   do not promote a post-hoc history.

## Commands

```
python examples/BSimReservoirPlanPocketOscSI_D0/run_d0.py
python examples/BSimReservoirPlanPocketOscSI_D0/check_d0.py
python -m unittest examples.BSimReservoirPlanPocketOscSI_D0.tests.test_d0_equations
```
