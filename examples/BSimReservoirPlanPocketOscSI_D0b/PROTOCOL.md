# D0b protocol — bulk delay-DDE, predeclared AHL kick, Figure 4

**Gate:** D0-reopen (`D0b`). Not D1.  
**Prior standing:** `examples/PocketDish/D0_BULK_ORACLE_STANDING.md` remains **FAIL**.  
This file does not overwrite that FAIL. Diagnosis
`examples/PocketDish/D0_FAIL_DIAGNOSIS.md` is cited, not a pass artifact.

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, before any D0b identity-grid run.

Equations: [`EQUATIONS.md`](EQUATIONS.md). Integrator and period checker
are copies of `examples/BSimReservoirPlanPocketOscSI_D0/` (sound D0
oracle). PocketOsc-SI is **not** patched.

No NARMA, GFP maps, Java, spatial \(N=200\), D1g, \(\mu_{\mathrm{bus}}\),
`QS_KMLA`, HybridDish, packed rods, ACs, or parameter retune after seeing
periods. Do not hunt ICs. Do not add \(\mu\) after occupancy. Do not
convert SI-scaled time.

## Scientific question

At TAKEN parameters, \(d=0.5\), bulk \(D_1=0\), SI Hill \(P\):

Does a **predeclared AHL-kick IVP** produce sustained LuxI oscillations
on a **predeclared \(\mu\) line**, and does period increase with \(\mu\)
on that line in the Fig. 4b sense?

Hypothesis: **unknown**. Post-FAIL diagnosis suggested a finite OSC
window near \(\mu\approx 0.3\)–\(0.8\) and a possible **non-monotonic**
period. This reopen freezes the test so a non-monotonic result is
**identity FAIL**, not a silent pass.

Main-text Fig. 4b (PMC author manuscript caption):

> “Period of oscillations as a function of the flow rate \(\mu\) at cell
> density \(d=0.5\) (top panel).”

Fig. 4c caption:

> “Oscillations occur over a finite range of cell densities, and period
> increases with \(\mu\) after the bifurcation line is crossed.”

Main-text modeling (PMC):

> “The period grows with the external AHL flow rate (effective
> degradation) and the amplitude of the oscillations, in good agreement
> with the experiments (compare Fig. 4b with Figs. 3c and d).”

Main-text AHL kick (PMC; wave-array language, used here as bulk
justification because \(P\) depends on delayed \(H_i\), not on \(I\)):

> “A small AHL perturbation in the middle of the array initiates waves
> of LuxI concentration (Fig. 4c)”

Experimental comparison class (main text): channel flow 180–296 µm/min
gave periods 52–90 min. Absolute experimental minutes are **not** a
fudge target. Model periods are SI-scaled; the published **model class**
is about 48–88 scaled/min-comparable units as used in Fig. 4b. Identity
window used here is the documented margin \([40,110]\). Do not retune
\(\alpha\) to hit 55/90 experimental minutes.

## A. TAKEN table (do not retune)

Copied from D0. No change to \(\alpha,\tau,\gamma_A,\gamma_I,\gamma_H,
k_1,d\).

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

Time: SI-scaled. No `TIME_ADJ=60`. Do not convert after seeing periods.
\(\mu\) is the SI external-AHL decay rate. Not channel µm/min, not
\(\mu_{\mathrm{bus}}\).

## B. History / IC ensemble (frozen; `PREDECLARED_AHL_KICK`)

The SI does not state bulk \(A,I,H_i,H_e\) or \(H_i(t<0)\). D0’s primary
`SI_BASAL_PERTURB` (\(I=1\), \(H_i=0\), history 0) is the **wrong bulk
kick**: \(P\) depends on delayed \(H_i\), so history 0 forces \(P=\delta\)
on \([0,\tau]\). That arm latched ON for \(\mu\le 0.03\) and stayed OFF
for \(\mu\ge 0.04\).

Primary is an **AHL kick**, justified from the main-text “small AHL
perturbation”, not SI Fig. 6 protein-only \(I=1\).

Kick amplitude **chosen once now** from the diagnosis note, not from the
new grid:

- \(H_i=H_e=0.05\) occupied at diagnosis \(\mu=0.4\) (period 63.6).
- \(H_i=H_e=0.01\) was `NO_PERIOD` at the same \(\mu\).
- \(0.10\) and \(0.20\) also occupied; they are not the primary.

**One primary:** `AHL_KICK_005`. Do not choose the kick after seeing the
new grid. Do not hunt ICs after traces. Do not promote
`AHL_KICK_020` into the identity arm.

All arms: \(d=0.5\), constant delay history.

| id | \(A,I,H_i,H_e\) at \(t=0\) | \(H_i(t<0)\) | Role |
|---|---|---|---|
| `AHL_KICK_005` | \(0, 0, 0.05, 0.05\) | \(0.05\) (equals state) | **PRIMARY.** Main-text small AHL perturbation; 0.05 occupied in diagnosis at \(\mu=0.4\); 0.01 did not |
| `SI_BASAL_PERTURB` | \(0, 1, 0, 0\) | \(0\) | failed D0 primary; expect OFF except possibly tiny \(\mu\) |
| `AHL_KICK_001` | \(0, 0, 0.01, 0.01\) | \(0.01\) | subthreshold kick from diagnosis |
| `FAILED_PRIOR_IVP` | \(0, 0, 1, 1\) | \(1\) | PocketOsc-SI ENGINEERING start |
| `AHL_KICK_020` | \(0, 0, 0.20, 0.20\) | \(0.20\) | OPTIONAL basin-check; **not** the identity arm |

Control \(\mu\) (frozen now, all non-primary arms): \(\mu=0.40\).

Additional D0 comparison (not IC hunting): `FAILED_PRIOR_IVP` at
\(\mu=1.5\) on both \(P\) forms.

## C. \(\mu\) grid (frozen; do not add after occupancy)

Identity line sits inside a basin that can cycle, **not** in \(\mu=1\)–\(2\)
(that D0 line is the OFF well on the TAKEN table).

**Identity, \(d=0.5\) (scored):**

\(\mu = 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60\).

**OFF / edge extras (reported; not used for identity occupancy or
Spearman):**

\(\mu = 0, 0.10, 0.20, 0.28, 0.80, 1.00, 1.20, 1.50, 2.00\).

Keep \(\mu=1\)–\(2\) as extras so the failed D0 line is shown as
**outside** the window, not dropped.

Do not add points after occupancy. Do not convert \(\mu\) from channel
µm/min or \(\mu_{\mathrm{bus}}\).

## D. Horizon and period rule

Reuse the D0 period checker. No documented bug was found in
`period_check.py`. Peak rules are **not** changed. SI Data Analysis
(quoted in D0 PROTOCOL): “Peak-to-peak values were taken for all period
measurements.” That is already the D0 rule.

Keep SI-scaled time. No `TIME_ADJ=60`.

Frozen (copy of D0):

- \(t_{\mathrm{end}} = 1000\)
- discard \(t < 180\)
- `sample_dt = 0.5`
- readout = intracellular \(I\)
- min spacing 25, relative prominence 0.20, `min_peaks = 4`
- persistence and relative-amplitude occupancy gates from D0:
  1. \(\ge 4\) accepted \(I\) peaks after discard;
  2. at least one of those peaks lies in the last 30% of \([180,1000]\);
  3. \(P_{95}-P_{5}\) of \(I\) on the last 40% of the accepted window
     is at least \(0.40\) times that of the first 40%;
  4. relative amplitude \((P_{95}-P_{5})/\mathrm{median}(I)\ge 0.25\).
- period = mean inter-peak interval on the accepted window
- else `NO_PERIOD`

Integrator (copy of D0):

- Constant delay \(\tau=10\).
- Method of steps: scipy `solve_ivp` BDF on successive intervals of
  length \(\tau\), dense output as the delay interpolant.
- Delay lookup **raises** if \(t-\tau\) falls in a gap.
- Nominal tolerances: `rtol=1e-8`, `atol=1e-10`, `max_step=0.05`.
- Tightness: 2× tighter `rtol=5e-9`, `atol=5e-11` on a **predeclared
  subset** (below). OSC/NO_PERIOD must match; period relative change
  \(\le 0.02\) when both are OSC.
- No hidden clipping. `NEGATIVE_STATE` is not scored as OSC.

### Predeclared tightness subset (frozen now)

| ensemble | \(\mu\) | expected role |
|---|---:|---|
| `AHL_KICK_005` | 0.40 | identity interior |
| `AHL_KICK_005` | 1.20 | D0-line extra; expected OFF |
| `SI_BASAL_PERTURB` | 0.40 | wrong kick; expected `NO_PERIOD` |

## E. Identity gates (predeclared)

D0b **PASS** only if **ALL** of these hold on the **primary** arm
`AHL_KICK_005`:

1. At least 6 of 8 identity \(\mu\) points are `OSC`.
2. Periods on those OSC points lie in the published model class 40–110
   SI-scaled units (48–88 is the comparison core; do not fudge \(\alpha\)
   to hit 55/90 experimental minutes).
3. Spearman correlation of \((\mu, T)\) on OSC identity points is
   \(\ge 0.95\) **and** period **strictly increases** along the identity
   line (consecutive OSC identity points ordered by \(\mu\) satisfy
   \(T_{i+1}>T_i\)).
4. Frozen extras at \(\mu\ge 1.2\) are `NO_PERIOD` or collapsed OFF.
5. `SI_BASAL_PERTURB` at identity \(\mu=0.40\) is `NO_PERIOD`
   (wrong kick remains wrong).
6. Tighter rtol/atol 2× does not change OSC/NO_PERIOD calls on the
   predeclared tightness subset, and OSC periods shift by less than 2%.

Scoring if not PASS:

- If the primary arm does **not** occupy the identity line
  (\(<6/8\) OSC): **`FAIL`**. Do not hunt another kick. Do not start D1.
- If the primary arm oscillates (\(\ge 6/8\) OSC) but period does **not**
  increase with \(\mu\) (Spearman \(<0.95\) or not strictly increasing):
  **`FAIL_NO_IDENTITY`**. Report the measured trend. Do not retune. Do
  not start D1.
- Any other failed gate after occupancy and a passing trend: **`FAIL`**.
  Do not retune. Do not start D1.

A non-monotonic period on the frozen identity line is an identity
failure, even if the basin is occupied. Do not retune to force Fig. 4b.

## Required evidence

1. Frozen PROTOCOL and configs with `frozen_before_traces=true` **before**
   the identity-grid run.
2. Equation transcription (copied from D0; SI Hill \(P\)).
3. Primary-arm traces over the frozen \(\mu\) grid.
4. Period table vs \(\mu\) with OSC / NO_PERIOD.
5. Spearman and monotonicity on OSC identity points.
6. Frozen extras, including \(\mu=1\)–\(2\).
7. Control arms at \(\mu=0.40\), including `SI_BASAL_PERTURB`.
8. Tightness subset.
9. Fixtures: at least one OSC \(\mu\) and `SI_BASAL_PERTURB`.
10. Standing `examples/PocketDish/D0B_BULK_ORACLE_STANDING.md`.
    Do **not** rewrite `D0_BULK_ORACLE_STANDING.md`.

## Commands

Working directory: `examples/BSimReservoirPlanPocketOscSI_D0b`

```
python -m unittest tests.test_d0b_equations tests.test_d0b_protocol
python run_d0b.py
python check_d0b.py
```

From the repository root:

```
python -m unittest examples.BSimReservoirPlanPocketOscSI_D0b.tests.test_d0b_equations examples.BSimReservoirPlanPocketOscSI_D0b.tests.test_d0b_protocol
python examples/BSimReservoirPlanPocketOscSI_D0b/run_d0b.py
python examples/BSimReservoirPlanPocketOscSI_D0b/check_d0b.py
```

`--smoke` runs only primary \(\mu=0.40\). It must not change ICs,
\(\mu\), or parameters.

## What this protocol is not

Not a D0 PASS. Not a promotion of `diagnose_ahl_kick.py`. Not D1.
Not a retune of the TAKEN table. Not a hunt for a monotonic \(\mu\)
sub-window after traces.
