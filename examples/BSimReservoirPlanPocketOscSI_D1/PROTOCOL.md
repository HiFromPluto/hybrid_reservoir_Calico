# D1 protocol — Java circuit parity of Object B (NOT_FIG4B)

**Gate:** D1 Java circuit parity.  
**Organism:** `DaninoSI_OccupiedDF` (Object B). **NOT_FIG4B.**  
**Not ported:** Object A `Danino2010_Fig4b_bulk_twin` (FAIL / FAIL_NO_IDENTITY).

Claim freeze (must exist before this run):
[`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).

D0 remains FAIL. D0b remains FAIL_NO_IDENTITY. This protocol does not
reopen Fig. 4b identity, does not shrink \(\mu=0.32\)–\(0.40\), does not
retune \(\alpha,\tau,\gamma_A,\gamma_I,\gamma_H,k_1,d\), and does not
hunt ICs.

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after claim freeze, before any D1 Java integration.

No space, rods, ACs, NARMA, `TIME_ADJ=60`, \(\mu_{\mathrm{bus}}\), D1g,
`QS_KMLA`, or HybridDish.

## Scientific question

Does a Java constant-delay bulk DDE, at the TAKEN table with SI Hill
\(P\) and delay only in \(H_i\), match the frozen D0b Python fixtures of
Object B within predeclared tolerances?

This is **not** a Fig. 4b test. Periods are compared to D0b Python, not
to 55/90 experimental minutes.

## TAKEN table (copied; do not retune)

Same as D0b: \(C_A=1\), \(C_I=4\), \(\delta=10^{-3}\), \(\alpha=2500\),
\(\tau=10\), \(k=1\), \(k_1=0.1\), \(b=0.06\), \(\gamma_A=15\),
\(\gamma_I=24\), \(\gamma_H=0.01\), \(f=0.3\), \(g=0.01\), \(d_0=0.88\),
\(D=2.5\), \(d=0.5\), \(D_1=0\). SI-scaled time.

\(P=(\delta+\alpha H_\tau^2)/(1+k_1 H_\tau^2)\), \(H_\tau=H_i(t-\tau)\).

## Frozen fixtures (copy; do not regenerate to fit Java)

From `examples/BSimReservoirPlanPocketOscSI_D0b/results/fixtures/`:

| file | IVP | \(\mu\) | D0b flag | D0b period |
|---|---|---:|---|---:|
| `primary_mu_0p40.*` | `AHL_KICK_005` | 0.40 | OSC | 63.583 |
| `primary_osc.*` | `AHL_KICK_005` | 0.32 | OSC | 56.308 |
| `si_basal_perturb_mu_0p40.*` | `SI_BASAL_PERTURB` | 0.40 | NO_PERIOD | — |
| `failed_prior_mu_1p5.*` | `FAILED_PRIOR_IVP` | 1.50 | NO_PERIOD | — |

Also integrate primary `AHL_KICK_005` at \(\mu=1.5\) (D0b extra; no
separate npz). Must remain `NO_PERIOD` / OFF-class as in D0b.

## Integrator (numerical; not a biological retune)

- Constant delay \(\tau=10\). Delay only in \(H_i\).
- RK4 with fixed `rk_dt = 0.001` SI-scaled units.
- Linear interpolant of stored \(H_i\) for \(t-\tau\). For \(t-\tau<0\),
  frozen constant history.
- Sample grid copied from D0b: \(t_{\mathrm{end}}=1000\),
  `sample_dt=0.5`.
- Injected `BSimRandom` is accepted and **unused**. Circuit is
  deterministic.
- Deprecated `bsim.dde.BSimDdeSolver` is **not** used (known-bad
  history construction).

## Period rule

Copied from D0b, unchanged: discard 180, spacing 25, rel prominence
0.20, min_peaks 4, persist last 30% with ≥1 peak, amplitude persist
ratio 0.40, relative amplitude ≥ 0.25. Readout \(I\). NARMA-blind.

## Predeclared agreement gates

D1 **PASS** only if all hold, with **no** TAKEN-table change:

1. Claim freeze file exists and names Object B as NOT_FIG4B.
2. On the frozen sample grid, for every fixture above, every sample,
   every state \(A,I,H_i,H_e\):
   \[
   |x_{\mathrm{java}}-x_{\mathrm{py}}|
   \le \texttt{atol} + \texttt{rtol}\,|x_{\mathrm{py}}|
   \]
   with \(\texttt{rtol}=0.02\), \(\texttt{atol}=10^{-4}\).
3. Additionally, on \([0,\tau]\) of each fixture (delay is pure
   history): \(\texttt{rtol}=10^{-3}\), \(\texttt{atol}=10^{-6}\).
4. Java period checker: primary \(\mu=0.40\) is `OSC` and
   \(|T_{\mathrm{java}}-T_{\mathrm{py}}|/T_{\mathrm{py}}<0.02\)
   with \(T_{\mathrm{py}}=63.583\ldots\) from the fixture.
5. Java `SI_BASAL_PERTURB` at \(\mu=0.40\) is `NO_PERIOD`.
6. Java `failed_prior` at \(\mu=1.5\) and primary at \(\mu=1.5\) are
   `NO_PERIOD`.
7. Java primary \(\mu=0.32\) is `OSC` and period relative error
   \(<0.02\) versus fixture \(T=56.308\ldots\).

If Java cannot match without retuning biology: **D1=FAIL**. Stop.

## Commands

```
ant d1-parity
python examples/BSimReservoirPlanPocketOscSI_D1/check_d1.py
```

Standing: `examples/PocketDish/D1_JAVA_CIRCUIT_PARITY_STANDING.md`.
Must repeat NOT_FIG4B and that D0/D0b standings are unchanged.
