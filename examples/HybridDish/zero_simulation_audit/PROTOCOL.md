# Zero-simulation evidence audit protocol

Frozen before any audit score was calculated on 2026-08-16.

## Scope and prohibitions

This is an analysis-only package. Do not run BSim or Java. Do not modify
Java, simulation configurations, task inputs, `GATE_EVIDENCE.md`, closed
Overall lines, or existing result CSVs. Do not regenerate NARMA, waveform,
Mackey–Glass, or Lorenz inputs. Do not retune kinetics, ridge rules,
templates, thresholds, source parameters, or task definitions. Use every
existing seed; never select a best seed.

The audit may create files only in
`examples/HybridDish/zero_simulation_audit/` and
`examples/HybridDish/external_analysis/`. Existing project evidence remains
read-only.

## Closed readout protocol

Every scalar trainable readout uses the behavior of
`examples/HybridDish/matlab/ridge_closed.m` and the existing Python checkers:

- washout 40 windows;
- train 110 windows and test 50 windows;
- inner validation is the final 22 train rows: windows 128..149;
- inner standardization is fit only on windows 40..127;
- final standardization is fit on all windows 40..149;
- test windows are never used for standardization or model selection;
- zero-variance columns are retained as standardized zero columns;
- an unregularized intercept is included;
- lambda grid `{1e-6,1e-4,1e-2,1,1e2,1e4,1e6}`;
- exact-score ties choose the larger lambda;
- after validation, refit on all 110 training rows;
- NRMSE is RMSE divided by the population standard deviation of the
  corresponding target split.

Waveform classification preserves the frozen checker exactly: three
one-vs-rest ridge heads share one lambda, selected by highest validation
macro OVR AUC with larger-lambda tie-breaking. Test metrics are window macro
OVR AUC, window accuracy, per-class AUC, confusion, and block-pooled macro
AUC.

## Stage 0 — mandatory sanity

Before dependent analyses, reproduce within `1e-3`:

- NARMA driven F408 `0.9312 / 0.9228 / 0.9310`, field `1.0289`,
  Brownian mean approximately `1.1624`, silent `1.1622`, persistence
  `1.0068`, and train-intercept `1.1622`.
- Waveform driven macro AUC `0.7969 / 0.8287 / 0.8411`, field `0.8603`,
  and driven accuracy mean approximately `0.73`.
- BenchA Lorenz driven mean approximately `0.9895`, field approximately
  `0.095`; Mackey–Glass driven mean `0.2602`, field `0.0229`.

If one sanity family fails, stop analyses dependent on that family and
continue independent modules. Report the exact failure. Thresholds may not
be changed.

## Module A — frozen cell-free kinetic surrogate

Use the driven voxel files for NARMA seeds `111/222/333` and waveform seeds
`101/202/303`. Use `AHL_uM_*` samples in chronological order. Derive time as
`Window*300 + TimeInWindow_s`. Initial `R=L=0`; retain the 40-window washout.
For each of 200 voxels:

`H(C) = C^2 / (1.6^2 + C^2)`

`dR/dt = (H(C) - R) / 15`

`dL/dt = (R - L) / 1500`

Primary integration linearly interpolates stored AHL and uses internal
`dt <= 1 s`. Repeat at `dt=0.5 s` for convergence. Sensitivity uses
zero-order hold on the previous sample. Validate against a constant-C
analytical solution or an independent unit test.

Per seed, form normal per-window averages for:

1. `UNMASKED_RL` (400 features);
2. `OCCUPANCY_MASKED_RL` (400), setting R/L to zero at each exported sample
   wherever the driven `Den_* == 0`;
3. `UNMASKED_R` (200);
4. `UNMASKED_L` (200);
5. `MASKED_R` (200);
6. `MASKED_L` (200).

Do not add actual deaths. The biological comparator is driven `FRL`
(`Receiver_R_* + Lum_Mean_*`, 400 features); official F408 is historical
reference only.

NARMA rows report seed, family, integration method, feature count, lambda,
test NRMSE, and test R2. Waveform rows report seed, family, integration
method, feature count, lambda, window macro OVR AUC, accuracy, per-class AUC,
and block-pooled macro AUC. Comparisons are descriptive and create no new
gate.

Attribution language is fixed: field to unmasked estimates the contribution
of the assumed Hill/reporter kinetics; unmasked to masked estimates living
spatial occupancy/sampling; masked to actual FRL contains effects not
reproduced by the voxel cascade, including sub-voxel cell-position sampling
and population dynamics. A surrogate match does not mean nothing biological
contributes: R/L are modeled biological kinetics.

## Module B — NARMA direct and trivial baselines

Use frozen `narma10_target.csv` and the closed split:

1. `TRAIN_INTERCEPT`: constant training-target mean;
2. `PERSISTENCE`: teacher-forced `y[n] -> y[n+1]`;
3. `LINEAR_U_DELAY_10`: `[u[n],...,u[n-9]]`, zero padded;
4. `NARMA_INFORMED_INPUT`: the ten taps plus exactly `u[n]*u[n-9]`;
5. `TEST_MEAN_ORACLE_REFERENCE`: constant test-target mean, labelled
   “non-predictive oracle normalization reference; uses test labels.”

No additional terms may be added after scores are seen. Trainable baselines
use the closed ridge. Report feature count, lambda, train, validation and
test NRMSE, and test R2. Rank only legal baselines and state whether driven
mean NRMSE `0.9283` beats each.

## Module C — leakage, falsification, and feature audit

For every NARMA driven seed, evaluate official biology and field-only against
targets circularly shifted by offsets `20,40,60,80,100`. Every shift reruns
standardization, validation, lambda selection, and final refit; true-label
lambdas may not be reused. These are falsification controls, not alternative
tasks.

Only one NARMA input realization exists. Do not mislabel bacterial seeds as
independent inputs. Implement and test a state-trajectory/target-trajectory
mismatch function, but record the experiment as `PENDING` until independent
inputs exist. Circular shifts are the available misalignment null.

Audit explicitly with PASS/FAIL and code locations: validation is windows
128..149; test is excluded from lambda selection; inner standardization uses
only inner train; final standardization uses all 110 train rows; intercept is
unregularized; zero-variance columns cannot generate NaN; ties select larger
lambda. Reprint actual F408, FRL, R-only, and L-only; do not select a new
official family.

## Module D — Lorenz zero-simulation baselines

Use both frozen target files. Reconstruct 201 Lorenz states using checker
RK4/generator logic, verify stored `u`, `x`, `y_next`, declared hashes, and
then derive `x[n+1]`.

Analyze BenchA `delta_t_sample=0.02` and BenchA2
`delta_t_sample=1.0`. Targets are `AUTO_X=x[n+1]` and
`CROSS_Y=y[n+1]`.

- AUTO_X: train-intercept, persistence `x[n]`, and linear AR delay vectors
  with `m in {1,3,5,7,10}`.
- CROSS_Y: train-intercept, scalar linear `x[n]`, and x-delay vectors with
  `m in {3,5,7,10}`.
- Both: test-mean oracle reference, clearly non-predictive.

Report every development m. Select m on validation only, never test, and
report the selected test result. Include existing field-only and driven
scores where alignment permits. Recommend L1/L2/L3 without generating inputs
or running BSim.

## Module E — IPC affordability

Record current geometry: 200 total, 160 post-washout, 110 train, 50 test,
408 official features, observed linear MC approximately 1.22, and useful
lags approximately five.

Compare at least a broad Dambre-style basis and a narrow predeclared panel.
For each report target count, samples per target, feature/sample rank limit,
compute burden, estimator/null burden, and recommendation. Exact sample needs
are protocol-dependent; do not present `10^4–10^5` as a theorem. Current 200
windows are plainly insufficient for broad IPC.

Print exactly one primary line:

- `IPC_DECISION: KILL_BROAD_KEEP_MC`, or
- `IPC_DECISION: NARROW_PANEL_FEASIBLE`, or
- `IPC_DECISION: LONG_RANDOM_DRIVE_REQUIRED`.

## Module F — waveform moment diagnosis

Use frozen five-point templates and verify mean, variance/population SD,
power, min, max, and range. Expected means are sine `0.25`, square `0.29`,
triangle `0.21`. Do not assume moments explain classification before testing.

Use blocks as observations:

- washout blocks 0..7;
- train blocks 8..29;
- test blocks 30..39;
- inner train blocks correspond to windows 40..124;
- validation is whole blocks 25..29 (windows 125..149).

No block may be split. Fit frozen OVR ridge classifiers for `MEAN_ONLY`,
`POWER_ONLY`, `VARIANCE_ONLY`, `MOMENTS` (mean, variance, power, min, max,
range), and `RAW_TEMPLATE_5`. Report macro AUC, accuracy, per-class AUC,
confusion, and the ten-test-block limitation.

Audit the proposed eight-point sine/square/triangle moments and unequal
SD/power cue. Propose without generating production files: (a) recognizable
moment-aware waveforms with explicit moment baselines and (b)
equal-histogram temporal-order templates with matched scalar moments. Neither
is final without a frozen Waveform2 protocol.

## Module G — external provenance

Copy the three named external markdown files byte-for-byte, preserve their
filenames, and do not edit the copies. `MANIFEST.md` records original absolute
path, acquisition date, SHA-256, byte size, and label
“external, untrusted analysis—not project evidence”. Verify source and
destination hashes match. Do not vendor chat review unless it exists as a
source file.

## Interpretation and completion

AC/field nonlinearities are legitimate system-level hybrid computation.
Component baselines attribute computation; they do not subtract it from the
system claim. Test-mean is an oracle reference, never a predictive baseline.
Closed NARMA PASS, waveform FAIL, Lorenz FAIL, MG PASS, and sweep verdicts
remain unchanged. No audit result becomes a retrospective gate.

The integrated report must contain the executive verdict, sanity results,
kinetic attribution, NARMA baselines, leakage/null controls, Lorenz verdict,
IPC decision, waveform diagnosis, provenance link, and exact next simulation
recommendations ordered by evidence value. It ends with a short
“Zero-simulation audit completed” link to
`examples/HybridDish/HIGH_TIER_SUBMISSION_EVIDENCE_PLAN.md` without rewriting
that plan’s rationale or closed decisions.

Before completion verify: no BSim/Java command ran; no Java/config/input/gate
or existing result CSV changed; all seeds were used; all sanity checks are
printed; one command reruns the scripts; vendored hashes match sources; CSVs
and report agree; working files remain in the two authorized directories.
Stop after the report and plan link. Do not start independent NARMA,
Waveform2, Lorenz, IPC BSim, AC implementation, or Stage99.
