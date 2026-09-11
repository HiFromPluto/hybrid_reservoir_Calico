# E5.0 — biomarker patient protocol audit (no BSim)

Frozen: 2026-08-17, before any baseline score, lambda, or threshold.

This package is analysis and design-freeze only. It does not authorize
Java or BSim. It does not edit `GATE_EVIDENCE.md`. It does not retune
`K`, `n`, `tau_R`, `tau_L`, clamp, mortality, flow, or layout. It does
not reopen Track B. It does not start glucose, Danino, vesicles, C1,
Lorenz 202/303, Stage99, E4.1 222/333, or E4.2 retunes.

Track B remains Overall FAIL (synthetic 5-channel classification;
field/plumes classify the patient). This audit is a different object
and is not a rewrite of that FAIL.

The scientific question (audit, not a living claim):

> Can a HybridDish biomarker-patient protocol be frozen that is not
> already solved by a subject-legal map of the raw channels, and that
> has a valid subject split and a valid isolated/sequential design on
> this dish’s timescales?

If the answer is no, write `TASK_VOID` and stop. That is a successful
audit.

## Standing claims (E0–E4.2) — not reopened

Allowed claims remain those in
`examples/HybridDish/HIGH_TIER_SUBMISSION_EVIDENCE_PLAN.md` after
E0–E4.2. This package adds no living-dish score.

## Source and vendor

External tree (untrusted; port method and files, not scores):

`C:\Users\Ceylin\Downloads\PRC_Module8_Integrated_QS_Audit01\PRC_Module8_Integrated_QS_Audit01\reservoir_new\`

Vendored byte-for-byte, not regenerated:

`examples/HybridDish/external_analysis/synthetic_biomarker_BSim_inputs.mat`

copied from

`reservoir_new/data/synthetic_biomarker_BSim_inputs.mat`

Generator (quote-only; not executed):

`reservoir_new/matlab/MAT_biomarker_progression.m`

Batch construction (method only; 6 s windows are not HybridDish):

`reservoir_new/matlab/helpers/regenerate_batch_inputs.m`
`reservoir_new/biomarker_inputs/protocol.json`
`reservoir_new/biomarker_inputs/batch_structure.json`

Reset timescale notes:

`examples/HybridDish/design_space_preflight/E0_PREFLIGHT_REPORT.md`
section 8.

Do not copy `reservoir_new` NRMSE, AUC, IPC, or GR numbers.

## Module 0 — inventory (frozen questions)

Literal answers go in `results/BIOMARKER_INVENTORY.md`.

- `COHORT_STATUS`: `SYNTHETIC_COHORT` or `CLINICAL_ABSENT`
- N patients, T windows, B channels
- channel names and original units before min-max
- missingness
- how Y was generated (quote the MATLAB)
- whether any PhysioNet / MIMIC / eICU file exists in either tree
- Track B relationship: different object; FAIL stays
- HybridDish-illegal maps: glucose field, acid-as-lactate analog wire,
  5-AC Track B site table

This is not a clinical prediction task. Do not call `SYNTHETIC_COHORT`
a clinical dataset.

## Module 1 — subject-level split (frozen before any score)

Split unit: patient index. A window from patient `p` may appear in
only one of `{development, validation, confirmation}`.

Naive blocked cut, **not used** for scoring if the `.mat` is
class-blocked:

| Split | Patient indices (0-based) | Size |
|---|---|---|
| development | 0 .. 399 | 400 |
| validation | 400 .. 599 | 200 |
| confirmation | 600 .. 999 | 400 |

`MAT_biomarker_progression.m` assigns `Y` by construction in index
order (`i <= N/3` class 0, `i <= 2N/3` class 1, else class 2). If
that blocking is present in the vendored `.mat`, the naive cut is
**not stratified**. Do not reshuffle after seeing baselines.

**Frozen stratified rule, seed `20260817`, written before metrics:**

1. `rng = numpy.random.RandomState(20260817)`.
2. For class `c` in `{0, 1, 2}` in that order: take sorted patient
   indices with `Y==c`, then `rng.shuffle`.
3. Hamilton / largest-remainder allocation of 400 / 200 / 400 seats
   with class weights `n_c / 1000`. Truncate quotas toward zero.
   Remaining seats go to the largest fractional remainder; remainder
   ties break to the smaller class index. Confirmation receives the
   leftover indices in shuffle order after development then
   validation slices.
4. Confirmation is untouched: no architecture, lambda, or threshold
   selection on it.

Write `results/patient_split.json` with patient index lists and
SHA-256 of the concatenated comma-separated index strings (UTF-8),
one hash per split, lists in sorted ascending order for identity.
Also store within-class shuffle order for scout draws.

Do not rebalance after scores.

## Module 2 — direct-input baselines (no BSim)

Load `X`, `Xnorm`, `Y` from the vendored `.mat`. Features use
`Xnorm` as stored. All fits on development only. Lambda on
validation only. Confirmation scored once, last.

Ridge matches `examples/HybridDish/matlab/ridge_closed.m`:
population standardization fit on development rows; intercept
column unregularized; grid

`{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`.

Classification lambda: highest validation macro one-vs-rest AUC;
ties take the larger lambda. Forecast lambda: lowest validation
patient-pooled NRMSE; ties take the larger lambda. After lambda
selection, weights remain the development-only fit. Do not refit
on development+validation.

### Task A — patient-level 3-class (`Y`)

One row per patient. No window from another patient. Chance = 1/3.
Majority-class baseline required (development mode; smaller class
index on a count tie).

Features, each a 5-vector or scalar, no window leakage:

| Name | Definition |
|---|---|
| `MEAN_XNORM` | mean of `Xnorm` over T windows, 5-D |
| `LAST_XNORM` | last window of `Xnorm`, 5-D |
| `SLOPE_XNORM` | ordinary-least-squares slope of each channel vs window index `0..T-1`, 5-D |
| `MEAN_SLOPE` | scalar mean of those five slopes |

Classifiers:

- linear ridge one-vs-rest (primary)
- multinomial logistic, same lambda grid via `C = 1/lambda`,
  class-probability scores for AUC

Metrics on each split: accuracy and macro one-vs-rest AUC.
Confirmation is the void split.

### Task B — window-level horizon-1 forecast of channel 0

Channel 0 is CRP / `Xnorm[:,:,0]`. Predict `x_{t+1}` from causal
delays of all five channels at `t`. Delays never cross patient
bounds.

| Baseline | Features | Evaluation windows |
|---|---|---|
| persistence | `x_t → x_{t+1}` | `t = 0 .. T-2` |
| train-mean intercept | development mean of those targets | same as the compared model |
| linear delay-10 | 50-D: channels `0..4`, lags `0..9` (lag 0 = time `t`), channel-major | `t = 9 .. T-2` |

Persistence is also reported on the delay-10 window set as a
fairness diagnostic. The void uses the specified persistence
(`t = 0 .. T-2`) and the delay-10 model.

Metric: NRMSE = RMSE / population std of that patient’s targets,
then **patient-pooled**: mean of per-patient NRMSE. Do not concat
windows across patients for the gate. A concat NRMSE may be printed
only as a labelled diagnostic.

### Predeclared void rules

- `TASK_A_VOID` if confirmation accuracy ≥ 0.80 **or** macro AUC ≥ 0.90
  on `SLOPE_XNORM` or `MEAN_XNORM` alone (either allowed classifier).
- `TASK_B_VOID` if confirmation persistence NRMSE ≤ 0.35 **or**
  delay-10 NRMSE ≤ 0.35 (patient-pooled).
- If both void: overall `TASK_VOID`. Do not design a living BSim.
- If only one void: keep the non-void task and still do Module 3.
- Do not invent a new clinical target to escape a void.
- Do not use confirmation to choose the task.

## Module 3 — isolated vs sequential (paper freeze, no BSim)

HybridDish times, not `reservoir_new` 6 s windows:

- analysis window = 300 s
- `tau_L` = 1500 s
- warmup = 18000 s (same as the claim dish)
- E0 surrogate: `L` to 1% needs on the order of 25–30 zero windows
- 5 × 6 s reset is forbidden as a HybridDish protocol

Primary design: `ISOLATED_PATIENT`. Sequential batching is **not**
authorized until a later living reset study compares sequential vs
isolated on the **same hashed patients**.

Chemical map: do **not** freeze a 5-AC HybridDish layout here.
Constraints only:

- no glucose field
- no acid-as-lactate analog wire
- no Track B site table
- official living readout remains 408 (`R`, `L`, tagged deaths)
- any extra chemical needs its own field-only baseline

A later E5.1 prompt may pick a legal one-AHL or two-chemical map
only if a task survives.

Affordable isolated scout, if a task survives, is declared in
`results/RESET_ISOLATED_DESIGN.md`. Do not run it in this audit.

## Forbidden

- living or cell-free BSim
- copying `reservoir_new` NRMSE, AUC, IPC, GR
- glucose, vesicles, Danino, acid-as-wire
- reopening Track B FAIL
- using confirmation to choose the task
- a 100-patient sequential dish
- calling `SYNTHETIC_COHORT` a clinical dataset

## Deliverables

- `PROTOCOL.md` (this file)
- `results/BIOMARKER_INVENTORY.md`
- `results/patient_split.json`
- `results/DIRECT_INPUT_BASELINES.md`
- `results/direct_input_baselines.csv`
- `results/RESET_ISOLATED_DESIGN.md`
- `results/E5_0_VERDICT.md`

Checker: `python examples/HybridDish/biomarker_protocol_audit/run_e5_0_audit.py`
from the repository root. No Java.
