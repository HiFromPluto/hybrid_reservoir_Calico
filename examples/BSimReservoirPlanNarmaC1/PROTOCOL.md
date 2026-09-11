# C1 independent-input NARMA-10 protocol (frozen before NRMSE)

Independent-input NARMA-10 replication on the frozen HybridDish claim
dish. Clone of `examples/BSimReservoirPlanNarma10b/`. Change only the
`u` file each run points at. Mechanisms are not retuned. The unit of
generalization is the **independent input realization**, not the
bacterial seed.

Narma10b already showed driven F408 **0.928** vs Brownian/silent
**1.162** vs field **1.029** vs 10-tap of \(u\) **0.683** on **one**
`u` sequence (SHA-256 `d6c0cdfb…4c1e`), three bacterial seeds. C1 asks
whether the Brownian/silent/field result reproduces on new input
trajectories. It does **not** ask whether cells beat 0.683.

Zero-simulation audit is complete. Do not rerun it. Do not retune the
dish. Do not discard a sequence because it scores badly. Do not start
E5, washout twins, two-way, IPC, Stage99, WaveformS, or layout-flow.
`GATE_EVIDENCE.md` files are not edited. E5.1 stays `NOT_SCORED`.
E5.2 stays a separate occupancy freeze. Claim dish is not replaced.

## Scientific question

Does driven biology beat Brownian, silent, and field-only on
independent NARMA-10 trajectories without changing the dish,
kinetics, readout, or analysis rule?

If beating the legal 10-tap baseline is required for a headline, stop.
This package’s headline is reproducibility vs Brownian / silent / field.

## Dish and analysis (copied, not retuned)

- `1000 × 500 × 10` µm, `dt=0.05`, `FLOW_SPEED=0`, `NO_FLUX`
- Field `50×25×1`, readout `20×10` + `4×2`
- `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL `1.28e8` at `(500, 250, 5)` (CENTER)
- Acid `2e11` at `(300, 375, 5)`, held `0.5` in analysis windows, off
  in warmup
- Attractant AC0/AC2 stay at `u=0`
- Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, `warmup.ahl.input=0.5`
- 200 windows, 300 s, pulse 75 s, sample every 20 s
- Last sample `199;15;299.95`. CSV 200 / 3200 / 3200
- Official 408: `Receiver_R_*`, `Lum_Mean_*`, `Input_Driven_Death_*`
- Brownian ridge: `Den_*` 20×10 only
- Washout `0..39` (40), train `40..149` (110), test `150..199` (50)
- Inner lambda on windows `128..149` only. Then refit 110.
- Grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`, unregularized intercept,
  train-only standardization, ties → larger lambda
- Independent lambda per arm × trajectory × seed
- NRMSE = test RMSE / pop-std of that trajectory’s test labels

Recurrence (unchanged):

```
y[0] = 0
y[n+1] = 0.3 y[n] + 0.05 y[n] * sum_{i=0..9} y[n-i]
         + 1.5 u[n-9] u[n] + 0.1
```

with `y[k]=0` and `u[k]=0` for `k<0`. Window `n` predicts `y[n+1]`.

## Module 0 — hash `u` before any BSim

`u ~ Uniform[0, 0.5]` via `random.Random(seed)`, 200 values, 12
decimal places. Uniform`[0,1]` is forbidden.

| TrajID | `random.Random` seed | Provenance |
|---|---|---|
| 00 | 20260814 | Existing Narma10b. Hash must be `d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`. Do not rewrite the file. |
| 01 | 2026081501 | new |
| 02 | 2026081502 | new |
| 03 | 2026081503 | new |
| 04 | 2026081504 | new |
| 05 | 2026081505 | new |
| 06 | 2026081506 | new |
| 07 | 2026081507 | new |
| 08 | 2026081508 | new |
| 09 | 2026081509 | new |
| 10 | 2026081510 | new |

Write `input_ahl_narma200_trajXX.txt` and `narma10_target_trajXX.csv`
for 01–10. Hash every file in `results/U_FREEZE.md` **before** the
first new Java job. Abort if any new hash collides with traj 00 or
with another new traj.

**Interaction subset (frozen here, not after scores):** traj **01**
and **02**. Those two, and only those two new sequences, also get
bacterial seeds 222 and 333. Traj 00 already has 111/222/333 in
Narma10b.

No sequence may be dropped because NRMSE is bad. Do not add traj 11.

## Module 1 — living BSim

Smoke: traj 01, seed 111, print CENTER, `FLOW=0`, 200 windows, finite
AHL, last sample `199;15;299.95`. That smoke **is** the production
traj 01 / seed 111 run.

Then:

| Run | n | Notes |
|---|---|---|
| driven traj 01–10, seed 111 | 10 | primary layer; traj 01 is the smoke |
| driven traj 01–02, seeds 222 and 333 | 4 | interaction subset |
| traj 00 driven 111/222/333 | 0 new | score existing Narma10b CSVs |

**Brownian and silent.** Reuse Narma10b CSVs for seeds 111, 222, 333.
Score each frozen state trajectory against every new `y`. Label:
“reused state trajectories, independently scored targets.”

The checker must prove:

1. last sample `199;15;299.95` and 200/3200/3200
2. Brownian features are `Den_*` only (no `AHL`, `R`, `L`)
3. silent sources were off in the reused run
4. the reused feature matrices are byte-identical to Narma10b

If any of those fail, do **not** invent per-input Brownian jobs.
Run one new Brownian and one new silent per bacterial seed that you
actually use, then re-score those states against every `y`. Never one
Brownian per input realization.

Field-only AHL is extracted from **that trajectory’s** driven voxels.
Do not reuse traj 00 field for traj 01.

## Module 2 — scores (every trajectory, not only the mean)

For each traj \(j\) at seed 111, and for 00/01/02 also at 222/333:

\[
\Delta_{B,j} = \mathrm{NRMSE}_{Brownian,j} - \mathrm{NRMSE}_{Driven,j}
\]
\[
\Delta_{S,j} = \mathrm{NRMSE}_{Silent,j} - \mathrm{NRMSE}_{Driven,j}
\]
\[
\Delta_{F,j} = \mathrm{NRMSE}_{Field,j} - \mathrm{NRMSE}_{Driven,j}
\]

Positive favours biology. Report every row.

Primary layer uses seed **111** only (11 trajectories: 00–10).
Seeds 222/333 on {00,01,02} are nested replicates, not new inputs.

Also report, same splits, same lambda rule, **per trajectory**:

- train-mean intercept
- teacher-forced persistence
- legal 10-tap \([u_n,\ldots,u_{n-9}]\)
- NARMA-informed input: those 10 taps plus the one product \(u_n u_{n-9}\)
- occupancy-masked kinetic surrogate (\(K=1.6\), \(n=2\), \(\tau_R=15\),
  \(\tau_L=1500\), piecewise-constant AHL, occupancy from that driven
  run)
- unmasked kinetic surrogate (diagnostic)
- \(R\)-only and \(L\)-only (diagnostic; do not replace official 408)

### Statistics (primary layer, n=11 inputs)

Mean, median, s.e., and a paired bootstrap CI for each delta.
Bootstrap: 10000 resamples of the 11 paired deltas, seed `20260818`,
percentile 95% interval. Exact sign count (how many traj have
\(\Delta>0\)). Do not substitute a \(p\)-value for effect size. No
“non-overlapping error bars” shorthand.

### Leakage (analysis only)

On traj 00 (already audited) and traj 01 (predeclared, not chosen
after scores):

- circular target shifts 20/40/60 windows; full lambda re-selection
  on the shifted labels; do not reuse true-label lambda
- mismatch: F408 of traj 00 vs `y` of traj 01, and the reverse.
  Expect collapse toward intercept. Do not pick a prettier pair.

Train/test target mean and variance for every traj.

## Standing (predeclared)

| Label | Rule |
|---|---|
| System (primary) | seed-111 driven beats Brownian and silent on the sign count and mean \(\Delta_B,\Delta_S\) across traj 00–10 |
| Living-layer (diagnostic) | whether \(\Delta_F>0\) also reproduces; trajectory-dependent is allowed language |
| vs 10-tap / informed input | reported every traj; **not** a C1 kill gate; cells are not required to beat ~0.68 |
| Occupancy | mean \(R\) per traj; ALIVE if mean \(R \ge 0.05\) |
| Narma10b Overall | unchanged |
| Claim dish | not replaced |

Decision language (copy the matching sentence, do not invent a
stronger one):

- If biology beats Brownian/silent across independent inputs:
  “The frozen living state reproducibly carries NARMA-10 information
  beyond density and silent-dish controls across independent drives.”
- If \(\Delta_F>0\) also repeats:
  “On NARMA-10, the slow living-state readout adds weak predictive
  information beyond a linear voxel-AHL delay line.”
- If field-beating does not replicate but Brownian/silent does:
  “The living pathway is real, but its advantage over the carrier was
  trajectory-dependent.”

Do not retune after any of these.

## Compile

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanNarmaC1.java VoxelAnalyzer.java
```

## Runs

```
python examples/BSimReservoirPlanNarmaC1/generate_c1_inputs.py
python examples/BSimReservoirPlanNarmaC1/freeze_c1_sources.py
python examples/BSimReservoirPlanNarmaC1/run_c1_jobs.py --smoke
python examples/BSimReservoirPlanNarmaC1/run_c1_jobs.py
python examples/BSimReservoirPlanNarmaC1/check_c1.py
python examples/BSimReservoirPlanNarmaC1/plot_c1.py
```

Smoke uses `sim_config_c1_driven_traj01_seed111.properties` and is
kept as the production traj 01 / seed 111 CSV.

## Report

`results/U_FREEZE.md`, `results/C1_SCOUT.md`, `results/c1_narma.csv`,
`results/c1_deltas.csv`. Hash Java, configs, every `u` file, checker,
and result manifests.

## Stop list

- Do not retune kinetics, `Jmax`, layout, or flow
- Do not drop a `u` sequence after seeing NRMSE
- Do not add traj 11 after the freeze
- Do not claim cells beat the 10-tap baseline
- Do not rewrite Narma10b / Stage 6 Overall PASS
- Do not start E5, two-way, IPC, Stage99, WaveformS, or layout-flow
- Do not treat bacterial seeds as independent inputs
