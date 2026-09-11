# HybridDish

A BSim dish in which a one-way AHL source drives *E. coli* with a Hill
receiver \(R\) and a slow luminescence integrator \(L\). Acid is held
constant as a death channel. The ridge reads cell state, not vesicle
devices, glucose, or Danino oscillators.

This folder is the working hybrid (Plan Stage 6 biology). Zip and send
**this directory**. It is meant to sit at `examples/HybridDish` inside
the BSim tree. Status and purpose: `RESEARCH_REPORT.md`.

```
HybridDish/
  bsim/      Java simulator
  config/    run files
  input/     AHL, acid, NARMA target, optional waveform
  matlab/    load voxels → window-average → ridge
  output/    BSim writes run folders here
```

---

## Where to stand

All commands below assume the current directory is this folder:

```
…/bsim_clean/examples/HybridDish
```

Config paths (`input/…`, `output/…`) are relative to **here**, not to
`bsim/` or `config/`.

You need a JDK on `PATH` and a compiled BSim engine:

| Default | Meaning |
|---|---|
| `../../dist/build` | compiled `bsim.*` classes |
| `../../lib/core.jar` | Processing |
| `../../lib/vecmath.jar` | |
| `../../lib/objimport.jar` | |

If you unzipped HybridDish somewhere else:

```
set BSIM_ROOT=C:\path\to\bsim_clean
```

MATLAB is only for the readout after a run. It is not required to
compile or watch the GUI.

---

## What to run

### 1. Compile only

Windows:

```
compile_and_run.cmd
```

Linux / macOS:

```
chmod +x compile_and_run.sh
./compile_and_run.sh
```

Class files land in `HybridDish/*.class` under this folder. No simulation
starts until you pass a config.

On Linux / macOS, use `./compile_and_run.sh` and forward slashes in
paths (`config/smoke.properties`). The rest of this section uses the
Windows form.

### 2. Smoke test (simple / short)

Writes two short windows so you can check that Java starts, reads
`input/`, and creates `output/smoke/`. **Do not** treat the NRMSE as a
result. Kinetics and timing are not the production protocol.

```
compile_and_run.cmd config\smoke.properties
```

When it finishes you should have:

```
output/smoke/voxels.csv
output/smoke/window_summary.csv
output/smoke/results.csv
output/smoke/matlab_meta.txt
output/smoke/feature_contract.txt
output/smoke/run_status.txt
```

`run_status.txt` should show `summary_rows=2` and a handful of voxel
rows. That is enough to prove the I/O path. It is not 200 windows.

### 3. GUI / preview

Same smoke job, with the Processing drawer:

```
compile_and_run.cmd config\smoke.properties preview
```

`preview` overrides `headless=true`. The drawer uses the proven
`ortho` → `camera` → `perspective` sequence for this 1000×500×10 µm
box (eye pulled back by 500 µm, not the 10 µm thickness). The
production configs also accept `preview`, but a full NARMA-10 run is a
long headless job; use smoke if you only want to see cells and the AHL
pulse.

Close the window or let the short smoke clock finish.

### 4. Production NARMA-10 (the real experiment)

Warmup 18 000 s + 200 × 300 s windows at `dt=0.05`. This is hours, not
minutes. Output goes to `output/driven_seed101/`.

```
compile_and_run.cmd config\driven_seed101.properties
```

Controls (same seeds `101 / 202 / 303`):

| Config | Role |
|---|---|
| `config/driven_seed101.properties` | AHL carries NARMA-10; biology readout |
| `config/silent_seed101.properties` | Sources off |
| `config/brownian_seed101.properties` | Same pulses; ridge uses density only |
| `config/waveform_driven_seed101.properties` | Optional sine/square/triangle AHL |

### 5. MATLAB readout (after a production run)

From this same folder:

```matlab
addpath('matlab')
S = ridge_narma10('output/driven_seed101');
```

Official biology features (408): `Receiver_R_*` (window mean),
`Lum_Mean_*` (window mean), `Input_Driven_Death_*` (last sample).
Voxel AHL is a diagnostic only:

```matlab
Sf = ridge_narma10('output/driven_seed101', 'input/narma10_target.csv', ...
                   'Family', 'field');
```

Smoke-test CSVs are too short for this script (it expects 200 windows).

Details: `matlab/README.md`. Protocol numbers: `PROTOCOL.md`.

---

## Frozen kinetics

Do not retune these if you are reproducing the gated dish:

- Hill receiver: \(K = 1.6\) µM, \(n = 2\), \(\tau = 15\) s
- Luminescence: \(\alpha = \delta = 1/1500\) s⁻¹
- AHL source \(1.28\times10^8\), acid source \(2\times10^{11}\), \(K_\max = 0.002\)
- Acid input held at \(0.5\) every analysis window
- `FLOW_SPEED = 0`; attractant ACs stay silent

NARMA-10 split: washout 40, train 110, test 50. Lambda is chosen on
windows 128–149 only, then the readout is refit on all 110 train windows.

---

## How this dish was reached

The earlier generation (`reservoir_new`) had the right **architecture**
(a 1000×500×10 µm monolayer, windowed voxel export, a Brownian density
null, NARMA-10 as the task). The **operating point** did not: 6 s
windows against a 30 min generation, a QS switch whose field never
reached threshold, vesicle stores that did not track \(u\), and (in one
merge) unitless Monod glucose. Plan stages rebuilt one mechanism at a
time from the Stage 10 lineage, not by patching that merge.

| Stage | What changed | Gate |
|---|---|---|
| 1–2 | Affordable grid `50×25×1`, explicit `20×10` + `4×2` voxels, 300 s windows, 18 000 s warmup, population near the clamp | Protocol / plateau, not a computing claim |
| 3 | Danino 4-ODE as the dish receiver | **FAIL** — field tracks input; intracellular \(q\) does not (\(r(q,u)\approx 0\), occupancy ~0.001) |
| 3B | Replace that receiver with a fast Hill \(R\) on the same plume | **PASS** — \(R\) tracks input; \(K=1.6\) is trace-matched to this field, not a LuxR EC50 |
| 4 | Mixed-acid point source + input-tagged deaths | Mixed: \(r(\text{acid}, \text{Input\_Driven\_Deaths})\approx 0.89\); 2× fold on `Total_Deaths` **FAIL**. Ridge may use tagged deaths only |
| 5 | \(L += (R-L)\,dt/1500\) | **PASS** — \(L\) is slower than \(R\), not a duplicate |
| 6 | NARMA-10 on frozen 3B+4+5 biology | **PASS** vs Brownian and silent; biology also below a voxel-AHL delay line |
| 7 | Four-source SVD rank | **FAIL** — silent already spans the rank rule |
| 9 | Product-bit classification | **FAIL** — driven AUC ≈ field `AHL_uM_*` |
| Bench A | Lorenz / Mackey–Glass | Lorenz **FAIL**; Mackey–Glass beats Brownian, but field-only is much stronger |
| Track B | 5-AC “patients” | **FAIL** — chemical maps do the work |
| Waveform | Sine / square / triangle | Biology well above chance and Brownian; field-only still a bit better |

Citation boxes (not copied into this dish): Danino in a small bath **D1 FAIL**, **D1g PASS** once dilution/ICs are restored; Monod **M1–M3 PASS** in isolation. Glucose and Danino stay out of HybridDish.

The honest computing claim for **this package** is narrow: on NARMA-10,
the 408-feature biology readout beat a density null, a silent dish, and
a linear map of voxel AHL. Test NRMSE ≈ **0.93** (Stage 6 mean 0.9306;
replication seeds 111/222/333 mean 0.928). That is a weak predictor
(\(R^2\) on the order of 0.15). It is not a tuned echo-state net. Later
tasks mostly show that a map of the plume is enough, or better.

---

## Comparison (earlier dish vs this one)

Labels are conservative. “Stronger” here means the operating point is
alive and gated, not that the score is impressive.

| | `reservoir_new` (what ran) | HybridDish (this folder) |
|---|---|---|
| Domain | 1000×500×10 µm | Same |
| Windows | 6 s pulse 1.5 s | 300 s pulse 75 s |
| Warmup | 10 s (audit) / 600 s (on-disk) | 18 000 s (for \(L\)) |
| QS | Hysteretic threshold \(K\sim U[12,24]\) molecules/µm³; field ~0.05, switch off | Continuous Hill \(R\), \(K=1.6\) µM, \(n=2\), \(\tau=15\) s |
| Light | Wanted, not a live lux state in the on-disk run | \(L\) integrator, \(\tau=1500\) s |
| Sources | Vesicle ODE + store; no AHL-secreting AC | Point `addQuantity`; AHL at (500,250), acid held |
| Growth | Unitless Monod in that generation | Fixed 30 min rate, no glucose |
| Flow | On | Off |
| Readout | Field-inherited voxels (thousands–10⁵ columns) | Explicit 20×10 + 4×2; ridge **408** |
| NARMA-10 | Operating point could not support an honest test | Driven 0.93 vs Brownian 1.16 vs field-only ~1.01–1.03 |

What was kept from the old project: the dish idea, voxel export, a
Brownian null, and NARMA-10 as the thing that has to beat density.

What was not kept: vesicle ACs, glucose as a computing channel, Danino
in the 300 s dish, 6 s count features, and a QS switch below threshold.

---

## What this package is not

It is not Stage 11, `reservoir_new` vesicles, Monod glucose, or Danino
inside the reservoir. Those were tried, gated, and left in their own
folders. Do not merge them back to “improve” 0.93.
