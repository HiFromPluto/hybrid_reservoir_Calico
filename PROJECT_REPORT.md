# BSim Clean — concise project report

**What this repository is.** A working copy of [BSim](https://cellsimulationlabs.github.io/tools/bsim/) (an agent-based *E. coli* simulator) plus a staged rebuild of a millimetre-scale **hybrid bacterial reservoir**: living cells as the nonlinear state, a chemical source as the only drive, and a linear readout trained after the fact.

---

## One-sentence task (for a paper or coursework report)

The computational task is **NARMA-10**: a random input sequence \(u[n]\) is delivered as AHL pulses; a linear ridge trained on the living spatial state must predict a tenth-order nonlinear target \(y[n+1]\) that depends on the last ten outputs and on \(u[n]u[n-9]\), and must beat a Brownian density null (and a silent dish).

---

## What the simulator does

BSim treats each bacterium as an agent in a thin dish (here **1000 × 500 × 10 µm**). Cells grow, divide, die, run and tumble, and respond to diffusing chemicals. Population-level patterns (synchrony, chemotaxis, reporter maps) come from those local rules, not from a single ODE for the whole colony.

This tree keeps the original BSim engine and adds a **gated experimental programme**: each mechanism is copied forward, measured against a pre-declared number, and either kept or archived. Failed merges are left in place as evidence, not patched until they “look good.”

## What the hybrid dish is asking

The object is a **one-way hybrid reservoir**, not a trained physical neural net and not two-way device–cell feedback.

1. A point source injects **AHL** (the NARMA carrier). A second **acid** source is held constant so a death channel exists but is not a second free input.
2. Each cell runs a fast Hill receiver \(R\) and a slower luminescence integrator \(L\).
3. Time is biological: **30 min** doubling, **300 s** analysis windows, **75 s** pulses, **18 000 s** warmup so \(L\) is not still climbing.
4. After each window, a **voxel map** of \(R\), \(L\), and tagged acid deaths is exported. A **ridge regression** (linear readout only) is fit on a train split and scored on held-out windows.

**Question.** Can that living state carry a temporal task above a density null (Brownian particles in the same dish) and a silent dish, when the drive is encoded only in AHL?

**Not the question.** Whether the dish matches a trained echo-state network (those quote much lower NRMSE after training recurrent weights), or whether cells must beat a linear map of the AHL plume itself.

## Main result (NARMA-10)

On the claim dish, driven biology beat the nulls: test NRMSE about **0.93** versus about **1.16** for Brownian and silent arms (NRMSE = 1 is a mean predictor). That is a real but **weak** predictor — enough to say the living state is not just density — not an echo-state-network score. Other tasks in the tree (waveform, Mackey–Glass, Lorenz) are scored separately; some beat Brownian, some do not beat the chemical field.

Mechanisms that did **not** survive as the 300 s dish receiver were isolated rather than forced: a threshold / Danino quorum-sensing operating point that never occupies, vesicle-store artificial cells that do not track the input, and glucose/Monod as a computing channel.

## What a reader will find in the folders

| Area | Role |
|---|---|
| `src/` | BSim engine (agents, fields, physics) |
| `examples/BSimReservoirPlanStage6/` and HybridDish packs | Claim dish, NARMA-10 protocol, ridge, gate evidence |
| `examples/BSimDaninoD1*`, `BSimMonodM*` | Citation-box identities (QS bath, Monod/uptake) — not merged back into the dish |
| `examples/PocketDish/`, later Stage99 / lane packs | Next-device geometry (pocket + bus), not a retune of the claim dish |
| `STAGED_IMPLEMENTATION_PLAN.md`, `COMPARISON_RESERVOIR_NEW.md` | Why the rebuild exists and what was refused |

Build with Ant or an IDE (IntelliJ / Eclipse); Java examples take a `.properties` config. Evaluation scripts are Python. See the root `README.md` for the upstream BSim build notes.

## How to cite BSim itself

Matyjaszkiewicz et al., *ACS Synthetic Biology* (2017), doi:10.1021/acssynbio.7b00121; Gorochowski et al., *PLoS ONE* (2012), doi:10.1371/journal.pone.0042790.
