# C1 — independent-input NARMA-10 on the frozen claim dish

Track C confirmatory package. Frozen HybridDish / Narma10b claim dish.
`GATE_EVIDENCE.md` is not edited. Kinetics, layout, and flow are not
retuned. Claim dish stays CENTER / `FLOW=0`.

C1 asks whether driven biology beats Brownian, silent, and field-only
on **new independent NARMA-10 input trajectories**. Cells are not
required to beat the legal 10-tap baseline (~0.68). Bacterial seeds
are not independent inputs.

Traj 00 is the existing Narma10b `u` (SHA-256
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`)
and is scored in place. Traj 01–10 are new. Interaction subset
(frozen before scores): traj 01 and 02 also run seeds 222 and 333.

## Status

**System HOLD** on seed 111, traj 00–10 (n=11 independent inputs).
Mean ΔB `0.294` (9/11, bootstrap 95% CI `[0.155, 0.416]`), mean ΔS
`0.302` (9/11, CI `[0.152, 0.439]`). Occupancy **ALIVE** on every
primary trajectory. Traj 06 and 10 were not dropped.

Living-layer diagnostic: mean ΔF `0.051` (8/11); the ΔF CI includes 0.
10-tap remains better than cells on this task and is **not** a C1 kill
gate. Narma10b Overall is unchanged. Claim dish is not replaced.

See `PROTOCOL.md`, `results/U_FREEZE.md`, and `results/C1_SCOUT.md`.

Do not start E5, two-way, IPC, Stage99, WaveformS, or layout-flow
from this package.
