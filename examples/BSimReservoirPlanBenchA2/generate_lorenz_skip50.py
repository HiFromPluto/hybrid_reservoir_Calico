#!/usr/bin/env python3
"""Generate the frozen Track A2 subsampled Lorenz'63 AHL sequence and target.

Deterministic ODE. SKIP=50 RK4 steps between samples is frozen. Do not retune
sigma/rho/beta, dt, IC, transient, or SKIP after seeing NRMSE. Acid hold is
copied from Stage 6, not regenerated.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import numpy as np

SIGMA = 10.0
RHO = 28.0
BETA = 8.0 / 3.0
DT_L = 0.02
SKIP = 50
IC = (1.0, 1.0, 1.0)
TRANSIENT_STEPS = 5000
NUM_RECORDED = 201
NUM_WINDOWS = 200
U_HIGH = 0.5


def lorenz_deriv(state: np.ndarray) -> np.ndarray:
    x, y, z = state
    return np.array([
        SIGMA * (y - x),
        x * (RHO - z) - y,
        x * y - BETA * z,
    ], dtype=float)


def rk4_step(state: np.ndarray, dt: float) -> np.ndarray:
    k1 = lorenz_deriv(state)
    k2 = lorenz_deriv(state + 0.5 * dt * k1)
    k3 = lorenz_deriv(state + 0.5 * dt * k2)
    k4 = lorenz_deriv(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate() -> np.ndarray:
    """Sample i is the state after 5000 + i*SKIP RK4 steps from the IC."""
    state = np.array(IC, dtype=float)
    for _ in range(TRANSIENT_STEPS):
        state = rk4_step(state, DT_L)
    samples = np.empty((NUM_RECORDED, 3), dtype=float)
    samples[0] = state
    for index in range(1, NUM_RECORDED):
        for _ in range(SKIP):
            state = rk4_step(state, DT_L)
        samples[index] = state
    return samples


def affine_u(x: np.ndarray) -> tuple[np.ndarray, float, float]:
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if xmax == xmin:
        raise SystemExit("generator is wrong: xmax == xmin")
    u = U_HIGH * (x - xmin) / (xmax - xmin)
    if np.any(u < 0.0) or np.any(u > U_HIGH):
        raise SystemExit("u escaped [0, 0.5]")
    return u, xmin, xmax


def sha256_u(u: np.ndarray) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def main() -> None:
    directory = Path(__file__).resolve().parent
    stage6_acid = directory.parent / "BSimReservoirPlanStage6" / "input_acid_held05_200.txt"
    samples = integrate()
    x_drive = samples[:NUM_WINDOWS, 0]
    y_next = samples[1:NUM_WINDOWS + 1, 1]
    u, xmin, xmax = affine_u(x_drive)
    digest = sha256_u(u)

    ahl_path = directory / "input_ahl_lorenz_skip50.txt"
    acid_path = directory / "input_acid_held05_200.txt"
    target_path = directory / "lorenz_skip50_target.csv"

    ahl_path.write_text(
        "# Frozen Track A2 subsampled Lorenz'63 AHL input. Affine x[0..199] onto [0, 0.5].\n"
        "# sigma=10, rho=28, beta=8/3, RK4 dt=0.02, IC=(1,1,1), discard 5000 steps.\n"
        "# SKIP=50 RK4 steps between samples (dt_sample=1.0). Do not regenerate after NRMSE.\n"
        + "\n".join(f"{value:.12f}" for value in u)
        + "\n",
        encoding="utf-8",
    )
    shutil.copyfile(stage6_acid, acid_path)
    rows = ["n;u;x;y_next"]
    rows.extend(
        f"{index};{u[index]:.12f};{x_drive[index]:.12f};{y_next[index]:.12f}"
        for index in range(NUM_WINDOWS)
    )
    target_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"windows={NUM_WINDOWS}")
    print(f"sigma={SIGMA} rho={RHO} beta={BETA}")
    print(f"dt_L={DT_L} skip={SKIP} dt_sample={SKIP * DT_L} ic={IC}")
    print(f"transient_steps={TRANSIENT_STEPS}")
    print(f"xmin={xmin:.12f}")
    print(f"xmax={xmax:.12f}")
    print(f"u_min={float(np.min(u)):.12f} u_max={float(np.max(u)):.12f}")
    print(f"u_sha256={digest}")
    print(f"y_next_min={float(np.min(y_next)):.12f} y_next_max={float(np.max(y_next)):.12f}")
    print(f"ahl={ahl_path.name}")
    print(f"acid={acid_path.name}")
    print(f"target={target_path.name}")


if __name__ == "__main__":
    main()
