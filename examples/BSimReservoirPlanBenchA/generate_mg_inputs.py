#!/usr/bin/env python3
"""Generate the frozen Track A Mackey–Glass AHL sequence and target sidecar.

Discrete map used in RC papers (dt=1, tau=17). Deterministic. Do not retune
after seeing NRMSE. Acid hold is the Stage 6 copy already in this package.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

BETA = 0.2
GAMMA = 0.1
N_POW = 10
TAU = 17
X_NEG = 1.2
TRANSIENT_STEPS = 1000
NUM_RECORDED = 201
NUM_WINDOWS = 200
U_HIGH = 0.5


def mackey_glass_series(n_total: int) -> np.ndarray:
    """x[t+1] = x[t] + 0.2 x[t-17] / (1 + x[t-17]**10) - 0.1 x[t].

    For t-17 < 0, x = 1.2. x[0] is 1.2 (history includes the initial state).
    """
    x = np.empty(n_total, dtype=float)
    x[0] = X_NEG
    for t in range(n_total - 1):
        x_tau = x[t - TAU] if t - TAU >= 0 else X_NEG
        x[t + 1] = x[t] + BETA * x_tau / (1.0 + x_tau ** N_POW) - GAMMA * x[t]
    return x


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
    acid_path = directory / "input_acid_held05_200.txt"
    if not acid_path.exists():
        raise SystemExit("input_acid_held05_200.txt missing; run generate_lorenz_inputs.py first")

    series = mackey_glass_series(TRANSIENT_STEPS + NUM_RECORDED)
    recorded = series[TRANSIENT_STEPS:TRANSIENT_STEPS + NUM_RECORDED]
    x_drive = recorded[:NUM_WINDOWS]
    x_next = recorded[1:NUM_WINDOWS + 1]
    u, xmin, xmax = affine_u(x_drive)
    digest = sha256_u(u)

    ahl_path = directory / "input_ahl_mg200.txt"
    target_path = directory / "mg_target.csv"

    ahl_path.write_text(
        "# Frozen Track A Mackey-Glass AHL input. Affine x[0..199] onto [0, 0.5].\n"
        "# x[t+1] = x[t] + 0.2 x[t-17]/(1+x[t-17]**10) - 0.1 x[t], x=1.2 for t-17<0.\n"
        "# 1000 transient steps, then 201 recorded samples. Do not regenerate after NRMSE.\n"
        + "\n".join(f"{value:.12f}" for value in u)
        + "\n",
        encoding="utf-8",
    )
    rows = ["n;u;x;x_next"]
    rows.extend(
        f"{index};{u[index]:.12f};{x_drive[index]:.12f};{x_next[index]:.12f}"
        for index in range(NUM_WINDOWS)
    )
    target_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"windows={NUM_WINDOWS}")
    print(f"beta={BETA} gamma={GAMMA} n={N_POW} tau={TAU} x_neg={X_NEG}")
    print(f"transient_steps={TRANSIENT_STEPS}")
    print(f"xmin={xmin:.12f}")
    print(f"xmax={xmax:.12f}")
    print(f"u_min={float(np.min(u)):.12f} u_max={float(np.max(u)):.12f}")
    print(f"u_sha256={digest}")
    print(f"x_next_min={float(np.min(x_next)):.12f} x_next_max={float(np.max(x_next)):.12f}")
    print(f"ahl={ahl_path.name}")
    print(f"acid={acid_path.name}")
    print(f"target={target_path.name}")


if __name__ == "__main__":
    main()
