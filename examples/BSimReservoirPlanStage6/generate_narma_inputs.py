#!/usr/bin/env python3
"""Generate the frozen Stage 6 NARMA-10 AHL sequence and constant acid hold.

This script is the sequence provenance. Re-running it with the same seed must
reproduce the committed input files. Do not retune the seed after seeing NRMSE.
"""

import hashlib
import random
from pathlib import Path

SEQUENCE_SEED = 20260814
NUM_WINDOWS = 200
U_LOW = 0.0
U_HIGH = 0.5
ACID_HOLD = 0.5


def narma10(u):
    """y[n+1] = 0.3 y[n] + 0.05 y[n] sum_{i=0..9} y[n-i] + 1.5 u[n-9] u[n] + 0.1

    y[0] = 0. Target for window n is y[n+1], which depends on the u[n] just applied.
    """
    y = [0.0] * (len(u) + 1)
    for t, u_t in enumerate(u):
        acc = sum(y[t - i] if t - i >= 0 else 0.0 for i in range(10))
        u_lag = u[t - 9] if t >= 9 else 0.0
        y[t + 1] = 0.3 * y[t] + 0.05 * y[t] * acc + 1.5 * u_lag * u_t + 0.1
    return y


def main():
    directory = Path(__file__).resolve().parent
    rng = random.Random(SEQUENCE_SEED)
    u = [U_LOW + (U_HIGH - U_LOW) * rng.random() for _ in range(NUM_WINDOWS)]
    y = narma10(u)
    digest = hashlib.sha256(
        (",".join(f"{value:.12f}" for value in u)).encode("ascii")
    ).hexdigest()

    ahl_path = directory / "input_ahl_narma200.txt"
    acid_path = directory / "input_acid_held05_200.txt"
    target_path = directory / "narma10_target.csv"

    ahl_path.write_text(
        "# Frozen Stage 6 NARMA-10 AHL input. u ~ Uniform[0, 0.5], seed=20260814.\n"
        "# Standard NARMA-10 range; subset of [0,1]. Do not regenerate after seeing NRMSE.\n"
        + "\n".join(f"{value:.12f}" for value in u)
        + "\n",
        encoding="utf-8",
    )
    acid_path.write_text(
        "# Frozen Stage 6 acid hold. Every window is 0.5 (option a). Not a second NARMA channel.\n"
        + "\n".join(f"{ACID_HOLD:.2f}" for _ in u)
        + "\n",
        encoding="utf-8",
    )
    rows = ["window;u;y_next"]
    rows.extend(
        f"{index};{u[index]:.12f};{y[index + 1]:.12f}" for index in range(NUM_WINDOWS)
    )
    target_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"windows={NUM_WINDOWS}")
    print(f"sequence_seed={SEQUENCE_SEED}")
    print(f"u_range=[{U_LOW}, {U_HIGH}]")
    print(f"u_sha256={digest}")
    print(f"y_min={min(y[1:]):.6f} y_max={max(y[1:]):.6f}")
    print(f"ahl={ahl_path.name}")
    print(f"acid={acid_path.name}")
    print(f"target={target_path.name}")


if __name__ == "__main__":
    main()
