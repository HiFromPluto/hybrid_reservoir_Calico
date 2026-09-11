#!/usr/bin/env python3
"""Per-trajectory C1 NRMSE and delta plots. Keep every trajectory visible."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "results" / "c1_deltas.csv"
NARMA_CSV = HERE / "results" / "c1_narma.csv"
FIG_DIR = HERE / "results" / "figures"


def load_primary(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["seed"] == "111"]
    rows.sort(key=lambda row: row["traj"])
    return rows


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    delta_rows = load_primary(CSV_PATH)
    narma_rows = load_primary(NARMA_CSV)
    if len(delta_rows) != 11 or len(narma_rows) != 11:
        raise SystemExit("expected 11 seed-111 rows")
    traj = [row["traj"] for row in narma_rows]
    x = np.arange(len(traj))

    driven = np.array([float(row["driven_nrmse"]) for row in narma_rows])
    brown = np.array([float(row["brownian_nrmse"]) for row in narma_rows])
    silent = np.array([float(row["silent_nrmse"]) for row in narma_rows])
    field = np.array([float(row["field_nrmse"]) for row in narma_rows])
    tap = np.array([float(row["tap10_nrmse"]) for row in narma_rows])

    fig, ax = plt.subplots(figsize=(11, 5.5))
    for values, label, marker in (
        (driven, "driven F408", "o"),
        (brown, "Brownian Den", "s"),
        (silent, "silent F408", "D"),
        (field, "field AHL", "^"),
        (tap, "legal 10-tap u", "x"),
    ):
        ax.plot(x, values, marker=marker, linewidth=1.5, label=label)
    ax.set_xticks(x, traj)
    ax.set_xlabel("independent input trajectory")
    ax.set_ylabel("test NRMSE")
    ax.set_title("C1 NARMA-10 per trajectory (seed 111)")
    ax.legend(frameon=False, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "c1_nrmse_per_traj.png", dpi=150)
    plt.close(fig)

    db = np.array([float(row["delta_B"]) for row in delta_rows])
    ds = np.array([float(row["delta_S"]) for row in delta_rows])
    df = np.array([float(row["delta_F"]) for row in delta_rows])
    fig, ax = plt.subplots(figsize=(11, 5.5))
    width = 0.25
    ax.bar(x - width, db, width, label=r"$\Delta_B$ Brownian − driven")
    ax.bar(x, ds, width, label=r"$\Delta_S$ silent − driven")
    ax.bar(x + width, df, width, label=r"$\Delta_F$ field − driven")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x, traj)
    ax.set_xlabel("independent input trajectory")
    ax.set_ylabel("NRMSE delta (positive favours biology)")
    ax.set_title("C1 deltas per trajectory (seed 111)")
    ax.legend(frameon=False)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "c1_deltas_per_traj.png", dpi=150)
    plt.close(fig)
    print(f"wrote {FIG_DIR / 'c1_nrmse_per_traj.png'}")
    print(f"wrote {FIG_DIR / 'c1_deltas_per_traj.png'}")


if __name__ == "__main__":
    main()
