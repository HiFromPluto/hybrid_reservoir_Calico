#!/usr/bin/env python3
"""Post-FAIL diagnosis only. Does not retune. Does not reopen D0 as PASS.

Question: at TAKEN parameters and frozen primary IC, is there a bulk
oscillatory window between the mu=0 ON latch and the mu>=0.25 OFF state?
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from danino_si_dde import DaninoSIDDE, load_json
from period_check import extract_period, refuse_narma

refuse_narma()

OUT = HERE / "results" / "diagnosis_mu_window.csv"
protocol = load_json("protocol.json")
ensemble = load_json("ensemble.json")
primary = next(a for a in ensemble["arms"] if a["id"] == ensemble["primary_id"])
kick = next(a for a in ensemble["arms"] if a["id"] == "MODERATE_AHL")
model = DaninoSIDDE()

MUS = [
    0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10,
    0.12, 0.15, 0.18, 0.20, 0.22, 0.25,
]
DS = [0.5]


def run_arm(arm: dict, d: float, mu: float, t_end: float) -> dict:
    y0 = np.array([arm["A"], arm["I"], arm["H_i"], arm["H_e"]], dtype=float)
    hist = float(arm["H_i_history"])
    traj = model.integrate(d, mu, y0, hist, t_end, float(protocol["sample_dt"]))
    t, y = traj["t"], traj["Y"]
    lux = y[:, 1]
    peaks = extract_period(t, lux, protocol)
    post = t >= float(protocol["t_discard"]) - 1e-12
    return {
        "arm": arm["id"],
        "d": d,
        "mu": mu,
        "t_end": t_end,
        "flag": peaks["flag"],
        "period": peaks["period"],
        "n_peaks": peaks["n_peaks"],
        "mean_I": peaks["mean_I"],
        "I_max": float(np.max(lux)),
        "I_max_t": float(t[int(np.argmax(lux))]),
        "I_end": float(lux[-1]),
        "Hi_end": float(y[-1, 2]),
        "post_p5": float(np.percentile(lux[post], 5)) if np.any(post) else float("nan"),
        "post_p95": float(np.percentile(lux[post], 95)) if np.any(post) else float("nan"),
        "rel_amplitude": peaks.get("rel_amplitude", float("nan")),
        "amp_persist": peaks.get("amp_persist", float("nan")),
    }


rows = []
print("arm d mu t_end flag period n_peaks I_max I_end mean_I", flush=True)
for d in DS:
    for mu in MUS:
        for arm in (primary, kick):
            row = run_arm(arm, d, mu, float(protocol["t_end"]))
            rows.append(row)
            per = row["period"]
            per_s = f"{per:.2f}" if per == per else "nan"
            print(
                f"{row['arm']:16} {d:.2f} {mu:5.2f} {row['t_end']:.0f} "
                f"{row['flag']:10} {per_s:>7} {row['n_peaks']:3} "
                f"{row['I_max']:.3e} {row['I_end']:.3e} {row['mean_I']:.3e}",
                flush=True,
            )

# Longer horizon at a few interior points, primary only.
for mu in (0.05, 0.10, 0.15):
    row = run_arm(primary, 0.5, mu, 4000.0)
    rows.append(row)
    per = row["period"]
    per_s = f"{per:.2f}" if per == per else "nan"
    print(
        f"{row['arm']:16} 0.50 {mu:5.2f} {row['t_end']:.0f} "
        f"{row['flag']:10} {per_s:>7} {row['n_peaks']:3} "
        f"{row['I_max']:.3e} {row['I_end']:.3e} {row['mean_I']:.3e}",
        flush=True,
    )

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
print(f"wrote {OUT}", flush=True)
print("DIAGNOSIS ONLY. D0 remains FAIL. No retune.", flush=True)
