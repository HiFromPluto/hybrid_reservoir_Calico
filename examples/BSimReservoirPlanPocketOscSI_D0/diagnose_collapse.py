#!/usr/bin/env python3
"""Diagnosis: moderate-AHL kick from mu=0.28 to 1.5. No retune."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from danino_si_dde import DaninoSIDDE, load_json
from period_check import extract_period, refuse_narma

refuse_narma()
protocol = load_json("protocol.json")
model = DaninoSIDDE()
y0 = np.array([0.0, 1.0, 0.2, 0.2])
print("mu flag n_peaks period mean_I I_end p5 p95 rel_amp persist", flush=True)
for mu in (0.28, 0.32, 0.36, 0.40, 0.50, 0.60, 0.80, 1.00, 1.20, 1.50):
    traj = model.integrate(0.5, mu, y0, 0.2, 1000.0, float(protocol["sample_dt"]))
    lux = traj["Y"][:, 1]
    p = extract_period(traj["t"], lux, protocol)
    post = traj["t"] >= float(protocol["t_discard"]) - 1e-12
    per = p["period"]
    ps = f"{per:.2f}" if per == per else "nan"
    print(
        f"{mu:4.2f} {p['flag']:10} {p['n_peaks']:3} {ps:>7} "
        f"{p['mean_I']:.3e} {lux[-1]:.3e} "
        f"{np.percentile(lux[post], 5):.3e} {np.percentile(lux[post], 95):.3e} "
        f"{p.get('rel_amplitude', float('nan')):.3e} "
        f"{p.get('amp_persist', float('nan')):.3e}",
        flush=True,
    )
