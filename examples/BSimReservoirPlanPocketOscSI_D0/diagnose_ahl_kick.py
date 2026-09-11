#!/usr/bin/env python3
"""Diagnosis: AHL-kick ICs at mu=0.4, d=0.5. No retune."""

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
print("IC flag n_peaks period I_max I_end rel_amp persist", flush=True)
ics = [
    ("I1_H0", np.array([0.0, 1.0, 0.0, 0.0]), 0.0),
    ("I0_H0.05", np.array([0.0, 0.0, 0.05, 0.05]), 0.05),
    ("I0_H0.1", np.array([0.0, 0.0, 0.10, 0.10]), 0.10),
    ("I0_H0.2", np.array([0.0, 0.0, 0.20, 0.20]), 0.20),
    ("I1_H0.2", np.array([0.0, 1.0, 0.20, 0.20]), 0.20),
    ("I0_H0.01", np.array([0.0, 0.0, 0.01, 0.01]), 0.01),
]
for name, y0, hist in ics:
    traj = model.integrate(0.5, 0.4, y0, hist, 1000.0, float(protocol["sample_dt"]))
    lux = traj["Y"][:, 1]
    p = extract_period(traj["t"], lux, protocol)
    per = p["period"]
    ps = f"{per:.2f}" if per == per else "nan"
    print(
        f"{name:10} {p['flag']:10} {p['n_peaks']:3} {ps:>7} "
        f"{np.max(lux):.3e} {lux[-1]:.3e} "
        f"{p.get('rel_amplitude', float('nan')):.3e} "
        f"{p.get('amp_persist', float('nan')):.3e}",
        flush=True,
    )
