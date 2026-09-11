#!/usr/bin/env python3
"""PocketOsc-SI bulk delay-DDE. TAKEN SI table. Not D1g. No NARMA.

Equations copied from examples/PocketDish/DANINO_SI_MODEL.md.
ICs ENGINEERING (SI silent): A=I=0, Hi=He=1, Hi(t<0)=1.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
CONFIGS = HERE / "configs"


def load_json(name: str) -> dict:
    with (CONFIGS / name).open(encoding="utf-8") as handle:
        return json.load(handle)


PARAMS = load_json("params.json")
ARMS = {row["id"]: row for row in load_json("arms.json")["arms"]}
COVER_IDS = tuple(load_json("arms.json")["cover_ids"])

C_A = float(PARAMS["C_A"])
C_I = float(PARAMS["C_I"])
DELTA = float(PARAMS["delta"])
ALPHA = float(PARAMS["alpha"])
TAU = float(PARAMS["tau"])
K = float(PARAMS["k"])
K1 = float(PARAMS["k1"])
B = float(PARAMS["b"])
GAMMA_A = float(PARAMS["gamma_A"])
GAMMA_I = float(PARAMS["gamma_I"])
GAMMA_H = float(PARAMS["gamma_H"])
F = float(PARAMS["f"])
G = float(PARAMS["g"])
D0 = float(PARAMS["d0"])
D_MEM = float(PARAMS["D_membrane"])
T_END = float(PARAMS["t_end"])
T_DISCARD = float(PARAMS["t_discard"])
SAMPLE_DT = float(PARAMS["sample_dt"])
PEAK_SPACING = float(PARAMS["peak_min_spacing"])
REL_PROM = float(PARAMS["peak_rel_prominence"])
MIN_PEAKS = int(PARAMS["min_peaks"])
RTOL = float(PARAMS["integrator"]["rtol"])
ATOL = float(PARAMS["integrator"]["atol"])
MAX_STEP = float(PARAMS["integrator"]["max_step"])

ICS = PARAMS["ics"]
Y0 = np.array([float(ICS["A"]), float(ICS["I"]), float(ICS["H_i"]), float(ICS["H_e"])], dtype=float)
HI_HISTORY = float(ICS["H_i_history"])


def print_ics() -> None:
    print(
        "ICS_ENGINEERING: "
        f"A={ICS['A']} I={ICS['I']} Hi={ICS['H_i']} He={ICS['H_e']} "
        f"history_Hi(t<0)={ICS['H_i_history']}",
        flush=True,
    )


def production(h_tau: float) -> float:
    h2 = h_tau * h_tau
    return DELTA + ALPHA * h2 / (1.0 + K1 * h2)


def rhs(y: np.ndarray, h_tau: float, d: float, mu: float) -> np.ndarray:
    a, i, h_i, h_e = np.maximum(y, 0.0)
    p = production(max(h_tau, 0.0))
    dens = 1.0 - (d / D0) ** 4
    denom_f = 1.0 + F * (a + i)
    da = C_A * dens * p - GAMMA_A * a / denom_f
    di = C_I * dens * p - GAMMA_I * i / denom_f
    dhi = B * i / (1.0 + K * i) - (GAMMA_H * a * h_i) / (1.0 + G * a) + D_MEM * (h_e - h_i)
    dhe = -(d / (1.0 - d)) * D_MEM * (h_e - h_i) - mu * h_e
    return np.array([da, di, dhi, dhe], dtype=float)


def integrate_bulk(d: float, mu: float, t_end: float = T_END, sample_dt: float = SAMPLE_DT) -> dict:
    """Constant-delay DDE by method of steps. Bulk: D1=0."""
    segments: list[tuple[float, float, object]] = []

    def hi_delayed(t: float) -> float:
        s = t - TAU
        if s <= 0.0:
            return HI_HISTORY
        for t0, t1, sol in segments:
            if t0 - 1e-14 <= s <= t1 + 1e-14:
                ss = min(max(s, t0), t1)
                return float(max(sol.sol(ss)[2], 0.0))
        raise RuntimeError(f"delay lookup failed at t={t}")

    def fun(t: float, y: np.ndarray) -> np.ndarray:
        return rhs(y, hi_delayed(t), d, mu)

    t_left = 0.0
    y_left = Y0.copy()
    while t_left < t_end - 1e-15:
        t_right = min(t_left + TAU, t_end)
        sol = solve_ivp(
            fun,
            (t_left, t_right),
            y_left,
            method="BDF",
            rtol=RTOL,
            atol=ATOL,
            dense_output=True,
            max_step=MAX_STEP,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        segments.append((t_left, t_right, sol))
        y_left = sol.y[:, -1]
        t_left = t_right

    t_eval = np.arange(0.0, t_end + 0.5 * sample_dt, sample_dt)
    t_eval[-1] = t_end
    y_out = np.zeros((len(t_eval), 4), dtype=float)
    for i, t in enumerate(t_eval):
        if t <= 0.0:
            y_out[i] = Y0
            continue
        found = False
        for t0, t1, sol in segments:
            if t0 - 1e-14 <= t <= t1 + 1e-14:
                y_out[i] = sol.sol(min(max(t, t0), t1))
                found = True
                break
        if not found:
            raise RuntimeError(f"sample lookup failed at t={t}")
    y_out = np.maximum(y_out, 0.0)
    return {"t": t_eval, "Y": y_out, "d": float(d), "mu": float(mu), "D1": 0.0}


def i_peaks(t: np.ndarray, lux_i: np.ndarray) -> dict:
    mask = t >= T_DISCARD - 1e-12
    t_w = t[mask]
    y_w = lux_i[mask]
    if len(y_w) < 3:
        return {"n_peaks": 0, "period": float("nan"), "peak_times": []}
    amp = float(np.max(y_w) - np.min(y_w))
    if not math.isfinite(amp) or amp <= 0.0:
        return {"n_peaks": 0, "period": float("nan"), "peak_times": []}
    prominence = REL_PROM * amp
    distance = max(1, int(math.ceil(PEAK_SPACING / SAMPLE_DT)))
    idx, _ = find_peaks(y_w, prominence=prominence, distance=distance)
    times = [float(t_w[j]) for j in idx]
    n = len(times)
    if n < MIN_PEAKS:
        period = float("nan")
    else:
        period = float(np.mean(np.diff(np.array(times))))
    return {"n_peaks": n, "period": period, "peak_times": times}


def post_means(t: np.ndarray, y: np.ndarray) -> dict:
    mask = t >= T_DISCARD - 1e-12
    z = y[mask]
    return {
        "mean_A": float(np.mean(z[:, 0])),
        "mean_I": float(np.mean(z[:, 1])),
        "mean_Hi": float(np.mean(z[:, 2])),
        "mean_He": float(np.mean(z[:, 3])),
    }
