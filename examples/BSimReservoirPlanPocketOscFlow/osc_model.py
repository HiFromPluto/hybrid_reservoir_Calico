#!/usr/bin/env python3
"""PocketOsc-Flow 0-D mean-field. D50 ODE + frozen mu_bus(v). No NARMA."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
CONFIGS = HERE / "configs"

TIME_ADJ = 60.0
QS_DELTA1 = 0.8487 / TIME_ADJ
QS_DELTA2 = 0.0234 / TIME_ADJ
QS_G = 0.0412
QS_KP2 = 9.0 / TIME_ADJ
QS_KR1OFF = 6e-6 / TIME_ADJ
QS_KR1ON = 5.99e-5 / TIME_ADJ
QS_KCAT_AIIA = 2631.4 / TIME_ADJ
QS_T_A = 0.00276 / TIME_ADJ
QS_T_LA = 0.024 / TIME_ADJ
QS_A0LI = 7.785e-6 / TIME_ADJ
QS_A0AA = 6.183e-6 / TIME_ADJ
QS_KPLI = 0.9 / TIME_ADJ
QS_KPAA = 0.9 / TIME_ADJ
QS_KMLA = 0.01
QS_KMAA = 1200.0
QS_LTOT = 15.0
QS_N = 2.0
CELL_WALL_DIFF = 3.0 / TIME_ADJ
MU_GROWTH = math.log(2.0) / 1800.0
K_DANINO = 2.76e-3 / 60.0

PRIMARY_RTOL, PRIMARY_ATOL = 1e-6, 1e-9


def load_json(name: str) -> dict:
    with (CONFIGS / name).open(encoding="utf-8") as handle:
        return json.load(handle)


DISH = load_json("dish.json")
ARMS = {row["id"]: row for row in load_json("arms.json")["arms"]}

CONV = float(DISH["molecules_per_um3_per_uM"])
CELL_VOL = float(DISH["cell_vol_um3"])
D0 = float(DISH["D_um2_s"])
FOOTPRINT = float(DISH["footprint_um"][0]) * float(DISH["footprint_um"][1])
LX = float(DISH["footprint_um"][0])
LY = float(DISH["footprint_um"][1])
T_END = float(DISH["t_end_s"])
T_DISCARD = float(DISH["t_discard_s"])
OCC_WINDOW = float(DISH["occupancy_window_s"])
PEAK_SPACING = float(DISH["peak_min_spacing_s"])
PEAK_PROM = float(DISH["peak_prominence_uM"])
D50_DURATION = 14400.0


def hill_qs(x: np.ndarray | float) -> np.ndarray | float:
    x = np.maximum(x, 0.0)
    kn = QS_KMLA ** QS_N
    cn = np.power(x, QS_N)
    return cn / (kn + cn + 1e-30)


def occupancy_flag(mean_h: float) -> str:
    if not math.isfinite(mean_h):
        return "NA"
    if mean_h > 0.95:
        return "SATURATED"
    if mean_h >= 0.05:
        return "ALIVE"
    return "DEAD"


def packing(d: float, height_um: float) -> dict:
    d = float(d)
    if not (0.0 < d < 1.0):
        raise ValueError(f"d must be in (0,1), got {d}")
    v_trap = FOOTPRINT * float(height_um)
    ratio = d / (1.0 - d)
    a_open = LX * float(height_um)
    v_ext = (1.0 - d) * v_trap
    return {
        "d": d,
        "d_over_one_minus_d": ratio,
        "height_um": float(height_um),
        "V_trap_um3": v_trap,
        "V_ext_um3": v_ext,
        "A_open_um2": a_open,
        "N_implied_engineering": d * v_trap / CELL_VOL,
    }


def mu_bus(v_um_s: float, pack: dict) -> float:
    """ENGINEERING two-resistance leak. mu(0)=0. Not fitted to 55/90 min."""
    v = float(v_um_s)
    if v <= 0.0:
        return 0.0
    return pack["A_open_um2"] / pack["V_ext_um3"] * (v * D0) / (v * LY + D0)


def danino_rhs(y: np.ndarray, ahl_ext: float) -> np.ndarray:
    y = np.maximum(y, 0.0)
    extra = y[1] - ahl_ext
    hill_la = hill_qs(y[3])
    dy = np.zeros(4)
    dy[0] = (
        QS_A0LI
        + QS_KPLI * hill_la
        - (QS_DELTA1 * y[0]) / (QS_G * (y[0] + y[2]) + 1.0)
        - MU_GROWTH * y[0]
    )
    dy[1] = (
        QS_KP2 * y[0]
        - QS_KR1ON * (QS_LTOT - y[3]) * y[1]
        + QS_KR1OFF * y[3]
        - (QS_KCAT_AIIA * y[2] * y[1]) / (QS_KMAA + y[1])
        - QS_T_A * y[1]
        - CELL_WALL_DIFF * extra
    )
    dy[2] = (
        QS_A0AA
        + QS_KPAA * hill_la
        - (QS_DELTA2 * y[2]) / (QS_G * (y[0] + y[2]) + 1.0)
        - MU_GROWTH * y[2]
    )
    dy[3] = (
        QS_KR1ON * (QS_LTOT - y[3]) * y[1]
        - QS_KR1OFF * y[3]
        - QS_T_LA * y[3]
        - MU_GROWTH * y[3]
    )
    return dy


def integrate_arm(d: float, height_um: float, v_um_s: float, duration: float, sample_dt: float = 60.0) -> dict:
    pack = packing(d, height_um)
    ratio = pack["d_over_one_minus_d"]
    mu = mu_bus(v_um_s, pack)
    k = K_DANINO

    def fun(_t: float, state: np.ndarray) -> np.ndarray:
        c = max(float(state[0]), 0.0)
        y = np.maximum(state[1:], 0.0)
        dy = danino_rhs(y, c)
        dc = ratio * CELL_WALL_DIFF * (y[1] - c) - k * c - mu * c
        return np.concatenate(([dc], dy))

    t_eval = np.arange(0.0, duration + 0.5 * sample_dt, sample_dt)
    t_eval[-1] = duration
    sol = solve_ivp(
        fun,
        (0.0, duration),
        np.zeros(5),
        method="BDF",
        rtol=PRIMARY_RTOL,
        atol=PRIMARY_ATOL,
        t_eval=t_eval,
        max_step=30.0,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return {
        "t": sol.t,
        "C": sol.y[0],
        "Y": sol.y[1:].T,
        "sol": sol,
        "pack": pack,
        "k_decay": k,
        "mu_bus": mu,
        "v_um_s": float(v_um_s),
    }


def la_peaks(t: np.ndarray, la: np.ndarray) -> dict:
    mask = t >= T_DISCARD
    t_w = t[mask]
    y_w = la[mask]
    if len(y_w) < 3:
        return {"n_peaks": 0, "period_s": float("nan"), "peak_times_s": []}
    cand = []
    for i in range(1, len(y_w) - 1):
        if y_w[i] > y_w[i - 1] and y_w[i] >= y_w[i + 1] and y_w[i] >= PEAK_PROM:
            if not cand or (t_w[i] - cand[-1]) >= PEAK_SPACING:
                cand.append(float(t_w[i]))
    n = len(cand)
    if n < 3:
        period = float("nan")
    else:
        period = float(np.mean(np.diff(np.array(cand))))
    return {"n_peaks": n, "period_s": period, "peak_times_s": cand}


def field_budget(out: dict, duration: float, dt: float = 0.25) -> dict:
    pack = out["pack"]
    v_ext = pack["V_ext_um3"]
    ratio = pack["d_over_one_minus_d"]
    k = out["k_decay"]
    mu = out["mu_bus"]
    t = np.arange(0.0, duration + 0.5 * dt, dt)
    t[-1] = duration
    z = out["sol"].sol(t)
    c = np.maximum(z[0], 0.0)
    hi = np.maximum(z[2], 0.0)
    mass = c * v_ext * CONV
    flux_mol = v_ext * CONV * ratio * CELL_WALL_DIFF * (hi - c)
    decay_rate = k * mass
    leak_rate = mu * mass
    membrane = float(np.trapezoid(flux_mol, t))
    decay = float(np.trapezoid(decay_rate, t))
    leak = float(np.trapezoid(leak_rate, t))
    initial = float(mass[0])
    remaining = float(mass[-1])
    residual = initial + membrane - remaining - decay - leak
    dominant = max(abs(initial), abs(membrane), abs(remaining), abs(decay), abs(leak), 1.0)
    frac = abs(residual) / dominant
    c_nn = bool(np.min(c) >= -1e-12)
    mass_ok = frac <= 0.01 and c_nn
    return {
        "initial_molecules": initial,
        "membrane_molecules": membrane,
        "remaining_molecules": remaining,
        "decay_molecules": decay,
        "leak_molecules": leak,
        "residual_molecules": residual,
        "residual_frac": frac,
        "C_nonnegative": c_nn,
        "mass_ok": mass_ok,
        "mass_class": "LE1PCT" if mass_ok else "FAIL",
    }
