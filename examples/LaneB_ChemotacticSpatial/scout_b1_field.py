#!/usr/bin/env python3
"""B1 field-only scout. No cells. Pulse t_off=60 s then J=0.

Predicts kappa_field(t) on the N0 millimetre dish (100x50x1, D=800, k=0,
J_max=2e4, four A0 sources). Sets T_hold = first sample after t_off with
|kappa_field| < 0.15, else FAIL at the 400 s cap.

Does not raise D. Does not add hydrolysis. Not B0. Not NARMA.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL_JSON = HERE / "configs" / "protocol_b1.json"
PROTOCOL_MD = HERE / "PROTOCOL_PAINT_HOLD.md"

LX, LY, LZ = 1000.0, 500.0, 10.0
NX, NY, NZ = 100, 50, 1
DX, DY, DZ = LX / NX, LY / NY, LZ / NZ
D = 800.0
J_MAX = 20000.0
T_OFF = 60.0
DT = 0.02
SAMPLE_DT = 0.10
CAP = 400.0
KAPPA_HOLD = 0.15
EPS = 1e-12
ACS = (
    (200.0, 150.0, "att"),
    (200.0, 350.0, "att"),
    (800.0, 150.0, "rep"),
    (800.0, 350.0, "rep"),
)


def pde_i(x: float) -> int:
    i = int(x / DX)
    return max(0, min(NX - 1, i))


def pde_j(y: float) -> int:
    j = int(y / DY)
    return max(0, min(NY - 1, j))


def diffuse_2d(q: np.ndarray, h: float) -> np.ndarray:
    """Conservative no-flux FTCS on molecules/voxel (NZ=1, z faces no-flux)."""
    # Uniform voxels: F = D (q_i - q_{i+1}) / dx^2 molecules/s (C = q/V).
    fx = D * (q[:-1, :] - q[1:, :]) / (DX * DX)
    fy = D * (q[:, :-1] - q[:, 1:]) / (DY * DY)
    out = q.copy()
    out[:-1, :] -= h * fx
    out[1:, :] += h * fx
    out[:, :-1] -= h * fy
    out[:, 1:] += h * fy
    return out


def kappa_field(att: np.ndarray, rep: np.ndarray) -> float:
    mid = NX // 2
    s_l = att[:mid, :].sum() + rep[mid:, :].sum()
    s_r = att[mid:, :].sum() + rep[:mid, :].sum()
    return float((s_l - s_r) / (s_l + s_r + EPS))


def mean_abs_l(att: np.ndarray, rep: np.ndarray) -> float:
    vol = DX * DY * DZ
    return float(np.mean(np.abs(att / vol - rep / vol)))


def overlap(t_start: float, dt: float, t_on: float, t_off: float) -> float:
    lo = max(t_start, t_on)
    hi = min(t_start + dt, t_off)
    return max(0.0, hi - lo)


def deposit(att: np.ndarray, rep: np.ndarray, t_start: float, dt: float) -> None:
    dm = J_MAX * overlap(t_start, dt, 0.0, T_OFF)
    if dm <= 0.0:
        return
    for x, y, payload in ACS:
        i, j = pde_i(x), pde_j(y)
        if payload == "att":
            att[i, j] += dm
        else:
            rep[i, j] += dm


def subcycle_h(duration: float) -> float:
    gate = D * duration * (1.0 / DX**2 + 1.0 / DY**2 + 1.0 / DZ**2)
    n = max(1, int(np.ceil(gate / 0.5)))
    return duration / n


def advance(q: np.ndarray, duration: float) -> np.ndarray:
    h = subcycle_h(duration)
    n = int(round(duration / h))
    for _ in range(n):
        q = diffuse_2d(q, h)
    return q


def main() -> int:
    refuse = " ".join(sys.argv).lower()
    if any(k in refuse for k in ("narma", "charc", "ipc")):
        print("scout refuses NARMA/CHARC/IPC")
        return 2

    att = np.zeros((NX, NY), dtype=np.float64)
    rep = np.zeros((NX, NY), dtype=np.float64)
    times: list[float] = []
    kappas: list[float] = []
    abs_l: list[float] = []

    n_steps = int(round(CAP / DT))
    times.append(0.0)
    kappas.append(kappa_field(att, rep))
    abs_l.append(mean_abs_l(att, rep))

    half = 0.5 * DT
    sample_stride = int(round(SAMPLE_DT / DT))
    for step in range(n_steps):
        t_start = step * DT
        att = advance(att, half)
        rep = advance(rep, half)
        deposit(att, rep, t_start, DT)
        att = advance(att, half)
        rep = advance(rep, half)
        t_end = (step + 1) * DT
        if (step + 1) % sample_stride == 0:
            times.append(t_end)
            kappas.append(kappa_field(att, rep))
            abs_l.append(mean_abs_l(att, rep))

    t_arr = np.asarray(times)
    k_arr = np.asarray(kappas)
    t_hold = None
    k_at_hold = None
    for t, k in zip(t_arr, k_arr):
        if t <= T_OFF + 1e-12:
            continue
        if abs(k) < KAPPA_HOLD:
            t_hold = float(t)
            k_at_hold = float(k)
            break

    k_at_off = float(k_arr[np.argmin(np.abs(t_arr - T_OFF))])
    pulse = (t_arr > 0.0) & (t_arr <= T_OFF + 1e-12)
    mean_l_pulse = float(np.mean(np.asarray(abs_l)[pulse])) if np.any(pulse) else float("nan")
    k_at_cap = float(k_arr[-1])

    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / "b1_field_scout.csv"
    with csv_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("t,kappa_field,mean_abs_L\n")
        for t, k, ell in zip(times, kappas, abs_l):
            handle.write(f"{t:.6g},{k:.8g},{ell:.8g}\n")

    scout_ok = t_hold is not None
    summary = {
        "object": "LANE_B_CHEMOTACTIC_SPATIAL",
        "gate": "B1_PAINT_HOLD",
        "device": "LANE_B_MM_DISH_FLOW0",
        "scout": "field_only",
        "t_off_s": T_OFF,
        "observe_cap_s": CAP,
        "kappa_field_hold_max": KAPPA_HOLD,
        "kappa_field_at_t_off": k_at_off,
        "kappa_field_at_cap": k_at_cap,
        "mean_abs_L_during_pulse": mean_l_pulse,
        "T_hold_s": t_hold,
        "kappa_field_at_T_hold": k_at_hold,
        "scout_pass": scout_ok,
        "hydrolysis_added": False,
        "D_raised": False,
        "narma": False,
        "cells": False,
    }
    (RESULTS / "b1_field_scout.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    proto = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if scout_ok:
        proto["T_hold_s"] = t_hold
        proto["T_hold_status"] = "SCOUT_SET_BEFORE_MOTILE_TRACES"
        proto["kappa_field_scout_at_T_hold"] = k_at_hold
        proto["kappa_field_scout_at_t_off"] = k_at_off
    else:
        proto["T_hold_s"] = None
        proto["T_hold_status"] = "SCOUT_FAIL_CAP_400"
        proto["kappa_field_scout_at_cap"] = k_at_cap
    PROTOCOL_JSON.write_text(json.dumps(proto, indent=2) + "\n", encoding="utf-8")

    print(
        f"B1 field scout t_off={T_OFF:g}s kappa_field(t_off)={k_at_off:.6g} "
        f"mean_|L|_pulse={mean_l_pulse:.6g}"
    )
    if scout_ok:
        print(
            f"T_hold={t_hold:.6g} s  kappa_field(T_hold)={k_at_hold:.6g}  "
            f"SCOUT=PASS  (wrote protocol_b1.json)"
        )
        return 0
    print(
        f"T_hold unset. |kappa_field| never < {KAPPA_HOLD} by {CAP:g} s "
        f"(kappa_field(cap)={k_at_cap:.6g}). SCOUT=FAIL. Stop. "
        "Do not add hydrolysis. Do not raise D. Do not start Java paint arms."
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
