#!/usr/bin/env python3
"""B3 field-only scout. No cells. Distributed slab current, then J=0.

Solves J_slab so median |L(x+ell)-L(x)| at t=25 s is >= 0.3 (ell=17.2 um).
FAIL if mean |L| at 25 s > 100. Sets T_hold after pulse-off. Cap 40 s.
Does not drop D. Does not add hydrolysis. Not B0/B1/B2. Not NARMA.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL_JSON = HERE / "configs" / "protocol_b3.json"

LX, LY, LZ = 1000.0, 500.0, 10.0
NX, NY, NZ = 100, 50, 1
DX, DY, DZ = LX / NX, LY / NY, LZ / NZ
VOL = DX * DY * DZ
D = 800.0
K = 0.0
DT = 0.02
SAMPLE_DT = 0.10
SLAB = 100.0
ELL = 17.2
DELTA_L_TARGET = 0.3
MEAN_L_CAP = 100.0
T_LIVE = 25.0
KAPPA_HOLD = 0.15
T_HOLD_CAP = 40.0
EPS = 1e-12
X_SAMPLE = 1.0


def even_att_slab(x: float) -> bool:
    return (np.floor(x / SLAB) % 2.0) == 0.0


def even_mask() -> np.ndarray:
    xs = (np.arange(NX) + 0.5) * DX
    return np.array([even_att_slab(x) for x in xs], dtype=bool)


EVEN = even_mask()


def ligand_centers(att: np.ndarray, rep: np.ndarray) -> np.ndarray:
    return (att[:, 0] - rep[:, 0]) / VOL


def interp_l(l_cent: np.ndarray, x: float) -> float:
    xc = (np.arange(NX) + 0.5) * DX
    return float(np.interp(x, xc, l_cent, left=l_cent[0], right=l_cent[-1]))


def median_delta_l(att: np.ndarray, rep: np.ndarray, ell: float = ELL) -> float:
    l_cent = ligand_centers(att, rep)
    xs = np.arange(0.0, LX - ell + 1e-12, X_SAMPLE)
    d = np.empty(xs.size, dtype=np.float64)
    for n, x in enumerate(xs):
        d[n] = abs(interp_l(l_cent, x + ell) - interp_l(l_cent, x))
    return float(np.median(d))


def mean_abs_l(att: np.ndarray, rep: np.ndarray) -> float:
    return float(np.mean(np.abs(att / VOL - rep / VOL)))


def kappa_lambda(att: np.ndarray, rep: np.ndarray) -> float:
    s_att = 0.0
    s_rep = 0.0
    for i in range(NX):
        a = float(att[i, :].sum())
        r = float(rep[i, :].sum())
        if EVEN[i]:
            s_att += a
            s_rep += r
        else:
            s_att += r
            s_rep += a
    return float((s_att - s_rep) / (s_att + s_rep + EPS))


def diffuse_2d(q: np.ndarray, h: float) -> np.ndarray:
    fx = D * (q[:-1, :] - q[1:, :]) / (DX * DX)
    fy = D * (q[:, :-1] - q[:, 1:]) / (DY * DY)
    out = q.copy()
    out[:-1, :] -= h * fx
    out[1:, :] += h * fx
    out[:, :-1] -= h * fy
    out[:, 1:] += h * fy
    return out


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


def overlap(t_start: float, dt: float, t_on: float, t_off: float) -> float:
    lo = max(t_start, t_on)
    hi = min(t_start + dt, t_off)
    return max(0.0, hi - lo)


def deposit(att: np.ndarray, rep: np.ndarray, j_slab: float, t_start: float, dt: float) -> None:
    live = overlap(t_start, dt, 0.0, T_LIVE)
    if live <= 0.0 or j_slab == 0.0:
        return
    dm = j_slab * VOL * live
    att[EVEN, :] += dm
    rep[~EVEN, :] += dm


def integrate(j_slab: float, t_end: float) -> tuple[list[float], list[float], list[float], list[float]]:
    att = np.zeros((NX, NY), dtype=np.float64)
    rep = np.zeros((NX, NY), dtype=np.float64)
    times: list[float] = []
    kappas: list[float] = []
    deltas: list[float] = []
    abs_l: list[float] = []
    n_steps = int(round(t_end / DT))
    sample_stride = int(round(SAMPLE_DT / DT))
    times.append(0.0)
    kappas.append(kappa_lambda(att, rep))
    deltas.append(median_delta_l(att, rep))
    abs_l.append(mean_abs_l(att, rep))
    half = 0.5 * DT
    for step in range(n_steps):
        t_start = step * DT
        att = advance(att, half)
        rep = advance(rep, half)
        deposit(att, rep, j_slab, t_start, DT)
        att = advance(att, half)
        rep = advance(rep, half)
        t_end_step = (step + 1) * DT
        if (step + 1) % sample_stride == 0:
            times.append(t_end_step)
            kappas.append(kappa_lambda(att, rep))
            deltas.append(median_delta_l(att, rep))
            abs_l.append(mean_abs_l(att, rep))
    return times, kappas, deltas, abs_l


def main() -> int:
    refuse = " ".join(sys.argv).lower()
    if any(k in refuse for k in ("narma", "charc", "ipc")):
        print("scout refuses NARMA/CHARC/IPC")
        return 2
    if K != 0.0 or D != 800.0:
        print("scout identity is D=800 k=0. Stop.")
        return 2

    times_u, _, deltas_u, abs_u = integrate(1.0, T_LIVE)
    t_u = np.asarray(times_u)
    i25 = int(np.argmin(np.abs(t_u - T_LIVE)))
    d25_unit = float(np.asarray(deltas_u)[i25])
    if d25_unit <= 0.0 or not np.isfinite(d25_unit):
        print(f"SCOUT=FAIL unit median DeltaL(t=25)={d25_unit}. Name the blocker. Do not drop D.")
        return 1
    j_slab = DELTA_L_TARGET / d25_unit
    mean_l_unit = float(np.asarray(abs_u)[i25])
    mean_l_pred = mean_l_unit * j_slab
    if mean_l_pred > MEAN_L_CAP:
        print(
            f"SCOUT=FAIL mean_|L|(25)={mean_l_pred:.6g} > {MEAN_L_CAP:g} "
            f"to hit DeltaL>=0.3 (J_slab={j_slab:.6g}). Do not drop D. Do not start Java."
        )
        return 1

    times, kappas, deltas, abs_l = integrate(j_slab, T_HOLD_CAP)
    t_arr = np.asarray(times)
    k_arr = np.asarray(kappas)
    d_arr = np.asarray(deltas)
    l_arr = np.asarray(abs_l)
    i25 = int(np.argmin(np.abs(t_arr - T_LIVE)))
    delta_t25 = float(d_arr[i25])
    mean_l_25 = float(l_arr[i25])
    k_at_off = float(k_arr[i25])

    if delta_t25 < DELTA_L_TARGET - 1e-9:
        bump = DELTA_L_TARGET / delta_t25
        j_slab *= bump
        mean_l_25 *= bump
        if mean_l_25 > MEAN_L_CAP:
            print(
                f"SCOUT=FAIL mean_|L|(25)={mean_l_25:.6g} > {MEAN_L_CAP:g} after bump. "
                "Do not drop D. Do not start Java."
            )
            return 1
        times, kappas, deltas, abs_l = integrate(j_slab, T_HOLD_CAP)
        t_arr = np.asarray(times)
        k_arr = np.asarray(kappas)
        d_arr = np.asarray(deltas)
        l_arr = np.asarray(abs_l)
        i25 = int(np.argmin(np.abs(t_arr - T_LIVE)))
        delta_t25 = float(d_arr[i25])
        mean_l_25 = float(l_arr[i25])
        k_at_off = float(k_arr[i25])

    t_hold = None
    k_at_hold = None
    for t, k in zip(t_arr, k_arr):
        if t <= T_LIVE + 1e-12:
            continue
        if abs(k) < KAPPA_HOLD:
            t_hold = float(t)
            k_at_hold = float(k)
            break

    k_at_cap = float(k_arr[-1])
    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / "b3_field_scout.csv"
    with csv_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("t,kappa_lambda_field,median_abs_dL_ell,mean_abs_L\n")
        for t, k, ell, ellabs in zip(times, kappas, deltas, abs_l):
            handle.write(f"{t:.6g},{k:.8g},{ell:.8g},{ellabs:.8g}\n")

    grad_ok = delta_t25 >= DELTA_L_TARGET - 1e-9 and mean_l_25 <= MEAN_L_CAP
    hold_ok = t_hold is not None
    scout_ok = grad_ok and hold_ok

    summary = {
        "object": "LANE_B_CHEMOTACTIC_SPATIAL",
        "gate": "B3_SUSTAINED_STRIPE",
        "device": "LANE_B_MM_DISH_FLOW0",
        "scout": "field_only",
        "phase_even_slabs": "ATTRACTANT",
        "lambda_um": 200.0,
        "slab_um": SLAB,
        "ell_run_um": ELL,
        "D": D,
        "k": K,
        "T_live_s": T_LIVE,
        "J_slab": j_slab,
        "delta_L_ell_t25": delta_t25,
        "delta_L_unit_t25": d25_unit,
        "mean_abs_L_t25": mean_l_25,
        "t_occ_s": T_LIVE,
        "T_hold_s": t_hold,
        "kappa_field_at_t_occ": k_at_off,
        "kappa_field_at_T_hold": k_at_hold,
        "kappa_field_at_cap": k_at_cap,
        "observe_cap_s": T_HOLD_CAP,
        "scout_pass": scout_ok,
        "hydrolysis_added": False,
        "D_dropped": False,
        "narma": False,
        "cells": False,
    }
    (RESULTS / "b3_field_scout.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    proto = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    proto["J_slab"] = j_slab
    proto["delta_L_ell_t25"] = delta_t25
    proto["mean_abs_L_t25"] = mean_l_25
    proto["t_occ_s"] = T_LIVE
    if scout_ok:
        proto["T_hold_s"] = t_hold
        proto["T_hold_status"] = "SCOUT_SET_BEFORE_MOTILE_TRACES"
        proto["kappa_field_scout_at_t_occ"] = k_at_off
        proto["kappa_field_scout_at_T_hold"] = k_at_hold
    else:
        proto["T_hold_s"] = t_hold
        proto["T_hold_status"] = "SCOUT_FAIL"
        proto["kappa_field_scout_at_cap"] = k_at_cap
        proto["kappa_field_scout_at_t_occ"] = k_at_off
    PROTOCOL_JSON.write_text(json.dumps(proto, indent=2) + "\n", encoding="utf-8")

    print(
        f"B3 field scout J_slab={j_slab:.6g} DeltaL(ell,t=25)={delta_t25:.6g} "
        f"mean_|L|(25)={mean_l_25:.6g} kappa_field(25)={k_at_off:.6g}"
    )
    if scout_ok:
        print(
            f"T_hold={t_hold:.6g} s  kappa_field(T_hold)={k_at_hold:.6g}  "
            f"SCOUT=PASS  (wrote protocol_b3.json)"
        )
        return 0
    if not grad_ok:
        print(
            f"SCOUT=FAIL gradient. DeltaL(t=25)={delta_t25:.6g} mean_|L|={mean_l_25:.6g}. "
            "Do not drop D. Do not start Java."
        )
        return 1
    print(
        f"T_hold unset. |kappa_lambda| never < {KAPPA_HOLD} by {T_HOLD_CAP:g} s "
        f"(kappa_field(cap)={k_at_cap:.6g}). SCOUT=FAIL. Stop. "
        "Do not add hydrolysis. Do not drop D. Do not start Java."
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
