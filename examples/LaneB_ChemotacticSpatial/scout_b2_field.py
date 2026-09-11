#!/usr/bin/env python3
"""B2 field-only scout. No cells. Square-wave IC, then J=0.

Solves L_slab so median |L(x+ell)-L(x)| at t=1 s is >= 0.3 (ell=17.2 um).
Sets t_occ and T_hold from stripe kappa_lambda(t). Cap L_slab=100.
Does not drop D. Does not add hydrolysis. Not B0. Not B1. Not NARMA.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL_JSON = HERE / "configs" / "protocol_b2.json"

LX, LY, LZ = 1000.0, 500.0, 10.0
NX, NY, NZ = 100, 50, 1
DX, DY, DZ = LX / NX, LY / NY, LZ / NZ
D = 800.0
K = 0.0
DT = 0.02
SAMPLE_DT = 0.10
SLAB = 100.0
ELL = 17.2
DELTA_L_TARGET = 0.3
L_SLAB_CAP = 100.0
T_GRAD = 1.0
T_OCC_AFTER = 1.0
KAPPA_OCC = 0.40
T_OCC_FALLBACK = 2.0
T_OCC_MAX = 15.0
KAPPA_HOLD = 0.15
T_HOLD_CAP = 40.0
EPS = 1e-12
X_SAMPLE = 1.0


def even_att_slab(x: float) -> bool:
    return (np.floor(x / SLAB) % 2.0) == 0.0


def square_wave(amp: float) -> tuple[np.ndarray, np.ndarray]:
    att = np.zeros((NX, NY), dtype=np.float64)
    rep = np.zeros((NX, NY), dtype=np.float64)
    vol = DX * DY * DZ
    for i in range(NX):
        x_c = (i + 0.5) * DX
        if even_att_slab(x_c):
            att[i, :] = amp * vol
        else:
            rep[i, :] = amp * vol
    return att, rep


def ligand_centers(att: np.ndarray, rep: np.ndarray) -> np.ndarray:
    vol = DX * DY * DZ
    # Uniform in y: use column mean (identical columns after IC + no-flux y).
    return (att[:, 0] - rep[:, 0]) / vol


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


def kappa_lambda(att: np.ndarray, rep: np.ndarray) -> float:
    s_att = 0.0
    s_rep = 0.0
    for i in range(NX):
        x_c = (i + 0.5) * DX
        a = float(att[i, :].sum())
        r = float(rep[i, :].sum())
        if even_att_slab(x_c):
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


def integrate(amp: float, t_end: float, record: bool) -> tuple[np.ndarray, np.ndarray, list[float], list[float], list[float]]:
    att, rep = square_wave(amp)
    times: list[float] = []
    kappas: list[float] = []
    deltas: list[float] = []
    n_steps = int(round(t_end / DT))
    sample_stride = int(round(SAMPLE_DT / DT))
    if record:
        times.append(0.0)
        kappas.append(kappa_lambda(att, rep))
        deltas.append(median_delta_l(att, rep))
    half = 0.5 * DT
    for step in range(n_steps):
        att = advance(att, half)
        rep = advance(rep, half)
        att = advance(att, half)
        rep = advance(rep, half)
        t_end_step = (step + 1) * DT
        if record and (step + 1) % sample_stride == 0:
            times.append(t_end_step)
            kappas.append(kappa_lambda(att, rep))
            deltas.append(median_delta_l(att, rep))
    return att, rep, times, kappas, deltas


def main() -> int:
    refuse = " ".join(sys.argv).lower()
    if any(k in refuse for k in ("narma", "charc", "ipc")):
        print("scout refuses NARMA/CHARC/IPC")
        return 2
    if K != 0.0 or D != 800.0:
        print("scout identity is D=800 k=0. Stop.")
        return 2

    att1, rep1, _, _, _ = integrate(1.0, T_GRAD, record=False)
    d1_unit = median_delta_l(att1, rep1)
    if d1_unit <= 0.0 or not np.isfinite(d1_unit):
        print(f"SCOUT=FAIL unit median DeltaL(t=1)={d1_unit}. D=800 is the blocker.")
        return 1
    l_slab = DELTA_L_TARGET / d1_unit
    if l_slab > L_SLAB_CAP:
        print(
            f"SCOUT=FAIL L_slab={l_slab:.6g} > {L_SLAB_CAP:g} to hit DeltaL>=0.3 "
            f"(unit DeltaL={d1_unit:.6g}). D=800 is the blocker. Do not drop D. "
            "Do not start Java."
        )
        return 1

    att, rep, times, kappas, deltas = integrate(l_slab, T_HOLD_CAP, record=True)
    t_arr = np.asarray(times)
    k_arr = np.asarray(kappas)
    d_arr = np.asarray(deltas)

    i1 = int(np.argmin(np.abs(t_arr - T_GRAD)))
    delta_t1 = float(d_arr[i1])
    if delta_t1 < DELTA_L_TARGET:
        # Linear solve should hit exactly; floating leftover is a fail-safe bump.
        bump = DELTA_L_TARGET / delta_t1
        l_slab *= bump
        if l_slab > L_SLAB_CAP:
            print(
                f"SCOUT=FAIL L_slab={l_slab:.6g} > {L_SLAB_CAP:g} after bump. "
                "D=800 is the blocker. Do not start Java."
            )
            return 1
        att, rep, times, kappas, deltas = integrate(l_slab, T_HOLD_CAP, record=True)
        t_arr = np.asarray(times)
        k_arr = np.asarray(kappas)
        d_arr = np.asarray(deltas)
        i1 = int(np.argmin(np.abs(t_arr - T_GRAD)))
        delta_t1 = float(d_arr[i1])

    t_occ = None
    k_at_occ = None
    for t, k in zip(t_arr, k_arr):
        if t <= T_OCC_AFTER + 1e-12:
            continue
        if abs(k) >= KAPPA_OCC:
            t_occ = float(t)
            k_at_occ = float(k)
            break
    if t_occ is None:
        t_occ = T_OCC_FALLBACK
        k_at_occ = float(k_arr[np.argmin(np.abs(t_arr - t_occ))])

    t_hold = None
    k_at_hold = None
    for t, k in zip(t_arr, k_arr):
        if t <= t_occ + 1e-12:
            continue
        if abs(k) < KAPPA_HOLD:
            t_hold = float(t)
            k_at_hold = float(k)
            break

    k_at_cap = float(k_arr[-1])
    d_at_occ = float(d_arr[np.argmin(np.abs(t_arr - t_occ))])

    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / "b2_field_scout.csv"
    with csv_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("t,kappa_lambda_field,median_abs_dL_ell\n")
        for t, k, ell in zip(times, kappas, deltas):
            handle.write(f"{t:.6g},{k:.8g},{ell:.8g}\n")

    occ_ok = t_occ is not None and t_occ <= T_OCC_MAX + 1e-12
    hold_ok = t_hold is not None
    grad_ok = delta_t1 >= DELTA_L_TARGET - 1e-9 and l_slab <= L_SLAB_CAP
    scout_ok = occ_ok and hold_ok and grad_ok

    summary = {
        "object": "LANE_B_CHEMOTACTIC_SPATIAL",
        "gate": "B2_STRIPE_OCCUPY",
        "device": "LANE_B_MM_DISH_FLOW0",
        "scout": "field_only",
        "phase_even_slabs": "ATTRACTANT",
        "lambda_um": 200.0,
        "slab_um": SLAB,
        "ell_run_um": ELL,
        "D": D,
        "k": K,
        "L_slab": l_slab,
        "delta_L_ell_t1": delta_t1,
        "delta_L_unit_t1": d1_unit,
        "t_occ_s": t_occ,
        "T_hold_s": t_hold,
        "kappa_field_at_t_occ": k_at_occ,
        "kappa_field_at_T_hold": k_at_hold,
        "kappa_field_at_cap": k_at_cap,
        "delta_L_ell_at_t_occ": d_at_occ,
        "observe_cap_s": T_HOLD_CAP,
        "scout_pass": scout_ok,
        "hydrolysis_added": False,
        "D_dropped": False,
        "narma": False,
        "cells": False,
    }
    (RESULTS / "b2_field_scout.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    proto = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    proto["L_slab"] = l_slab
    proto["delta_L_ell_t1"] = delta_t1
    proto["t_occ_s"] = t_occ
    if scout_ok:
        proto["T_hold_s"] = t_hold
        proto["T_hold_status"] = "SCOUT_SET_BEFORE_MOTILE_TRACES"
        proto["kappa_field_scout_at_t_occ"] = k_at_occ
        proto["kappa_field_scout_at_T_hold"] = k_at_hold
    else:
        proto["T_hold_s"] = t_hold
        proto["T_hold_status"] = "SCOUT_FAIL"
        proto["kappa_field_scout_at_cap"] = k_at_cap
        proto["kappa_field_scout_at_t_occ"] = k_at_occ
    PROTOCOL_JSON.write_text(json.dumps(proto, indent=2) + "\n", encoding="utf-8")

    print(
        f"B2 field scout L_slab={l_slab:.6g} DeltaL(ell,t=1)={delta_t1:.6g} "
        f"(unit={d1_unit:.6g}) t_occ={t_occ:.6g}s kappa_field(t_occ)={k_at_occ:.6g}"
    )
    if scout_ok:
        print(
            f"T_hold={t_hold:.6g} s  kappa_field(T_hold)={k_at_hold:.6g}  "
            f"SCOUT=PASS  (wrote protocol_b2.json)"
        )
        return 0
    if not grad_ok:
        print(
            f"SCOUT=FAIL gradient. DeltaL(t=1)={delta_t1:.6g} L_slab={l_slab:.6g}. "
            "D=800 is the blocker. Do not drop D. Do not start Java."
        )
        return 1
    if not occ_ok:
        print(
            f"SCOUT=FAIL t_occ={t_occ} > {T_OCC_MAX:g} s. Stop. Do not start Java."
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
