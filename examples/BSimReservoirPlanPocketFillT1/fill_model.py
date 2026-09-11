#!/usr/bin/env python3
"""PocketFill-T1 reduced OPEN_BUS_3 field + optional mean-field Danino ODE.

No NARMA. D1g QS_* copied; cells write the field (D1g did not).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix, lil_matrix

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
MU = math.log(2.0) / 1800.0
K_DANINO = 2.76e-3 / 60.0

PRIMARY_RTOL, PRIMARY_ATOL = 1e-6, 1e-9


def load_json(name: str) -> dict:
    with (CONFIGS / name).open(encoding="utf-8") as handle:
        return json.load(handle)


DISH = load_json("dish.json")
ARMS = {row["id"]: row for row in load_json("arms.json")["arms"]}

D0 = float(DISH["D_um2_s"])
K_HYBRID = float(DISH["k_hybriddish_s"])
K_HILL = float(DISH["K_hybriddish_uM"])
N_HILL = float(DISH["n_hill"])
TAU_R = float(DISH["tau_R_s"])
CONV = float(DISH["molecules_per_um3_per_uM"])
J_MAX = float(DISH["J_max_molecules_s"])
DX = float(DISH["dx_um"])
DY = float(DISH["dy_um"])
DEPTH = float(DISH["pocket_um"][2])
VOXEL = DX * DY * DEPTH
POCKET_LX = float(DISH["pocket_um"][0])
POCKET_LY = float(DISH["pocket_um"][1])
N_PACK = int(DISH["N_pack"])
CELL_VOL = float(DISH["cell_vol_um3"])


def hill(c: np.ndarray | float, k: float = K_HILL, n: float = N_HILL) -> np.ndarray | float:
    c = np.maximum(c, 0.0)
    kn = k ** n
    cn = np.power(c, n)
    return cn / (kn + cn)


def occupancy_flag(mean_h: float) -> str:
    if not math.isfinite(mean_h):
        return "NA"
    if mean_h > 0.95:
        return "SATURATED"
    if mean_h >= 0.05:
        return "ALIVE"
    return "DEAD"


def danino_rhs(y: np.ndarray, ahl_ext: float) -> np.ndarray:
    y = np.maximum(y, 0.0)
    extra = y[1] - ahl_ext
    hill_la = (y[3] ** QS_N) / (QS_KMLA ** QS_N + y[3] ** QS_N + 1e-30)
    dy = np.zeros(4)
    dy[0] = (
        QS_A0LI
        + QS_KPLI * hill_la
        - (QS_DELTA1 * y[0]) / (QS_G * (y[0] + y[2]) + 1.0)
        - MU * y[0]
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
        - MU * y[2]
    )
    dy[3] = QS_KR1ON * (QS_LTOT - y[3]) * y[1] - QS_KR1OFF * y[3] - QS_T_LA * y[3] - MU * y[3]
    return dy


@dataclass
class Geometry:
    nx: int
    ny: int
    active: np.ndarray
    pocket: np.ndarray
    bus: np.ndarray
    idx: np.ndarray
    n: int
    pocket_idx: np.ndarray
    M: csr_matrix
    source_b: np.ndarray
    k_decay: float


def _empty_grid(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.zeros((nx, ny), dtype=bool),
        np.zeros((nx, ny), dtype=bool),
        np.zeros((nx, ny), dtype=bool),
    )


def _index_active(active: np.ndarray) -> tuple[np.ndarray, int]:
    idx = np.full(active.shape, -1, dtype=int)
    n = int(active.sum())
    idx[active] = np.arange(n)
    return idx, n


def _d_map(pocket: np.ndarray, factor: float) -> np.ndarray:
    d = np.full(pocket.shape, D0, dtype=float)
    d[pocket] = D0 * factor
    return d


def _assemble_diffusion(M: lil_matrix, geom_idx: np.ndarray, active: np.ndarray, dmap: np.ndarray) -> None:
    nx, ny = active.shape
    for i in range(nx - 1):
        for j in range(ny):
            if not (active[i, j] and active[i + 1, j]):
                continue
            di, dj = dmap[i, j], dmap[i + 1, j]
            d_face = 0.0 if (di + dj) == 0 else 2.0 * di * dj / (di + dj)
            coef = d_face / (DX * DX)
            p, q = int(geom_idx[i, j]), int(geom_idx[i + 1, j])
            M[p, q] += coef
            M[p, p] -= coef
            M[q, p] += coef
            M[q, q] -= coef
    for i in range(nx):
        for j in range(ny - 1):
            if not (active[i, j] and active[i, j + 1]):
                continue
            di, dj = dmap[i, j], dmap[i, j + 1]
            d_face = 0.0 if (di + dj) == 0 else 2.0 * di * dj / (di + dj)
            coef = d_face / (DY * DY)
            p, q = int(geom_idx[i, j]), int(geom_idx[i, j + 1])
            M[p, q] += coef
            M[p, p] -= coef
            M[q, p] += coef
            M[q, q] -= coef


def _add_decay(M: lil_matrix, n: int, k: float) -> None:
    for p in range(n):
        M[p, p] -= k


def _add_upwind_x(M: lil_matrix, geom_idx: np.ndarray, mask: np.ndarray, velocity: float) -> None:
    if velocity == 0:
        return
    nx, ny = mask.shape
    coef = velocity / DX
    for i in range(nx):
        for j in range(ny):
            if not mask[i, j]:
                continue
            p = int(geom_idx[i, j])
            M[p, p] -= coef
            if i > 0 and mask[i - 1, j]:
                M[p, int(geom_idx[i - 1, j])] += coef


def make_bus_geometry(k_decay: float, d_pocket_factor: float) -> Geometry:
    nx = int(round(float(DISH["bus_length_um"]) / DX))
    ny = int(round((POCKET_LY + float(DISH["bus_width_um"])) / DY))
    active, pocket, bus = _empty_grid(nx, ny)
    x0p = float(DISH["pocket_x0_um"])
    y0p = float(DISH["pocket_y0_um"])
    x0b = float(DISH["bus_x0_um"])
    y0b = float(DISH["bus_y0_um"])
    for i in range(nx):
        x = (i + 0.5) * DX
        for j in range(ny):
            y = (j + 0.5) * DY
            in_p = (x0p <= x < x0p + POCKET_LX) and (y0p <= y < y0p + POCKET_LY)
            in_b = (x0b <= x < x0b + float(DISH["bus_length_um"])) and (
                y0b <= y < y0b + float(DISH["bus_width_um"])
            )
            pocket[i, j] = in_p
            bus[i, j] = in_b
            active[i, j] = in_p or in_b
    idx, n = _index_active(active)
    M = lil_matrix((n, n), dtype=float)
    _assemble_diffusion(M, idx, active, _d_map(pocket, d_pocket_factor))
    _add_decay(M, n, k_decay)
    _add_upwind_x(M, idx, bus, 3.0)
    sx = int((x0p + float(DISH["source_pocket_um"][0])) / DX)
    sy = int((y0p + float(DISH["source_pocket_um"][1])) / DY)
    source_b = np.zeros(n)
    source_b[int(idx[sx, sy])] = J_MAX / (CONV * VOXEL)
    pocket_idx = idx[pocket]
    return Geometry(
        nx=nx,
        ny=ny,
        active=active,
        pocket=pocket,
        bus=bus,
        idx=idx,
        n=n,
        pocket_idx=pocket_idx[pocket_idx >= 0],
        M=M.tocsr(),
        source_b=source_b,
        k_decay=k_decay,
    )


def unpack(geom: Geometry, c_active: np.ndarray) -> np.ndarray:
    field = np.full((geom.nx, geom.ny), np.nan)
    field[geom.active] = c_active
    return field


def pocket_mean(geom: Geometry, c_active: np.ndarray) -> float:
    return float(np.mean(c_active[geom.pocket_idx]))


def integrate_point(geom: Geometry, u: float, duration: float, sample_dt: float = 20.0) -> dict:
    def fun(_t: float, state: np.ndarray) -> np.ndarray:
        return geom.M @ state + u * geom.source_b

    t_eval = np.arange(0.0, duration + 0.5 * sample_dt, sample_dt)
    t_eval[-1] = duration
    sol = solve_ivp(
        fun,
        (0.0, duration),
        np.zeros(geom.n),
        method="BDF",
        jac=lambda _t, _y: geom.M,
        rtol=PRIMARY_RTOL,
        atol=PRIMARY_ATOL,
        t_eval=t_eval,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    c = sol.y.T
    r = np.zeros_like(c)
    for i in range(len(sol.t) - 1):
        dt = float(sol.t[i + 1] - sol.t[i])
        h = hill(c[i + 1])
        r[i + 1] = h + (r[i] - h) * math.exp(-dt / TAU_R)
    return {"t": sol.t, "C": c, "R": r}


def integrate_volumetric(geom: Geometry, duration: float, sample_dt: float = 60.0) -> dict:
    n = geom.n
    n_p = len(geom.pocket_idx)
    share = 1.0 / (n_p * CONV * VOXEL)

    def fun(_t: float, state: np.ndarray) -> np.ndarray:
        c = state[:n]
        y = np.maximum(state[n:], 0.0)
        c_mean = float(np.mean(c[geom.pocket_idx]))
        dy = danino_rhs(y, c_mean)
        flux_mol = N_PACK * CELL_VOL * CONV * CELL_WALL_DIFF * (y[1] - c_mean)
        dc = geom.M @ c
        dc[geom.pocket_idx] += flux_mol * share
        return np.concatenate([dc, dy])

    t_eval = np.arange(0.0, duration + 0.5 * sample_dt, sample_dt)
    t_eval[-1] = duration
    z0 = np.zeros(n + 4)
    sol = solve_ivp(
        fun,
        (0.0, duration),
        z0,
        method="BDF",
        rtol=PRIMARY_RTOL,
        atol=PRIMARY_ATOL,
        t_eval=t_eval,
        max_step=30.0,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    c = sol.y[:n].T
    y = sol.y[n:].T
    return {"t": sol.t, "C": c, "Y": y}
