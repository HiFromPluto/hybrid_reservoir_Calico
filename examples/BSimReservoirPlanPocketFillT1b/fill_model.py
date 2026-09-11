#!/usr/bin/env python3
"""PocketFill-T1b reduced OPEN_BUS_3 field + mean-field Danino ODE.

Clone of T1 reduced field. QS_* and geometry copied verbatim.
No NARMA. Cells write the field. Seed arms are predeclared.
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
ARMS_DOC = load_json("arms.json")
ARMS = {row["id"]: row for row in ARMS_DOC["arms"]}
BATH_UM = float(ARMS_DOC["BATH_UM"])

D0 = float(DISH["D_um2_s"])
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


def hill(c: np.ndarray | float, k: float = QS_KMLA, n: float = QS_N) -> np.ndarray | float:
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


def mass_molecules(geom: Geometry, c_active: np.ndarray) -> float:
    return float(np.sum(c_active) * VOXEL * CONV)


def integrate_volumetric(
    geom: Geometry,
    duration: float,
    sample_dt: float = 60.0,
    seed: str = "none",
    bath_um: float = BATH_UM,
) -> dict:
    n = geom.n
    n_p = len(geom.pocket_idx)
    share = 1.0 / (n_p * CONV * VOXEL)
    hold_bath = seed == "bath"

    z0 = np.zeros(n + 4)
    if seed in ("pulse_pocket", "bath"):
        z0[geom.pocket_idx] = bath_um

    def fun(_t: float, state: np.ndarray) -> np.ndarray:
        c = state[:n].copy()
        y = np.maximum(state[n:], 0.0)
        if hold_bath:
            c[geom.pocket_idx] = bath_um
            c_mean = bath_um
        else:
            c_mean = float(np.mean(c[geom.pocket_idx]))
        dy = danino_rhs(y, c_mean)
        flux_mol = N_PACK * CELL_VOL * CONV * CELL_WALL_DIFF * (y[1] - c_mean)
        dc = geom.M @ c
        dc[geom.pocket_idx] += flux_mol * share
        if hold_bath:
            dc[geom.pocket_idx] = 0.0
        return np.concatenate([dc, dy])

    t_eval = np.arange(0.0, duration + 0.5 * sample_dt, sample_dt)
    t_eval[-1] = duration
    sol = solve_ivp(
        fun,
        (0.0, duration),
        z0,
        method="BDF",
        rtol=PRIMARY_RTOL,
        atol=PRIMARY_ATOL,
        t_eval=t_eval,
        max_step=30.0,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    c = sol.y[:n].T
    y = sol.y[n:].T
    if hold_bath:
        c[:, geom.pocket_idx] = bath_um
    return {"t": sol.t, "C": c, "Y": y, "sol": sol, "seed": seed}


def field_budget(geom: Geometry, sol, seed: str, duration: float, dt: float = 0.25) -> dict:
    n = geom.n
    n_p = len(geom.pocket_idx)
    share = 1.0 / (n_p * CONV * VOXEL)
    hold_bath = seed == "bath"
    t = np.arange(0.0, duration + 0.5 * dt, dt)
    t[-1] = duration
    z = sol.sol(t)
    c = z[:n].T
    y = z[n:].T
    if hold_bath:
        c[:, geom.pocket_idx] = BATH_UM
    mass = np.sum(c, axis=1) * VOXEL * CONV
    c_mean = np.mean(c[:, geom.pocket_idx], axis=1)
    flux_mol = N_PACK * CELL_VOL * CONV * CELL_WALL_DIFF * (y[:, 1] - c_mean)
    decay_rate = geom.k_decay * mass
    dc_m = (geom.M @ c.T).T
    leak_rate = -(np.sum(dc_m, axis=1) * VOXEL * CONV + decay_rate)
    if hold_bath:
        dc_p = dc_m[:, geom.pocket_idx] + flux_mol[:, None] * share
        bath_rate = -np.sum(dc_p, axis=1) * VOXEL * CONV
    else:
        bath_rate = np.zeros(len(t))
    membrane = float(np.trapezoid(flux_mol, t))
    decay = float(np.trapezoid(decay_rate, t))
    leak = float(np.trapezoid(leak_rate, t))
    bath = float(np.trapezoid(bath_rate, t))
    initial = float(mass[0])
    remaining = float(mass[-1])
    residual = initial + membrane + bath - remaining - decay - leak
    dominant = max(
        abs(initial),
        abs(membrane),
        abs(bath),
        abs(remaining),
        abs(decay),
        abs(leak),
        1.0,
    )
    frac = abs(residual) / dominant
    c_nn = bool(np.min(c) >= -1e-12)
    mass_ok = frac <= 0.01 and c_nn
    return {
        "initial_molecules": initial,
        "membrane_molecules": membrane,
        "bath_molecules": bath,
        "remaining_molecules": remaining,
        "decay_molecules": decay,
        "leak_molecules": leak,
        "residual_molecules": residual,
        "residual_frac": frac,
        "C_nonnegative": c_nn,
        "mass_ok": mass_ok,
        "mass_class": "LE1PCT" if mass_ok else "FAIL",
    }
