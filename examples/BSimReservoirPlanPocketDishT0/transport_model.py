#!/usr/bin/env python3
"""Reduced 2-D PocketDish-A transport + Hill/L surrogate.

E0.2-style cell-centred finite volume. No cells, no NARMA, no ridge.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix, lil_matrix


HERE = Path(__file__).resolve().parent
CONFIGS = HERE / "configs"

PRIMARY_RTOL, PRIMARY_ATOL = 1e-7, 1e-10
TIGHT_RTOL, TIGHT_ATOL = 1e-9, 1e-12


def load_json(name: str) -> dict:
    with (CONFIGS / name).open(encoding="utf-8") as handle:
        return json.load(handle)


DISH = load_json("dish.json")
ARMS = {row["id"]: row for row in load_json("arms.json")["arms"]}
DRIVES = {row["id"]: row for row in load_json("drives.json")["drives"]}

D = float(DISH["D_um2_s"])
K_DECAY = float(DISH["k_s"])
K_HILL = float(DISH["K_uM"])
N_HILL = float(DISH["n_hill"])
TAU_R = float(DISH["tau_R_s"])
TAU_L = float(DISH["tau_L_s"])
CONV = float(DISH["molecules_per_um3_per_uM"])
J_MAX = float(DISH["J_max_molecules_s"])
DX = float(DISH["dx_um"])
DY = float(DISH["dy_um"])
DEPTH = float(DISH["pocket_um"][2])
VOXEL = DX * DY * DEPTH
POCKET_LX, POCKET_LY = float(DISH["pocket_um"][0]), float(DISH["pocket_um"][1])
POCKET_V = POCKET_LX * POCKET_LY * DEPTH


def hill(c: np.ndarray) -> np.ndarray:
    c = np.maximum(c, 0.0)
    kn = K_HILL ** N_HILL
    cn = np.power(c, N_HILL)
    return cn / (kn + cn)


def occupancy_flag(mean_r: float) -> str:
    if not math.isfinite(mean_r):
        return "NA"
    if mean_r > 0.95:
        return "SATURATED"
    if mean_r >= 0.05:
        return "ALIVE"
    return "DEAD"


def closed_form_css(u: float) -> float:
    return J_MAX * u / (K_DECAY * POCKET_V * CONV)


def transport_scales(arm_id: str) -> dict[str, float]:
    length = POCKET_LX
    l2d = length * length / D
    da = K_DECAY * l2d
    spec = ARMS[arm_id]
    if spec["id"] == "THROUGH_8":
        v = float(spec["pocket_flow_um_s"])
        pe = v * length / D
        transit = length / v
    elif spec["id"] == "OPEN_BUS_3":
        v = float(spec["bus_flow_um_s"])
        pe = v * float(DISH["bus_width_um"]) / D
        transit = float(DISH["bus_length_um"]) / v
    else:
        pe = 0.0
        transit = float("inf")
    return {
        "L_um": length,
        "L2_over_D_s": l2d,
        "Damkohler": da,
        "Pe": pe,
        "Transit_s": transit,
        "tau_AHL_s": 1.0 / K_DECAY,
    }


@dataclass
class Geometry:
    arm_id: str
    nx: int
    ny: int
    active: np.ndarray
    pocket: np.ndarray
    bus: np.ndarray
    idx: np.ndarray
    n: int
    M: csr_matrix
    source_b: np.ndarray
    source_ij: tuple[int, int]
    pocket_flow: float
    bus_flow: float
    origin_x: float
    origin_y: float
    interface_pairs: list[tuple[int, int, int, int]] = field(default_factory=list)

    @property
    def x_centers(self) -> np.ndarray:
        return self.origin_x + (np.arange(self.nx) + 0.5) * DX

    @property
    def y_centers(self) -> np.ndarray:
        return self.origin_y + (np.arange(self.ny) + 0.5) * DY


def _empty_grid(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    active = np.zeros((nx, ny), dtype=bool)
    pocket = np.zeros((nx, ny), dtype=bool)
    bus = np.zeros((nx, ny), dtype=bool)
    return active, pocket, bus


def _index_active(active: np.ndarray) -> tuple[np.ndarray, int]:
    idx = np.full(active.shape, -1, dtype=int)
    n = int(active.sum())
    idx[active] = np.arange(n)
    return idx, n


def _assemble_diffusion(M: lil_matrix, geom_idx: np.ndarray, active: np.ndarray) -> None:
    nx, ny = active.shape
    coef_x = D / (DX * DX)
    coef_y = D / (DY * DY)
    for i in range(nx - 1):
        for j in range(ny):
            if not (active[i, j] and active[i + 1, j]):
                continue
            p, q = int(geom_idx[i, j]), int(geom_idx[i + 1, j])
            M[p, q] += coef_x
            M[p, p] -= coef_x
            M[q, p] += coef_x
            M[q, q] -= coef_x
    for i in range(nx):
        for j in range(ny - 1):
            if not (active[i, j] and active[i, j + 1]):
                continue
            p, q = int(geom_idx[i, j]), int(geom_idx[i, j + 1])
            M[p, q] += coef_y
            M[p, p] -= coef_y
            M[q, p] += coef_y
            M[q, q] -= coef_y


def _add_decay(M: lil_matrix, n: int) -> None:
    for p in range(n):
        M[p, p] -= K_DECAY


def _add_absorbing_y1(M: lil_matrix, geom_idx: np.ndarray, active: np.ndarray) -> None:
    nx, ny = active.shape
    coef = 2.0 * D / (DY * DY)
    j = ny - 1
    for i in range(nx):
        if not active[i, j]:
            continue
        p = int(geom_idx[i, j])
        M[p, p] -= coef


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


def _add_upwind_y(M: lil_matrix, geom_idx: np.ndarray, mask: np.ndarray, velocity: float) -> None:
    if velocity == 0:
        return
    nx, ny = mask.shape
    coef = velocity / DY
    for i in range(nx):
        for j in range(ny):
            if not mask[i, j]:
                continue
            p = int(geom_idx[i, j])
            M[p, p] -= coef
            if j > 0 and mask[i, j - 1]:
                M[p, int(geom_idx[i, j - 1])] += coef


def _source_b(idx: np.ndarray, ij: tuple[int, int], n: int) -> np.ndarray:
    b = np.zeros(n, dtype=float)
    p = int(idx[ij])
    if p < 0:
        raise RuntimeError(f"source cell {ij} is inactive")
    b[p] = J_MAX / (CONV * VOXEL)
    return b


def _interface_pairs(pocket: np.ndarray, bus: np.ndarray) -> list[tuple[int, int, int, int]]:
    nx, ny = pocket.shape
    pairs = []
    for i in range(nx):
        for j in range(ny - 1):
            if pocket[i, j] and bus[i, j + 1]:
                pairs.append((i, j, i, j + 1))
            if bus[i, j] and pocket[i, j + 1]:
                pairs.append((i, j + 1, i, j))
    return pairs


def make_geometry(arm_id: str) -> Geometry:
    spec = ARMS[arm_id]
    if spec["bus"]:
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
        origin_x, origin_y = 0.0, 0.0
        sx = int((x0p + float(DISH["source_pocket_um"][0])) / DX)
        sy = int((y0p + float(DISH["source_pocket_um"][1])) / DY)
    else:
        nx = int(round(POCKET_LX / DX))
        ny = int(round(POCKET_LY / DY))
        active, pocket, bus = _empty_grid(nx, ny)
        active[:, :] = True
        pocket[:, :] = True
        origin_x, origin_y = 0.0, 0.0
        sx = int(float(DISH["source_pocket_um"][0]) / DX)
        sy = int(float(DISH["source_pocket_um"][1]) / DY)

    idx, n = _index_active(active)
    M = lil_matrix((n, n), dtype=float)
    _assemble_diffusion(M, idx, active)
    _add_decay(M, n)
    if spec["boundary"] == "ABSORBING_OPEN_Y":
        _add_absorbing_y1(M, idx, active)
    if spec["boundary"] == "OUTFLOW_Y":
        _add_upwind_y(M, idx, pocket, float(spec["pocket_flow_um_s"]))
    if spec["boundary"] == "BUS_COUPLED":
        _add_upwind_x(M, idx, bus, float(spec["bus_flow_um_s"]))
    M_csr = M.tocsr()
    source_ij = (sx, sy)
    if not active[sx, sy]:
        raise RuntimeError(f"{arm_id}: source {source_ij} inactive")
    return Geometry(
        arm_id=arm_id,
        nx=nx,
        ny=ny,
        active=active,
        pocket=pocket,
        bus=bus,
        idx=idx,
        n=n,
        M=M_csr,
        source_b=_source_b(idx, source_ij, n),
        source_ij=source_ij,
        pocket_flow=float(spec["pocket_flow_um_s"]),
        bus_flow=float(spec["bus_flow_um_s"]),
        origin_x=origin_x,
        origin_y=origin_y,
        interface_pairs=_interface_pairs(pocket, bus),
    )


def unpack(geom: Geometry, c_active: np.ndarray) -> np.ndarray:
    field = np.full((geom.nx, geom.ny), np.nan, dtype=float)
    field[geom.active] = c_active
    return field


def pocket_values(geom: Geometry, c_active: np.ndarray) -> np.ndarray:
    field = unpack(geom, c_active)
    return field[geom.pocket]


def mass_molecules(geom: Geometry, c_active: np.ndarray) -> float:
    return float(np.sum(c_active) * VOXEL * CONV)


def injected_molecules(u_on: float, t_on: float) -> float:
    return J_MAX * u_on * t_on


def boundary_rate(geom: Geometry, c_active: np.ndarray) -> float:
    dcdt = geom.M @ c_active
    decay = K_DECAY * float(np.sum(c_active)) * VOXEL * CONV
    return float(-np.sum(dcdt) * VOXEL * CONV - decay)


def interface_rate(geom: Geometry, c_active: np.ndarray) -> float:
    if not geom.interface_pairs:
        return 0.0
    field = unpack(geom, c_active)
    rate = 0.0
    area = DX * DEPTH
    for i_p, j_p, i_b, j_b in geom.interface_pairs:
        cp, cb = field[i_p, j_p], field[i_b, j_b]
        rate += D * (cp - cb) / DY * area * CONV
    return float(rate)


def _u_constant(value: float) -> Callable[[float], float]:
    return lambda _t: value


def integrate_concentration(
    geom: Geometry,
    u_fn: Callable[[float], float],
    duration: float,
    c0: np.ndarray | None = None,
    *,
    sample_dt: float = 20.0,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
    extra_times: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    if duration <= 0:
        c0 = np.zeros(geom.n) if c0 is None else c0.copy()
        t = np.array([0.0])
        return {"t": t, "C": c0.reshape(1, -1)}
    c0 = np.zeros(geom.n, dtype=float) if c0 is None else np.asarray(c0, dtype=float).copy()

    def fun(time: float, state: np.ndarray) -> np.ndarray:
        return geom.M @ state + u_fn(time) * geom.source_b

    times = [0.0]
    t = sample_dt
    while t < duration - 1e-12:
        times.append(t)
        t += sample_dt
    times.append(duration)
    if extra_times is not None:
        times.extend(float(x) for x in extra_times if 0.0 < x < duration)
    t_eval = np.unique(np.asarray(times, dtype=float))
    sol = solve_ivp(
        fun,
        (0.0, duration),
        c0,
        method="BDF",
        jac=lambda _t, _y: geom.M,
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"{geom.arm_id}: BDF failed: {sol.message}")
    return {"t": sol.t.copy(), "C": sol.y.T.copy(), "sol": sol}


def advance_reporter(c_samples: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_t, n = c_samples.shape
    r = np.zeros((n_t, n), dtype=float)
    lum = np.zeros((n_t, n), dtype=float)
    for i in range(n_t - 1):
        dt = float(t[i + 1] - t[i])
        h = hill(c_samples[i + 1])
        decay_r = math.exp(-dt / TAU_R)
        r[i + 1] = h + (r[i] - h) * decay_r
        r_mid = 0.5 * (r[i] + r[i + 1])
        decay_l = math.exp(-dt / TAU_L)
        lum[i + 1] = r_mid + (lum[i] - r_mid) * decay_l
    return r, lum


def _resample_segment(sol, duration: float, t0: float, dt: float) -> tuple[np.ndarray, np.ndarray]:
    if duration <= 0:
        return np.array([t0]), sol.sol(np.array([0.0])).T
    n = max(2, int(round(duration / dt)) + 1)
    local = np.linspace(0.0, duration, n)
    return local + t0, sol.sol(local).T


def simulate_drive(
    geom: Geometry,
    drive_id: str,
    *,
    duration_on: float | None = None,
    duration_off: float | None = None,
    sample_dt: float = 1.0,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
) -> dict:
    drive = DRIVES[drive_id]
    t_on = float(drive["t_on_s"] if duration_on is None else duration_on)
    t_off = float(drive["t_off_s"] if duration_off is None else duration_off)
    u_on = float(drive["u_hold"])
    first = integrate_concentration(
        geom, _u_constant(u_on), t_on, sample_dt=min(20.0, max(t_on, 1.0)), rtol=rtol, atol=atol
    )
    t, c = _resample_segment(first["sol"], t_on, 0.0, sample_dt)
    if t_off > 0:
        rest = integrate_concentration(
            geom,
            _u_constant(0.0),
            t_off,
            c0=first["sol"].sol(np.array([t_on])).ravel(),
            sample_dt=min(20.0, max(t_off, 1.0)),
            rtol=rtol,
            atol=atol,
        )
        t_rest, c_rest = _resample_segment(rest["sol"], t_off, t_on, sample_dt)
        t = np.concatenate([t, t_rest[1:]])
        c = np.vstack([c, c_rest[1:]])
    r, lum = advance_reporter(c, t)
    return {
        "t": t,
        "C": c,
        "R": r,
        "L": lum,
        "u_on": u_on,
        "t_on": t_on,
        "t_off": t_off,
        "injected": injected_molecules(u_on, t_on),
    }


def u_trace(t: np.ndarray, t_on: float, u_on: float) -> np.ndarray:
    return np.where(t <= t_on + 1e-12, u_on, 0.0)


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 2 or float(np.std(a)) == 0.0 or float(np.std(b)) == 0.0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def mass_budget(initial: float, injected: float, remaining: float, decay: float, boundary: float) -> tuple[float, float]:
    residual = initial + injected - remaining - decay - boundary
    denom = initial + injected
    fraction = abs(residual) / denom if denom > 0 else abs(residual)
    return residual, fraction


def integrate_budget(geom: Geometry, t: np.ndarray, c: np.ndarray) -> tuple[float, float, float]:
    decay_rate = K_DECAY * np.sum(c, axis=1) * VOXEL * CONV
    bound_rate = np.array([boundary_rate(geom, row) for row in c])
    iface_rate = np.array([interface_rate(geom, row) for row in c])
    decay = float(np.trapezoid(decay_rate, t))
    boundary = float(np.trapezoid(bound_rate, t))
    interface = float(np.trapezoid(iface_rate, t))
    return decay, boundary, interface


def corner_interior(geom: Geometry, field: np.ndarray) -> dict[str, float]:
    pocket = geom.pocket
    ii, jj = np.where(pocket)
    i0, i1 = int(ii.min()), int(ii.max())
    j0, j1 = int(jj.min()), int(jj.max())
    corners = {
        "SW": (i0, j0),
        "SE": (i1, j0),
        "NW": (i0, j1),
        "NE": (i1, j1),
    }
    interior = np.zeros_like(pocket)
    interior[i0 + 2 : i1 - 1, j0 + 2 : j1 - 1] = True
    interior &= pocket
    interior_vals = field[interior]
    interior_mean = float(np.nanmean(interior_vals)) if interior_vals.size else float("nan")

    def at(key: str) -> float:
        i, j = corners[key]
        return float(field[i, j])

    closed = 0.5 * (at("SW") + at("SE"))
    all4 = 0.25 * (at("SW") + at("SE") + at("NW") + at("NE"))
    ratio = (lambda num: num / interior_mean if interior_mean and interior_mean > 0 else float("inf"))
    wall = np.zeros_like(pocket)
    wall[i0, j0 : j1 + 1] = True
    wall[i1, j0 : j1 + 1] = True
    wall[i0 : i1 + 1, j0] = True
    if geom.arm_id == "CLOSED_NOFLUX":
        wall[i0 : i1 + 1, j1] = True
    wall &= pocket
    wall &= ~interior
    wall_mean = float(np.nanmean(field[wall])) if np.any(wall) else float("nan")
    return {
        "C_SW": at("SW"),
        "C_SE": at("SE"),
        "C_NW": at("NW"),
        "C_NE": at("NE"),
        "C_interior": interior_mean,
        "CornerInteriorRatio_closed_end": ratio(closed),
        "CornerInteriorRatio_all_four": ratio(all4),
        "WallInteriorRatio": wall_mean / interior_mean if interior_mean > 0 else float("inf"),
    }


def plateau_ok(t: np.ndarray, mean_r: np.ndarray, last_fraction: float, rel_tol: float) -> bool:
    t1 = float(t[-1])
    t0 = t1 * (1.0 - last_fraction)
    tail = mean_r[t >= t0]
    if tail.size == 0:
        return False
    terminal = float(mean_r[-1])
    ref = float(np.mean(tail))
    if abs(ref) < 1e-15:
        return abs(terminal) < rel_tol
    return abs(terminal - ref) / abs(ref) <= rel_tol
