#!/usr/bin/env python3
"""Reduced 2-D PocketNeck-T2 transport + Hill/L surrogate.

Clone of T0 transport_model.py (cell-centred FV, sparse stencil, BDF).
Neck voxels are inserted *between* pocket and bus; pocket volume is not
eaten. Occupancy means use pocket voxels only. No cells, no NARMA.
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

T0_OPEN_BUS_3_MEAN_R = 0.0318


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
BUS_LENGTH = float(DISH["bus_length_um"])
BUS_WIDTH = float(DISH["bus_width_um"])
POCKET_X0 = float(DISH["pocket_x0_um"])
POCKET_Y0 = float(DISH["pocket_y0_um"])
BUS_X0 = float(DISH["bus_x0_um"])
FLUSH_BUS_Y0 = float(DISH["flush_bus_y0_um"])
NECKED_BUS_Y0 = float(DISH["necked_bus_y0_um"])
NECK_Y0 = float(DISH["neck_y0_um"])
NECK_LENGTH = float(DISH["neck_length_um"])


def apply_grid(dx: float, dy: float | None = None) -> None:
    """Override grid spacing (CONVERGENCE extra only). Restores via caller."""
    global DX, DY, VOXEL
    DX = float(dx)
    DY = float(dx if dy is None else dy)
    VOXEL = DX * DY * DEPTH


def restore_primary_grid() -> None:
    apply_grid(float(DISH["dx_um"]), float(DISH["dy_um"]))


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


def t0_flush_match(mean_r: float) -> bool:
    target = float(DISH.get("T0_OPEN_BUS_3_mean_R", T0_OPEN_BUS_3_MEAN_R))
    tol = float(DISH.get("T0_OPEN_BUS_3_mean_R_rel_tol", 0.10))
    if not math.isfinite(mean_r) or target <= 0:
        return False
    return abs(mean_r - target) / target <= tol


def transport_scales(arm_id: str) -> dict[str, float]:
    length = POCKET_LX
    l2d = length * length / D
    da = K_DECAY * l2d
    spec = ARMS[arm_id]
    v = float(spec["bus_flow_um_s"])
    pe = v * BUS_WIDTH / D if v else 0.0
    transit = BUS_LENGTH / v if v else float("inf")
    return {
        "L_um": length,
        "L2_over_D_s": l2d,
        "Damkohler": da,
        "Pe": pe,
        "Transit_s": transit,
        "tau_AHL_s": 1.0 / K_DECAY,
        "W_um": float(spec["W_um"]),
        "L_n_um": float(spec["L_n_um"]),
    }


@dataclass
class Geometry:
    arm_id: str
    nx: int
    ny: int
    active: np.ndarray
    pocket: np.ndarray
    neck: np.ndarray
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
    dx: float
    dy: float
    voxel: float
    interface_pairs: list[tuple[int, int, int, int]] = field(default_factory=list)

    @property
    def x_centers(self) -> np.ndarray:
        return self.origin_x + (np.arange(self.nx) + 0.5) * self.dx

    @property
    def y_centers(self) -> np.ndarray:
        return self.origin_y + (np.arange(self.ny) + 0.5) * self.dy


def _empty_grid(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    active = np.zeros((nx, ny), dtype=bool)
    pocket = np.zeros((nx, ny), dtype=bool)
    neck = np.zeros((nx, ny), dtype=bool)
    bus = np.zeros((nx, ny), dtype=bool)
    return active, pocket, neck, bus


def _index_active(active: np.ndarray) -> tuple[np.ndarray, int]:
    idx = np.full(active.shape, -1, dtype=int)
    n = int(active.sum())
    idx[active] = np.arange(n)
    return idx, n


def _assemble_diffusion(M: lil_matrix, geom_idx: np.ndarray, active: np.ndarray, dx: float, dy: float) -> None:
    nx, ny = active.shape
    coef_x = D / (dx * dx)
    coef_y = D / (dy * dy)
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


def _add_upwind_x(M: lil_matrix, geom_idx: np.ndarray, mask: np.ndarray, velocity: float, dx: float) -> None:
    if velocity == 0:
        return
    nx, ny = mask.shape
    coef = velocity / dx
    for i in range(nx):
        for j in range(ny):
            if not mask[i, j]:
                continue
            p = int(geom_idx[i, j])
            M[p, p] -= coef
            if i > 0 and mask[i - 1, j]:
                M[p, int(geom_idx[i - 1, j])] += coef


def _add_upwind_y(M: lil_matrix, geom_idx: np.ndarray, mask: np.ndarray, velocity: float, dy: float) -> None:
    if velocity == 0:
        return
    nx, ny = mask.shape
    coef = velocity / dy
    for i in range(nx):
        for j in range(ny):
            if not mask[i, j]:
                continue
            p = int(geom_idx[i, j])
            M[p, p] -= coef
            if j > 0 and mask[i, j - 1]:
                M[p, int(geom_idx[i, j - 1])] += coef


def _source_b(idx: np.ndarray, ij: tuple[int, int], n: int, voxel: float) -> np.ndarray:
    b = np.zeros(n, dtype=float)
    p = int(idx[ij])
    if p < 0:
        raise RuntimeError(f"source cell {ij} is inactive")
    b[p] = J_MAX / (CONV * voxel)
    return b


def _interface_pairs(upstream: np.ndarray, downstream: np.ndarray) -> list[tuple[int, int, int, int]]:
    nx, ny = upstream.shape
    pairs = []
    for i in range(nx):
        for j in range(ny - 1):
            if upstream[i, j] and downstream[i, j + 1]:
                pairs.append((i, j, i, j + 1))
            if downstream[i, j] and upstream[i, j + 1]:
                pairs.append((i, j + 1, i, j))
    return pairs


def _neck_x_window(width_um: float) -> tuple[float, float]:
    x_mid = POCKET_X0 + 0.5 * POCKET_LX
    return x_mid - 0.5 * width_um, x_mid + 0.5 * width_um


def make_geometry(arm_id: str) -> Geometry:
    spec = ARMS[arm_id]
    dx, dy, voxel = DX, DY, VOXEL
    necked = bool(spec["neck"])
    if necked:
        domain_y = POCKET_LY + float(spec["L_n_um"]) + BUS_WIDTH
        bus_y0 = NECKED_BUS_Y0
    else:
        domain_y = POCKET_LY + BUS_WIDTH
        bus_y0 = FLUSH_BUS_Y0
    nx = int(round(BUS_LENGTH / dx))
    ny = int(round(domain_y / dy))
    active, pocket, neck, bus = _empty_grid(nx, ny)
    x0p, y0p = POCKET_X0, POCKET_Y0
    x0b = BUS_X0
    w_um = float(spec["W_um"])
    ln_um = float(spec["L_n_um"])
    neck_x0, neck_x1 = _neck_x_window(w_um)
    for i in range(nx):
        x = (i + 0.5) * dx
        for j in range(ny):
            y = (j + 0.5) * dy
            in_p = (x0p <= x < x0p + POCKET_LX) and (y0p <= y < y0p + POCKET_LY)
            in_b = (x0b <= x < x0b + BUS_LENGTH) and (bus_y0 <= y < bus_y0 + BUS_WIDTH)
            in_n = False
            if necked:
                in_n = (neck_x0 <= x < neck_x1) and (NECK_Y0 <= y < NECK_Y0 + ln_um)
            if in_p and in_n:
                raise RuntimeError(f"{arm_id}: neck stole a pocket voxel at ({i},{j})")
            pocket[i, j] = in_p
            neck[i, j] = in_n
            bus[i, j] = in_b
            active[i, j] = in_p or in_n or in_b
    origin_x, origin_y = 0.0, 0.0
    sx = int((x0p + float(DISH["source_pocket_um"][0])) / dx)
    sy = int((y0p + float(DISH["source_pocket_um"][1])) / dy)

    idx, n = _index_active(active)
    M = lil_matrix((n, n), dtype=float)
    _assemble_diffusion(M, idx, active, dx, dy)
    _add_decay(M, n)
    if spec["boundary"] == "BUS_COUPLED":
        _add_upwind_x(M, idx, bus, float(spec["bus_flow_um_s"]), dx)
    M_csr = M.tocsr()
    source_ij = (sx, sy)
    if not active[sx, sy] or not pocket[sx, sy]:
        raise RuntimeError(f"{arm_id}: source {source_ij} not in pocket")
    downstream = neck if necked else bus
    return Geometry(
        arm_id=arm_id,
        nx=nx,
        ny=ny,
        active=active,
        pocket=pocket,
        neck=neck,
        bus=bus,
        idx=idx,
        n=n,
        M=M_csr,
        source_b=_source_b(idx, source_ij, n, voxel),
        source_ij=source_ij,
        pocket_flow=float(spec["pocket_flow_um_s"]),
        bus_flow=float(spec["bus_flow_um_s"]),
        origin_x=origin_x,
        origin_y=origin_y,
        dx=dx,
        dy=dy,
        voxel=voxel,
        interface_pairs=_interface_pairs(pocket, downstream),
    )


def unpack(geom: Geometry, c_active: np.ndarray) -> np.ndarray:
    field = np.full((geom.nx, geom.ny), np.nan, dtype=float)
    field[geom.active] = c_active
    return field


def pocket_values(geom: Geometry, c_active: np.ndarray) -> np.ndarray:
    field = unpack(geom, c_active)
    return field[geom.pocket]


def mass_molecules(geom: Geometry, c_active: np.ndarray) -> float:
    return float(np.sum(c_active) * geom.voxel * CONV)


def injected_molecules(u_on: float, t_on: float) -> float:
    return J_MAX * u_on * t_on


def boundary_rate(geom: Geometry, c_active: np.ndarray) -> float:
    dcdt = geom.M @ c_active
    decay = K_DECAY * float(np.sum(c_active)) * geom.voxel * CONV
    return float(-np.sum(dcdt) * geom.voxel * CONV - decay)


def interface_rate(geom: Geometry, c_active: np.ndarray) -> float:
    if not geom.interface_pairs:
        return 0.0
    field = unpack(geom, c_active)
    rate = 0.0
    area = geom.dx * DEPTH
    for i_p, j_p, i_d, j_d in geom.interface_pairs:
        cp, cd = field[i_p, j_p], field[i_d, j_d]
        rate += D * (cp - cd) / geom.dy * area * CONV
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
    sample_dt: float = 20.0,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
) -> dict:
    drive = DRIVES[drive_id]
    t_on = float(drive["t_on_s"] if duration_on is None else duration_on)
    t_off = float(drive["t_off_s"] if duration_off is None else duration_off)
    u_on = float(drive["u_hold"])
    if sample_dt > 60.0 + 1e-12:
        raise RuntimeError(f"sample_dt={sample_dt} s exceeds the T2 cap of 60 s")
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


def mass_budget(initial: float, injected: float, remaining: float, decay: float, boundary: float) -> tuple[float, float]:
    residual = initial + injected - remaining - decay - boundary
    denom = initial + injected
    fraction = abs(residual) / denom if denom > 0 else abs(residual)
    return residual, fraction


def residual_vs_dominant(residual: float, injected: float, decay: float, boundary: float, interface: float) -> tuple[float, float]:
    dominant = max(abs(injected), abs(decay), abs(boundary), abs(interface), 1.0)
    return dominant, abs(residual) / dominant


def integrate_budget(geom: Geometry, t: np.ndarray, c: np.ndarray) -> tuple[float, float, float]:
    decay_rate = K_DECAY * np.sum(c, axis=1) * geom.voxel * CONV
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


def voxel_report(geom: Geometry) -> dict[str, object]:
    spec = ARMS[geom.arm_id]
    w = float(spec["W_um"])
    ln = float(spec["L_n_um"])
    w_over = w / geom.dx
    ln_over = ln / geom.dy if ln > 0 else float("nan")
    necked = bool(spec["neck"])
    if necked:
        w_int = abs(w_over - round(w_over)) < 1e-9
        ln_int = abs(ln_over - round(ln_over)) < 1e-9
    else:
        w_int = True
        ln_int = True
    n_pocket = int(geom.pocket.sum())
    expected_pocket = int(round(POCKET_LX / geom.dx) * round(POCKET_LY / geom.dy))
    n_neck = int(geom.neck.sum())
    expected_neck = int(round(w_over) * round(ln_over)) if necked else 0
    return {
        "arm": geom.arm_id,
        "dx_um": geom.dx,
        "W_um": w,
        "L_n_um": ln,
        "W_over_dx": w_over,
        "Ln_over_dx": ln_over,
        "W_over_dx_integer": w_int,
        "Ln_over_dx_integer": ln_int,
        "n_pocket": n_pocket,
        "n_neck": n_neck,
        "n_bus": int(geom.bus.sum()),
        "n_active": geom.n,
        "n_interface": len(geom.interface_pairs),
        "expected_pocket": expected_pocket,
        "expected_neck": expected_neck,
        "pocket_volume_ok": n_pocket == expected_pocket,
        "neck_count_ok": n_neck == expected_neck,
        "source_in_pocket": bool(geom.pocket[geom.source_ij]),
    }


def print_voxel_check(geom: Geometry) -> dict[str, object]:
    row = voxel_report(geom)
    print(
        f"VOXEL {row['arm']}: W/dx={row['W_over_dx']:.4g} "
        f"Ln/dx={row['Ln_over_dx'] if math.isfinite(float(row['Ln_over_dx'])) else 'NA'} "
        f"W_int={row['W_over_dx_integer']} Ln_int={row['Ln_over_dx_integer']} "
        f"n_pocket={row['n_pocket']} (expect {row['expected_pocket']}) "
        f"n_neck={row['n_neck']} (expect {row['expected_neck']}) "
        f"n_bus={row['n_bus']} n_iface={row['n_interface']}"
    )
    if not row["pocket_volume_ok"]:
        raise RuntimeError(f"{geom.arm_id}: pocket voxel count changed; closed-box V would drift")
    if not row["neck_count_ok"]:
        raise RuntimeError(f"{geom.arm_id}: neck voxel count {row['n_neck']} != {row['expected_neck']}")
    if ARMS[geom.arm_id]["neck"] and not (row["W_over_dx_integer"] and row["Ln_over_dx_integer"]):
        raise RuntimeError(f"{geom.arm_id}: W/dx or Ln/dx is not an integer")
    return row
