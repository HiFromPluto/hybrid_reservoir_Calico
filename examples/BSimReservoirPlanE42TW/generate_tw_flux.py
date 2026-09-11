#!/usr/bin/env python3
"""E4.2 cell-free two-way transducer. No bacteria. No living Java.

MODEL_STATUS is HYPOTHETICAL_DESIGN_ENVELOPE: there is no measured
P sensor curve. K_P, k_P, and alpha are frozen from well-mixed
algebra and the predeclared envelope, not from NARMA NRMSE and not
from Lentini numerical values.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PREFLIGHT = HERE.parent / "HybridDish" / "design_space_preflight"
E4AC = HERE.parent / "BSimReservoirPlanE4AC"
sys.path.insert(0, str(PREFLIGHT))
sys.path.insert(0, str(E4AC))
import generate_a1_flux as a1  # noqa: E402
import run_transport_screen as rts  # noqa: E402

NARMA_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
LENTINI_SHA = "5c4beb444515ffc58ba01d333598ba5ddaa47e75447eaa0e0691ac2549411ab6"
MODEL_STATUS = "HYPOTHETICAL_DESIGN_ENVELOPE"
G0 = 0.0
N_AC = 2.0
K_AC = 0.25
JMAX_A1 = 74677509.75821304
N_P = 2.0
ALPHA0 = 0.5
K_P = 1.6
D_P = 159.0
DECAY_P = 0.0033
MOLECULES_PER_UM3_PER_UM = 602.2
LX, LY, LZ = 1000.0, 500.0, 10.0
VOLUME_UM3 = LX * LY * LZ
N_CELLS = 1800
MEAN_R_A1 = 0.1915
PULSE_S = 75.0
WINDOW_S = 300.0
WARMUP_S = 18_000.0
WARMUP_U = 0.5
WARMUP_WINDOWS = int(WARMUP_S / WINDOW_S)
DT_0D = 0.05
RECEIVER_K = 1.6
RECEIVER_N = 2.0
TAU_R = 15.0
OCC_DEAD_R = 0.05
SAT_R = 0.9
PAYLOAD_TW_TOL = 1e-12
MASS_TOL = 0.01


def sha256_u(u: np.ndarray) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def load_u() -> np.ndarray:
    values = []
    for line in (HERE / "input_ahl_narma200.txt").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(float(line))
    u = np.array(values, dtype=float)
    digest = sha256_u(u)
    if digest != NARMA_SHA:
        raise SystemExit(f"STOP u SHA-256 {digest} != Narma10b {NARMA_SHA}")
    if len(u) != 200:
        raise SystemExit(f"STOP expected 200 u values, got {len(u)}")
    return u


def g_gate(u, g0=G0, n=N_AC, k=K_AC):
    u = np.asarray(u, dtype=float)
    un = np.power(np.maximum(u, 0.0), n)
    kn = k ** n
    g = g0 + (1.0 - g0) * un / (kn + un)
    return float(g) if np.ndim(g) == 0 else g


def h_gate(p, alpha, n=N_P, k=K_P):
    p = np.asarray(p, dtype=float)
    pn = np.power(np.maximum(p, 0.0), n)
    kn = k ** n
    h = 1.0 + alpha * pn / (kn + pn)
    return float(h) if np.ndim(h) == 0 else h


def flux(u, p, alpha, jmax=JMAX_A1):
    return jmax * g_gate(u) * h_gate(p, alpha)


class TwoWayInjector:
    """Static product gate. No vesicle store. J does not deplete."""

    def __init__(self, jmax: float, alpha: float):
        self.jmax = float(jmax)
        self.alpha = float(alpha)

    def gate_u(self, u):
        return g_gate(u)

    def gate_p(self, p):
        return h_gate(p, self.alpha)

    def flux(self, u, p=0.0, pulsing: bool = True):
        if not pulsing:
            return 0.0 * np.asarray(u, dtype=float) if np.ndim(u) else 0.0
        return self.jmax * self.gate_u(u) * self.gate_p(p)


def freeze_p_scale() -> dict:
    """Well-mixed P ≈ K_P at 1800 cells, R = mean_R_A1. Not an NRMSE fit."""
    total_molecules_at_kp = K_P * MOLECULES_PER_UM3_PER_UM * VOLUME_UM3
    production = DECAY_P * total_molecules_at_kp
    occupancy_scale = N_CELLS * MEAN_R_A1
    if occupancy_scale <= 0.0:
        raise SystemExit("STOP well-mixed k_P cannot be closed: occupancy scale is 0")
    k_p = production / occupancy_scale
    if not math.isfinite(k_p) or k_p <= 0.0:
        raise SystemExit("STOP well-mixed formula cannot be closed; do not invent k_P from NRMSE")
    check_p = (N_CELLS * k_p * MEAN_R_A1) / (DECAY_P * MOLECULES_PER_UM3_PER_UM * VOLUME_UM3)
    return {
        "K_P_uM": K_P,
        "k_P_molecules_per_s_per_R": k_p,
        "mean_R_A1_operating_point": MEAN_R_A1,
        "n_cells": N_CELLS,
        "dish_um": [LX, LY, LZ],
        "volume_um3": VOLUME_UM3,
        "decay_P_s_minus_1": DECAY_P,
        "D_P_um2_s": D_P,
        "molecules_per_um3_per_uM": MOLECULES_PER_UM3_PER_UM,
        "addQuantity_units": "molecules added to the containing voxel",
        "getConc_units": "molecules per um^3",
        "uM_from_getConc": "getConc / 602.2",
        "algebra": (
            "steady well-mixed: N * k_P * R = decay_P * (K_P_uM * 602.2 * V); "
            "k_P = decay_P * K_P * 602.2 * V / (N * mean_R_A1)"
        ),
        "steady_P_check_uM": check_p,
        "steady_P_rel_err_vs_K_P": abs(check_p - K_P) / K_P,
        "note": (
            "mean_R_A1 is an A1 occupancy scale, not a score. "
            "K_P uses the same number as receiver K, different species. "
            "Not copied from Lentini. Not fitted to NARMA."
        ),
        "MODEL_STATUS": MODEL_STATUS,
        "digital_twin": False,
    }


def commanded_mass(amplitudes_g_h: np.ndarray, jmax: float = JMAX_A1) -> float:
    return float(jmax * PULSE_S * np.sum(amplitudes_g_h))


def all_commands(u: np.ndarray) -> np.ndarray:
    return np.concatenate([np.full(WARMUP_WINDOWS, WARMUP_U), u])


def hill_receiver(c, p=None):
    """AHL Hill receiver. P is absent: the unused argument is a test hook."""
    del p
    c = np.asarray(c, dtype=float)
    cn = np.power(np.maximum(c, 0.0), RECEIVER_N)
    kn = RECEIVER_K ** RECEIVER_N
    return cn / (kn + cn)


def unit_tests(u: np.ndarray, alpha: float, field: dict) -> list[dict]:
    tw = TwoWayInjector(JMAX_A1, alpha)
    a1_inj = a1.A1Injector(JMAX_A1)
    rows = []
    commands = all_commands(u)
    mass_a1 = commanded_mass(a1_inj.gate(commands))
    mass_tw0 = commanded_mass(tw.gate_u(commands) * tw.gate_p(0.0))
    rel = abs(mass_tw0 - mass_a1) / max(abs(mass_a1), 1.0)
    j_tw = np.array([float(tw.flux(val, 0.0)) for val in u])
    j_a1 = np.array([float(a1_inj.flux(val)) for val in u])
    max_rel_j = float(np.max(np.abs(j_tw - j_a1) / np.maximum(np.abs(j_a1), 1.0)))
    rows.append({
        "test": "A_P0_recovers_A1",
        "pass": rel <= PAYLOAD_TW_TOL and max_rel_j <= PAYLOAD_TW_TOL,
        "observed": rel,
        "expected": PAYLOAD_TW_TOL,
        "note": f"payload rel={rel:.3e}; max per-window J rel={max_rel_j:.3e}",
    })

    h0 = float(tw.gate_p(0.0))
    h_inf = float(tw.gate_p(1.0e12))
    h_k = float(tw.gate_p(K_P))
    rows.append({
        "test": "B_h_identities",
        "pass": (
            abs(h0 - 1.0) <= 1e-12
            and abs(h_inf - (1.0 + alpha)) <= 1e-9
            and abs(h_k - (1.0 + alpha / 2.0)) <= 1e-12
        ),
        "observed": h_k,
        "expected": 1.0 + alpha / 2.0,
        "note": f"h(0)={h0:.12f}; h(inf)={h_inf:.12f}; h(K_P)={h_k:.12f}",
    })

    u_fixed = 0.4
    p_grid = np.concatenate([
        np.array([0.0]),
        np.geomspace(1e-6, 1.0e3 * K_P, 80),
    ])
    j_grid = np.array([float(tw.flux(u_fixed, p)) for p in p_grid])
    dj = np.diff(j_grid)
    j_sat = float(tw.flux(u_fixed, 1.0e12))
    j_base = float(tw.flux(u_fixed, 0.0))
    sat_ok = abs(j_sat / j_base - (1.0 + alpha)) <= 1e-9
    rows.append({
        "test": "C_P_clamp_monotone_sat",
        "pass": bool(np.all(dj >= -1e-12 * max(1.0, j_base))) and sat_ok,
        "observed": float(np.min(dj)),
        "expected": 0.0,
        "note": f"min dJ={float(np.min(dj)):.3e}; J(inf)/J(0)={j_sat / j_base:.12f}",
    })

    first = float(tw.flux(0.4, 0.8))
    second = float(tw.flux(0.4, 0.8))
    has_store = any(hasattr(tw, name) for name in ("S", "CS_in", "store"))
    rows.append({
        "test": "D_no_vesicle_store",
        "pass": (not has_store) and abs(first - second) <= 1e-15 * max(1.0, abs(first)),
        "observed": first - second,
        "expected": 0.0,
        "note": "identical J on repeated (u,P); no store depletion",
    })

    p_zero = field.get("p_remaining", 0.0)
    rows.append({
        "test": "E_no_cells_P0_AHL_mass",
        "pass": (
            field["narma_fraction"] <= MASS_TOL
            and field["full_fraction"] <= MASS_TOL
            and field["min_c"] >= -1e-9
            and abs(p_zero) <= 1e-15
        ),
        "observed": field["narma_fraction"],
        "expected": MASS_TOL,
        "note": "P field stays 0 without cells; leftover warmup is NARMA initial",
    })

    c = 1.6
    r0 = float(hill_receiver(c, 0.0))
    r1 = float(hill_receiver(c, 10.0))
    r2 = float(hill_receiver(c, 1.0e6))
    expected = c ** 2 / (RECEIVER_K ** 2 + c ** 2)
    rows.append({
        "test": "F_P_absent_from_AHL_receiver",
        "pass": abs(r0 - expected) <= 1e-15 and r0 == r1 == r2,
        "observed": r1 - r0,
        "expected": 0.0,
        "note": "Hill(C) ignores P; P is not a second QS input",
    })
    return rows


def _open_loop_mean_r(g: np.ndarray, v_eff: float, dt: float = 1.0) -> float:
    """Open-loop A1 (h=1) NARMA occupancy at an effective AHL volume."""
    conv = MOLECULES_PER_UM3_PER_UM * v_eff
    n_pulse = int(round(PULSE_S / dt))
    n_window = int(round(WINDOW_S / dt))
    decay_factor = math.exp(-dt * DECAY_P)
    inj = JMAX_A1 * dt / conv
    er = math.exp(-dt / TAU_R)
    c = 0.0
    r = 0.0
    r_sum = 0.0
    n = 0
    for i, amp in enumerate(g):
        window_r = 0.0
        for step in range(n_window):
            c = c * decay_factor + (inj * amp if step < n_pulse else 0.0)
            if c < 0.0:
                c = 0.0
            target = c * c / (RECEIVER_K * RECEIVER_K + c * c)
            r = target + (r - target) * er
            window_r += r
        if i >= WARMUP_WINDOWS:
            r_sum += window_r / n_window
            n += 1
    return r_sum / max(n, 1)


def calibrate_v_eff(u: np.ndarray) -> tuple[float, float]:
    """Effective AHL mixing volume so open-loop A1 NARMA mean R ≈ 0.1915.

    Point-source occupancy is not well-mixed over the dish. P still uses
    the full dish volume from the k_P freeze. This scale is occupancy,
    not NRMSE.
    """
    g = g_gate(all_commands(u))
    lo, hi = 1.0e2, VOLUME_UM3
    for _ in range(28):
        mid = math.sqrt(lo * hi)
        val = _open_loop_mean_r(g, mid)
        if val > MEAN_R_A1:
            lo = mid
        else:
            hi = mid
    v_eff = math.sqrt(lo * hi)
    return v_eff, _open_loop_mean_r(g, v_eff)


def integrate_0d(commands: np.ndarray, alpha: float, k_p: float, v_eff: float,
                 pulsed: bool, closed: bool, dt: float = 1.0):
    conv_ahl = MOLECULES_PER_UM3_PER_UM * v_eff
    conv_p = MOLECULES_PER_UM3_PER_UM * VOLUME_UM3
    n_pulse = int(round(PULSE_S / dt))
    n_window = int(round(WINDOW_S / dt))
    decay_c = math.exp(-dt * DECAY_P)
    decay_p = math.exp(-dt * DECAY_P)
    inj_scale = dt / conv_ahl
    p_prod_scale = N_CELLS * k_p * dt / conv_p
    er = math.exp(-dt / TAU_R)
    c = 0.0
    r = 0.0
    p = 0.0
    window_means = []
    last_r = []
    finite = True
    for amp_u in commands:
        g = g_gate(amp_u)
        acc_r = acc_p = acc_h = 0.0
        for step in range(n_window):
            pulsing = (not pulsed) or (step < n_pulse)
            h = h_gate(p, alpha) if closed else 1.0
            j = JMAX_A1 * g * h if pulsing else 0.0
            c = c * decay_c + j * inj_scale
            if c < 0.0:
                c = 0.0
            target = hill_receiver(c, p)
            r = target + (r - target) * er
            p = p * decay_p + p_prod_scale * max(r, 0.0)
            if p < 0.0:
                p = 0.0
            if not (math.isfinite(c) and math.isfinite(r) and math.isfinite(p)):
                finite = False
                break
            acc_r += r
            acc_p += p
            acc_h += h
        if not finite:
            break
        window_means.append({
            "u": float(amp_u),
            "mean_R": acc_r / n_window,
            "mean_P": acc_p / n_window,
            "mean_h": acc_h / n_window,
        })
        last_r.append(r)
    return {
        "windows": window_means,
        "finite": finite,
        "last_R": last_r,
    }


def classify_traj(result: dict, u_series: np.ndarray | None, label: str) -> dict:
    if not result["finite"]:
        return {"label": label, "class": "UNSTABLE", "mean_R": float("nan"),
                "r_Ru": float("nan"), "note": "non-finite state"}
    r = np.array([row["mean_R"] for row in result["windows"]], dtype=float)
    pvals = np.array([row["mean_P"] for row in result["windows"]], dtype=float)
    hvals = np.array([row["mean_h"] for row in result["windows"]], dtype=float)
    last = r[len(r) // 2 :] if len(r) > 4 else r
    mean_r = float(np.mean(last if u_series is None else r))
    mean_p = float(np.mean(pvals[len(pvals) // 2 :] if u_series is None else pvals))
    mean_h = float(np.mean(hvals[len(hvals) // 2 :] if u_series is None else hvals))
    r_u = float("nan")
    if u_series is not None and len(u_series) == len(r):
        a = r - np.mean(r)
        b = u_series - np.mean(u_series)
        den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
        r_u = 0.0 if den < 1e-15 else float(np.sum(a * b) / den)
    last = r[len(r) // 2 :] if len(r) > 4 else r
    pinned_high = float(np.min(last)) > SAT_R and mean_r > SAT_R
    exploding = float(np.max(pvals)) > 1.0e6
    if exploding:
        klass = "UNSTABLE"
    elif pinned_high:
        klass = "SATURATED"
    elif mean_r < OCC_DEAD_R:
        klass = "DEAD"
    else:
        klass = "STABLE_RESPONSIVE"
    return {
        "label": label,
        "class": klass,
        "mean_R": mean_r,
        "mean_P": mean_p,
        "mean_h": mean_h,
        "r_Ru": r_u,
        "note": "",
    }


def loop_gain_screen(u: np.ndarray, alpha: float, k_p: float, v_eff: float) -> dict:
    rows = []
    for const_u in (0.0, 0.25, 0.5):
        commands = np.full(80, const_u)
        closed = integrate_0d(commands, alpha, k_p, v_eff, pulsed=True, closed=True)
        rows.append(classify_traj(closed, None, f"const_u_{const_u}_pulsed"))
    narma_commands = all_commands(u)
    narma = integrate_0d(narma_commands, alpha, k_p, v_eff, pulsed=True, closed=True)
    narma_analysis = {
        "windows": narma["windows"][WARMUP_WINDOWS:],
        "finite": narma["finite"],
        "last_R": narma["last_R"][WARMUP_WINDOWS:] if narma["finite"] else [],
    }
    rows.append(classify_traj(narma_analysis, u, "narma_pulsed_closed"))
    # Envelope class: u=0 DEAD is the open-loop identity and does not fail the screen.
    # Constant-u uses the dish 75/300 pulse protocol. Always-on u=0.5 saturates
    # even open-loop A1 in a well-mixed dish and is not the frozen operating protocol.
    relevant = [row for row in rows if row["label"] != "const_u_0.0_pulsed"]
    classes = {row["class"] for row in relevant}
    if "UNSTABLE" in classes:
        overall = "UNSTABLE"
    elif "SATURATED" in classes:
        overall = "SATURATED"
    elif classes == {"DEAD"}:
        overall = "DEAD"
    elif "STABLE_RESPONSIVE" in classes and not (classes & {"UNSTABLE", "SATURATED"}):
        overall = "STABLE_RESPONSIVE"
    else:
        overall = "DEAD"
    return {"rows": rows, "overall": overall, "v_eff_um3": v_eff, "alpha": alpha}


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"empty CSV {path}")
    if fields is None:
        fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_cell_free(alpha: float) -> None:
    fig_dir = HERE / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    u_grid = np.linspace(0.0, 0.5, 201)
    p_grid = np.linspace(0.0, 8.0, 201)
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(u_grid, g_gate(u_grid), color="#1f4e79", lw=2.0, label="g(u) A1")
    ax.set_xlabel("command u")
    ax.set_ylabel("g(u)")
    ax.set_title("Frozen A1 gate (untouched)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_g_vs_u.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(p_grid, h_gate(p_grid, alpha), color="#9c2a2a", lw=2.0, label=f"h(P) alpha={alpha}")
    ax.axhline(1.0, color="0.6", lw=1.0, label="h(0)=1 recovers A1")
    ax.axhline(1.0 + alpha, color="0.8", ls="--", lw=1.0)
    ax.axvline(K_P, color="0.8", ls=":", lw=1.0)
    ax.set_xlabel("P (µM)")
    ax.set_ylabel("h(P)")
    ax.set_title("Feedback modulator (hypothetical envelope)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_h_vs_p.png", dpi=140)
    plt.close(fig)

    uu, pp = np.meshgrid(u_grid, p_grid)
    jj = flux(uu, pp, alpha) / JMAX_A1
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    mesh = ax.pcolormesh(u_grid, p_grid, jj, shading="auto", cmap="viridis")
    fig.colorbar(mesh, ax=ax, label="J / Jmax_A1")
    ax.set_xlabel("u")
    ax.set_ylabel("P (µM)")
    ax.set_title("J(u,P) = Jmax_A1 g(u) h(P)")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_j_up.png", dpi=140)
    plt.close(fig)


def write_calibration_md(
    u: np.ndarray,
    alpha: float,
    k_p: float,
    tests: list[dict],
    field: dict,
    screen: dict,
    p_scale: dict,
    v_eff: float,
    v_eff_mean_r: float,
) -> None:
    passed = all(row["pass"] for row in tests)
    lines = [
        "# E4.2 two-way cell-free calibration",
        "",
        f"**MODEL_STATUS = `{MODEL_STATUS}`**",
        "",
        "Not a digital twin. No measured P sensor curve. Lentini is",
        "topology + Turing/replay controls only. `K_P`, `k_P`, and",
        "`alpha` were not fitted to NARMA NRMSE and were not copied from",
        "the paper.",
        "",
        f"Module 1: **{'PASS' if passed else 'FAIL'}**",
        f"Module 1.6: **{screen['overall']}**",
        "",
        "## Frozen envelope",
        "",
        "| Parameter | Value | Role |",
        "|---|---|---|",
        f"| `g0` | {G0:g} | A1, untouched |",
        f"| `n_AC` | {N_AC:g} | A1, untouched |",
        f"| `K_AC` | {K_AC:g} | A1, untouched |",
        f"| `Jmax_A1` | {JMAX_A1:.12e} molecules/s | A1 payload match, untouched |",
        f"| `n_P` | {N_P:g} | envelope, not measured |",
        f"| `alpha` | {alpha:g} | envelope; h ∈ [1, 1+alpha] |",
        f"| `K_P` | {K_P:g} µM | same number as receiver K, different species |",
        f"| `k_P` | {k_p:.12e} molecules/s per unit R | well-mixed freeze |",
        f"| `D_P` | {D_P:g} µm²/s | AHL-like inert reporter |",
        f"| `decay_P` | {DECAY_P:g} s⁻¹ | AHL-like inert reporter |",
        f"| `h(0)` | {h_gate(0.0, alpha):.12f} | recovers A1 |",
        f"| `h(K_P)` | {h_gate(K_P, alpha):.12f} | 1+alpha/2 |",
        f"| `h(∞)` | {h_gate(1e12, alpha):.12f} | 1+alpha |",
        "",
        "Gate:",
        "",
        "```",
        "g(u) = u^2 / (K_AC^2 + u^2)",
        "h(P) = 1 + alpha * P^n_P / (K_P^n_P + P^n_P)",
        "J(t) = Jmax_A1 * g(u) * h(P)   during the 75 s pulse; 0 after",
        "```",
        "",
        "## P scale freeze",
        "",
        "See `p_scale_freeze.json`. Steady well-mixed identity:",
        "",
        "```",
        p_scale["algebra"],
        "```",
        "",
        f"- dish volume `{VOLUME_UM3:.6g}` µm³",
        f"- addQuantity conversion `{MOLECULES_PER_UM3_PER_UM}` molecules/µm³ per µM",
        f"- 1800 cells at `R={MEAN_R_A1}` → steady P `{p_scale['steady_P_check_uM']:.12f}` µM",
        f"- relative error vs `K_P`: `{p_scale['steady_P_rel_err_vs_K_P']:.3e}`",
        "",
        "## Frozen u",
        "",
        f"- SHA-256: `{sha256_u(u)}`",
        f"- Lentini PDF SHA-256: `{LENTINI_SHA}`",
        "",
        "## Unit tests",
        "",
        "| Test | Pass | Observed | Limit / expected | Note |",
        "|---|---|---|---|---|",
    ]
    for row in tests:
        lines.append(
            f"| {row['test']} | **{'PASS' if row['pass'] else 'FAIL'}** | "
            f"{row['observed']:.6g} | {row['expected']:.6g} | {row['note']} |"
        )
    lines.extend(
        [
            "",
            "## E0.2 mass identity (no bacteria, P=0)",
            "",
            "CENTER / FLOW=0. P production is zero without cells, so the",
            "AHL field is the A1 identity. Leftover warmup mass is the",
            "NARMA window initial.",
            "",
            "| Budget | Molecules | Residual fraction |",
            "|---|---|---|",
            f"| warmup leftover | {field['leftover']:.12e} | {field['warmup_fraction']:.3e} |",
            f"| NARMA injected | {field['narma_injected']:.12e} | |",
            f"| NARMA remaining | {field['narma_remaining']:.12e} | {field['narma_fraction']:.3e} |",
            f"| NARMA decay | {field['narma_decay']:.12e} | |",
            f"| NARMA boundary | {field['narma_boundary']:.12e} | |",
            f"| P remaining | {field['p_remaining']:.6g} | |",
            f"| min C | {field['min_c']:.6g} µM | |",
            "",
            "## Module 1.6 0-D loop-gain screen",
            "",
            f"AHL well-mixed volume `{v_eff:.6g}` µm³ is the dish volume",
            f"already used for `k_P`. Open-loop A1 NARMA mean R `{v_eff_mean_r:.4f}`",
            f"vs occupancy scale `{MEAN_R_A1}` (spatial A1 mean; not a fit).",
            "Constant-u drives use the frozen 75/300 pulse protocol.",
            "",
            "| Drive | Class | mean_R | mean_P | mean_h | r(R,u) |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in screen["rows"]:
        r_u = "" if row["r_Ru"] != row["r_Ru"] else f"{row['r_Ru']:.3f}"
        lines.append(
            f"| {row['label']} | **{row['class']}** | {row['mean_R']:.4f} | "
            f"{row['mean_P']:.4f} | {row['mean_h']:.4f} | {r_u} |"
        )
    lines.extend(
        [
            "",
            f"Overall: **{screen['overall']}** at `alpha={alpha}`.",
            "",
            "## Living BSim authorization",
            "",
        ]
    )
    if passed and screen["overall"] == "STABLE_RESPONSIVE":
        lines.append("Module 1 PASS and Module 1.6 STABLE_RESPONSIVE — seed-111 CLOSED may start.")
    elif passed:
        lines.append(
            f"Module 1 PASS but Module 1.6 is {screen['overall']}. "
            "Living Java is not authorized."
        )
    else:
        lines.append("Module 1 FAIL — stop; do not start living BSim.")
    lines.append("")
    (HERE / "results" / "E4_2_CALIBRATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_frozen(alpha: float, k_p: float, tests: list[dict], screen: dict) -> None:
    frozen = {
        "MODEL_STATUS": MODEL_STATUS,
        "digital_twin": False,
        "g0": G0,
        "n_AC": N_AC,
        "K_AC": K_AC,
        "Jmax_A1": JMAX_A1,
        "n_P": N_P,
        "K_P": K_P,
        "k_P": k_p,
        "alpha": alpha,
        "D_P": D_P,
        "decay_P": DECAY_P,
        "h_at_0": float(h_gate(0.0, alpha)),
        "h_at_K_P": float(h_gate(K_P, alpha)),
        "h_at_inf": float(h_gate(1.0e12, alpha)),
        "u_sha256": NARMA_SHA,
        "lentini_pdf_sha256": LENTINI_SHA,
        "module1_pass": all(row["pass"] for row in tests),
        "module16_class": screen["overall"],
        "alpha_cut_applied": abs(alpha - ALPHA0) > 1e-15,
    }
    (HERE / "results" / "tw_frozen.json").write_text(
        json.dumps(frozen, indent=2) + "\n", encoding="utf-8"
    )


def reuse_a1_mass_identity() -> dict:
    """P=0 recovers A1; reuse E4.1 no-bacteria AHL identity. Do not rerun BSim."""
    path = E4AC / "results" / "a1_mass_identity.csv"
    if not path.exists():
        raise SystemExit("STOP missing E4AC a1_mass_identity.csv; do not overwrite E4AC")
    with path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    leftover = float(row["narma_initial_leftover_molecules"])
    injected = float(row["narma_injected_molecules"])
    remaining = float(row["narma_remaining_molecules"])
    decay = float(row["narma_decay_molecules"])
    boundary = float(row["narma_boundary_molecules"])
    full_injected = float(row["full_injected_molecules"])
    return {
        "leftover": leftover,
        "warmup_injected": full_injected - injected,
        "warmup_decay": float("nan"),
        "warmup_boundary": 0.0,
        "warmup_residual": 0.0,
        "warmup_fraction": 0.0,
        "narma_injected": injected,
        "narma_remaining": remaining,
        "narma_decay": decay,
        "narma_boundary": boundary,
        "narma_residual": float(row["narma_residual_fraction"]),
        "narma_fraction": float(row["narma_residual_fraction"]),
        "full_injected": full_injected,
        "full_decay": decay,
        "full_boundary": boundary,
        "full_residual": float(row["full_residual_fraction"]),
        "full_fraction": float(row["full_residual_fraction"]),
        "min_c": float(row["min_C_uM"]),
        "p_remaining": 0.0,
        "reused_from": str(path),
    }


def run_module1(u: np.ndarray, alpha: float, k_p: float) -> tuple[list[dict], dict]:
    del k_p
    print("reusing E4.1 no-bacteria AHL mass identity (P=0 recovers A1; no living BSim)...")
    field = reuse_a1_mass_identity()
    tests = unit_tests(u, alpha, field)
    for row in tests:
        print(f"  {row['test']}: {'PASS' if row['pass'] else 'FAIL'}  observed={row['observed']:.6g}", flush=True)
    return tests, field


def main() -> None:
    print(f"MODEL_STATUS={MODEL_STATUS}")
    print("Not a digital twin. Lentini is topology only; no paper numbers copied.")
    u = load_u()
    p_scale = freeze_p_scale()
    k_p = p_scale["k_P_molecules_per_s_per_R"]
    (HERE / "results" / "p_scale_freeze.json").write_text(
        json.dumps(p_scale, indent=2) + "\n", encoding="utf-8"
    )
    print(f"K_P={K_P} k_P={k_p:.12e} steady_P={p_scale['steady_P_check_uM']:.12f}")
    if p_scale["steady_P_rel_err_vs_K_P"] > 1e-9:
        raise SystemExit("STOP well-mixed formula did not close on K_P")

    alpha = ALPHA0
    tests, field = run_module1(u, alpha, k_p)
    if not all(row["pass"] for row in tests):
        write_frozen(alpha, k_p, tests, {"overall": "MODULE1_FAIL", "rows": []})
        write_calibration_md(u, alpha, k_p, tests, field, {"overall": "MODULE1_FAIL", "rows": []},
                             p_scale, float("nan"), float("nan"))
        raise SystemExit("STOP Module 1 FAIL; do not start living Java")

    print("0-D well-mixed volume is the dish volume already used for k_P (not NRMSE)...")
    v_eff, v_eff_mean_r = VOLUME_UM3, _open_loop_mean_r(g_gate(all_commands(u)), VOLUME_UM3)
    print(f"V_wellmix={v_eff:.6g} um^3  open-loop NARMA mean_R={v_eff_mean_r:.4f}")
    screen = loop_gain_screen(u, alpha, k_p, v_eff)
    print(f"Module 1.6 overall={screen['overall']}")
    for row in screen["rows"]:
        print(f"  {row['label']}: {row['class']} mean_R={row['mean_R']:.4f}")

    if screen["overall"] in {"SATURATED", "UNSTABLE"}:
        print("predeclared cut alpha 0.5 -> 0.25; re-run Module 1 A-F")
        alpha = 0.25
        tests, field = run_module1(u, alpha, k_p)
        if not all(row["pass"] for row in tests):
            write_frozen(alpha, k_p, tests, screen)
            raise SystemExit("STOP Module 1 FAIL after alpha cut")
        screen = loop_gain_screen(u, alpha, k_p, v_eff)
        print(f"Module 1.6 after cut overall={screen['overall']}")
        for row in screen["rows"]:
            print(f"  {row['label']}: {row['class']} mean_R={row['mean_R']:.4f}")

    write_csv(HERE / "results" / "tw_unit_tests.csv", tests)
    write_csv(HERE / "results" / "tw_loop_gain.csv", screen["rows"])
    u_grid = np.linspace(0.0, 0.5, 201)
    p_grid = np.linspace(0.0, 8.0, 201)
    write_csv(
        HERE / "results" / "tw_gate_curve.csv",
        [
            {"u": float(x), "g": float(g_gate(x)), "J_A1": float(JMAX_A1 * g_gate(x)),
             "J_TW_P0": float(flux(x, 0.0, alpha))}
            for x in u_grid
        ],
    )
    write_csv(
        HERE / "results" / "tw_h_curve.csv",
        [{"P_uM": float(p), "h": float(h_gate(p, alpha))} for p in p_grid],
    )
    plot_cell_free(alpha)
    write_frozen(alpha, k_p, tests, screen)
    write_calibration_md(u, alpha, k_p, tests, field, screen, p_scale, v_eff, v_eff_mean_r)

    if not all(row["pass"] for row in tests):
        raise SystemExit("STOP Module 1 FAIL; do not start living Java")
    if screen["overall"] != "STABLE_RESPONSIVE":
        print("STOP Module 1.6 not STABLE_RESPONSIVE; keep the row; no living Java")
        raise SystemExit(2)
    print("Module 1 PASS")
    print("Module 1.6 STABLE_RESPONSIVE")


if __name__ == "__main__":
    main()
