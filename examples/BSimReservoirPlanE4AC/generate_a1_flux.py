#!/usr/bin/env python3
"""E4.1 cell-free A1 transducer. No bacteria. No Java.

MODEL_STATUS is HYPOTHETICAL_DESIGN_ENVELOPE: there is no measured
u→flux curve. g0, n_AC, and K_AC are frozen engineering envelope
values. Jmax_A1 is chosen only to match A0 commanded mass. Nothing
here is fitted to NARMA NRMSE.
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
sys.path.insert(0, str(PREFLIGHT))
import run_transport_screen as rts  # noqa: E402

NARMA_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
MODEL_STATUS = "HYPOTHETICAL_DESIGN_ENVELOPE"
JMAX_A0 = 128_000_000.0
G0 = 0.0
N_AC = 2.0
K_AC = 0.25
PULSE_S = 75.0
WINDOW_S = 300.0
WARMUP_S = 18_000.0
WARMUP_U = 0.5
WARMUP_WINDOWS = int(WARMUP_S / WINDOW_S)
MASS_TOL = 0.01

NEAR_IJ = (25, 12)  # source voxel, (500, 250) µm on 20 µm grid
MID_IJ = (37, 12)   # +240 µm in x
FAR_IJ = (48, 12)   # near +x wall


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
    if np.any(u < 0.0) or np.any(u > 0.5):
        raise SystemExit("STOP u escaped [0, 0.5]")
    return u


class A1Injector:
    """Static Hill gate. No vesicle store. J does not deplete."""

    def __init__(self, jmax: float, g0: float = G0, n: float = N_AC, k: float = K_AC):
        self.jmax = float(jmax)
        self.g0 = float(g0)
        self.n = float(n)
        self.k = float(k)

    def gate(self, u: np.ndarray | float) -> np.ndarray | float:
        u = np.asarray(u, dtype=float)
        un = np.power(np.maximum(u, 0.0), self.n)
        kn = self.k ** self.n
        g = self.g0 + (1.0 - self.g0) * un / (kn + un)
        return float(g) if g.ndim == 0 else g

    def flux(self, u: np.ndarray | float, pulsing: bool = True) -> np.ndarray | float:
        if not pulsing:
            return 0.0 * np.asarray(u, dtype=float) if np.ndim(u) else 0.0
        return self.jmax * self.gate(u)


def commanded_mass(jmax: float, amplitudes: np.ndarray) -> float:
    return float(jmax * PULSE_S * np.sum(amplitudes))


def all_commands(u: np.ndarray) -> np.ndarray:
    return np.concatenate([np.full(WARMUP_WINDOWS, WARMUP_U), u])


def choose_jmax_a1(u: np.ndarray) -> tuple[float, dict[str, float]]:
    commands = all_commands(u)
    a0 = A1Injector(JMAX_A0)
    probe = A1Injector(1.0)
    mass_a0 = commanded_mass(JMAX_A0, commands)
    g_sum = float(np.sum(probe.gate(commands)))
    jmax_a1 = mass_a0 / (PULSE_S * g_sum)
    a1 = A1Injector(jmax_a1)
    mass_a1 = commanded_mass(jmax_a1, a1.gate(commands))
    rel = abs(mass_a1 - mass_a0) / mass_a0
    return jmax_a1, {
        "mass_a0": mass_a0,
        "mass_a1": mass_a1,
        "relative_error": rel,
        "warmup_u_sum": float(WARMUP_WINDOWS * WARMUP_U),
        "narma_u_sum": float(np.sum(u)),
        "warmup_g_sum": float(WARMUP_WINDOWS * a1.gate(WARMUP_U)),
        "narma_g_sum": float(np.sum(a1.gate(u))),
        "g_at_0": float(a1.gate(0.0)),
        "g_at_0p5": float(a1.gate(0.5)),
        "g_at_inf": float(a1.gate(1.0e9)),
    }


def pulse_timeseries(u: np.ndarray, injector: A1Injector, include_warmup: bool = True):
    commands = all_commands(u) if include_warmup else u
    start = 0.0 if include_warmup else WARMUP_S
    times = []
    fluxes = []
    gates = []
    commands_out = []
    cumulative = []
    mass = 0.0
    t = start
    dt_plot = 1.0
    for amplitude in commands:
        for elapsed in np.arange(0.0, WINDOW_S, dt_plot):
            pulsing = elapsed < PULSE_S
            j = injector.flux(amplitude, pulsing=pulsing)
            g = injector.gate(amplitude) if pulsing else 0.0
            times.append(t + elapsed)
            fluxes.append(float(j))
            gates.append(float(g))
            commands_out.append(float(amplitude))
            mass += float(j) * dt_plot
            cumulative.append(mass)
        t += WINDOW_S
    return {
        "t": np.array(times),
        "j": np.array(fluxes),
        "g": np.array(gates),
        "u": np.array(commands_out),
        "mass": np.array(cumulative),
    }


def voxel_index(i: int, j: int) -> int:
    return i * rts.NY + j


def integrate_field(effective_amplitudes: np.ndarray, collect: bool):
    condition = rts.make_conditions()[0]
    assert condition.condition_id == "PRI_CENTER_F0p0"
    state = rts.initial_state()
    warmup = effective_amplitudes[:WARMUP_WINDOWS]
    narma = effective_amplitudes[WARMUP_WINDOWS:]
    for amplitude in warmup:
        state, _ = rts.simulate_window(state, condition, float(amplitude), collect=False)
    leftover = rts.state_mass_molecules(state)
    warmup_decay = float(state[-2])
    warmup_boundary = float(state[-1])
    warmup_injected = rts.commanded_injected_mass(warmup)
    warmup_residual, warmup_fraction = rts.mass_budget(
        0.0, warmup_injected, leftover, warmup_decay, warmup_boundary
    )
    state = rts.reset_budget(state)
    samples = []
    for amplitude in narma:
        state, window_samples = rts.simulate_window(
            state, condition, float(amplitude), collect=collect
        )
        if collect and window_samples is not None:
            samples.append(window_samples)
    remaining = rts.state_mass_molecules(state)
    narma_injected = rts.commanded_injected_mass(narma)
    narma_decay = float(state[-2])
    narma_boundary = float(state[-1])
    narma_residual, narma_fraction = rts.mass_budget(
        leftover, narma_injected, remaining, narma_decay, narma_boundary
    )
    full_injected = warmup_injected + narma_injected
    full_decay = warmup_decay + narma_decay
    full_boundary = warmup_boundary + narma_boundary
    full_residual, full_fraction = rts.mass_budget(
        0.0, full_injected, remaining, full_decay, full_boundary
    )
    sample_stack = np.vstack(samples) if samples else None
    return {
        "leftover": leftover,
        "warmup_injected": warmup_injected,
        "warmup_decay": warmup_decay,
        "warmup_boundary": warmup_boundary,
        "warmup_residual": warmup_residual,
        "warmup_fraction": warmup_fraction,
        "narma_injected": narma_injected,
        "narma_remaining": remaining,
        "narma_decay": narma_decay,
        "narma_boundary": narma_boundary,
        "narma_residual": narma_residual,
        "narma_fraction": narma_fraction,
        "full_injected": full_injected,
        "full_decay": full_decay,
        "full_boundary": full_boundary,
        "full_residual": full_residual,
        "full_fraction": full_fraction,
        "min_c": float(np.min(state[: rts.N])),
        "samples": sample_stack,
        "final_state": state,
    }


def unit_tests(u: np.ndarray, jmax_a1: float, stats: dict[str, float], field: dict) -> list[dict]:
    a1 = A1Injector(jmax_a1)
    rows = []

    j0 = float(a1.flux(0.0))
    expected0 = jmax_a1 * G0
    rows.append({
        "test": "A_u0_equals_Jmax_g0",
        "pass": abs(j0 - expected0) <= 1e-12 * max(1.0, abs(jmax_a1)),
        "observed": j0,
        "expected": expected0,
        "note": "g0=0 so J(0)=0",
    })

    g05 = float(a1.gate(0.5))
    g_inf = float(a1.gate(1.0e9))
    expected05 = (0.5 ** N_AC) / (K_AC ** N_AC + 0.5 ** N_AC)
    rows.append({
        "test": "B_g_at_0p5_and_infinity",
        "pass": abs(g05 - expected05) <= 1e-12 and abs(g_inf - 1.0) <= 1e-9,
        "observed": g05,
        "expected": expected05,
        "note": f"g(0.5)={g05:.12f}; g(inf)={g_inf:.12f} (documented, not a fit)",
    })

    rows.append({
        "test": "C_narma_window_mass_identity",
        "pass": field["narma_fraction"] <= MASS_TOL
        and field["full_fraction"] <= MASS_TOL
        and field["min_c"] >= -1e-9,
        "observed": field["narma_fraction"],
        "expected": MASS_TOL,
        "note": "leftover warmup as NARMA initial; remaining+decay+boundary=initial+injected",
    })

    rows.append({
        "test": "D_payload_match_vs_A0",
        "pass": stats["relative_error"] <= MASS_TOL,
        "observed": stats["relative_error"],
        "expected": MASS_TOL,
        "note": "commanded mass warmup+200 windows, same u, including warmup 0.5",
    })

    first = float(a1.flux(0.4))
    second = float(a1.flux(0.4))
    has_store = hasattr(a1, "S") or hasattr(a1, "CS_in") or hasattr(a1, "store")
    rows.append({
        "test": "E_no_vesicle_store",
        "pass": (not has_store) and abs(first - second) <= 1e-15 * max(1.0, abs(first)),
        "observed": first - second,
        "expected": 0.0,
        "note": "identical J on a repeated u; no CS_in / store depletion",
    })
    return rows


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"empty CSV {path}")
    if fields is None:
        fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_cell_free(u: np.ndarray, jmax_a1: float, series_a0, series_a1, probes) -> None:
    fig_dir = HERE / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    a1 = A1Injector(jmax_a1)
    u_grid = np.linspace(0.0, 0.5, 201)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(u_grid, u_grid, color="0.6", lw=1.2, label="A0  g=u")
    ax.plot(u_grid, a1.gate(u_grid), color="#1f4e79", lw=2.0, label="A1  Hill gate")
    ax.axvline(0.5, color="0.8", ls="--", lw=0.8)
    ax.set_xlabel("command u")
    ax.set_ylabel("gate g(u)")
    ax.set_title("Cell-free A1 gate (hypothetical envelope)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_g_vs_u.pdf")
    fig.savefig(fig_dir / "fig_g_vs_u.png", dpi=140)
    plt.close(fig)

    t0 = series_a0["t"] / 3600.0
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.plot(t0, series_a0["j"] / 1e8, color="0.55", lw=0.8, label="A0 J")
    ax.plot(t0, series_a1["j"] / 1e8, color="#1f4e79", lw=0.8, label="A1 J")
    ax.set_xlabel("time (h, warmup + NARMA)")
    ax.set_ylabel("flux (10^8 molecules/s)")
    ax.set_title("Commanded flux, same u file")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_j_vs_time.pdf")
    fig.savefig(fig_dir / "fig_j_vs_time.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(t0, series_a0["mass"] / 1e11, color="0.55", lw=1.6, label="A0")
    ax.plot(t0, series_a1["mass"] / 1e11, color="#1f4e79", lw=1.6, label="A1")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("cumulative commanded mass (10^11 molecules)")
    ax.set_title("Payload match (commanded, not occupancy)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_cumulative_mass.pdf")
    fig.savefig(fig_dir / "fig_cumulative_mass.png", dpi=140)
    plt.close(fig)

    if probes is not None:
        t_narma = np.arange(len(probes["near"])) * (WINDOW_S / 16.0) / 3600.0
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
        ax.plot(t_narma, probes["near"], label="near (source)", lw=1.2)
        ax.plot(t_narma, probes["mid"], label="mid +240 µm", lw=1.2)
        ax.plot(t_narma, probes["far"], label="far near +x wall", lw=1.2)
        ax.set_xlabel("NARMA time (h)")
        ax.set_ylabel("AHL (µM)")
        ax.set_title("A1 cell-free plume probes (E0.2 integrator)")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(fig_dir / "fig_probes.pdf")
        fig.savefig(fig_dir / "fig_probes.png", dpi=140)
        plt.close(fig)


def write_calibration_md(
    u: np.ndarray,
    jmax_a1: float,
    stats: dict[str, float],
    tests: list[dict],
    field: dict,
) -> None:
    passed = all(row["pass"] for row in tests)
    lines = [
        "# E4.1 A1 cell-free calibration",
        "",
        f"**MODEL_STATUS = `{MODEL_STATUS}`**",
        "",
        "Not a digital twin. No measured `u→flux` curve exists in the",
        "repository or project lab notes. See `AC_SPEC_INVENTORY.md`.",
        "`K_AC`, `n_AC`, and `g0` were not fitted to NARMA NRMSE.",
        "`Jmax_A1` was chosen only to match A0 commanded mass.",
        "",
        f"Module 1: **{'PASS' if passed else 'FAIL'}**",
        "",
        "## Frozen envelope",
        "",
        "| Parameter | Value | Role |",
        "|---|---|---|",
        f"| `g0` | {G0:g} | envelope, not measured |",
        f"| `n_AC` | {N_AC:g} | envelope, not measured |",
        f"| `K_AC` | {K_AC:g} | envelope, not measured |",
        f"| `Jmax_A0` | {JMAX_A0:.6e} molecules/s | frozen HybridDish / Narma10b |",
        f"| `Jmax_A1` | {jmax_a1:.12e} molecules/s | payload match only |",
        f"| `g(0)` | {stats['g_at_0']:.12f} | unit test A |",
        f"| `g(0.5)` | {stats['g_at_0p5']:.12f} | unit test B |",
        f"| `g(∞)` | {stats['g_at_inf']:.12f} | documented saturation |",
        "",
        "Gate:",
        "",
        "```",
        "g(u) = g0 + (1-g0) * u^n / (K_AC^n + u^n)",
        "J(t) = Jmax_A1 * g(u(t))   during the 75 s pulse; 0 after",
        "```",
        "",
        "## Frozen u",
        "",
        f"- file: `input_ahl_narma200.txt` (copied from Narma10b, not regenerated)",
        f"- SHA-256: `{sha256_u(u)}`",
        f"- warmup command 0.5 for {WARMUP_WINDOWS} windows, then 200 NARMA windows",
        "",
        "## Payload match",
        "",
        "| Quantity | Molecules |",
        "|---|---|",
        f"| A0 commanded (warmup+NARMA) | {stats['mass_a0']:.12e} |",
        f"| A1 commanded (warmup+NARMA) | {stats['mass_a1']:.12e} |",
        f"| relative error | {stats['relative_error']:.3e} |",
        f"| A0 warmup Σu | {stats['warmup_u_sum']:.12f} |",
        f"| A0 NARMA Σu | {stats['narma_u_sum']:.12f} |",
        f"| A1 warmup Σg | {stats['warmup_g_sum']:.12f} |",
        f"| A1 NARMA Σg | {stats['narma_g_sum']:.12f} |",
        "",
        "Relative error ≤ 1% is required before living BSim. `K_AC` was",
        "not changed after this match.",
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
            "## E0.2 mass identity (A1 no-bacteria field)",
            "",
            "CENTER / FLOW=0 / NO_FLUX. Leftover warmup mass is the NARMA",
            "window initial. Bacteria are absent.",
            "",
            "| Budget | Molecules | Residual fraction |",
            "|---|---|---|",
            f"| warmup leftover | {field['leftover']:.12e} | {field['warmup_fraction']:.3e} |",
            f"| NARMA injected | {field['narma_injected']:.12e} | |",
            f"| NARMA remaining | {field['narma_remaining']:.12e} | {field['narma_fraction']:.3e} |",
            f"| NARMA decay | {field['narma_decay']:.12e} | |",
            f"| NARMA boundary | {field['narma_boundary']:.12e} | |",
            f"| full-horizon injected | {field['full_injected']:.12e} | {field['full_fraction']:.3e} |",
            f"| min C | {field['min_c']:.6g} µM | |",
            "",
            "No vesicle store. Two-way communication was not started.",
            "",
            "## Living BSim authorization",
            "",
            "Living Java is authorized only if every unit test PASSes.",
            f"This run: **{'PASS — seed-111 A1 driven and A1 Brownian may start' if passed else 'FAIL — stop; do not start living BSim'}**.",
            "",
        ]
    )
    (HERE / "results" / "E4_AC_CALIBRATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    print(f"MODEL_STATUS={MODEL_STATUS}")
    print("Not a digital twin. Envelope g0, n_AC, K_AC frozen; Jmax_A1 from mass match only.")
    u = load_u()
    jmax_a1, stats = choose_jmax_a1(u)
    print(f"g(0)={stats['g_at_0']:.12f} g(0.5)={stats['g_at_0p5']:.12f} g(inf)={stats['g_at_inf']:.12f}")
    print(f"Jmax_A1={jmax_a1:.12e} payload_rel_err={stats['relative_error']:.3e}")
    if stats["relative_error"] > MASS_TOL:
        raise SystemExit("STOP payload cannot be matched at frozen K_AC; do not retune K_AC to NRMSE")

    a0 = A1Injector(JMAX_A0)
    a1 = A1Injector(jmax_a1)
    commands = all_commands(u)
    effective_a1 = a1.gate(commands) * (jmax_a1 / JMAX_A0)
    print("integrating A1 no-bacteria field (E0.2 CENTER/FLOW=0)...")
    field = integrate_field(effective_a1, collect=True)
    tests = unit_tests(u, jmax_a1, stats, field)
    for row in tests:
        print(f"  {row['test']}: {'PASS' if row['pass'] else 'FAIL'}  observed={row['observed']:.6g}")

    series_a0 = pulse_timeseries(u, a0)
    series_a1 = pulse_timeseries(u, a1)
    probes = None
    samples = field["samples"]
    if samples is not None:
        c = samples[:, : rts.N]
        probes = {
            "near": c[:, voxel_index(*NEAR_IJ)],
            "mid": c[:, voxel_index(*MID_IJ)],
            "far": c[:, voxel_index(*FAR_IJ)],
        }
        write_csv(
            HERE / "results" / "a1_probes.csv",
            [
                {
                    "sample": i,
                    "AHL_near_uM": float(probes["near"][i]),
                    "AHL_mid_uM": float(probes["mid"][i]),
                    "AHL_far_uM": float(probes["far"][i]),
                }
                for i in range(len(probes["near"]))
            ],
        )

    write_csv(
        HERE / "results" / "a1_payload_match.csv",
        [
            {
                "MODEL_STATUS": MODEL_STATUS,
                "g0": G0,
                "n_AC": N_AC,
                "K_AC": K_AC,
                "Jmax_A0_molecules_per_s": JMAX_A0,
                "Jmax_A1_molecules_per_s": jmax_a1,
                "g_at_0": stats["g_at_0"],
                "g_at_0p5": stats["g_at_0p5"],
                "g_at_inf": stats["g_at_inf"],
                "A0_commanded_mass_molecules": stats["mass_a0"],
                "A1_commanded_mass_molecules": stats["mass_a1"],
                "relative_error": stats["relative_error"],
                "u_sha256": sha256_u(u),
                "digital_twin": False,
            }
        ],
    )
    write_csv(HERE / "results" / "a1_unit_tests.csv", tests)
    write_csv(
        HERE / "results" / "a1_mass_identity.csv",
        [
            {
                "narma_initial_leftover_molecules": field["leftover"],
                "narma_injected_molecules": field["narma_injected"],
                "narma_remaining_molecules": field["narma_remaining"],
                "narma_decay_molecules": field["narma_decay"],
                "narma_boundary_molecules": field["narma_boundary"],
                "narma_residual_fraction": field["narma_fraction"],
                "full_injected_molecules": field["full_injected"],
                "full_residual_fraction": field["full_fraction"],
                "min_C_uM": field["min_c"],
            }
        ],
    )
    u_grid = np.linspace(0.0, 0.5, 201)
    write_csv(
        HERE / "results" / "a1_gate_curve.csv",
        [
            {"u": float(x), "g_A0": float(x), "g_A1": float(a1.gate(x)), "J_A0": float(a0.flux(x)), "J_A1": float(a1.flux(x))}
            for x in u_grid
        ],
    )
    # Compact pulse table: one row per window (warmup + NARMA), not 1 s samples.
    write_csv(
        HERE / "results" / "a1_window_flux.csv",
        [
            {
                "window": i - WARMUP_WINDOWS,
                "phase": "warmup" if i < WARMUP_WINDOWS else "narma",
                "u": float(commands[i]),
                "g_A1": float(a1.gate(commands[i])),
                "J_A0_molecules_per_s": float(a0.flux(commands[i])),
                "J_A1_molecules_per_s": float(a1.flux(commands[i])),
                "mass_A0_window_molecules": float(a0.flux(commands[i]) * PULSE_S),
                "mass_A1_window_molecules": float(a1.flux(commands[i]) * PULSE_S),
            }
            for i in range(len(commands))
        ],
    )
    frozen = {
        "MODEL_STATUS": MODEL_STATUS,
        "digital_twin": False,
        "g0": G0,
        "n_AC": N_AC,
        "K_AC": K_AC,
        "Jmax_A0": JMAX_A0,
        "Jmax_A1": jmax_a1,
        "g_at_0": stats["g_at_0"],
        "g_at_0p5": stats["g_at_0p5"],
        "u_sha256": sha256_u(u),
        "payload_relative_error": stats["relative_error"],
        "module1_pass": all(row["pass"] for row in tests),
    }
    (HERE / "results" / "a1_frozen.json").write_text(json.dumps(frozen, indent=2) + "\n", encoding="utf-8")
    plot_cell_free(u, jmax_a1, series_a0, series_a1, probes)
    write_calibration_md(u, jmax_a1, stats, tests, field)
    if not all(row["pass"] for row in tests):
        raise SystemExit("STOP Module 1 FAIL; do not start living Java")
    print("Module 1 PASS")


if __name__ == "__main__":
    main()
