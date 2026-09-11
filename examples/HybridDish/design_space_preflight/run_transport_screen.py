#!/usr/bin/env python3
"""Analysis-only HybridDish transport/kinetic preflight.

This script never invokes Java or BSim. It reads frozen voxel traces and
writes the E0.2 validation/screening products beside this file.
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.integrate import solve_ivp


HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parents[1]
NARMA = EXAMPLES / "BSimReservoirPlanNarma10b"

NX, NY = 50, 25
RX, RY = 20, 10
LX, LY, DEPTH = 1000.0, 500.0, 10.0
DX, DY = LX / NX, LY / NY
VOXEL_VOLUME = DX * DY * DEPTH
MOLECULES_PER_UM3_PER_UM = 602.2
N = NX * NY
K_HILL, TAU_R, TAU_L = 1.6, 15.0, 1500.0
SOURCE_RATE = 128_000_000.0
SOURCE_CONC_RATE = SOURCE_RATE / (MOLECULES_PER_UM3_PER_UM * VOXEL_VOLUME)
WARMUP_S, WINDOW_S, PULSE_S = 18_000.0, 300.0, 75.0
SAMPLE_OFFSETS = np.array([*np.arange(0.05, 280.06, 20.0), 300.0])
PRIMARY_RTOL, PRIMARY_ATOL, PRIMARY_MAX_STEP = 1e-7, 1e-10, 2.0
TIGHT_RTOL, TIGHT_ATOL, TIGHT_MAX_STEP = 1e-9, 1e-12, 0.5
ROBIN_H = 0.159

POSITIONS = (
    ("CENTER", 500.0, 250.0),
    ("UPSTREAM_CENTER", 200.0, 250.0),
    ("DOWNSTREAM_CENTER", 800.0, 250.0),
    ("UPSTREAM_OFFAXIS", 200.0, 150.0),
    ("CROSSSTREAM_OFFAXIS", 500.0, 350.0),
)
FLOWS = (0.0, 0.25, 0.5, 0.666667, 1.0, 2.0, 3.333333, 8.0)

RANDOM_DRIVE = (
    0.956808, 0.643339, 0.339706, 0.685943, 0.512018, 0.550269, 0.695872, 0.905832,
    0.643954, 0.136880, 0.869236, 0.147003, 0.900775, 0.448573, 0.077669, 0.642089,
    0.771637, 0.962792, 0.490324, 0.100637, 0.293225, 0.010041, 0.700897, 0.085965,
    0.073583, 0.745776, 0.256190, 0.066602, 0.580787, 0.173086, 0.039567, 0.527662,
    0.611706, 0.933062, 0.401168, 0.772280, 0.484610, 0.291416, 0.560397, 0.145653,
    0.898381, 0.487822, 0.333184, 0.418505, 0.136513, 0.144141, 0.885746, 0.682077,
    0.367768, 0.583934, 0.260984, 0.026392, 0.439821, 0.782357, 0.641776, 0.698631,
    0.328100, 0.625604, 0.053712, 0.777457, 0.908097, 0.539165, 0.945477, 0.703801,
)
RANDOM_SHA256 = "5919f8fe318308163b1d798967fd8f89e99ca5eaa268af0b8a62e30e0b67307f"

DRIVES = {
    "SINGLE_PULSE": (1.0,) + (0.0,) * 11,
    "STEP_ON_OFF": (0.0, 0.0) + (1.0,) * 4 + (0.0,) * 6,
    "ALTERNATING_BINARY": (1.0, 0.0) * 6,
    "FIXED_RANDOM": RANDOM_DRIVE,
}

CONTROLS = (
    "driven living cells; silent cells; Brownian density; field-only; "
    "unmasked kinetic surrogate; occupancy-masked surrogate; direct-input baseline"
)


@dataclass(frozen=True)
class Condition:
    condition_id: str
    family: str
    position_name: str
    x_um: float
    y_um: float
    flow_um_s: float
    diffusion_um2_s: float
    decay_s_1: float
    boundary: str
    range_status: str
    notes: str


def canonical_random_bytes() -> bytes:
    return ("".join(f"{value:.6f}\n" for value in RANDOM_DRIVE)).encode("ascii")


def verify_random_drive() -> None:
    got = hashlib.sha256(canonical_random_bytes()).hexdigest()
    if got != RANDOM_SHA256:
        raise RuntimeError(f"fixed random drive hash mismatch: {got}")


def make_conditions() -> list[Condition]:
    conditions: list[Condition] = []
    for name, x, y in POSITIONS:
        for flow in FLOWS:
            token = str(flow).replace(".", "p")
            conditions.append(
                Condition(
                    f"PRI_{name}_F{token}",
                    "PRIMARY_FACTORIAL",
                    name,
                    x,
                    y,
                    flow,
                    159.0,
                    0.0033,
                    "NO_FLUX" if flow == 0 else "OUTFLOW",
                    "IMPLEMENTED_OR_DERIVED_FLOW",
                    "Fixed-payload one-source position by predeclared flow.",
                )
            )
    for flow in (0.0, 0.666667):
        for boundary in ("ABSORBING", "ROBIN_LEAKY"):
            token = str(flow).replace(".", "p")
            conditions.append(
                Condition(
                    f"BND_CENTER_F{token}_{boundary}",
                    "BOUNDARY_SENSITIVITY",
                    "CENTER",
                    500.0,
                    250.0,
                    flow,
                    159.0,
                    0.0033,
                    boundary,
                    "HYPOTHETICAL_DESIGN_ENVELOPE",
                    "Boundary sensitivity; not a measured permeability.",
                )
            )
    for condition_id, diffusion, decay, note in (
        ("HYP_CENTER_D80", 80.0, 0.0033, "Hypothetical low-D one-factor arm."),
        ("HYP_CENTER_D240", 240.0, 0.0033, "Hypothetical high-D one-factor arm."),
        ("HYP_CENTER_K0p00165", 159.0, 0.00165, "Hypothetical half-decay one-factor arm."),
        ("HYP_CENTER_K0p0066", 159.0, 0.0066, "Hypothetical double-decay one-factor arm."),
    ):
        conditions.append(
            Condition(
                condition_id,
                "HYPOTHETICAL_TRANSPORT",
                "CENTER",
                500.0,
                250.0,
                0.0,
                diffusion,
                decay,
                "NO_FLUX",
                "HYPOTHETICAL_DESIGN_ENVELOPE",
                note,
            )
        )
    return conditions


def write_conditions(conditions: Iterable[Condition]) -> None:
    fields = [
        "ConditionID", "Family", "SourcePosition", "SourceX_um", "SourceY_um",
        "Flow_um_s", "TransitTime_s", "Diffusion_um2_s", "Decay_s_1",
        "DecayTime_s", "ChemicalBoundary", "PayloadRule", "RangeStatus", "Notes",
    ]
    with (HERE / "transport_screen_conditions.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for c in conditions:
            writer.writerow(
                {
                    "ConditionID": c.condition_id,
                    "Family": c.family,
                    "SourcePosition": c.position_name,
                    "SourceX_um": f"{c.x_um:.6g}",
                    "SourceY_um": f"{c.y_um:.6g}",
                    "Flow_um_s": f"{c.flow_um_s:.6g}",
                    "TransitTime_s": "inf" if c.flow_um_s == 0 else f"{LX / c.flow_um_s:.9g}",
                    "Diffusion_um2_s": f"{c.diffusion_um2_s:.9g}",
                    "Decay_s_1": f"{c.decay_s_1:.9g}",
                    "DecayTime_s": f"{1 / c.decay_s_1:.9g}",
                    "ChemicalBoundary": c.boundary,
                    "PayloadRule": "fixed source rate x command x 75 s",
                    "RangeStatus": c.range_status,
                    "Notes": c.notes,
                }
            )


def source_index(condition: Condition) -> tuple[int, int]:
    return min(NX - 1, int(condition.x_um / DX)), min(NY - 1, int(condition.y_um / DY))


def initial_state() -> np.ndarray:
    return np.zeros(3 * N + 2, dtype=float)


def state_mass_molecules(state: np.ndarray) -> float:
    return float(np.sum(state[:N]) * VOXEL_VOLUME * MOLECULES_PER_UM3_PER_UM)


def make_rhs(condition: Condition, amplitude: float):
    diffusion = condition.diffusion_um2_s
    decay = condition.decay_s_1
    flow = condition.flow_um_s
    boundary = condition.boundary
    sx, sy = source_index(condition)
    face_x = DY * DEPTH
    face_y = DX * DEPTH

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        c = state[:N].reshape(NX, NY)
        r = state[N : 2 * N].reshape(NX, NY)
        lum = state[2 * N : 3 * N].reshape(NX, NY)
        dc = np.zeros_like(c)

        delta_x = c[1:, :] - c[:-1, :]
        pair_x = diffusion * delta_x / (DX * DX)
        dc[:-1, :] += pair_x
        dc[1:, :] -= pair_x
        delta_y = c[:, 1:] - c[:, :-1]
        pair_y = diffusion * delta_y / (DY * DY)
        dc[:, :-1] += pair_y
        dc[:, 1:] -= pair_y

        boundary_loss_um3_um_s = 0.0
        if boundary == "ABSORBING":
            dc[0, :] -= 2 * diffusion * c[0, :] / (DX * DX)
            dc[-1, :] -= 2 * diffusion * c[-1, :] / (DX * DX)
            dc[:, 0] -= 2 * diffusion * c[:, 0] / (DY * DY)
            dc[:, -1] -= 2 * diffusion * c[:, -1] / (DY * DY)
            boundary_loss_um3_um_s += float(
                np.sum(2 * diffusion * c[0, :] / DX) * face_x
                + np.sum(2 * diffusion * c[-1, :] / DX) * face_x
                + np.sum(2 * diffusion * c[:, 0] / DY) * face_y
                + np.sum(2 * diffusion * c[:, -1] / DY) * face_y
            )
        elif boundary == "ROBIN_LEAKY":
            dc[0, :] -= ROBIN_H * c[0, :] / DX
            dc[-1, :] -= ROBIN_H * c[-1, :] / DX
            dc[:, 0] -= ROBIN_H * c[:, 0] / DY
            dc[:, -1] -= ROBIN_H * c[:, -1] / DY
            boundary_loss_um3_um_s += float(
                ROBIN_H
                * (
                    (np.sum(c[0, :]) + np.sum(c[-1, :])) * face_x
                    + (np.sum(c[:, 0]) + np.sum(c[:, -1])) * face_y
                )
            )

        if flow > 0:
            dc[0, :] -= flow * c[0, :] / DX
            dc[1:, :] -= flow * (c[1:, :] - c[:-1, :]) / DX
            boundary_loss_um3_um_s += float(flow * np.sum(c[-1, :]) * face_x)

        dc -= decay * c
        dc[sx, sy] += SOURCE_CONC_RATE * amplitude
        h = c * c / (K_HILL * K_HILL + c * c)
        dr = (h - r) / TAU_R
        dl = (r - lum) / TAU_L

        derivative = np.empty_like(state)
        derivative[:N] = dc.ravel()
        derivative[N : 2 * N] = dr.ravel()
        derivative[2 * N : 3 * N] = dl.ravel()
        derivative[-2] = decay * float(np.sum(c)) * VOXEL_VOLUME * MOLECULES_PER_UM3_PER_UM
        derivative[-1] = boundary_loss_um3_um_s * MOLECULES_PER_UM3_PER_UM
        return derivative

    return rhs


def integrate_segment(
    state: np.ndarray,
    condition: Condition,
    duration: float,
    amplitude: float,
    sample_times: np.ndarray | None = None,
    *,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
    max_step: float = PRIMARY_MAX_STEP,
) -> tuple[np.ndarray, np.ndarray | None]:
    if duration <= 0:
        return state.copy(), None
    sol = solve_ivp(
        make_rhs(condition, amplitude),
        (0.0, duration),
        state,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        dense_output=sample_times is not None and len(sample_times) > 0,
    )
    if not sol.success:
        raise RuntimeError(f"{condition.condition_id}: solver failed: {sol.message}")
    sampled = None
    if sample_times is not None and len(sample_times) > 0:
        sampled = sol.sol(sample_times).T
    return sol.y[:, -1].copy(), sampled


def simulate_window(
    state: np.ndarray,
    condition: Condition,
    amplitude: float,
    *,
    collect: bool,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
    max_step: float = PRIMARY_MAX_STEP,
) -> tuple[np.ndarray, np.ndarray | None]:
    pulse_samples = SAMPLE_OFFSETS[SAMPLE_OFFSETS <= PULSE_S] if collect else None
    rest_offsets = SAMPLE_OFFSETS[SAMPLE_OFFSETS > PULSE_S] if collect else None
    state, first = integrate_segment(
        state, condition, PULSE_S, amplitude, pulse_samples,
        rtol=rtol, atol=atol, max_step=max_step,
    )
    shifted = None if rest_offsets is None else rest_offsets - PULSE_S
    state, second = integrate_segment(
        state, condition, WINDOW_S - PULSE_S, 0.0, shifted,
        rtol=rtol, atol=atol, max_step=max_step,
    )
    if not collect:
        return state, None
    parts = [part for part in (first, second) if part is not None]
    return state, np.vstack(parts)


def warm_state(
    condition: Condition,
    *,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
    max_step: float = PRIMARY_MAX_STEP,
) -> np.ndarray:
    state = initial_state()
    for _ in range(int(WARMUP_S / WINDOW_S)):
        state, _ = simulate_window(
            state, condition, 0.5, collect=False, rtol=rtol, atol=atol, max_step=max_step
        )
    return state


def reset_budget(state: np.ndarray) -> np.ndarray:
    result = state.copy()
    result[-2:] = 0.0
    return result


def commanded_injected_mass(amplitudes: Iterable[float]) -> float:
    return SOURCE_RATE * PULSE_S * float(sum(amplitudes))


def mass_budget(
    initial: float,
    injected: float,
    remaining: float,
    decay: float,
    boundary: float,
) -> tuple[float, float]:
    residual = initial + injected - remaining - decay - boundary
    denominator = initial + injected
    fraction = abs(residual) / denominator if denominator > 0 else abs(residual)
    return residual, fraction


def simulate_drive(
    warm: np.ndarray,
    condition: Condition,
    drive: Iterable[float],
    *,
    clear_budget: bool = True,
    collect: bool = True,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
    max_step: float = PRIMARY_MAX_STEP,
) -> tuple[np.ndarray, np.ndarray | None]:
    state = reset_budget(warm) if clear_budget else warm.copy()
    samples: list[np.ndarray] = []
    for amplitude in drive:
        state, window_samples = simulate_window(
            state, condition, amplitude, collect=collect,
            rtol=rtol, atol=atol, max_step=max_step,
        )
        if collect:
            assert window_samples is not None
            samples.append(window_samples)
    return state, (np.vstack(samples) if collect else None)


def readout_indices() -> np.ndarray:
    indices = []
    for x in range(RX):
        for y in range(RY):
            px, py = (x + 0.5) * LX / RX, (y + 0.5) * LY / RY
            indices.append(min(NX - 1, int(px / DX)) * NY + min(NY - 1, int(py / DY)))
    return np.asarray(indices, dtype=int)


READOUT_INDICES = readout_indices()


def readout_state(samples: np.ndarray, offset: int) -> np.ndarray:
    return samples[:, offset + READOUT_INDICES]


def source_relative_probe_indices(condition: Condition) -> tuple[int, int, int]:
    readout_x = np.array([(i + 0.5) * LX / RX for i in range(RX)])
    readout_y = np.array([(i + 0.5) * LY / RY for i in range(RY)])
    y = int(np.argmin(np.abs(readout_y - condition.y_um)))
    targets_x = (
        condition.x_um,
        condition.x_um + 0.5 * (LX - condition.x_um),
        LX - LX / (2 * RX),
    )
    result = []
    for target_x in targets_x:
        x = int(np.argmin(np.abs(readout_x - target_x)))
        result.append(int(READOUT_INDICES[x * RY + y]))
    return tuple(result)


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def nrmse_range(pred: np.ndarray, obs: np.ndarray) -> float:
    scale = float(np.max(obs) - np.min(obs))
    return float(np.sqrt(np.mean((pred - obs) ** 2)) / scale) if scale > 0 else float("inf")


def relative_error(pred: float, obs: float) -> float:
    return abs(pred - obs) / abs(obs) if obs != 0 else float("inf")


def lag1(values: np.ndarray) -> float:
    return safe_corr(values[:-1], values[1:]) if len(values) > 2 else float("nan")


def memory_diagnostic(features: np.ndarray, drive: np.ndarray) -> tuple[float, list[float]]:
    scores: list[float] = []
    for delay in range(1, 11):
        rows = np.arange(16, 64)
        valid = rows[rows >= delay]
        train_rows = valid[valid <= 47]
        test_rows = valid[valid >= 48]
        x_train, x_test = features[train_rows], features[test_rows]
        y_train, y_test = drive[train_rows - delay], drive[test_rows - delay]
        mean = x_train.mean(axis=0)
        std = x_train.std(axis=0)
        std[std == 0] = 1.0
        x_train = (x_train - mean) / std
        x_test = (x_test - mean) / std
        x_train = np.column_stack([np.ones(len(x_train)), x_train])
        x_test = np.column_stack([np.ones(len(x_test)), x_test])
        beta = np.linalg.lstsq(x_train, y_train, rcond=None)[0]
        prediction = x_test @ beta
        denominator = float(np.sum((y_test - y_test.mean()) ** 2))
        score = 1.0 - float(np.sum((y_test - prediction) ** 2)) / denominator if denominator > 0 else float("nan")
        scores.append(score)
    return float(sum(max(0.0, score) for score in scores if np.isfinite(score))), scores


def fmt_trace(values: np.ndarray) -> str:
    return "|".join(f"{value:.8g}" for value in values)


def drive_metrics(
    condition: Condition,
    drive_name: str,
    drive: tuple[float, ...],
    initial_mass: float,
    final_state: np.ndarray,
    samples: np.ndarray,
) -> dict[str, object]:
    c = samples[:, :N]
    r = samples[:, N : 2 * N]
    lum = samples[:, 2 * N : 3 * N]
    h = c * c / (K_HILL * K_HILL + c * c)
    wall_mask = np.zeros((NX, NY), dtype=bool)
    wall_mask[[0, -1], :] = True
    wall_mask[:, [0, -1]] = True
    wall = wall_mask.ravel()
    wall_mean = float(np.mean(c[:, wall]))
    interior_mean = float(np.mean(c[:, ~wall]))
    wall_ratio = wall_mean / interior_mean if interior_mean > 0 else float("inf")
    wall_mass_fraction = float(np.sum(c[:, wall]) / np.sum(c)) if np.sum(c) > 0 else 0.0
    probe_ids = source_relative_probe_indices(condition)
    per_window = len(SAMPLE_OFFSETS)
    end_rows = np.arange(per_window - 1, len(samples), per_window)
    mean_c_end = np.mean(c[end_rows], axis=1)
    mean_h_end = np.mean(h[end_rows], axis=1)
    mean_r_end = np.mean(r[end_rows], axis=1)
    mean_l_end = np.mean(lum[end_rows], axis=1)
    near_end, middle_end, far_end = (c[end_rows, probe] for probe in probe_ids)
    injected = commanded_injected_mass(drive)
    remaining = state_mass_molecules(final_state)
    decay_loss = float(final_state[-2])
    boundary_loss = float(final_state[-1])
    residual, mass_rel = mass_budget(initial_mass, injected, remaining, decay_loss, boundary_loss)
    memory_sum, memory_scores = (float("nan"), [])
    if drive_name == "FIXED_RANDOM":
        features = np.column_stack(
            [mean_c_end, mean_h_end, mean_r_end, mean_l_end, near_end, middle_end, far_end]
        )
        memory_sum, memory_scores = memory_diagnostic(features, np.asarray(drive))
    min_c = float(np.min(c))
    numerical_failure = min_c < -1e-9 or mass_rel > 0.01 or not np.all(np.isfinite(samples))
    return {
        "ConditionID": condition.condition_id,
        "DecisionGrade": "PENDING_VALIDATION",
        "Drive": drive_name,
        "Windows": len(drive),
        "InitialMass_molecules": initial_mass,
        "InjectedMass_molecules": injected,
        "RemainingMass_molecules": remaining,
        "DecayLoss_molecules": decay_loss,
        "BoundaryLoss_molecules": boundary_loss,
        "MassResidual_molecules": residual,
        "MassResidualFraction": mass_rel,
        "MinAHL_uM": min_c,
        "PeakAHL_uM": float(np.max(c)),
        "MeanAHL_uM": float(np.mean(c)),
        "PlumeCoverage_H_ge_0p10": float(np.mean(h >= 0.10)),
        "WallInteriorRatio": wall_ratio,
        "WallMassFraction": wall_mass_fraction,
        "MeanH": float(np.mean(h)),
        "VarH": float(np.var(h)),
        "MeanR": float(np.mean(r)),
        "VarR": float(np.var(r)),
        "MeanL": float(np.mean(lum)),
        "VarL": float(np.var(lum)),
        "WeakFraction_H_lt_0p05": float(np.mean(h < 0.05)),
        "SaturatedFraction_H_gt_0p90": float(np.mean(h > 0.90)),
        "Lag1MeanC": lag1(mean_c_end),
        "Lag1MeanR": lag1(mean_r_end),
        "Lag1MeanL": lag1(mean_l_end),
        "NearProbeWindowEnd_uM": fmt_trace(near_end),
        "MiddleProbeWindowEnd_uM": fmt_trace(middle_end),
        "FarProbeWindowEnd_uM": fmt_trace(far_end),
        "MemoryR2Delay1to10": fmt_trace(np.asarray(memory_scores)) if memory_scores else "",
        "LinearMemorySumPositiveR2": memory_sum,
        "NumericalFailure": numerical_failure,
    }


def first_crossing(times: np.ndarray, residual: np.ndarray, threshold: float) -> float:
    peak_index = int(np.argmax(residual))
    candidates = np.flatnonzero(residual[peak_index:] <= threshold)
    if len(candidates) == 0:
        return float("nan")
    index = peak_index + int(candidates[0])
    if index == peak_index or index == 0:
        return float(times[index])
    t0, t1 = times[index - 1], times[index]
    y0, y1 = residual[index - 1], residual[index]
    if y1 == y0:
        return float(t1)
    return float(t0 + (threshold - y0) * (t1 - t0) / (y1 - y0))


def simulate_reset(
    warm: np.ndarray,
    condition: Condition,
    *,
    rtol: float = PRIMARY_RTOL,
    atol: float = PRIMARY_ATOL,
    max_step: float = PRIMARY_MAX_STEP,
) -> tuple[list[dict[str, object]], dict[str, float]]:
    state = reset_budget(warm)
    state, _ = integrate_segment(
        state, condition, PULSE_S, 1.0, None, rtol=rtol, atol=atol, max_step=max_step
    )
    times: list[np.ndarray] = []
    norms_c: list[np.ndarray] = []
    norms_r: list[np.ndarray] = []
    norms_l: list[np.ndarray] = []
    end_indices: list[int] = []
    elapsed = 0.0

    def consume(duration: float) -> None:
        nonlocal state, elapsed
        local_times = np.arange(1.0, duration + 0.1, 1.0)
        state, sampled = integrate_segment(
            state, condition, duration, 0.0, local_times,
            rtol=rtol, atol=atol, max_step=max_step,
        )
        assert sampled is not None
        times.append(elapsed + local_times)
        norms_c.append(np.linalg.norm(sampled[:, :N], axis=1))
        norms_r.append(np.linalg.norm(sampled[:, N : 2 * N], axis=1))
        norms_l.append(np.linalg.norm(sampled[:, 2 * N : 3 * N], axis=1))
        elapsed += duration

    consume(WINDOW_S - PULSE_S)
    for _ in range(40):
        consume(WINDOW_S)
        end_indices.append(sum(len(chunk) for chunk in times) - 1)

    all_times = np.concatenate(times)
    raw = {
        "C": np.concatenate(norms_c),
        "R": np.concatenate(norms_r),
        "L": np.concatenate(norms_l),
    }
    residuals = {
        key: values / np.max(values) if np.max(values) > 0 else np.zeros_like(values)
        for key, values in raw.items()
    }
    rows: list[dict[str, object]] = []
    for zero_window, index in enumerate(end_indices, start=1):
        rows.append(
            {
                "ConditionID": condition.condition_id,
                "DecisionGrade": "PENDING_VALIDATION",
                "ZeroWindow": zero_window,
                "TimeSinceShutoff_s": all_times[index],
                "CResidual": residuals["C"][index],
                "RResidual": residuals["R"][index],
                "LResidual": residuals["L"][index],
            }
        )
    summary: dict[str, float] = {}
    for key in ("C", "R", "L"):
        for threshold, label in ((0.37, "37"), (0.14, "14"), (0.05, "5"), (0.01, "1")):
            summary[f"{key}_t{label}_s"] = first_crossing(all_times, residuals[key], threshold)
    summary["CResidual30"] = float(rows[29]["CResidual"])
    summary["RResidual30"] = float(rows[29]["RResidual"])
    summary["LResidual30"] = float(rows[29]["LResidual"])
    return rows, summary


def screen_worker(condition: Condition) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    warm = warm_state(condition)
    initial_mass = state_mass_molecules(warm)
    result_rows: list[dict[str, object]] = []
    for drive_name, drive in DRIVES.items():
        final_state, samples = simulate_drive(warm, condition, drive)
        result_rows.append(drive_metrics(condition, drive_name, drive, initial_mass, final_state, samples))
    reset_rows, reset_summary = simulate_reset(warm, condition)
    for row in result_rows:
        row.update(reset_summary)
    random_row = next(row for row in result_rows if row["Drive"] == "FIXED_RANDOM")
    exclusions: list[str] = []
    if any(bool(row["NumericalFailure"]) for row in result_rows):
        exclusions.append("NUMERICAL_OR_MASS_BALANCE")
    if float(random_row["MeanR"]) < 0.02 or float(random_row["WeakFraction_H_lt_0p05"]) > 0.90:
        exclusions.append("RECEIVER_EFFECTIVELY_DEAD")
    if float(random_row["MeanH"]) > 0.90 or float(random_row["SaturatedFraction_H_gt_0p90"]) > 0.90:
        exclusions.append("RECEIVER_NEARLY_ALWAYS_SATURATED")
    if float(random_row["WallInteriorRatio"]) > 3 and float(random_row["WallMassFraction"]) > 0.5:
        exclusions.append("SEVERE_WALL_HOTSPOT")
    if max(reset_summary["CResidual30"], reset_summary["RResidual30"], reset_summary["LResidual30"]) > 0.01:
        exclusions.append("RESET_INCOMPATIBLE_30_WINDOWS")
    if math.isclose(condition.flow_um_s, 8.0):
        exclusions.append("KNOWN_DEAD_FLOW_ANCHOR")
    exclusion_text = "|".join(exclusions)
    for row in result_rows:
        row["Exclusions"] = exclusion_text
        row["EligibleForPromotion"] = not exclusions
    return result_rows, reset_rows


def load_input(path: Path) -> tuple[float, ...]:
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(float(line))
    return tuple(values)


def load_validation_seed(seed: int) -> dict[str, np.ndarray | str]:
    path = NARMA / "results" / f"narma10b_driven_seed{seed}" / "voxels.csv"
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    base_names = ["Window", "Sample", "TimeInWindow_s"]
    ahl_names = [f"AHL_uM_{i}" for i in range(200)]
    den_names = [f"Den_{i}" for i in range(200)]
    r_names = [f"Receiver_R_{i}" for i in range(200)]
    names = base_names + ahl_names + den_names + r_names
    usecols = [header.index(name) for name in names]
    data = np.loadtxt(path, delimiter=";", skiprows=1, usecols=usecols)
    ahl = data[:, 3:203]
    den = data[:, 203:403]
    receiver = data[:, 403:603]
    digest = hashlib.sha256(np.ascontiguousarray(ahl).tobytes()).hexdigest()
    absolute_time = data[:, 0] * WINDOW_S + data[:, 2] + 0.05
    return {
        "time": absolute_time,
        "ahl": ahl,
        "den": den,
        "receiver": receiver,
        "hash": digest,
    }


def convergence_check(condition: Condition) -> tuple[bool, dict[str, float]]:
    warm_primary = warm_state(condition)
    primary_state, primary_samples = simulate_drive(warm_primary, condition, DRIVES["SINGLE_PULSE"])
    primary_reset, primary_reset_summary = simulate_reset(warm_primary, condition)
    del primary_reset

    warm_tight = warm_state(
        condition, rtol=TIGHT_RTOL, atol=TIGHT_ATOL, max_step=TIGHT_MAX_STEP
    )
    tight_state, tight_samples = simulate_drive(
        warm_tight,
        condition,
        DRIVES["SINGLE_PULSE"],
        rtol=TIGHT_RTOL,
        atol=TIGHT_ATOL,
        max_step=TIGHT_MAX_STEP,
    )
    tight_reset, tight_reset_summary = simulate_reset(
        warm_tight, condition, rtol=TIGHT_RTOL, atol=TIGHT_ATOL, max_step=TIGHT_MAX_STEP
    )
    del tight_reset
    metrics: dict[str, float] = {}
    for name, offset in (("C", 0), ("R", N), ("L", 2 * N)):
        primary = primary_samples[:, offset : offset + N]
        tight = tight_samples[:, offset : offset + N]
        for statistic, function in (("mean", np.mean), ("peak", np.max)):
            p_value, t_value = float(function(primary)), float(function(tight))
            metrics[f"{name}_{statistic}_relative_difference"] = relative_error(p_value, t_value)
    for loss_name, index in (("decay", -2), ("boundary", -1)):
        metrics[f"{loss_name}_relative_difference"] = relative_error(
            float(primary_state[index]), float(tight_state[index])
        ) if float(tight_state[index]) != 0 else abs(float(primary_state[index]))
    reset_differences = []
    for key, tight_value in tight_reset_summary.items():
        if "_t" not in key or not np.isfinite(tight_value):
            continue
        difference = abs(primary_reset_summary[key] - tight_value)
        metrics[f"{key}_absolute_difference_s"] = difference
        reset_differences.append(difference <= max(2.0, 0.01 * tight_value))
    relative_metrics = [
        value for key, value in metrics.items() if key.endswith("relative_difference")
    ]
    passed = all(value <= 0.01 for value in relative_metrics) and all(reset_differences)
    metrics["passed"] = float(passed)
    return passed, metrics


def validate_model() -> tuple[bool, list[dict[str, object]], dict[str, float], str]:
    condition = make_conditions()[0]
    assert condition.condition_id == "PRI_CENTER_F0p0"
    seeds = {seed: load_validation_seed(seed) for seed in (111, 222, 333)}
    hashes = [str(seeds[seed]["hash"]) for seed in seeds]
    fields_identical = len(set(hashes)) == 1
    drive = load_input(NARMA / "input_ahl_narma200.txt")
    warm = warm_state(condition)
    leftover = state_mass_molecules(warm)
    warmup_injected = commanded_injected_mass((0.5,) * int(WARMUP_S / WINDOW_S))
    warmup_decay = float(warm[-2])
    warmup_boundary = float(warm[-1])
    warmup_residual, warmup_fraction = mass_budget(
        0.0, warmup_injected, leftover, warmup_decay, warmup_boundary
    )
    final_state, samples = simulate_drive(warm, condition, drive)
    assert samples is not None
    model_ahl = readout_state(samples, 0)
    model_r = readout_state(samples, N)
    model_mean_ahl = np.mean(model_ahl, axis=1)
    convergence_pass, convergence = convergence_check(condition)
    rows: list[dict[str, object]] = []
    probe_readout_ids = (105, 155, 195)
    narma_injected = commanded_injected_mass(drive)
    narma_remaining = state_mass_molecules(final_state)
    narma_decay = float(final_state[-2])
    narma_boundary = float(final_state[-1])
    narma_residual, narma_fraction = mass_budget(
        leftover, narma_injected, narma_remaining, narma_decay, narma_boundary
    )
    full_injected = warmup_injected + narma_injected
    full_decay = warmup_decay + narma_decay
    full_boundary = warmup_boundary + narma_boundary
    full_residual, full_fraction = mass_budget(
        0.0, full_injected, narma_remaining, full_decay, full_boundary
    )
    mass_ok = narma_fraction <= 0.01 and full_fraction <= 0.01
    all_pass = fields_identical and convergence_pass and mass_ok
    for seed, loaded in seeds.items():
        actual_ahl = np.asarray(loaded["ahl"])
        den = np.asarray(loaded["den"])
        actual_r = np.asarray(loaded["receiver"])
        actual_mean_ahl = np.mean(actual_ahl, axis=1)
        observed_mean_r_trace = np.sum(actual_r * den, axis=1) / np.maximum(1.0, np.sum(den, axis=1))
        predicted_mean_r_trace = np.sum(model_r * den, axis=1) / np.maximum(1.0, np.sum(den, axis=1))
        observed_mean_r = float(np.mean(observed_mean_r_trace))
        predicted_mean_r = float(np.mean(predicted_mean_r_trace))
        time = np.asarray(loaded["time"])
        actual_integral = float(np.trapezoid(actual_mean_ahl, time))
        model_integral = float(np.trapezoid(model_mean_ahl, time))
        row: dict[str, object] = {
            "Seed": seed,
            "AHLTraceSHA256": loaded["hash"],
            "AHLFieldsIdenticalAcrossSeeds": fields_identical,
            "DomainMeanCorrelation": safe_corr(model_mean_ahl, actual_mean_ahl),
            "DomainMeanNRMSE": nrmse_range(model_mean_ahl, actual_mean_ahl),
            "PulseIntegratedAHLObserved_uM_s": actual_integral,
            "PulseIntegratedAHLPredicted_uM_s": model_integral,
            "PulseIntegratedAHLRelativeError": relative_error(model_integral, actual_integral),
            "MeanRObserved": observed_mean_r,
            "MeanRPredicted": predicted_mean_r,
            "MeanRRelativeError": relative_error(predicted_mean_r, observed_mean_r),
            "MassResidualFraction": narma_fraction,
            "NarmaWindowInitialMass_molecules": leftover,
            "NarmaWindowInjectedMass_molecules": narma_injected,
            "NarmaWindowRemainingMass_molecules": narma_remaining,
            "NarmaWindowDecayLoss_molecules": narma_decay,
            "NarmaWindowBoundaryLoss_molecules": narma_boundary,
            "NarmaWindowMassResidual_molecules": narma_residual,
            "FullHorizonInitialMass_molecules": 0.0,
            "FullHorizonInjectedMass_molecules": full_injected,
            "FullHorizonRemainingMass_molecules": narma_remaining,
            "FullHorizonDecayLoss_molecules": full_decay,
            "FullHorizonBoundaryLoss_molecules": full_boundary,
            "FullHorizonMassResidual_molecules": full_residual,
            "FullHorizonMassResidualFraction": full_fraction,
            "WarmupMassResidualFraction": warmup_fraction,
            "ConvergencePass": convergence_pass,
        }
        for label, readout_id in zip(("Near", "Middle", "Far"), probe_readout_ids):
            row[f"{label}Correlation"] = safe_corr(model_ahl[:, readout_id], actual_ahl[:, readout_id])
            row[f"{label}NRMSE"] = nrmse_range(model_ahl[:, readout_id], actual_ahl[:, readout_id])
        row_pass = (
            float(row["DomainMeanCorrelation"]) >= 0.98
            and float(row["DomainMeanNRMSE"]) <= 0.10
            and float(row["PulseIntegratedAHLRelativeError"]) <= 0.10
            and float(row["MeanRRelativeError"]) <= 0.10
            and narma_fraction <= 0.01
            and full_fraction <= 0.01
            and convergence_pass
        )
        row["Pass"] = row_pass
        all_pass = all_pass and row_pass
        rows.append(row)
    status = "VALIDATED_FOR_SCREENING" if all_pass else "NOT_VALIDATED"
    convergence.update(
        {
            "narma_window_initial_mass_molecules": leftover,
            "narma_window_injected_mass_molecules": narma_injected,
            "narma_window_remaining_mass_molecules": narma_remaining,
            "narma_window_decay_loss_molecules": narma_decay,
            "narma_window_boundary_loss_molecules": narma_boundary,
            "narma_window_mass_residual_molecules": narma_residual,
            "narma_window_mass_residual_fraction": narma_fraction,
            "full_horizon_injected_mass_molecules": full_injected,
            "full_horizon_remaining_mass_molecules": narma_remaining,
            "full_horizon_decay_loss_molecules": full_decay,
            "full_horizon_boundary_loss_molecules": full_boundary,
            "full_horizon_mass_residual_molecules": full_residual,
            "full_horizon_mass_residual_fraction": full_fraction,
            "warmup_injected_mass_molecules": warmup_injected,
            "warmup_decay_loss_molecules": warmup_decay,
            "warmup_mass_residual_molecules": warmup_residual,
            "warmup_mass_residual_fraction": warmup_fraction,
        }
    )
    return all_pass, rows, convergence, status


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_value(value: object) -> str:
    if isinstance(value, (float, np.floating)):
        return "NA" if not np.isfinite(value) else f"{float(value):.6g}"
    return str(value)


def write_validation_report(
    status: str,
    rows: list[dict[str, object]],
    convergence: dict[str, float],
) -> None:
    lines = [
        "# Reduced transport-model validation",
        "",
        "Frozen criteria are in `PROTOCOL.md`; no transport parameter was fitted.",
        "",
        f"`TRANSPORT_MODEL_STATUS: {status}`",
        "",
        "The three driven NARMA AHL files are deterministic duplicates when their",
        "reported AHL arrays have the same SHA-256. They are retained as three",
        "provenance rows, not counted as three independent field realizations.",
        "",
        "## Frozen-gate results",
        "",
        "| Seed | AHL r | AHL NRMSE | integral rel. error | mean-R rel. error | mass residual | pass |",
        "|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in rows:
        lines.append(
            "| {Seed} | {DomainMeanCorrelation:.5f} | {DomainMeanNRMSE:.5f} | "
            "{PulseIntegratedAHLRelativeError:.5f} | {MeanRRelativeError:.5f} | "
            "{MassResidualFraction:.3e} | {Pass} |".format(**row)
        )
    lines.extend(
        [
            "",
            "Mean R is compared after weighting the fixed-grid surrogate R by each",
            "stored sample's cell density. This tests the kinetic approximation without",
            "silently converting empty voxels or bacterial removal into chemical loss.",
            "",
            "## Representative probes",
            "",
            "| Seed | near r / NRMSE | middle r / NRMSE | far r / NRMSE |",
            "|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['Seed']} | {row['NearCorrelation']:.5f} / {row['NearNRMSE']:.5f} | "
            f"{row['MiddleCorrelation']:.5f} / {row['MiddleNRMSE']:.5f} | "
            f"{row['FarCorrelation']:.5f} / {row['FarNRMSE']:.5f} |"
        )
    lines.extend(
        [
            "",
            "Probe coordinates are source-near `(525,275)`, middle `(775,275)`,",
            "and far `(975,275) um` on the stored 20x10 readout.",
            "",
            "Mass residual in the table is the NARMA-window identity:",
            "leftover warmup mass as `remaining(t0)`, then NARMA injection,",
            "decay, and boundary only. A separate full-horizon identity starts",
            "from `C=0` before warmup. Historical 2026-08-17 rows that mixed",
            "those windows (`mass residual 0.378`) are retained in",
            "`transport_validation_failed_2026-08-17.csv` and",
            "`TRANSPORT_MODEL_VALIDATION_failed_2026-08-17.md`.",
            "",
            "## Numerical convergence and mass balance",
            "",
        ]
    )
    for key in sorted(convergence):
        lines.append(f"- `{key}`: {markdown_value(convergence[key])}")
    if status == "NOT_VALIDATED":
        lines.extend(
            [
                "",
                "Ranking and shortlisting are prohibited. Any sensitivity CSV is raw,",
                "non-decision-grade diagnostic output. The next builder must repair the",
                "surrogate/observation mapping and repeat this frozen validation without",
                "fitting D, decay, source strength, or boundary coefficients to pass.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "All frozen gates passed. The surrogate is authorized only for this",
                "physical screening role; it is not living-cell or BSim evidence.",
            ]
        )
    (HERE / "TRANSPORT_MODEL_VALIDATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def classify_regimes(random_rows: list[dict[str, object]]) -> dict[str, str]:
    eligible = [row for row in random_rows if bool(row["EligibleForPromotion"])]
    if not eligible:
        return {}
    values = np.asarray([float(row["MeanR"]) for row in eligible])
    low, high = np.quantile(values, [1 / 3, 2 / 3])
    regimes = {}
    for row in random_rows:
        value = float(row["MeanR"])
        regimes[str(row["ConditionID"])] = (
            "WEAK" if value <= low else "INTERMEDIATE" if value <= high else "STRONG"
        )
    return regimes


def build_shortlist(
    validation_pass: bool,
    conditions: list[Condition],
    screen_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    if not validation_pass:
        return []
    by_id = {c.condition_id: c for c in conditions}
    random_rows = [row for row in screen_rows if row["Drive"] == "FIXED_RANDOM"]
    random_by_id = {str(row["ConditionID"]): row for row in random_rows}
    regimes = classify_regimes(random_rows)
    selected: list[str] = []
    claim_id = "PRI_CENTER_F0p0"
    if bool(random_by_id[claim_id]["EligibleForPromotion"]):
        selected.append(claim_id)

    target_flows = (0.0, 0.666667, 3.333333)
    target_positions = ("CENTER", "UPSTREAM_CENTER", "DOWNSTREAM_CENTER")
    for position in target_positions:
        for flow in target_flows:
            candidates = [
                c for c in conditions
                if c.family == "PRIMARY_FACTORIAL"
                and c.position_name == position
                and math.isclose(c.flow_um_s, flow)
                and bool(random_by_id[c.condition_id]["EligibleForPromotion"])
            ]
            if candidates and candidates[0].condition_id not in selected:
                selected.append(candidates[0].condition_id)

    eligible_remaining = [
        c for c in conditions
        if c.family == "PRIMARY_FACTORIAL"
        and c.flow_um_s < 8
        and bool(random_by_id[c.condition_id]["EligibleForPromotion"])
        and c.condition_id not in selected
    ]
    for desired in ("WEAK", "INTERMEDIATE", "STRONG"):
        for c in eligible_remaining:
            if regimes.get(c.condition_id) == desired and c.condition_id not in selected:
                selected.append(c.condition_id)
                break
        if len(selected) >= 11:
            break
    selected = selected[:11]

    rows: list[dict[str, object]] = []
    for condition_id in selected:
        c = by_id[condition_id]
        rows.append(
            {
                "ConditionID": condition_id,
                "SourcePosition": f"{c.position_name} ({c.x_um:g},{c.y_um:g},5) um",
                "Flow": f"{c.flow_um_s:g} um/s",
                "Diffusion": f"{c.diffusion_um2_s:g} um^2/s",
                "Decay": f"{c.decay_s_1:g} 1/s",
                "ChemicalBoundary": c.boundary,
                "PayloadRule": "PRIMARY_FIXED_COMMANDED_AHL_PAYLOAD",
                "ReasonIncluded": (
                    "Frozen claim condition" if condition_id == claim_id
                    else "Predeclared layout-by-flow factorial survivor"
                ),
                "ExpectedTransportRegime": regimes.get(condition_id, "REFERENCE"),
                "CalibrationStatus": (
                    "FROZEN_CLAIM_REFERENCE" if condition_id == claim_id
                    else "EXPLORATORY_ENGINEERING_NOT_BIOLOGICALLY_CALIBRATED"
                ),
                "RequiredControls": CONTROLS,
            }
        )
    anchor = next(
        c for c in conditions
        if c.family == "PRIMARY_FACTORIAL" and c.position_name == "CENTER" and c.flow_um_s == 8.0
    )
    rows.append(
        {
            "ConditionID": "HISTORICAL_CENTER_F8_ANCHOR",
            "SourcePosition": f"CENTER ({anchor.x_um:g},{anchor.y_um:g},5) um",
            "Flow": "8 um/s",
            "Diffusion": "159 um^2/s",
            "Decay": "0.0033 1/s",
            "ChemicalBoundary": "OUTFLOW",
            "PayloadRule": "HISTORICAL_FIXED_PAYLOAD_ANCHOR_NO_NEW_RUN_REQUIRED",
            "ReasonIncluded": "Known occupancy-destroying flow anchor; excluded from promotion",
            "ExpectedTransportRegime": "DEAD_WASHOUT_ANCHOR",
            "CalibrationStatus": "HISTORICAL_ANCHOR_NOT_NEW_LIVING_RUN",
            "RequiredControls": "Existing historical evidence; rerun only if implementation comparability fails",
        }
    )
    return rows[:12]


def write_shortlist(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "ConditionID", "SourcePosition", "Flow", "Diffusion", "Decay",
        "ChemicalBoundary", "PayloadRule", "ReasonIncluded",
        "ExpectedTransportRegime", "CalibrationStatus", "RequiredControls",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_builder_prompt(status: str, shortlist: list[dict[str, object]]) -> None:
    lines = ["# E0.3 frozen builder prompt", ""]
    if status != "VALIDATED_FOR_SCREENING":
        lines.extend(
            [
                "Status: repair/validation only.",
                "",
                "`E0.3_BSIM_AUTHORIZATION: WITHHELD`",
                "",
                "Do not run an E0.3 living factorial. Repair the reduced surrogate or",
                "its observation mapping, rerun the already frozen validation criteria,",
                "and diagnose the discrepancy without fitting D, decay, source strength,",
                "or boundary coefficients to pass. Retain all failed outputs and hashes.",
                "Only a subsequent frozen prompt after `VALIDATED_FOR_SCREENING` may",
                "authorize living BSim.",
            ]
        )
    else:
        lines.extend(
            [
                "Status: development-only E0.3 request. This file does not execute BSim.",
                "",
                "`E0.3_BSIM_AUTHORIZATION: REQUESTED_FOR_SEPARATE_BUILDER`",
                "",
                "Build and run only the exact table below. The 8 um/s row is a",
                "historical anchor and does not require a new run unless the builder",
                "proves the old advection implementation is not comparable.",
                "",
                "| ConditionID | Source | Flow | D | decay | boundary | payload |",
                "|---|---|---|---|---|---|---|",
            ]
        )
        for row in shortlist:
            lines.append(
                f"| {row['ConditionID']} | {row['SourcePosition']} | {row['Flow']} | "
                f"{row['Diffusion']} | {row['Decay']} | {row['ChemicalBoundary']} | "
                f"{row['PayloadRule']} |"
            )
        lines.extend(
            [
                "",
                "## Frozen execution requirements",
                "",
                "- Development-only: no result is confirmatory and no historical gate changes.",
                "- Keep the active domain, field/readout grids, 18,000 s warmup, 300 s",
                "  windows, 75 s pulse, K=1.6 uM, n=2, tau_R=15 s, tau_L=1500 s,",
                "  growth, clamp, acid source, and mortality unchanged.",
                "- Primary comparisons hold total commanded AHL payload fixed. Any",
                "  matched-occupancy arm must be separately named and may not replace it.",
                "- For every promoted living condition run: driven living cells, silent",
                "  cells, Brownian density, field-only, unmasked kinetic surrogate,",
                "  occupancy-masked surrogate, and direct-input baseline.",
                "- Smoke first: exact printed geometry/source/flow/boundary, source inside",
                "  bounds, finite non-negative fields, expected row counts, and mass",
                "  accounting within 1%. A smoke pass does not count as evidence.",
                "- Sanity: reproduce the frozen claim condition before interpreting new",
                "  conditions. Preserve the 8 um/s occupancy-dead behavior as an anchor.",
                "- Occupancy gates: exclude mean R<0.02 or >90% weak samples; exclude",
                "  mean H>0.90 or >90% saturated samples.",
                "- Viability gates: report population, births, clamp/acid/OOB deaths and",
                "  spatial occupancy; do not tune clamp or mortality to rescue a condition.",
                "- Report plume coverage, wall/interior ratio, boundary loss, reset",
                "  residuals, population retention, memory diagnostics, and complete",
                "  task/control scores. Do not select on NRMSE/AUC alone.",
                "- Hash Java, configs, inputs, condition manifest, checker, and every",
                "  result manifest before analysis. Record the random/task hashes.",
                "- No best-seed selection. Use all predeclared seeds and retain every",
                "  completed, dead, failed, and excluded result.",
                "- Keep the complete development surface. Freeze at most a small",
                "  predeclared number of designs for untouched confirmation on new",
                "  inputs/device or cell draws; do not reuse development data as",
                "  confirmation.",
                f"- Random-drive SHA-256 used by the screen (not a task hash): `{RANDOM_SHA256}`.",
                "- Do not execute this prompt from the preflight package. A later builder",
                "  may run living BSim only from this frozen table.",
            ]
        )
    (HERE / "E0_3_FROZEN_BUILDER_PROMPT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _unit_test_condition(**overrides: object) -> Condition:
    values = {
        "condition_id": "UNIT_TEST",
        "family": "UNIT_TEST",
        "position_name": "CENTER",
        "x_um": 500.0,
        "y_um": 250.0,
        "flow_um_s": 0.0,
        "diffusion_um2_s": 159.0,
        "decay_s_1": 0.0033,
        "boundary": "NO_FLUX",
        "range_status": "UNIT_TEST",
        "notes": "mass-budget unit test",
    }
    values.update(overrides)
    return Condition(**values)  # type: ignore[arg-type]


def _assert_rel(name: str, value: float, ceiling: float) -> None:
    if not np.isfinite(value) or value > ceiling:
        raise RuntimeError(f"mass-budget unit test {name} failed: {value} > {ceiling}")


def run_mass_budget_unit_tests() -> None:
    """Accounting identities A-E. Does not fit D, k, source, or boundaries."""
    condition_a = _unit_test_condition(decay_s_1=0.0)
    state = initial_state()
    field = state[:N].reshape(NX, NY)
    field[20:30, 8:16] = 1.25
    state[:N] = field.ravel()
    mass0 = state_mass_molecules(state)
    state, _ = integrate_segment(state, condition_a, 800.0, 0.0)
    _assert_rel("A", abs(state_mass_molecules(state) - mass0) / mass0, 1e-12)
    _assert_rel("A_decay_counter", abs(float(state[-2])), 1e-6)
    _assert_rel("A_boundary_counter", abs(float(state[-1])), 1e-6)

    k = 0.0033
    duration = 400.0
    condition_b = _unit_test_condition(decay_s_1=k)
    state = initial_state()
    state[:N] = 2.0
    mass0 = state_mass_molecules(state)
    state, _ = integrate_segment(state, condition_b, duration, 0.0)
    expected = mass0 * math.exp(-k * duration)
    _assert_rel("B", abs(state_mass_molecules(state) - expected) / expected, 1e-12)

    condition_c = _unit_test_condition()
    state = initial_state()
    state, _ = integrate_segment(state, condition_c, PULSE_S, 1.0)
    state, _ = integrate_segment(state, condition_c, 5.0 / k, 0.0)
    injected = commanded_injected_mass((1.0,))
    residual, fraction = mass_budget(
        0.0,
        injected,
        state_mass_molecules(state),
        float(state[-2]),
        float(state[-1]),
    )
    _assert_rel("C", fraction, 1e-12)

    condition_d = _unit_test_condition()
    state = initial_state()
    for _ in range(2):
        state, _ = simulate_window(state, condition_d, 0.5, collect=False)
    leftover = state_mass_molecules(state)
    if leftover <= 0:
        raise RuntimeError("mass-budget unit test D failed: leftover warmup mass is not positive")
    drive = (1.0, 0.0, 0.5)
    final, _ = simulate_drive(state, condition_d, drive, clear_budget=True, collect=False)
    injected = commanded_injected_mass(drive)
    residual, fraction = mass_budget(
        leftover,
        injected,
        state_mass_molecules(final),
        float(final[-2]),
        float(final[-1]),
    )
    _assert_rel("D", fraction, 1e-12)
    wrong_injected = commanded_injected_mass((0.5, 0.5) + drive)
    _, wrong_fraction = mass_budget(
        0.0,
        wrong_injected,
        state_mass_molecules(final),
        float(final[-2]),
        float(final[-1]),
    )
    if wrong_fraction <= 0.01:
        raise RuntimeError(
            "mass-budget unit test D failed: mixed-horizon identity unexpectedly closed"
        )

    condition_e = _unit_test_condition()
    warmup = (0.5, 0.5)
    state = initial_state()
    for amplitude in warmup:
        state, _ = simulate_window(state, condition_e, amplitude, collect=False)
    warmup_decay = float(state[-2])
    if warmup_decay <= 0:
        raise RuntimeError("mass-budget unit test E failed: warmup decay was not accumulated")
    final, _ = simulate_drive(state, condition_e, drive, clear_budget=False, collect=False)
    injected = commanded_injected_mass(warmup + drive)
    residual, fraction = mass_budget(
        0.0,
        injected,
        state_mass_molecules(final),
        float(final[-2]),
        float(final[-1]),
    )
    _assert_rel("E", fraction, 1e-12)
    if float(final[-2]) <= warmup_decay:
        raise RuntimeError("mass-budget unit test E failed: full-horizon decay omitted the drive")
    print("mass-budget unit tests A-E passed", flush=True)


def main() -> None:
    verify_random_drive()
    run_mass_budget_unit_tests()
    conditions = make_conditions()
    write_conditions(conditions)

    validation_pass, validation_rows, convergence, status = validate_model()
    write_csv(HERE / "transport_validation.csv", validation_rows)
    write_validation_report(status, validation_rows, convergence)

    workers = min(4, max(1, os.cpu_count() or 1))
    screen_rows: list[dict[str, object]] = []
    reset_rows: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {executor.submit(screen_worker, condition): condition for condition in conditions}
        for future in as_completed(future_map):
            condition = future_map[future]
            result, reset = future.result()
            screen_rows.extend(result)
            reset_rows.extend(reset)
            print(f"completed {condition.condition_id}", flush=True)

    condition_order = {condition.condition_id: i for i, condition in enumerate(conditions)}
    drive_order = {name: i for i, name in enumerate(DRIVES)}
    screen_rows.sort(key=lambda row: (condition_order[str(row["ConditionID"])], drive_order[str(row["Drive"])]))
    reset_rows.sort(key=lambda row: (condition_order[str(row["ConditionID"])], int(row["ZeroWindow"])))
    decision_grade = "DECISION_GRADE_PHYSICAL_SCREEN" if validation_pass else "NON_DECISION_GRADE"
    for row in screen_rows:
        row["DecisionGrade"] = decision_grade
    for row in reset_rows:
        row["DecisionGrade"] = decision_grade
    write_csv(HERE / "transport_screen_results.csv", screen_rows)
    write_csv(HERE / "reset_screen_results.csv", reset_rows)

    shortlist = build_shortlist(validation_pass, conditions, screen_rows)
    write_shortlist(HERE / "E0_3_SHORTLIST.csv", shortlist)
    write_builder_prompt(status, shortlist)
    print(f"TRANSPORT_MODEL_STATUS: {status}")
    print(f"random_drive_sha256={RANDOM_SHA256}")
    print(f"conditions={len(conditions)} shortlist_rows={len(shortlist)}")


if __name__ == "__main__":
    main()

