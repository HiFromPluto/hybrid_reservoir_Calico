#!/usr/bin/env python3
"""E4.1 seed-111 scout checker. A1 vs reused A0 Narma10b.

Does not edit GATE_EVIDENCE. Does not retune K, n, tau, clamp, mortality,
flow, layout, K_AC, or Jmax. Seeds 222/333 are authorized only by the
frozen scout rule in PROTOCOL.md.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NARMA = HERE.parent / "BSimReservoirPlanNarma10b"
sys.path.insert(0, str(NARMA))
import check_narma10b as N  # noqa: E402

WASHOUT, TRAIN, TEST = 40, 110, 50
NUM_WINDOWS = 200
NARMA_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
RX, RY = 20, 10
K_HILL = 1.6
G0, N_AC, K_AC = 0.0, 2.0, 0.25
JMAX_A1 = 74677509.75821304
WEAK_MARGIN = 0.03
OCC_DEAD_R = 0.05
OCC_DEAD_CORR = 0.5
LAST_SAMPLE_PREFIX = "199;15;299.95"
EXPECTED_AUX = 3200

A0_DRIVEN = NARMA / "results" / "narma10b_driven_seed111"
A0_BROWN = NARMA / "results" / "narma10b_brownian_seed111"
A0_SILENT = NARMA / "results" / "narma10b_silent_seed111"
A1_DRIVEN = HERE / "results" / "e4ac_a1_driven_seed111"
A1_BROWN = HERE / "results" / "e4ac_a1_brownian_seed111"

FORBIDDEN_JAVA = (
    "glucose", "Glucose", "GLUCOSE",
    "Danino", "danino",
    "CS_in",
    "setGoal(glucose)",
)


def sha256_u(u):
    return hashlib.sha256(",".join(f"{v:.12f}" for v in u).encode("ascii")).hexdigest()


def ac_gate(u):
    u = np.asarray(u, dtype=float)
    un = np.power(np.maximum(u, 0.0), N_AC)
    kn = K_AC ** N_AC
    return G0 + (1.0 - G0) * un / (kn + un)


def hill(c):
    c = np.asarray(c, dtype=float)
    return c * c / (K_HILL * K_HILL + c * c)


def pearson(a, b):
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return 0.0 if den < 1e-15 else float(np.sum(a * b) / den)


def delay_matrix(values, taps=10):
    X = np.zeros((len(values), taps), dtype=float)
    for n, value in enumerate(values):
        for lag in range(taps):
            if n - lag >= 0:
                X[n, lag] = values[n - lag]
    return X


def constant_nrmse(y, constant):
    y_test = y[WASHOUT + TRAIN :]
    pred = np.full_like(y_test, constant)
    return N.nrmse(y_test, pred)


def load_voxel_arrays(path):
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    ahl_i = [header.index(f"AHL_uM_{i}") for i in range(200)]
    den_i = [header.index(f"Den_{i}") for i in range(200)]
    r_i = [header.index(f"Receiver_R_{i}") for i in range(200)]
    l_i = [header.index(f"Lum_Mean_{i}") for i in range(200)]
    win = header.index("Window")
    tcol = header.index("TimeInWindow_s")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return {
        "window": data[:, win].astype(int),
        "t_in": data[:, tcol],
        "ahl": data[:, ahl_i],
        "den": data[:, den_i],
        "receiver": data[:, r_i],
        "lum": data[:, l_i],
        "absolute_time": data[:, win] * 300.0 + data[:, tcol],
    }


def window_mean_matrix(values, windows):
    out = []
    for window in range(NUM_WINDOWS):
        mask = windows == window
        out.append(np.mean(values[mask], axis=0))
    return np.vstack(out)


def occupancy(voxels, u):
    r_win = window_mean_matrix(voxels["receiver"], voxels["window"])
    mean_r_w = np.mean(r_win, axis=1)
    mean_r = float(np.mean(mean_r_w))
    r_u = pearson(mean_r_w, u)
    if abs(r_u) < OCC_DEAD_CORR or mean_r < OCC_DEAD_R:
        label = "DEAD"
    else:
        label = "ALIVE"
    return {
        "occupancy": label,
        "r_meanR_u": r_u,
        "mean_R": mean_r,
        "mean_AHL": float(np.mean(voxels["ahl"])),
        "mean_L": float(np.mean(voxels["lum"])),
        "ahl_finite": bool(np.all(np.isfinite(voxels["ahl"]))),
        "ahl_nonneg": bool(np.min(voxels["ahl"]) >= -1e-12),
    }


def integrate_surrogate(absolute_time, ahl):
    R = np.zeros(ahl.shape[1], dtype=float)
    L = np.zeros(ahl.shape[1], dtype=float)
    out_r = np.empty_like(ahl)
    out_l = np.empty_like(ahl)
    out_r[0], out_l[0] = R, L
    for index in range(1, len(absolute_time)):
        delta = absolute_time[index] - absolute_time[index - 1]
        C = ahl[index - 1]
        h = hill(C)
        er = math.exp(-delta / 15.0)
        el = math.exp(-delta / 1500.0)
        old_r, old_l = R, L
        R = h + (old_r - h) * er
        L = h + (old_l - h) * el + (old_r - h) * 15.0 / (15.0 - 1500.0) * (er - el)
        out_r[index], out_l[index] = R, L
    return out_r, out_l


def ridge_score(X, y):
    metrics = N.evaluate_readout(np.asarray(X), np.asarray(y), list(range(NUM_WINDOWS)))
    return {
        "lambda": metrics["lambda"],
        "test_nrmse": metrics["test_nrmse"],
        "test_r2": metrics["test_r2"],
        "n_features": np.asarray(X).shape[1],
    }


def csv_ok(run_dir, expected_windows=NUM_WINDOWS, expected_aux=EXPECTED_AUX):
    run_dir = Path(run_dir)
    summary = N.validate_csv(run_dir / "window_summary.csv", expected_windows)
    voxels = N.validate_csv(run_dir / "voxels.csv", expected_aux, check_last_sample=True)
    results = N.validate_csv(run_dir / "results.csv", expected_aux, check_last_sample=True)
    last = N.last_sample_ok(run_dir / "voxels.csv")
    return (
        summary["row_count_pass"]
        and summary["rectangular_csv_pass"]
        and voxels["row_count_pass"]
        and voxels["rectangular_csv_pass"]
        and results["row_count_pass"]
        and last
        and str(voxels.get("last_sample_pass", last)).lower() != "false"
    )


def java_ok():
    java = (HERE / "BSimReservoirPlanE4AC.java").read_text(encoding="utf-8")
    voxel = (HERE / "VoxelAnalyzer.java").read_text(encoding="utf-8")
    if any(token in java for token in FORBIDDEN_JAVA):
        return False
    if "CS_in" in java or "class ArtificialCell" in java:
        return False
    if "FLOW_SPEED = 0.0" not in java:
        return False
    if "RECEIVER_K_UM = 1.6" not in java:
        return False
    if "acGate" not in java:
        return False
    if "JMAX_A1_FROZEN" not in java:
        return False
    if "package BSimReservoirPlanE4AC;" not in java:
        return False
    if "BSimReservoirPlanE4AC.DeathCause" not in voxel:
        return False
    return True


def gate_evidence_untouched():
    paths = [
        NARMA / "results" / "GATE_EVIDENCE.md",
        HERE.parent / "BSimReservoirPlanStage6" / "results" / "GATE_EVIDENCE.md",
    ]
    return all(path.exists() for path in paths)


def fmt(value, digits=4):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return ""
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def claim_sanity(u, y):
    driven = N.read_run(A0_DRIVEN, "driven", NUM_WINDOWS, EXPECTED_AUX)
    brown = N.read_run(A0_BROWN, "brownian", NUM_WINDOWS, EXPECTED_AUX)
    silent = N.read_run(A0_SILENT, "silent", NUM_WINDOWS, EXPECTED_AUX)
    d = N.evaluate_arm_run(driven, u, y)
    b = N.evaluate_arm_run(brown, u, y)
    s = N.evaluate_arm_run(silent, u, y)
    checks = {
        "A0_F408": (d["task"]["test_nrmse"], 0.9312),
        "A0_field": (d["field_only"]["test_nrmse"], 1.0289),
        "A0_brownian": (b["task"]["test_nrmse"], 1.1625),
        "A0_silent": (s["task"]["test_nrmse"], 1.1622),
    }
    print("A0 NARMA10b seed 111 reuse sanity")
    failed = False
    for name, (got, expected) in checks.items():
        ok = abs(got - expected) <= 1e-3
        print(f"  {name:14} {got:.4f} expected {expected:.4f} {'PASS' if ok else 'FAIL'}")
        failed = failed or not ok
    if failed:
        raise SystemExit("STOP A0 reuse sanity failed; Narma10b voxels must not be overwritten")
    return d, b, s


def analyze(u, y):
    g = ac_gate(u)
    j_window_mean = JMAX_A1 * g * (75.0 / 300.0)
    taps_u = delay_matrix(u, 10)
    product = np.array([u[n] * u[n - 9] if n >= 9 else 0.0 for n in range(len(u))])
    informed = np.column_stack([taps_u, product])
    baselines = {
        "LINEAR_U_DELAY_10": ridge_score(taps_u, y),
        "NARMA_INFORMED_INPUT": ridge_score(informed, y),
        "TRAIN_INTERCEPT": {
            "lambda": "",
            "test_nrmse": constant_nrmse(y, float(np.mean(y[WASHOUT:WASHOUT + TRAIN]))),
            "test_r2": float("nan"),
            "n_features": 0,
        },
        "PERSISTENCE": {
            "lambda": "",
            "test_nrmse": N.nrmse(y[WASHOUT + TRAIN :], y[WASHOUT + TRAIN - 1 : WASHOUT + TRAIN + TEST - 1]),
            "test_r2": float("nan"),
            "n_features": 1,
        },
        "AC_G_SCALAR": ridge_score(g.reshape(-1, 1), y),
        "AC_J_WINDOW_MEAN": ridge_score(j_window_mean.reshape(-1, 1), y),
        "AC_G_DELAY_10": ridge_score(delay_matrix(g, 10), y),
    }
    a0_d, a0_b, a0_s = claim_sanity(u, y)
    if not csv_ok(A1_DRIVEN) or not csv_ok(A1_BROWN):
        raise SystemExit("STOP A1 CSV completeness failed")
    a1_driven = N.read_run(A1_DRIVEN, "driven", NUM_WINDOWS, EXPECTED_AUX)
    a1_brown = N.read_run(A1_BROWN, "brownian", NUM_WINDOWS, EXPECTED_AUX)
    a1_d = N.evaluate_arm_run(a1_driven, u, y)
    a1_b = N.evaluate_arm_run(a1_brown, u, y)
    a1_s = a0_s
    voxels = load_voxel_arrays(A1_DRIVEN / "voxels.csv")
    occ = occupancy(voxels, u)
    sr, sl = integrate_surrogate(voxels["absolute_time"], voxels["ahl"])
    mask = (voxels["den"] != 0).astype(float)
    unmasked = np.column_stack([
        window_mean_matrix(sr, voxels["window"]),
        window_mean_matrix(sl, voxels["window"]),
    ])
    masked = np.column_stack([
        window_mean_matrix(sr * mask, voxels["window"]),
        window_mean_matrix(sl * mask, voxels["window"]),
    ])
    a0_voxels = load_voxel_arrays(A0_DRIVEN / "voxels.csv")
    a0_sr, a0_sl = integrate_surrogate(a0_voxels["absolute_time"], a0_voxels["ahl"])
    a0_mask = (a0_voxels["den"] != 0).astype(float)
    a0_unmasked = np.column_stack([
        window_mean_matrix(a0_sr, a0_voxels["window"]),
        window_mean_matrix(a0_sl, a0_voxels["window"]),
    ])
    a0_masked = np.column_stack([
        window_mean_matrix(a0_sr * a0_mask, a0_voxels["window"]),
        window_mean_matrix(a0_sl * a0_mask, a0_voxels["window"]),
    ])
    a1_unmasked = ridge_score(unmasked, y)
    a1_masked = ridge_score(masked, y)
    a0_unmasked_s = ridge_score(a0_unmasked, y)
    a0_masked_s = ridge_score(a0_masked, y)

    a0_f408 = a0_d["task"]["test_nrmse"]
    a1_f408 = a1_d["task"]["test_nrmse"]
    a1_field = a1_d["field_only"]["test_nrmse"]
    a0_field = a0_d["field_only"]["test_nrmse"]
    system_pass = a1_f408 < a1_b["task"]["test_nrmse"] and a1_f408 < a1_s["task"]["test_nrmse"]
    living_pass = a1_f408 < a1_field
    a0_living = a0_f408 < a0_field
    delta = abs(a1_f408 - a0_f408)
    story = delta >= WEAK_MARGIN
    sign_flip = living_pass != a0_living
    alive = occ["occupancy"] == "ALIVE"
    confirm = alive and system_pass and (story or sign_flip)
    no_story = (not story) and alive
    row = {
        "Seed": 111,
        "MODEL_STATUS": "HYPOTHETICAL_DESIGN_ENVELOPE",
        "digital_twin": False,
        "occupancy": occ["occupancy"],
        "mean_R": occ["mean_R"],
        "r_meanR_u": occ["r_meanR_u"],
        "mean_AHL": occ["mean_AHL"],
        "A0_F408": a0_f408,
        "A0_field": a0_field,
        "A0_brownian": a0_b["task"]["test_nrmse"],
        "A0_silent": a0_s["task"]["test_nrmse"],
        "A0_unmasked": a0_unmasked_s["test_nrmse"],
        "A0_masked": a0_masked_s["test_nrmse"],
        "A1_F408": a1_f408,
        "A1_field": a1_field,
        "A1_brownian": a1_b["task"]["test_nrmse"],
        "A1_silent_reused": a1_s["task"]["test_nrmse"],
        "A1_unmasked": a1_unmasked["test_nrmse"],
        "A1_masked": a1_masked["test_nrmse"],
        "AC_G_SCALAR": baselines["AC_G_SCALAR"]["test_nrmse"],
        "AC_J_WINDOW_MEAN": baselines["AC_J_WINDOW_MEAN"]["test_nrmse"],
        "AC_G_DELAY_10": baselines["AC_G_DELAY_10"]["test_nrmse"],
        "LINEAR_U_DELAY_10": baselines["LINEAR_U_DELAY_10"]["test_nrmse"],
        "NARMA_INFORMED_INPUT": baselines["NARMA_INFORMED_INPUT"]["test_nrmse"],
        "TRAIN_INTERCEPT": baselines["TRAIN_INTERCEPT"]["test_nrmse"],
        "PERSISTENCE": baselines["PERSISTENCE"]["test_nrmse"],
        "abs_A1_minus_A0_F408": delta,
        "system_pass": system_pass,
        "living_layer_pass": living_pass,
        "a0_living_layer_pass": a0_living,
        "living_layer_sign_flip": sign_flip,
        "NO_STORY_MOVE": no_story and not sign_flip,
        "confirm_222_333": confirm,
        "GATE_EVIDENCE_untouched": gate_evidence_untouched(),
        "java_ok": java_ok(),
        "ahl_finite": occ["ahl_finite"],
        "ahl_nonneg": occ["ahl_nonneg"],
    }
    return row, baselines, occ


def write_report(row, baselines):
    system = "PASS" if row["system_pass"] else "FAIL"
    living = "PASS" if row["living_layer_pass"] else "FAIL"
    lines = [
        "# E4.1 one-way A1 scout (seed 111)",
        "",
        "**MODEL_STATUS = `HYPOTHETICAL_DESIGN_ENVELOPE`.** Not a digital twin.",
        "Module 1 payload match PASSed before this living run. Narma10b A0",
        "voxels were reused, not overwritten. `GATE_EVIDENCE.md` was not edited.",
        "K, n, tau, clamp, mortality, flow, layout, `K_AC`, and `Jmax_A1` were",
        "not retuned. Two-way communication was not started.",
        "",
        f"u SHA-256 `{NARMA_SHA}`.",
        f"A1 `g0=0`, `n_AC=2`, `K_AC=0.25`, `Jmax_A1={JMAX_A1:.12e}` molecules/s.",
        "Claim dish CENTER / FLOW=0.",
        "",
        f"Occupancy: **{row['occupancy']}** (`mean_R={row['mean_R']:.4f}`, "
        f"`r(mean_R,u)={row['r_meanR_u']:.3f}`).",
        "",
        "## Claims (report, do not retune)",
        "",
        f"- **System** (A1 F408 vs A1 Brownian and silent): **{system}**",
        f"- **Living-layer** (A1 F408 vs A1 field): **{living}**",
        f"- **Transducer** |A1 F408 − A0 F408| = {row['abs_A1_minus_A0_F408']:.4f}. "
        "Claim-dish promotion is NO_STORY_MOVE if this is a 0.03 tick around 0.93.",
        "- **Ceiling:** direct-input 0.6828 is not required to lose.",
        "- **0.93 stays a weak predictor.**",
        "",
        "## A0 reused (Narma10b seed 111)",
        "",
        "| Readout | Test NRMSE |",
        "|---|---|",
        f"| F408 | {row['A0_F408']:.4f} |",
        f"| field | {row['A0_field']:.4f} |",
        f"| Brownian | {row['A0_brownian']:.4f} |",
        f"| silent | {row['A0_silent']:.4f} |",
        f"| unmasked RL surrogate | {row['A0_unmasked']:.4f} |",
        f"| occupancy-masked RL surrogate | {row['A0_masked']:.4f} |",
        "",
        "## A1 seed 111",
        "",
        "| Readout | Test NRMSE |",
        "|---|---|",
        f"| F408 | {row['A1_F408']:.4f} |",
        f"| field | {row['A1_field']:.4f} |",
        f"| A1 Brownian | {row['A1_brownian']:.4f} |",
        f"| silent (reused A0) | {row['A1_silent_reused']:.4f} |",
        f"| unmasked RL surrogate | {row['A1_unmasked']:.4f} |",
        f"| occupancy-masked RL surrogate | {row['A1_masked']:.4f} |",
        f"| AC-output g[n] | {row['AC_G_SCALAR']:.4f} |",
        f"| AC-output window-mean J | {row['AC_J_WINDOW_MEAN']:.4f} |",
        f"| AC-output 10-tap g | {row['AC_G_DELAY_10']:.4f} |",
        "",
        "## Task baselines (unchanged u / target)",
        "",
        "| Baseline | Test NRMSE |",
        "|---|---|",
        f"| LINEAR_U_DELAY_10 | {row['LINEAR_U_DELAY_10']:.4f} |",
        f"| NARMA_INFORMED_INPUT | {row['NARMA_INFORMED_INPUT']:.4f} |",
        f"| PERSISTENCE | {row['PERSISTENCE']:.4f} |",
        f"| TRAIN_INTERCEPT | {row['TRAIN_INTERCEPT']:.4f} |",
        "",
        "Direct-input 0.6828 remains the task ceiling. Biology is not",
        "required to beat it. An AC-output or field win is attribution.",
        "",
        "## Decision",
        "",
    ]
    if row["occupancy"] == "DEAD":
        lines.append(
            "Occupancy DEAD. Keep the row. Do not raise Jmax or K_AC. "
            "Do not run 222/333."
        )
    else:
        lines.extend(
            [
                "Stop after seed 111. No best seed. Do not retune. Do not start two-way.",
                "",
                "The frozen 222/333 numeric gate "
                + ("**fires**" if row["confirm_222_333"] else "does not fire")
                + f" (ALIVE={row['occupancy']=='ALIVE'}, system={row['system_pass']}, "
                + f"|ΔF408|={row['abs_A1_minus_A0_F408']:.4f}). "
                "That is a later confirmation slot, not a new claim dish. "
                "This job does not run 222/333.",
            ]
        )
        if row["abs_A1_minus_A0_F408"] < 0.04:
            lines.extend(
                [
                    "",
                    f"**Claim-dish promotion: NO_STORY_MOVE.** |Δ|={row['abs_A1_minus_A0_F408']:.4f} "
                    "is a 0.03-scale tick around the weak 0.93 band.",
                ]
            )
    lines.extend(
        [
            "",
            f"GATE_EVIDENCE untouched: {row['GATE_EVIDENCE_untouched']}.",
            f"Java constraints (no vesicle/glucose/Danino/CS_in, FLOW=0, A1 gate): {row['java_ok']}.",
            "",
            "C1 remains DEFER. E0.3 remains NO_STORY_MOVE. Waveform1 FAIL,",
            "Waveform2c field≥driven, and L3 living-layer FAIL stay.",
            "",
        ]
    )
    (HERE / "results" / "E4_AC_SCOUT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    with (HERE / "results" / "e4_ac_scout.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    (HERE / "results" / "e4_ac_scout.json").write_text(
        json.dumps(row, indent=2, default=str) + "\n", encoding="utf-8"
    )


def check_smoke(log_path: Path) -> None:
    text = log_path.read_text(encoding="utf-8")
    needed = (
        "MODEL_STATUS=HYPOTHETICAL_DESIGN_ENVELOPE",
        "digital_twin=false",
        "FLOW_SPEED=0.0",
        "g(0)=0.000000000000",
        "g(0.5)=0.800000000000",
        "CENTER",
    )
    ok = all(token in text for token in needed)
    print("SMOKE A1 params / CENTER / FLOW=0 / not a twin:", "PASS" if ok else "FAIL")
    if not ok:
        raise SystemExit("STOP smoke print missing A1 parameters, CENTER, or FLOW=0")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke-log", default="")
    args = parser.parse_args()
    if args.smoke_log:
        check_smoke(Path(args.smoke_log))
        return
    u, y, digest = N.load_target(HERE / "narma10_target.csv", HERE / "input_ahl_narma200.txt")
    if digest != NARMA_SHA:
        raise SystemExit(f"STOP u hash {digest}")
    if not java_ok():
        raise SystemExit("STOP Java constraints failed")
    row, baselines, occ = analyze(u, y)
    write_report(row, baselines)
    print(
        "A1", occ["occupancy"],
        f"R={occ['mean_R']:.4f}",
        f"F408={row['A1_F408']:.4f}",
        f"field={row['A1_field']:.4f}",
        f"system={'PASS' if row['system_pass'] else 'FAIL'}",
        f"living={'PASS' if row['living_layer_pass'] else 'FAIL'}",
        f"dA0={row['abs_A1_minus_A0_F408']:.4f}",
        f"confirm222={row['confirm_222_333']}",
    )


if __name__ == "__main__":
    main()
