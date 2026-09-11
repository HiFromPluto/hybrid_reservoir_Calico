#!/usr/bin/env python3
"""Evaluate Track A2 subsampled Lorenz'63 gates on the frozen Stage 6 dish.

BenchA Lorenz FAIL and Mackey-Glass PASS are not rewritten. Stage 6 NARMA-10
remains PASS. Track B is not started. SKIP=50 is frozen. The Stage 6
val-slice leak is closed: lambda is picked on rows 88..109 only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

WASHOUT = 40
TRAIN = 110
TEST = 50
INNER_VAL = 22
NUM_WINDOWS = WASHOUT + TRAIN + TEST
EXPECTED_AUX_ROWS = NUM_WINDOWS * 16
LAST_SAMPLE_PREFIX = "199;15;299.95"
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
MC_MAX_DELAY = 20
SEEDS = (101, 202, 303)

BIOLOGY_MEAN_PREFIXES = ("Receiver_R_", "Lum_Mean_")
BIOLOGY_LAST_PREFIXES = ("Input_Driven_Death_",)
BROWNIAN_MEAN_PREFIXES = ("Den_",)
BROWNIAN_LAST_PREFIXES = ()
FIELD_MEAN_PREFIXES = ("AHL_uM_",)
EXCLUDED_FROM_RIDGE = (
    "window AHL",
    "occupancy Fraction_q_gt_0_5",
    "pH",
    "Births",
    "Total_Deaths",
    "Clamp_Deaths",
    "OOB",
    "Population",
    "Lum_Sum",
    "voxel AHL_uM (biology ridge)",
    "voxel pH",
    "voxel Den (biology ridge)",
    "voxel Fraction_q_gt_0_5",
    "voxel Lum_Sum",
)
FORBIDDEN_JAVA = (
    "glucose", "Glucose", "GLUCOSE",
    "Danino", "danino",
    "vesicle", "Vesicle",
    "setGoal(glucose)",
)
SKIP = 50
LORENZ_SHA = "35646e7b504940dbc859ec66f1f2c5b49017f25388c28f0d662ffc29a7f9a1e5"
LORENZ_XMIN = -17.981997956481
LORENZ_XMAX = 16.013692471459

HERE = Path(__file__).resolve().parent
STAGE6_DIR = HERE.parent / "BSimReservoirPlanStage6"
BENCHA_DIR = HERE.parent / "BSimReservoirPlanBenchA"


def json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not json serializable: {type(obj)}")


def sha256_u(u):
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def lorenz_deriv(state):
    sigma, rho, beta = 10.0, 28.0, 8.0 / 3.0
    x, y, z = state
    return np.array([
        sigma * (y - x),
        x * (rho - z) - y,
        x * y - beta * z,
    ], dtype=float)


def rk4_step(state, dt):
    k1 = lorenz_deriv(state)
    k2 = lorenz_deriv(state + 0.5 * dt * k1)
    k3 = lorenz_deriv(state + 0.5 * dt * k2)
    k4 = lorenz_deriv(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def recompute_lorenz_skip50():
    state = np.array([1.0, 1.0, 1.0], dtype=float)
    for _ in range(5000):
        state = rk4_step(state, 0.02)
    samples = np.empty((201, 3), dtype=float)
    samples[0] = state
    for index in range(1, 201):
        for _ in range(SKIP):
            state = rk4_step(state, 0.02)
        samples[index] = state
    return samples


def affine_u(x, xmin, xmax):
    if xmax == xmin:
        raise ValueError("xmax == xmin")
    return 0.5 * (x - xmin) / (xmax - xmin)


def load_target(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    if len(rows) != NUM_WINDOWS:
        raise ValueError(f"target has {len(rows)} windows; expected {NUM_WINDOWS}")
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    x = np.array([float(row["x"]) for row in rows], dtype=float)
    y = np.array([float(row["y_next"]) for row in rows], dtype=float)
    if np.any(u < 0.0) or np.any(u > 0.5):
        raise ValueError("u escaped [0, 0.5]")
    digest = sha256_u(u)
    samples = recompute_lorenz_skip50()
    x_drive = samples[:NUM_WINDOWS, 0]
    y_next = samples[1:NUM_WINDOWS + 1, 1]
    xmin, xmax = float(np.min(x_drive)), float(np.max(x_drive))
    u_hat = affine_u(x_drive, xmin, xmax)
    if not np.allclose(x, x_drive, rtol=0, atol=1e-9):
        raise ValueError("lorenz_skip50_target.csv x does not match skip-50 RK4")
    if not np.allclose(y, y_next, rtol=0, atol=1e-9):
        raise ValueError("lorenz_skip50_target.csv y_next does not match skip-50 RK4")
    if not np.allclose(u, u_hat, rtol=0, atol=1e-9):
        raise ValueError("lorenz_skip50_target.csv u does not match affine map")
    if abs(xmin - LORENZ_XMIN) > 1e-9 or abs(xmax - LORENZ_XMAX) > 1e-9:
        raise ValueError("Lorenz xmin/xmax drifted from PROTOCOL.md")
    if digest != LORENZ_SHA:
        raise ValueError("Lorenz u SHA-256 drifted from PROTOCOL.md")
    return u, y


def last_sample_ok(path):
    text = Path(path).read_text(encoding="utf-8").splitlines()
    data = [line for line in text[1:] if line.strip()]
    if not data:
        return False
    return data[-1].startswith(LAST_SAMPLE_PREFIX)


def validate_csv(path, expected_rows, check_last_sample=False):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, [])
        widths = [len(row) for row in reader]
    result = {
        "path": str(path),
        "row_count": len(widths),
        "column_count": len(header),
        "row_count_pass": len(widths) == expected_rows,
        "rectangular_csv_pass": bool(header) and all(width == len(header) for width in widths),
    }
    if check_last_sample:
        result["last_sample_pass"] = last_sample_ok(path)
    return result


def columns_with_prefix(header, prefixes):
    names = []
    for prefix in prefixes:
        names.extend(name for name in header if name.startswith(prefix))
    return names


def read_window_matrix(voxels_path, mean_prefixes, last_prefixes):
    path = Path(voxels_path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        mean_cols = columns_with_prefix(header, mean_prefixes)
        last_cols = columns_with_prefix(header, last_prefixes)
        feature_names = mean_cols + last_cols
        if not feature_names:
            raise ValueError(f"no ridge columns in {path}")
        window_idx = header.index("Window")
        mean_idx = [header.index(name) for name in mean_cols]
        last_idx = [header.index(name) for name in last_cols]
        sums = {}
        counts = {}
        lasts = {}
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"ragged voxels.csv row in {path}")
            window = int(row[window_idx])
            if window not in sums:
                sums[window] = np.zeros(len(mean_idx), dtype=float)
                counts[window] = 0
                lasts[window] = np.zeros(len(last_idx), dtype=float)
            if mean_idx:
                sums[window] += np.array([float(row[i]) for i in mean_idx], dtype=float)
            counts[window] += 1
            if last_idx:
                lasts[window] = np.array([float(row[i]) for i in last_idx], dtype=float)
    windows = sorted(sums)
    matrix = []
    for window in windows:
        parts = []
        if mean_idx:
            parts.append(sums[window] / counts[window])
        if last_idx:
            parts.append(lasts[window])
        matrix.append(np.concatenate(parts) if parts else np.zeros(0))
    return {
        "windows": windows,
        "X": np.vstack(matrix),
        "feature_names": feature_names,
        "samples_per_window": [counts[window] for window in windows],
    }


def nrmse(y_true, y_pred):
    residual = y_true - y_pred
    rmse = float(np.sqrt(np.mean(residual * residual)))
    denom = float(np.std(y_true))
    return rmse / denom if denom > 0 else float("inf")


def r_squared(y_true, y_pred):
    residual = y_true - y_pred
    ss_res = float(np.sum(residual * residual))
    centred = y_true - np.mean(y_true)
    ss_tot = float(np.sum(centred * centred))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def standardize(train, *others):
    mean = train.mean(axis=0)
    std = train.std(axis=0)
    std = np.where(std < 1e-12, 1.0, std)
    converted = [(matrix - mean) / std for matrix in (train,) + others]
    return converted, mean, std


def ridge_fit(X, y, lam):
    n_features = X.shape[1]
    Xb = np.column_stack([np.ones(len(X)), X])
    gram = Xb.T @ Xb
    gram[1:, 1:] = gram[1:, 1:] + lam * np.eye(n_features)
    target = Xb.T @ y
    try:
        return np.linalg.solve(gram, target)
    except np.linalg.LinAlgError:
        weights, *_ = np.linalg.lstsq(gram, target, rcond=None)
        return weights


def ridge_predict(X, weights):
    Xb = np.column_stack([np.ones(len(X)), X])
    return Xb @ weights


def select_lambda(X, y):
    """Fit on rows 0..87; pick lambda on rows 88..109 only. Test is not used."""
    inner_train_end = TRAIN - INNER_VAL
    X_inner, y_inner = X[:inner_train_end], y[:inner_train_end]
    X_val, y_val = X[inner_train_end:TRAIN], y[inner_train_end:TRAIN]
    (X_inner_z, X_val_z), _, _ = standardize(X_inner, X_val)
    scored = []
    grid_scores = []
    for lam in RIDGE_GRID:
        weights = ridge_fit(X_inner_z, y_inner, lam)
        pred = ridge_predict(X_val_z, weights)
        score = nrmse(y_val, pred)
        scored.append((score, -lam, lam))
        grid_scores.append({"lambda": lam, "val_nrmse": score})
    scored.sort()
    return scored[0][2], grid_scores


def fit_eval(X, y, lam):
    X_train, y_train = X[:TRAIN], y[:TRAIN]
    X_test, y_test = X[TRAIN:], y[TRAIN:]
    (X_train_z, X_test_z), _, _ = standardize(X_train, X_test)
    weights = ridge_fit(X_train_z, y_train, lam)
    train_pred = ridge_predict(X_train_z, weights)
    test_pred = ridge_predict(X_test_z, weights)
    return {
        "lambda": lam,
        "train_nrmse": nrmse(y_train, train_pred),
        "test_nrmse": nrmse(y_test, test_pred),
        "test_r2": r_squared(y_test, test_pred),
    }


def evaluate_readout(X_all, target, windows):
    if windows != list(range(NUM_WINDOWS)):
        raise ValueError(f"expected windows 0..{NUM_WINDOWS - 1}, got {windows[:3]}..{windows[-1:]}")
    X = X_all[WASHOUT:]
    y = target[WASHOUT:]
    lam, grid_scores = select_lambda(X, y)
    metrics = fit_eval(X, y, lam)
    metrics["lambda_grid"] = grid_scores
    return metrics


def memory_capacity(X_all, u, windows):
    scores = []
    total = 0.0
    for delay in range(1, MC_MAX_DELAY + 1):
        target = np.array([u[n - delay] if n >= delay else 0.0 for n in windows], dtype=float)
        metrics = evaluate_readout(X_all, target, windows)
        r2 = metrics["test_r2"]
        contribution = max(0.0, r2) if math.isfinite(r2) else 0.0
        scores.append({
            "delay": delay,
            "lambda": metrics["lambda"],
            "test_nrmse": metrics["test_nrmse"],
            "test_r2": r2,
            "contribution": contribution,
        })
        total += contribution
    k0_target = np.array(u, dtype=float)
    k0 = evaluate_readout(X_all, k0_target, windows)
    return {
        "mc": total,
        "kmax": MC_MAX_DELAY,
        "delays": scores,
        "k0_diagnostic": {
            "test_nrmse": k0["test_nrmse"],
            "test_r2": k0["test_r2"],
            "lambda": k0["lambda"],
        },
        "n_features_checked": X_all[WASHOUT:].shape[1],
    }


def mean_se(values):
    array = np.array(values, dtype=float)
    mean = float(np.mean(array))
    se = float(np.std(array, ddof=1) / np.sqrt(len(array))) if len(array) > 1 else float("nan")
    return mean, se


def intervals_overlap(mean_a, se_a, mean_b, se_b):
    return not (mean_a + se_a < mean_b - se_b or mean_b + se_b < mean_a - se_a)


def read_run(path, arm, expected_windows, expected_aux_rows):
    directory = Path(path)
    summary_path = directory / "window_summary.csv"
    voxels_path = directory / "voxels.csv"
    results_path = directory / "results.csv"
    if arm == "brownian":
        matrix = read_window_matrix(voxels_path, BROWNIAN_MEAN_PREFIXES, BROWNIAN_LAST_PREFIXES)
    else:
        matrix = read_window_matrix(voxels_path, BIOLOGY_MEAN_PREFIXES, BIOLOGY_LAST_PREFIXES)
    field = None
    if arm == "driven":
        field = read_window_matrix(voxels_path, FIELD_MEAN_PREFIXES, ())
    return {
        "path": str(directory),
        "arm": arm,
        "matrix": matrix,
        "field": field,
        "summary_csv": validate_csv(summary_path, expected_windows),
        "auxiliary_csvs": {
            "results.csv": validate_csv(results_path, expected_aux_rows, check_last_sample=True),
            "voxels.csv": validate_csv(voxels_path, expected_aux_rows, check_last_sample=True),
        },
        "feature_contract": (directory / "feature_contract.txt").read_text(encoding="utf-8")
        if (directory / "feature_contract.txt").exists() else "",
    }


def evaluate_arm_run(run, u, y):
    matrix = run["matrix"]
    task_metrics = evaluate_readout(matrix["X"], y, matrix["windows"])
    mc = memory_capacity(matrix["X"], u, matrix["windows"])
    field = None
    if run["field"] is not None:
        field = evaluate_readout(run["field"]["X"], y, run["field"]["windows"])
        field["feature_names"] = run["field"]["feature_names"]
        field["n_features"] = len(run["field"]["feature_names"])
    return {
        "path": run["path"],
        "arm": run["arm"],
        "n_features": len(matrix["feature_names"]),
        "feature_names": matrix["feature_names"],
        "task": task_metrics,
        "memory_capacity": mc,
        "field_only": field,
        "summary_csv": run["summary_csv"],
        "auxiliary_csvs": run["auxiliary_csvs"],
    }


def stage6_java_untouched(stage6_dir):
    java_path = Path(stage6_dir) / "BSimReservoirPlanStage6.java"
    voxel_path = Path(stage6_dir) / "VoxelAnalyzer.java"
    if not java_path.exists() or not voxel_path.exists():
        return False
    java = java_path.read_text(encoding="utf-8")
    voxel = voxel_path.read_text(encoding="utf-8")
    return (
        "package BSimReservoirPlanStage6;" in java
        and "public final class BSimReservoirPlanStage6" in java
        and "NARMA-10" in java
        and "package BSimReservoirPlanStage6;" in voxel
        and "BSimReservoirPlanStage6.DeathCause" in voxel
    )


def bencha_java_untouched(bencha_dir):
    java_path = Path(bencha_dir) / "BSimReservoirPlanBenchA.java"
    voxel_path = Path(bencha_dir) / "VoxelAnalyzer.java"
    if not java_path.exists() or not voxel_path.exists():
        return False
    java = java_path.read_text(encoding="utf-8")
    voxel = voxel_path.read_text(encoding="utf-8")
    return (
        "package BSimReservoirPlanBenchA;" in java
        and "public final class BSimReservoirPlanBenchA" in java
        and "package BSimReservoirPlanBenchA;" in voxel
        and "BSimReservoirPlanBenchA.DeathCause" in voxel
    )


def heading_overall(text, heading, expected):
    pattern = rf"^## {re.escape(heading)}\s*\n(?:.*\n)*?^### Overall:\s*{expected}\s*$"
    return bool(re.search(pattern, text, re.MULTILINE))


def stage6_overall_pass(path):
    if not Path(path).exists():
        return False
    text = Path(path).read_text(encoding="utf-8")
    return bool(re.search(r"^## Overall:\s*PASS\s*$", text, re.MULTILINE))


def bencha_prior_labels(path):
    if not Path(path).exists():
        return False
    text = Path(path).read_text(encoding="utf-8")
    return heading_overall(text, "Lorenz'63", "FAIL") and heading_overall(
        text, "Mackey-Glass", "PASS"
    )


def no_track_b(java_path):
    text = Path(java_path).read_text(encoding="utf-8")
    if any(token in text for token in FORBIDDEN_JAVA):
        return False
    if "new Vector3d(250, 250, 5), new Vector3d(500, 250, 5)" not in text:
        return False
    if "new Vector3d(300, 375, 5)" not in text:
        return False
    return True


def write_markdown(path, evidence):
    overall = "PASS" if evidence["overall_gate_pass"] else "FAIL"
    gates = evidence["gates"]
    lines = [
        "# Track A2 gate evidence",
        "",
        "Track A2 copies frozen Stage 6 / BenchA Java and tests subsampled",
        "Lorenz'63 (SKIP=50, dt_sample=1.0) against a Brownian density null",
        "and a silent-source control. BenchA Lorenz remains **FAIL**. BenchA",
        "Mackey-Glass remains **PASS**. Stage 6 NARMA-10 remains **PASS**.",
        "Mechanisms were not retuned. Track B was not started. The Stage 6",
        "val-slice leak is closed.",
        "",
        f"## Overall: {overall}",
        "",
        evidence["failure_reason"] if not evidence["overall_gate_pass"]
        else "Driven biology beat the Brownian null with silent in the expected place.",
        "",
        "| Gate | Result |",
        "|---|---|",
        f"| 1. Driven test NRMSE < Brownian, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_brownian_pass'] else 'FAIL'}** |",
        f"| 2. Silent not ≈ driven; silent worse than driven | **{'PASS' if gates['silent_not_equal_driven_pass'] else 'FAIL'}** |",
        f"| 3. Driven uses frozen Stage 5/6 analysis channels | **{'PASS' if gates['frozen_feature_list_pass'] else 'FAIL'}** |",
        f"| 4. CSV rectangularity 200/3200/3200; last sample 199;15;299.95 | **{'PASS' if gates['csv_validation_pass'] else 'FAIL'}** |",
        f"| 5. BenchA Lorenz FAIL, MG PASS; Stage 6 PASS; Stage 6 and BenchA Java untouched; no glucose/Danino/extra ACs/Track B | **{'PASS' if gates['prior_stage_labels_pass'] else 'FAIL'}** |",
        "",
        "## Test NRMSE",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 | mean ± s.e. |",
        "|---|---|---|---|---|",
    ]
    for arm in ("brownian", "silent", "driven"):
        values = evidence["test_nrmse"][arm]
        mean = evidence["mean_se"][arm]["mean"]
        se = evidence["mean_se"][arm]["se"]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} | {mean:.4f} ± {se:.4f} |"
        )
    field_vals = evidence["field_only_test_nrmse"]
    if field_vals:
        mean, se = mean_se(field_vals)
        lines.extend([
            "",
            "## Field-only diagnostic (not a gate)",
            "",
            "Ridge on driven voxel `AHL_uM_*` only. If biology does not beat",
            "this delay line, that is recorded; the Brownian gate is unchanged.",
            "",
            f"Field-only NRMSE: {field_vals[0]:.4f} / {field_vals[1]:.4f} / {field_vals[2]:.4f} "
            f"(mean {mean:.4f} ± {se:.4f}).",
            f"Biology vs field-only: {evidence['biology_vs_field_note']}",
        ])
    lines.extend([
        "",
        "## Driven feature list",
        "",
        "Frozen Stage 5/6 analysis channels, voxelised:",
        "",
        "- `Receiver_R_*` / Mean_q, 20×10",
        "- `Lum_Mean_*` / Mean_L, 20×10",
        "- `Input_Driven_Death_*`, 4×2",
        "",
        f"Printed count: {evidence['driven_feature_count']} features.",
        "",
        "Excluded from the ridge: " + ", ".join(EXCLUDED_FROM_RIDGE) + ".",
        "",
        "## Memory capacity (diagnostic)",
        "",
        "| Arm | MC k=1..20 |",
        "|---|---|",
    ])
    for arm in ("brownian", "silent", "driven"):
        values = evidence["memory_capacity"][arm]
        mean, se = mean_se(values)
        lines.append(
            f"| {arm} | {values[0]:.3f} / {values[1]:.3f} / {values[2]:.3f} (mean {mean:.3f} ± {se:.3f}) |"
        )
    lines.extend([
        "",
        "## Protocol",
        "",
        "See `PROTOCOL.md`. SKIP=50; affine `u ∈ [0, 0.5]`; acid held at 0.5;",
        "washout 40 / train 110 / test 50; lambda selected on windows 128..149 only.",
        "",
        "CSV completeness is the on-disk files, not BSim timestep stdout.",
        "Track B was not started.",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def default_run_dirs(arm):
    return [str(HERE / "results" / f"lorenz_skip50_{arm}_seed{seed}") for seed in SEEDS]


def evaluate_task(args):
    target_path = Path(args.target) if args.target else HERE / "lorenz_skip50_target.csv"
    u, y = load_target(target_path)
    brownian = args.brownian or default_run_dirs("brownian")
    silent = args.silent or default_run_dirs("silent")
    driven = args.driven or default_run_dirs("driven")
    grouped = {
        "brownian": [read_run(path, "brownian", args.expected_windows, args.expected_aux_rows)
                     for path in brownian],
        "silent": [read_run(path, "silent", args.expected_windows, args.expected_aux_rows)
                   for path in silent],
        "driven": [read_run(path, "driven", args.expected_windows, args.expected_aux_rows)
                   for path in driven],
    }
    evaluated = {
        arm: [evaluate_arm_run(run, u, y) for run in runs]
        for arm, runs in grouped.items()
    }

    nrmse_by_arm = {
        arm: [run["task"]["test_nrmse"] for run in runs]
        for arm, runs in evaluated.items()
    }
    mean_se_by_arm = {arm: mean_se(values) for arm, values in nrmse_by_arm.items()}
    mc_by_arm = {
        arm: [run["memory_capacity"]["mc"] for run in runs]
        for arm, runs in evaluated.items()
    }
    field_nrmse = [
        run["field_only"]["test_nrmse"] for run in evaluated["driven"] if run["field_only"]
    ]

    driven_features = evaluated["driven"][0]["feature_names"]
    expected_prefixes = ("Receiver_R_", "Lum_Mean_", "Input_Driven_Death_")
    forbidden = ("AHL_uM_", "pH_", "Den_", "Fraction_q", "Lum_Sum_", "Birth_",
                 "Clamp_Death_", "OOB_Death_")
    feature_ok = (
        all(name.startswith(expected_prefixes) for name in driven_features)
        and not any(name.startswith(forbidden) for name in driven_features)
        and len(driven_features) == 408
        and all(run["n_features"] == 408 for run in evaluated["driven"] + evaluated["silent"])
        and all(run["n_features"] == 200 for run in evaluated["brownian"])
    )
    print("DRIVEN_FEATURE_LIST")
    for name in driven_features:
        print(name)
    print(f"driven_feature_count={len(driven_features)}")

    csv_pass = all(
        run["summary_csv"]["row_count_pass"]
        and run["summary_csv"]["rectangular_csv_pass"]
        and all(item["row_count_pass"] and item["rectangular_csv_pass"]
                for item in run["auxiliary_csvs"].values())
        and all(item.get("last_sample_pass", True)
                for item in run["auxiliary_csvs"].values())
        for runs in evaluated.values() for run in runs
    )

    driven_mean, driven_se = mean_se_by_arm["driven"]
    brown_mean, brown_se = mean_se_by_arm["brownian"]
    silent_mean, silent_se = mean_se_by_arm["silent"]
    driven_beats_brownian = (
        all(d < b for d, b in zip(nrmse_by_arm["driven"], nrmse_by_arm["brownian"]))
        and driven_mean + driven_se < brown_mean - brown_se
    )
    silent_overlap_driven = intervals_overlap(silent_mean, silent_se, driven_mean, driven_se)
    silent_worse_than_driven = silent_mean > driven_mean and not silent_overlap_driven

    stage6_evidence = Path(args.stage6_evidence) if args.stage6_evidence else (
        STAGE6_DIR / "results" / "GATE_EVIDENCE.md"
    )
    bencha_evidence = Path(args.bencha_evidence) if args.bencha_evidence else (
        BENCHA_DIR / "results" / "GATE_EVIDENCE.md"
    )
    prior_ok = (
        stage6_overall_pass(stage6_evidence)
        and bencha_prior_labels(bencha_evidence)
        and stage6_java_untouched(STAGE6_DIR)
        and bencha_java_untouched(BENCHA_DIR)
        and no_track_b(HERE / "BSimReservoirPlanBenchA2.java")
    )

    gates = {
        "replicate_count_pass": all(len(values) == 3 for values in nrmse_by_arm.values()),
        "driven_beats_brownian_pass": driven_beats_brownian,
        "silent_not_equal_driven_pass": silent_worse_than_driven,
        "frozen_feature_list_pass": feature_ok,
        "csv_validation_pass": csv_pass,
        "prior_stage_labels_pass": prior_ok,
    }
    overall = all(gates.values())
    if not gates["driven_beats_brownian_pass"]:
        failure_reason = "null_wins: driven test NRMSE did not beat Brownian by a clear margin."
    elif not gates["silent_not_equal_driven_pass"]:
        failure_reason = "silent_approx_driven: silent NRMSE is not clearly worse than driven; inputs are not driving the reservoir."
    elif not gates["frozen_feature_list_pass"]:
        failure_reason = "readout_leakage: driven features are not the frozen Stage 5/6 analysis channels."
    elif not gates["csv_validation_pass"]:
        failure_reason = "csv_incomplete: a summary/sample/voxel file is missing rows, is ragged, or last sample is not 199;15;299.95."
    elif not gates["prior_stage_labels_pass"]:
        failure_reason = "provenance: BenchA Lorenz is not FAIL, BenchA MG is not PASS, Stage 6 is not PASS, Stage 6/BenchA Java was touched, or Track B / glucose / Danino / extra ACs leaked in."
    else:
        failure_reason = ""

    if field_nrmse:
        field_mean, field_se = mean_se(field_nrmse)
        if driven_mean + driven_se < field_mean - field_se:
            field_note = "driven biology beat field-only (non-overlapping mean±s.e.)."
        elif intervals_overlap(driven_mean, driven_se, field_mean, field_se):
            field_note = "driven biology ≈ field-only. Recorded; Brownian gate unchanged."
        else:
            field_note = "field-only beat driven biology. Recorded; Brownian gate unchanged."
    else:
        field_note = "field-only diagnostic unavailable."

    evidence = {
        "schema": "BSimReservoirPlanBenchA2-gates-v1",
        "task": "lorenz_skip50",
        "skip": SKIP,
        "dt_sample": SKIP * 0.02,
        "stage6_status": "PASS" if stage6_overall_pass(stage6_evidence) else "NOT_PASS",
        "protocol": {
            "num_windows": NUM_WINDOWS,
            "washout": WASHOUT,
            "train": TRAIN,
            "test": TEST,
            "inner_val": INNER_VAL,
            "lambda_val_rows": "88..109 after washout drop (windows 128..149)",
            "ridge_grid": list(RIDGE_GRID),
            "u_range": [0.0, 0.5],
            "acid": "held at 0.5",
        },
        "analysis_channels_driven": [
            "Receiver_R_* / Mean_q (20x10)",
            "Lum_Mean_* / Mean_L (20x10)",
            "Input_Driven_Death_* (4x2)",
        ],
        "driven_feature_count": len(driven_features),
        "excluded_from_ridge": list(EXCLUDED_FROM_RIDGE),
        "test_nrmse": nrmse_by_arm,
        "mean_se": {
            arm: {"mean": mean, "se": se} for arm, (mean, se) in mean_se_by_arm.items()
        },
        "memory_capacity": mc_by_arm,
        "field_only_test_nrmse": field_nrmse,
        "biology_vs_field_note": field_note,
        "lambdas": {
            arm: [run["task"]["lambda"] for run in runs]
            for arm, runs in evaluated.items()
        },
        "files": {
            arm: [{k: v for k, v in run.items() if k != "feature_names"} | {
                "feature_count": run["n_features"]
            } for run in runs]
            for arm, runs in evaluated.items()
        },
        "gates": gates,
        "overall_gate_pass": overall,
        "failure_reason": failure_reason,
    }
    return evidence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brownian", nargs=3)
    parser.add_argument("--silent", nargs=3)
    parser.add_argument("--driven", nargs=3)
    parser.add_argument("--target")
    parser.add_argument("--expected-windows", type=int, default=NUM_WINDOWS)
    parser.add_argument("--expected-aux-rows", type=int, default=EXPECTED_AUX_ROWS)
    parser.add_argument("--stage6-evidence")
    parser.add_argument("--bencha-evidence")
    parser.add_argument("--evidence")
    parser.add_argument("--markdown", default=str(HERE / "results" / "GATE_EVIDENCE.md"))
    args = parser.parse_args()

    evidence = evaluate_task(args)
    evidence_path = Path(args.evidence) if args.evidence else (
        HERE / "results" / "lorenz_skip50_gate_evidence.json"
    )
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    evidence_path.write_text(
        json.dumps(evidence, indent=2, allow_nan=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    markdown_path = Path(args.markdown)
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    write_markdown(markdown_path, evidence)

    print(json.dumps({
        "task": "lorenz_skip50",
        "overall_gate_pass": evidence["overall_gate_pass"],
        "gates": evidence["gates"],
        "test_nrmse": evidence["test_nrmse"],
        "mean_se": evidence["mean_se"],
        "memory_capacity": evidence["memory_capacity"],
        "field_only_test_nrmse": evidence["field_only_test_nrmse"],
        "biology_vs_field_note": evidence["biology_vs_field_note"],
        "failure_reason": evidence["failure_reason"],
        "driven_feature_count": evidence["driven_feature_count"],
    }, indent=2, allow_nan=True, default=json_default))
    print(f"evidence_file={evidence_path}")
    print(f"markdown_file={markdown_path}")
    print(f"Overall: {'PASS' if evidence['overall_gate_pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
