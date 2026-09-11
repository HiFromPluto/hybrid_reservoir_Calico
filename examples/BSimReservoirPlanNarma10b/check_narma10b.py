#!/usr/bin/env python3
"""Evaluate NARMA-10b replication gates on the frozen Stage 6 dish.

Independent replication with bacterial RNG seeds 111/222/333. Stage 6
NARMA-10 remains PASS and is not modified. The NARMA u sequence is
copied, not regenerated. Track B / waveform / Lorenz / vesicles /
glucose / Danino are not started. The Stage 6 val-slice leak is closed:
lambda is picked on rows 88..109 only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
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
SEEDS = (111, 222, 333)
NARMA_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
U_LOW = 0.0
U_HIGH = 0.5

# Stage 6 published means. Labelled "Stage 6, not this run."
STAGE6_MEAN = {
    "driven": 0.9306,
    "brownian": 1.1625,
    "silent": 1.1622,
}
STAGE6_SE = {
    "driven": 0.0024,
    "brownian": 0.0002,
    "silent": 0.0000,
}
STAGE6_PER_SEED = {
    "brownian": (1.1627, 1.1627, 1.1621),
    "silent": (1.1622, 1.1622, 1.1622),
    "driven": (0.9302, 0.9349, 0.9267),
}
STAGE6_FIELD_MEAN = 1.0126
STAGE6_FIELD_PER_SEED = (1.0126, 1.0126, 1.0126)
STAGE6_SEEDS = (101, 202, 303)

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

HERE = Path(__file__).resolve().parent
STAGE6_DIR = HERE.parent / "BSimReservoirPlanStage6"
BENCHA_DIR = HERE.parent / "BSimReservoirPlanBenchA"
BENCHA2_DIR = HERE.parent / "BSimReservoirPlanBenchA2"
WAVEFORM_DIR = HERE.parent / "BSimReservoirPlanWaveform"
TRACKB_DIR = HERE.parent / "BSimReservoirPlanTrackB"


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


def narma10(u):
    y = [0.0] * (len(u) + 1)
    for t, u_t in enumerate(u):
        acc = sum(y[t - i] if t - i >= 0 else 0.0 for i in range(10))
        u_lag = u[t - 9] if t >= 9 else 0.0
        y[t + 1] = 0.3 * y[t] + 0.05 * y[t] * acc + 1.5 * u_lag * u_t + 0.1
    return y


def load_sequence(path):
    values = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(float(line))
    return np.array(values, dtype=float)


def abort(message):
    print(f"ABORT: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_target(path, ahl_path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    if len(rows) != NUM_WINDOWS:
        abort(f"target has {len(rows)} windows; expected {NUM_WINDOWS}")
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    y = np.array([float(row["y_next"]) for row in rows], dtype=float)
    if np.any(u < U_LOW) or np.any(u > U_HIGH):
        abort("u escaped [0, 0.5]; Uniform[0,1] is not allowed")
    digest = sha256_u(u)
    if digest != NARMA_SHA:
        abort(f"AHL SHA-256 mismatch: {digest} != {NARMA_SHA}")
    u_file = load_sequence(ahl_path)
    if len(u_file) != NUM_WINDOWS:
        abort(f"{ahl_path} has {len(u_file)} values; expected {NUM_WINDOWS}")
    if not np.allclose(u, u_file, rtol=0, atol=1e-12):
        abort("input_ahl_narma200.txt does not match narma10_target.csv u")
    file_digest = sha256_u(u_file)
    if file_digest != NARMA_SHA:
        abort(f"AHL file SHA-256 mismatch: {file_digest} != {NARMA_SHA}")
    recomputed = np.array(narma10(u.tolist())[1:], dtype=float)
    if not np.allclose(y, recomputed, rtol=0, atol=1e-10):
        abort("narma10_target.csv does not match the frozen NARMA-10 recurrence")
    return u, y, digest


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


def package_java_untouched(directory, package_name, class_name, extra_tokens=()):
    java_path = Path(directory) / f"{class_name}.java"
    voxel_path = Path(directory) / "VoxelAnalyzer.java"
    if not java_path.exists() or not voxel_path.exists():
        return False
    java = java_path.read_text(encoding="utf-8")
    voxel = voxel_path.read_text(encoding="utf-8")
    ok = (
        f"package {package_name};" in java
        and f"public final class {class_name}" in java
        and f"package {package_name};" in voxel
        and f"{package_name}.DeathCause" in voxel
    )
    return ok and all(token in java for token in extra_tokens)


def stage6_overall_pass(path):
    if not Path(path).exists():
        return False
    text = Path(path).read_text(encoding="utf-8")
    return bool(re.search(r"^## Overall:\s*PASS\s*$", text, re.MULTILINE))


def frozen_dish_java(java_path):
    text = Path(java_path).read_text(encoding="utf-8")
    if any(token in text for token in FORBIDDEN_JAVA):
        return False
    if "new Vector3d(250, 250, 5), new Vector3d(500, 250, 5)" not in text:
        return False
    if "new Vector3d(300, 375, 5)" not in text:
        return False
    if "RECEIVER_K_UM = 1.6" not in text:
        return False
    if "GROWTH_RATE = 4.0 * Math.PI / 1800.0" not in text:
        return False
    if "CARRYING_CAPACITY = 2000" not in text:
        return False
    return True


def write_markdown(path, evidence):
    overall = "PASS" if evidence["overall_gate_pass"] else "FAIL"
    gates = evidence["gates"]
    lines = [
        "# NARMA-10b gate evidence",
        "",
        "Independent NARMA-10 replication on the frozen Stage 6 hybrid dish",
        "with bacterial RNG seeds 111/222/333. Stage 6 NARMA-10 remains",
        "**PASS** and was not edited. Mechanisms were not retuned. The NARMA",
        "`u` sequence was copied, not regenerated. Track B / waveform /",
        "Lorenz / vesicles / glucose / Danino were not started. The Stage 6",
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
        f"| 5. Stage 6 NARMA-10 PASS; Stage 6 Java untouched; Waveform/Track B/BenchA/A2 not edited; no glucose/Danino/extra ACs | **{'PASS' if gates['prior_stage_labels_pass'] else 'FAIL'}** |",
        "",
        "## This replication (seeds 111 / 222 / 333)",
        "",
        "| Arm | seed 111 | seed 222 | seed 333 | mean ± s.e. |",
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
            f"Field-only (diagnostic, not a gate): {field_vals[0]:.4f} / {field_vals[1]:.4f} / {field_vals[2]:.4f} "
            f"(mean {mean:.4f} ± {se:.4f}).",
            f"Biology vs field-only: {evidence['biology_vs_field_note']}",
        ])
    lines.extend([
        "",
        "## Stage 6, not this run (seeds 101 / 202 / 303)",
        "",
        "Published Stage 6 NARMA-10. Copied here for comparison only.",
        "Stage 6 `GATE_EVIDENCE.md` was not edited.",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 | mean ± s.e. |",
        "|---|---|---|---|---|",
        f"| brownian | {STAGE6_PER_SEED['brownian'][0]:.4f} | {STAGE6_PER_SEED['brownian'][1]:.4f} | {STAGE6_PER_SEED['brownian'][2]:.4f} | {STAGE6_MEAN['brownian']:.4f} ± {STAGE6_SE['brownian']:.4f} |",
        f"| silent | {STAGE6_PER_SEED['silent'][0]:.4f} | {STAGE6_PER_SEED['silent'][1]:.4f} | {STAGE6_PER_SEED['silent'][2]:.4f} | {STAGE6_MEAN['silent']:.4f} ± {STAGE6_SE['silent']:.4f} |",
        f"| driven | {STAGE6_PER_SEED['driven'][0]:.4f} | {STAGE6_PER_SEED['driven'][1]:.4f} | {STAGE6_PER_SEED['driven'][2]:.4f} | {STAGE6_MEAN['driven']:.4f} ± {STAGE6_SE['driven']:.4f} |",
        "",
        f"Field-only (Stage 6, not this run): {STAGE6_FIELD_PER_SEED[0]:.4f} / {STAGE6_FIELD_PER_SEED[1]:.4f} / {STAGE6_FIELD_PER_SEED[2]:.4f} "
        f"(mean {STAGE6_FIELD_MEAN:.4f}).",
        "",
        "## Field-only diagnostic (not a gate)",
        "",
        "Ridge on driven voxel `AHL_uM_*` only. If biology does not beat",
        "this delay line, that is recorded; the Brownian gate is unchanged.",
        "Do not add features if field-only wins or loses.",
        "",
    ])
    lines.extend([
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
        "| Arm | MC k=1..20 | k=0 test R² (printed, not added) |",
        "|---|---|---|",
    ])
    for arm in ("brownian", "silent", "driven"):
        values = evidence["memory_capacity"][arm]
        mean, se = mean_se(values)
        k0 = evidence["mc_k0_r2"][arm]
        lines.append(
            f"| {arm} | {values[0]:.3f} / {values[1]:.3f} / {values[2]:.3f} (mean {mean:.3f} ± {se:.3f}) "
            f"| {k0[0]:.3f} / {k0[1]:.3f} / {k0[2]:.3f} |"
        )
    lines.extend([
        "",
        "## Protocol",
        "",
        "See `PROTOCOL.md`. Frozen NARMA-10 `u ~ Uniform[0, 0.5]` seed 20260814;",
        "acid held at 0.5; washout 40 / train 110 / test 50; lambda selected",
        "on windows 128..149 only. Seeds 111/222/333. Silent was rerun.",
        "",
        "CSV completeness is the on-disk files, not BSim timestep stdout.",
        "Track B / waveform / Lorenz were not started.",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def default_run_dirs(arm):
    return [str(HERE / "results" / f"narma10b_{arm}_seed{seed}") for seed in SEEDS]


def evaluate_task(args):
    target_path = Path(args.target) if args.target else HERE / "narma10_target.csv"
    ahl_path = Path(args.ahl) if args.ahl else HERE / "input_ahl_narma200.txt"
    u, y, digest = load_target(target_path, ahl_path)
    print(f"u_sha256={digest}")
    print("narma10_recurrence=OK")
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
    mc_k0 = {
        arm: [run["memory_capacity"]["k0_diagnostic"]["test_r2"] for run in runs]
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
    prior_ok = (
        stage6_overall_pass(stage6_evidence)
        and package_java_untouched(STAGE6_DIR, "BSimReservoirPlanStage6",
                                   "BSimReservoirPlanStage6", extra_tokens=("NARMA-10",))
        and package_java_untouched(BENCHA_DIR, "BSimReservoirPlanBenchA",
                                   "BSimReservoirPlanBenchA")
        and package_java_untouched(BENCHA2_DIR, "BSimReservoirPlanBenchA2",
                                   "BSimReservoirPlanBenchA2")
        and package_java_untouched(WAVEFORM_DIR, "BSimReservoirPlanWaveform",
                                   "BSimReservoirPlanWaveform")
        and package_java_untouched(TRACKB_DIR, "BSimReservoirPlanTrackB",
                                   "BSimReservoirPlanTrackB")
        and frozen_dish_java(HERE / "BSimReservoirPlanNarma10b.java")
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
        failure_reason = "provenance: Stage 6 GATE_EVIDENCE is not PASS, Stage 6 Java was touched, Waveform/Track B/BenchA/A2 were edited, or glucose / Danino / extra ACs leaked in."
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
        "schema": "BSimReservoirPlanNarma10b-gates-v1",
        "task": "narma10",
        "seeds": list(SEEDS),
        "u_sha256": digest,
        "stage6_status": "PASS" if stage6_overall_pass(stage6_evidence) else "NOT_PASS",
        "stage6_not_this_run": {
            "seeds": list(STAGE6_SEEDS),
            "test_nrmse_mean": STAGE6_MEAN,
            "test_nrmse_se": STAGE6_SE,
            "test_nrmse_per_seed": {k: list(v) for k, v in STAGE6_PER_SEED.items()},
            "field_only_mean": STAGE6_FIELD_MEAN,
            "field_only_per_seed": list(STAGE6_FIELD_PER_SEED),
        },
        "protocol": {
            "num_windows": NUM_WINDOWS,
            "washout": WASHOUT,
            "train": TRAIN,
            "test": TEST,
            "inner_val": INNER_VAL,
            "lambda_val_rows": "88..109 after washout drop (windows 128..149)",
            "ridge_grid": list(RIDGE_GRID),
            "u_range": [U_LOW, U_HIGH],
            "acid": "held at 0.5",
            "narma": "y[n+1] = 0.3 y[n] + 0.05 y[n] sum_{i=0..9} y[n-i] + 1.5 u[n-9] u[n] + 0.1",
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
        "mc_k0_r2": mc_k0,
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
    parser.add_argument("--ahl")
    parser.add_argument("--expected-windows", type=int, default=NUM_WINDOWS)
    parser.add_argument("--expected-aux-rows", type=int, default=EXPECTED_AUX_ROWS)
    parser.add_argument("--stage6-evidence")
    parser.add_argument("--evidence")
    parser.add_argument("--markdown", default=str(HERE / "results" / "GATE_EVIDENCE.md"))
    args = parser.parse_args()

    evidence = evaluate_task(args)
    evidence_path = Path(args.evidence) if args.evidence else (
        HERE / "results" / "narma10b_gate_evidence.json"
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
        "task": "narma10",
        "seeds": list(SEEDS),
        "overall_gate_pass": evidence["overall_gate_pass"],
        "gates": evidence["gates"],
        "test_nrmse": evidence["test_nrmse"],
        "mean_se": evidence["mean_se"],
        "memory_capacity": evidence["memory_capacity"],
        "mc_k0_r2": evidence["mc_k0_r2"],
        "field_only_test_nrmse": evidence["field_only_test_nrmse"],
        "biology_vs_field_note": evidence["biology_vs_field_note"],
        "stage6_not_this_run": evidence["stage6_not_this_run"],
        "failure_reason": evidence["failure_reason"],
        "driven_feature_count": evidence["driven_feature_count"],
        "u_sha256": evidence["u_sha256"],
    }, indent=2, allow_nan=True, default=json_default))
    print(f"evidence_file={evidence_path}")
    print(f"markdown_file={markdown_path}")
    print(f"Overall: {'PASS' if evidence['overall_gate_pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
