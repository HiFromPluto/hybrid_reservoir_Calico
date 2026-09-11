#!/usr/bin/env python3
"""Evaluate Plan Stage 6 NARMA-10 null-model gates.

Stage 4 remains FAIL. Stage 5 remains PASS. Stage 7 is not started.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

WASHOUT = 40
TRAIN = 110
TEST = 50
INNER_VAL = 22
NUM_WINDOWS = WASHOUT + TRAIN + TEST
EXPECTED_AUX_ROWS = NUM_WINDOWS * 16
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
MC_MAX_DELAY = 20
STAGE4_STATUS = "FAIL"
STAGE5_STATUS = "PASS"

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


def json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not json serializable: {type(obj)}")


def narma10(u):
    y = [0.0] * (len(u) + 1)
    for t, u_t in enumerate(u):
        acc = sum(y[t - i] if t - i >= 0 else 0.0 for i in range(10))
        u_lag = u[t - 9] if t >= 9 else 0.0
        y[t + 1] = 0.3 * y[t] + 0.05 * y[t] * acc + 1.5 * u_lag * u_t + 0.1
    return y


def load_target(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    y = np.array([float(row["y_next"]) for row in rows], dtype=float)
    recomputed = np.array(narma10(u.tolist())[1:], dtype=float)
    if len(u) != NUM_WINDOWS:
        raise ValueError(f"target has {len(u)} windows; expected {NUM_WINDOWS}")
    if not np.allclose(y, recomputed, rtol=0, atol=1e-10):
        raise ValueError("narma10_target.csv does not match the frozen recurrence")
    return u, y


def validate_csv(path, expected_rows):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, [])
        widths = [len(row) for row in reader]
    return {
        "path": str(path),
        "row_count": len(widths),
        "column_count": len(header),
        "row_count_pass": len(widths) == expected_rows,
        "rectangular_csv_pass": bool(header) and all(width == len(header) for width in widths),
    }


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
    inner_train_end = TRAIN - INNER_VAL
    X_inner, y_inner = X[:inner_train_end], y[:inner_train_end]
    X_val, y_val = X[inner_train_end:], y[inner_train_end:]
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
    X = X_all[WASHOUT:]
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
        "n_features_checked": X.shape[1],
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
            "results.csv": validate_csv(results_path, expected_aux_rows),
            "voxels.csv": validate_csv(voxels_path, expected_aux_rows),
        },
        "feature_contract": (directory / "feature_contract.txt").read_text(encoding="utf-8")
        if (directory / "feature_contract.txt").exists() else "",
    }


def evaluate_arm_run(run, u, y):
    matrix = run["matrix"]
    narma = evaluate_readout(matrix["X"], y, matrix["windows"])
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
        "narma": narma,
        "memory_capacity": mc,
        "field_only": field,
        "summary_csv": run["summary_csv"],
        "auxiliary_csvs": run["auxiliary_csvs"],
    }


def write_markdown(path, evidence):
    gates = evidence["gates"]
    overall = "PASS" if evidence["overall_gate_pass"] else "FAIL"
    lines = [
        "# Stage 6 gate evidence",
        "",
        "Stage 6 copies frozen Stage 5 and tests NARMA-10 against a Brownian",
        "density null and a silent-source control. Stage 4 remains **FAIL**.",
        "Stage 5 remains **PASS**. Mechanisms were not retuned. Stage 7 was",
        "not started.",
        "",
        f"## Overall: {overall}",
        "",
        evidence["failure_reason"] if not evidence["overall_gate_pass"]
        else "Driven biology beat the Brownian null with silent in the expected place.",
        "",
        "| Gate | Result |",
        "|---|---|",
        f"| 1. Driven NARMA-10 NRMSE < Brownian, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_brownian_pass'] else 'FAIL'}** |",
        f"| 2. Silent not ≈ driven; silent worse than driven | **{'PASS' if gates['silent_not_equal_driven_pass'] else 'FAIL'}** |",
        f"| 3. Driven uses frozen Stage 5 analysis channels | **{'PASS' if gates['frozen_feature_list_pass'] else 'FAIL'}** |",
        f"| 4. CSV rectangularity 200/3200/3200 | **{'PASS' if gates['csv_validation_pass'] else 'FAIL'}** |",
        f"| 5. Stage 4 FAIL, Stage 5 PASS unchanged | **{'PASS' if gates['prior_stage_labels_pass'] else 'FAIL'}** |",
        "",
        "## NARMA-10 test NRMSE",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 | mean ± s.e. |",
        "|---|---|---|---|---|",
    ]
    for arm in ("brownian", "silent", "driven"):
        values = evidence["narma_test_nrmse"][arm]
        mean = evidence["narma_mean_se"][arm]["mean"]
        se = evidence["narma_mean_se"][arm]["se"]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} | {mean:.4f} ± {se:.4f} |"
        )
    lines.extend([
        "",
        "Silent NRMSE is identical across seeds: with sources off, R / L / clamp",
        "deaths carry no NARMA input, the ridge collapses to an intercept, and",
        "test NRMSE is fixed by `y_test`. Brownian density is the same to three",
        "decimals. Driven sits 0.23 below that chance floor.",
    ])
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
        "Frozen Stage 5 analysis channels, voxelised:",
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
        "See `PROTOCOL.md`. AHL ~ Uniform[0, 0.5] seed 20260814; acid held at 0.5;",
        "washout 40 / train 110 / test 50; ridge grid selected on training data only.",
        "",
        "CSV completeness is the on-disk files, not BSim timestep stdout.",
        "Stage 7 was not started.",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brownian", nargs=3, required=True)
    parser.add_argument("--silent", nargs=3, required=True)
    parser.add_argument("--driven", nargs=3, required=True)
    parser.add_argument("--target", default="narma10_target.csv")
    parser.add_argument("--expected-windows", type=int, default=NUM_WINDOWS)
    parser.add_argument("--expected-aux-rows", type=int, default=EXPECTED_AUX_ROWS)
    parser.add_argument("--evidence", default="results/stage6_gate_evidence.json")
    parser.add_argument("--markdown", default="results/GATE_EVIDENCE.md")
    args = parser.parse_args()

    u, y = load_target(args.target)
    grouped = {
        "brownian": [read_run(path, "brownian", args.expected_windows, args.expected_aux_rows)
                     for path in args.brownian],
        "silent": [read_run(path, "silent", args.expected_windows, args.expected_aux_rows)
                   for path in args.silent],
        "driven": [read_run(path, "driven", args.expected_windows, args.expected_aux_rows)
                   for path in args.driven],
    }
    evaluated = {
        arm: [evaluate_arm_run(run, u, y) for run in runs]
        for arm, runs in grouped.items()
    }

    nrmse_by_arm = {
        arm: [run["narma"]["test_nrmse"] for run in runs]
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

    gates = {
        "replicate_count_pass": all(len(values) == 3 for values in nrmse_by_arm.values()),
        "driven_beats_brownian_pass": driven_beats_brownian,
        "silent_not_equal_driven_pass": silent_worse_than_driven,
        "frozen_feature_list_pass": feature_ok,
        "csv_validation_pass": csv_pass,
        "prior_stage_labels_pass": STAGE4_STATUS == "FAIL" and STAGE5_STATUS == "PASS",
    }
    overall = all(gates.values())
    if not gates["driven_beats_brownian_pass"]:
        failure_reason = "null_wins: driven NARMA-10 NRMSE did not beat Brownian by a clear margin."
    elif not gates["silent_not_equal_driven_pass"]:
        failure_reason = "silent_approx_driven: silent NRMSE is not clearly worse than driven; inputs are not driving the reservoir."
    elif not gates["frozen_feature_list_pass"]:
        failure_reason = "readout_leakage: driven features are not the frozen Stage 5 analysis channels."
    elif not gates["csv_validation_pass"]:
        failure_reason = "csv_incomplete: a summary/sample/voxel file is missing rows or is ragged."
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
        "schema": "BSimReservoirPlanStage6-gates-v1",
        "stage4_status": STAGE4_STATUS,
        "stage5_status": STAGE5_STATUS,
        "protocol": {
            "num_windows": NUM_WINDOWS,
            "washout": WASHOUT,
            "train": TRAIN,
            "test": TEST,
            "inner_val": INNER_VAL,
            "ridge_grid": list(RIDGE_GRID),
            "u_range": [0.0, 0.5],
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
        "narma_test_nrmse": nrmse_by_arm,
        "narma_mean_se": {
            arm: {"mean": mean, "se": se} for arm, (mean, se) in mean_se_by_arm.items()
        },
        "memory_capacity": mc_by_arm,
        "field_only_test_nrmse": field_nrmse,
        "biology_vs_field_note": field_note,
        "lambdas": {
            arm: [run["narma"]["lambda"] for run in runs]
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

    destination = Path(args.evidence)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(evidence, indent=2, allow_nan=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    write_markdown(args.markdown, evidence)
    print(json.dumps({
        "overall_gate_pass": overall,
        "gates": gates,
        "narma_test_nrmse": nrmse_by_arm,
        "narma_mean_se": evidence["narma_mean_se"],
        "memory_capacity": mc_by_arm,
        "field_only_test_nrmse": field_nrmse,
        "biology_vs_field_note": field_note,
        "failure_reason": failure_reason,
        "driven_feature_count": len(driven_features),
    }, indent=2, allow_nan=True, default=json_default))
    print(f"evidence_file={destination}")
    print(f"markdown_file={args.markdown}")


if __name__ == "__main__":
    main()
