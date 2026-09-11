#!/usr/bin/env python3
"""Evaluate Plan Stage 9 classification gates on frozen Stage 6 CSVs.

Stage 4 remains FAIL. Stage 5 remains PASS. Stage 6 remains PASS.
Stage 7 remains FAIL. Stage 8 is skipped. Do not modify Stage 6.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

WASHOUT = 40
TRAIN = 110
TEST = 50
INNER_VAL = 22
NUM_WINDOWS = WASHOUT + TRAIN + TEST
EXPECTED_AUX_ROWS = NUM_WINDOWS * 16
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
PRODUCT_THRESHOLD = 0.0625
AHL_SHA256 = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
LAST_SAMPLE = ("199", "15", "299.95")
STAGE4_STATUS = "FAIL"
STAGE5_STATUS = "PASS"
STAGE6_STATUS = "PASS"
STAGE7_STATUS = "FAIL"
STAGE8_STATUS = "SKIPPED"

BIOLOGY_MEAN_PREFIXES = ("Receiver_R_", "Lum_Mean_")
BIOLOGY_LAST_PREFIXES = ("Input_Driven_Death_",)
BROWNIAN_MEAN_PREFIXES = ("Den_",)
BROWNIAN_LAST_PREFIXES = ()
FIELD_MEAN_PREFIXES = ("AHL_uM_",)
EXCLUDED_FROM_RIDGE = (
    "window AHL",
    "occupancy",
    "pH",
    "Births",
    "Total_Deaths",
    "Clamp",
    "OOB",
    "Population",
    "Lum_Sum",
    "voxel pH / Den (biology arms)",
    "Att_*",
    "raw u[n] in the field baseline",
)


def json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not json serializable: {type(obj)}")


def load_ahl(path):
    values = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(float(line))
    digest = hashlib.sha256(
        ",".join(f"{value:.12f}" for value in values).encode()
    ).hexdigest()
    if digest != AHL_SHA256:
        raise ValueError(f"AHL SHA-256 {digest} does not match frozen {AHL_SHA256}")
    if len(values) != NUM_WINDOWS:
        raise ValueError(f"AHL sequence has {len(values)} windows; expected {NUM_WINDOWS}")
    return np.array(values, dtype=float)


def product_bit(u, threshold=PRODUCT_THRESHOLD):
    y = np.zeros(len(u), dtype=int)
    for n, u_n in enumerate(u):
        prev = u[n - 1] if n else 0.0
        y[n] = 1 if u_n * prev > threshold else 0
    return y


def class_balance(y, start, end, name):
    slice_y = y[start:end]
    positives = int(np.sum(slice_y))
    return {
        "split": name,
        "windows": f"{start}..{end - 1}",
        "n": int(len(slice_y)),
        "positives": positives,
        "fraction": float(positives / len(slice_y)) if len(slice_y) else float("nan"),
    }


def delay_line_matrix(u):
    X = np.zeros((len(u), 2), dtype=float)
    for n, u_n in enumerate(u):
        X[n, 0] = u_n
        X[n, 1] = u[n - 1] if n else 0.0
    return X


def validate_csv(path, expected_rows):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, [])
        rows = list(reader)
    widths = [len(row) for row in rows]
    last = rows[-1] if rows else []
    last_ok = (
        len(last) >= 3
        and last[0] == LAST_SAMPLE[0]
        and last[1] == LAST_SAMPLE[1]
        and last[2] == LAST_SAMPLE[2]
    )
    return {
        "path": str(path),
        "row_count": len(widths),
        "column_count": len(header),
        "row_count_pass": len(widths) == expected_rows,
        "rectangular_csv_pass": bool(header) and all(width == len(header) for width in widths),
        "last_sample": f"{last[0]};{last[1]};{last[2]}" if len(last) >= 3 else "",
        "last_sample_pass": last_ok if Path(path).name == "results.csv" else True,
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


def roc_auc(y_true, scores):
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    n_pos = int(np.sum(y_true == 1))
    n_neg = int(np.sum(y_true == 0))
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    sorted_scores = scores[order]
    ranks = np.empty(len(scores), dtype=float)
    i = 0
    n = len(scores)
    while i < n:
        j = i
        while j + 1 < n and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1
    pos_rank_sum = float(np.sum(ranks[y_true == 1]))
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def accuracy_at_half(y_true, scores):
    pred = (np.asarray(scores) >= 0.5).astype(int)
    return float(np.mean(pred == np.asarray(y_true)))


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
    """Fit 40..127, pick lambda on 128..149 only. Do not copy Stage 6's X[end:]."""
    inner_train_end = TRAIN - INNER_VAL
    X_inner, y_inner = X[:inner_train_end], y[:inner_train_end]
    X_val, y_val = X[inner_train_end:TRAIN], y[inner_train_end:TRAIN]
    (X_inner_z, X_val_z), _, _ = standardize(X_inner, X_val)
    scored = []
    grid_scores = []
    for lam in RIDGE_GRID:
        weights = ridge_fit(X_inner_z, y_inner, lam)
        pred = ridge_predict(X_val_z, weights)
        auc = roc_auc(y_val, pred)
        scored.append((-auc, -lam, lam, auc))
        grid_scores.append({"lambda": lam, "val_auc": auc})
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
        "train_auc": roc_auc(y_train, train_pred),
        "test_auc": roc_auc(y_test, test_pred),
        "test_accuracy_at_0_5": accuracy_at_half(y_test, test_pred),
    }


def evaluate_readout(X_all, target, windows):
    if windows != list(range(NUM_WINDOWS)):
        raise ValueError(
            f"expected windows 0..{NUM_WINDOWS - 1}, got {windows[:3]}..{windows[-1:]}"
        )
    X = X_all[WASHOUT:]
    y = target[WASHOUT:]
    lam, grid_scores = select_lambda(X, y)
    metrics = fit_eval(X, y, lam)
    metrics["lambda_grid"] = grid_scores
    return metrics


def mean_se(values):
    array = np.array(values, dtype=float)
    mean = float(np.mean(array))
    se = float(np.std(array, ddof=1) / np.sqrt(len(array))) if len(array) > 1 else float("nan")
    return mean, se


def intervals_overlap(mean_a, se_a, mean_b, se_b):
    return not (mean_a + se_a < mean_b - se_b or mean_b + se_b < mean_a - se_a)


def overlaps_chance(mean, se, chance=0.5):
    return mean - se <= chance <= mean + se


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
    }


def evaluate_arm_run(run, y):
    matrix = run["matrix"]
    cls = evaluate_readout(matrix["X"], y, matrix["windows"])
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
        "classification": cls,
        "field_only": field,
        "summary_csv": run["summary_csv"],
        "auxiliary_csvs": run["auxiliary_csvs"],
    }


def write_markdown(path, evidence):
    gates = evidence["gates"]
    overall = "PASS" if evidence["overall_gate_pass"] else "FAIL"
    lines = [
        "# Stage 9 gate evidence",
        "",
        "Stage 9 evaluates a frozen product-bit classification task on the",
        "existing Stage 6 CSVs. BSim was not rerun. Stage 4 remains **FAIL**.",
        "Stage 5 remains **PASS**. Stage 6 remains **PASS**. Stage 7 remains",
        "**FAIL**. Stage 8 is **skipped** (Danino already archived as Stage 3).",
        "Kinetics were not retuned. Two-way AC↔bacteria coupling was not started.",
        "",
        f"## Overall: {overall}",
        "",
        evidence["failure_reason"] if not evidence["overall_gate_pass"]
        else "Driven biology beat Brownian and the field-only linear baseline; silent sat at chance.",
        "",
        evidence["failure_detail"],
        "",
        "| Gate | Result |",
        "|---|---|",
        f"| 1. Driven test AUC > Brownian, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_brownian_pass'] else 'FAIL'}** |",
        f"| 2. Driven test AUC > field-only linear baseline, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_field_pass'] else 'FAIL'}** |",
        f"| 3. Silent not above chance by a clear margin; silent not ≈ driven | **{'PASS' if gates['silent_chance_pass'] else 'FAIL'}** |",
        f"| 4. Driven uses frozen 408 Stage 5/6 features | **{'PASS' if gates['frozen_feature_list_pass'] else 'FAIL'}** |",
        f"| 5. CSV rectangularity 200/3200/3200; last sample 199;15;299.95 | **{'PASS' if gates['csv_validation_pass'] else 'FAIL'}** |",
        f"| 6. Stage 4 FAIL, 5 PASS, 6 PASS, 7 FAIL; Stage 8 skipped | **{'PASS' if gates['prior_stage_labels_pass'] else 'FAIL'}** |",
        "",
        "## Frozen task",
        "",
        "`y[n] = 1` if `u[n] * u[n-1] > 0.0625` else `0`, with `u[-1] = 0`.",
        "AHL SHA-256 `d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.",
        "Threshold was not changed after seeing AUC.",
        "",
        "| Split | Windows | n | Positives | Fraction |",
        "|---|---|---|---|---|",
    ]
    for row in evidence["class_balance"]:
        lines.append(
            f"| {row['split']} | {row['windows']} | {row['n']} | "
            f"{row['positives']} | {row['fraction']:.3f} |"
        )
    lines.extend([
        "",
        "## Test AUC (primary)",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 | mean ± s.e. |",
        "|---|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "driven", "field"):
        values = evidence["test_auc"][arm]
        mean = evidence["auc_mean_se"][arm]["mean"]
        se = evidence["auc_mean_se"][arm]["se"]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} | {mean:.4f} ± {se:.4f} |"
        )
    lines.extend([
        "",
        "Per-seed driven vs field (gate 2). Driven must be higher in every seed.",
        "",
        "| Seed | Driven AUC | Field AUC | Driven > field |",
        "|---|---|---|---|",
    ])
    for seed, driven_auc, field_auc_seed in zip(
        (101, 202, 303),
        evidence["test_auc"]["driven"],
        evidence["test_auc"]["field"],
    ):
        lines.append(
            f"| {seed} | {driven_auc:.4f} | {field_auc_seed:.4f} | "
            f"{'yes' if driven_auc > field_auc_seed else 'no'} |"
        )
    lines.extend([
        "",
        "Accuracy at 0.5 is printed, not the gate.",
        "",
        "| Arm | seed 101 acc | seed 202 acc | seed 303 acc |",
        "|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "driven", "field"):
        values = evidence["test_accuracy"][arm]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} |"
        )
    delay = evidence["delay_line_diagnostic"]
    lines.extend([
        "",
        "## Delay-line diagnostic (not a gate)",
        "",
        "Ridge on `(u[n], u[n-1])` only. This knows the stimulus without a dish.",
        "It is not a substitute for the field-only baseline.",
        "",
        f"Delay-line test AUC `{delay['test_auc']:.4f}`, accuracy at 0.5 "
        f"`{delay['test_accuracy_at_0_5']:.4f}`, lambda `{delay['lambda']}`.",
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
        "Excluded from any ridge: " + ", ".join(EXCLUDED_FROM_RIDGE) + ".",
        "",
        "## Lambdas (independent per arm × seed)",
        "",
        "| Arm | 101 | 202 | 303 |",
        "|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "driven", "field"):
        values = evidence["lambdas"][arm]
        lines.append(f"| {arm} | {values[0]:g} | {values[1]:g} | {values[2]:g} |")
    lines.extend([
        "",
        "## Protocol",
        "",
        "See `PROTOCOL.md`. Classification on frozen Stage 6 CSVs. Washout 40 /",
        "train 110 / test 50. Lambda fit 40..127, pick on 128..149 only",
        "(`X[inner_train_end:TRAIN]`), refit 40..149. Stage 6's",
        "`X[inner_train_end:]` slice is not copied.",
        "",
        "CSV completeness is the on-disk Stage 6 files, last sample `199;15;299.95`.",
        "Stage 8 was skipped. Two-way coupling was not started.",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brownian", nargs=3, required=True)
    parser.add_argument("--silent", nargs=3, required=True)
    parser.add_argument("--driven", nargs=3, required=True)
    parser.add_argument("--ahl", default="../BSimReservoirPlanStage6/input_ahl_narma200.txt")
    parser.add_argument("--expected-windows", type=int, default=NUM_WINDOWS)
    parser.add_argument("--expected-aux-rows", type=int, default=EXPECTED_AUX_ROWS)
    parser.add_argument("--evidence", default="results/stage9_gate_evidence.json")
    parser.add_argument("--markdown", default="results/GATE_EVIDENCE.md")
    args = parser.parse_args()

    u = load_ahl(args.ahl)
    y = product_bit(u)
    balances = [
        class_balance(y, 0, WASHOUT, "washout"),
        class_balance(y, WASHOUT, WASHOUT + TRAIN, "train"),
        class_balance(y, WASHOUT + TRAIN, NUM_WINDOWS, "test"),
        class_balance(y, 0, NUM_WINDOWS, "all"),
    ]
    print("CLASS_BALANCE")
    for row in balances:
        print(
            f"{row['split']} windows={row['windows']} n={row['n']} "
            f"pos={row['positives']} frac={row['fraction']:.6f}"
        )

    grouped = {
        "brownian": [read_run(path, "brownian", args.expected_windows, args.expected_aux_rows)
                     for path in args.brownian],
        "silent": [read_run(path, "silent", args.expected_windows, args.expected_aux_rows)
                   for path in args.silent],
        "driven": [read_run(path, "driven", args.expected_windows, args.expected_aux_rows)
                   for path in args.driven],
    }
    evaluated = {
        arm: [evaluate_arm_run(run, y) for run in runs]
        for arm, runs in grouped.items()
    }

    auc_by_arm = {
        arm: [run["classification"]["test_auc"] for run in runs]
        for arm, runs in evaluated.items()
    }
    acc_by_arm = {
        arm: [run["classification"]["test_accuracy_at_0_5"] for run in runs]
        for arm, runs in evaluated.items()
    }
    field_auc = [run["field_only"]["test_auc"] for run in evaluated["driven"] if run["field_only"]]
    field_acc = [
        run["field_only"]["test_accuracy_at_0_5"] for run in evaluated["driven"] if run["field_only"]
    ]
    auc_by_arm["field"] = field_auc
    acc_by_arm["field"] = field_acc
    mean_se_by_arm = {arm: mean_se(values) for arm, values in auc_by_arm.items()}

    delay = evaluate_readout(delay_line_matrix(u), y, list(range(NUM_WINDOWS)))
    delay["feature_names"] = ["u[n]", "u[n-1]"]

    driven_features = evaluated["driven"][0]["feature_names"]
    expected_prefixes = ("Receiver_R_", "Lum_Mean_", "Input_Driven_Death_")
    forbidden = ("AHL_uM_", "pH_", "Den_", "Fraction_q", "Lum_Sum_", "Birth_",
                 "Clamp_Death_", "OOB_Death_", "Att_")
    feature_ok = (
        all(name.startswith(expected_prefixes) for name in driven_features)
        and not any(name.startswith(forbidden) for name in driven_features)
        and len(driven_features) == 408
        and all(run["n_features"] == 408 for run in evaluated["driven"] + evaluated["silent"])
        and all(run["n_features"] == 200 for run in evaluated["brownian"])
        and all(run["field_only"]["n_features"] == 200 for run in evaluated["driven"])
    )
    print("DRIVEN_FEATURE_LIST")
    for name in driven_features:
        print(name)
    print(f"driven_feature_count={len(driven_features)}")

    csv_pass = all(
        run["summary_csv"]["row_count_pass"]
        and run["summary_csv"]["rectangular_csv_pass"]
        and all(
            item["row_count_pass"]
            and item["rectangular_csv_pass"]
            and item["last_sample_pass"]
            for item in run["auxiliary_csvs"].values()
        )
        for runs in evaluated.values() for run in runs
    )

    driven_mean, driven_se = mean_se_by_arm["driven"]
    brown_mean, brown_se = mean_se_by_arm["brownian"]
    silent_mean, silent_se = mean_se_by_arm["silent"]
    field_mean, field_se = mean_se_by_arm["field"]

    driven_beats_brownian = (
        all(d > b for d, b in zip(auc_by_arm["driven"], auc_by_arm["brownian"]))
        and driven_mean - driven_se > brown_mean + brown_se
    )
    driven_beats_field = (
        all(d > f for d, f in zip(auc_by_arm["driven"], auc_by_arm["field"]))
        and driven_mean - driven_se > field_mean + field_se
    )
    silent_overlap_driven = intervals_overlap(silent_mean, silent_se, driven_mean, driven_se)
    silent_worse_than_driven = silent_mean + silent_se < driven_mean - driven_se
    silent_at_chance = overlaps_chance(silent_mean, silent_se)
    silent_chance_pass = (silent_at_chance or silent_worse_than_driven) and not silent_overlap_driven

    prior_ok = (
        STAGE4_STATUS == "FAIL"
        and STAGE5_STATUS == "PASS"
        and STAGE6_STATUS == "PASS"
        and STAGE7_STATUS == "FAIL"
        and STAGE8_STATUS == "SKIPPED"
    )

    gates = {
        "replicate_count_pass": all(len(values) == 3 for values in auc_by_arm.values()),
        "driven_beats_brownian_pass": driven_beats_brownian,
        "driven_beats_field_pass": driven_beats_field,
        "silent_chance_pass": silent_chance_pass,
        "frozen_feature_list_pass": feature_ok,
        "csv_validation_pass": csv_pass,
        "prior_stage_labels_pass": prior_ok,
    }
    overall = all(gates.values())
    if not gates["driven_beats_brownian_pass"]:
        failure_reason = (
            "brownian_wins: driven test AUC did not beat Brownian density in every "
            "seed with non-overlapping mean±s.e. The hybrid did not outperform "
            "passive occupancy."
        )
    elif not gates["driven_beats_field_pass"]:
        failure_reason = (
            "field_wins: driven test AUC did not beat the field-only AHL_uM_* "
            "baseline in every seed with non-overlapping mean±s.e. The AC plume "
            "is doing the work, not the hybrid medium."
        )
    elif not gates["silent_chance_pass"]:
        failure_reason = (
            "silent_approx_driven: silent test AUC is not at chance and overlaps "
            "driven. Inputs are not what the classifier is using."
        )
    elif not gates["frozen_feature_list_pass"]:
        failure_reason = (
            "readout_leakage: driven features are not the frozen Stage 5/6 408-list."
        )
    elif not gates["csv_validation_pass"]:
        failure_reason = (
            "csv_incomplete: a Stage 6 summary/sample/voxel file is missing rows, "
            "is ragged, or last sample is not 199;15;299.95."
        )
    elif not gates["prior_stage_labels_pass"]:
        failure_reason = "prior_stage_labels: Stage 4/5/6/7/8 labels were reopened."
    else:
        failure_reason = ""

    if not driven_beats_field:
        failure_detail = (
            "Gate 2 fails. Driven test AUC mean ± s.e. is "
            f"{driven_mean:.4f} ± {driven_se:.4f}; field-only AHL_uM_* is "
            f"{field_mean:.4f} ± {field_se:.4f}. "
            f"Per seed: driven {auc_by_arm['driven'][0]:.4f} / "
            f"{auc_by_arm['driven'][1]:.4f} / {auc_by_arm['driven'][2]:.4f} "
            f"vs field {auc_by_arm['field'][0]:.4f} / "
            f"{auc_by_arm['field'][1]:.4f} / {auc_by_arm['field'][2]:.4f}. "
            "The intervals overlap, and driven is not higher in every seed. "
            "The AC carrier plume already ranks the product bit; the hybrid "
            "408-feature readout does not add a clear margin. Kinetics, "
            "features, AC sites, vesicles, Danino, and two-way coupling were "
            "not changed after seeing this."
        )
    elif not driven_beats_brownian:
        failure_detail = (
            "Gate 1 fails. Driven did not beat Brownian density with "
            "non-overlapping mean±s.e. in every seed."
        )
    elif not silent_chance_pass:
        failure_detail = (
            "Gate 3 fails. Silent test AUC is not at chance and overlaps driven."
        )
    elif overall:
        failure_detail = (
            "No failure. Driven beat Brownian and the field baseline; silent "
            "sat at chance. Frozen 408-feature contract and Stage 6 CSVs used."
        )
    else:
        failure_detail = failure_reason

    evidence = {
        "schema": "BSimReservoirPlanStage9-gates-v1",
        "stage4_status": STAGE4_STATUS,
        "stage5_status": STAGE5_STATUS,
        "stage6_status": STAGE6_STATUS,
        "stage7_status": STAGE7_STATUS,
        "stage8_status": STAGE8_STATUS,
        "evaluated_on": "examples/BSimReservoirPlanStage6/results (no BSim rerun)",
        "protocol": {
            "num_windows": NUM_WINDOWS,
            "washout": WASHOUT,
            "train": TRAIN,
            "test": TEST,
            "inner_val": INNER_VAL,
            "lambda_validation_slice": "X[inner_train_end:TRAIN] after washout strip",
            "ridge_grid": list(RIDGE_GRID),
            "u_range": [0.0, 0.5],
            "ahl_sha256": AHL_SHA256,
            "task": "y[n]=1 if u[n]*u[n-1]>0.0625 else 0; u[-1]=0",
            "product_threshold": PRODUCT_THRESHOLD,
            "primary_metric": "test AUC",
        },
        "class_balance": balances,
        "analysis_channels_driven": [
            "Receiver_R_* / Mean_q (20x10)",
            "Lum_Mean_* / Mean_L (20x10)",
            "Input_Driven_Death_* (4x2)",
        ],
        "driven_feature_count": len(driven_features),
        "excluded_from_ridge": list(EXCLUDED_FROM_RIDGE),
        "test_auc": auc_by_arm,
        "test_accuracy": acc_by_arm,
        "auc_mean_se": {
            arm: {"mean": mean, "se": se} for arm, (mean, se) in mean_se_by_arm.items()
        },
        "delay_line_diagnostic": delay,
        "lambdas": {
            "brownian": [run["classification"]["lambda"] for run in evaluated["brownian"]],
            "silent": [run["classification"]["lambda"] for run in evaluated["silent"]],
            "driven": [run["classification"]["lambda"] for run in evaluated["driven"]],
            "field": [run["field_only"]["lambda"] for run in evaluated["driven"]],
            "delay_line": delay["lambda"],
        },
        "files": {
            arm: [{
                "path": run["path"],
                "n_features": run["n_features"],
                "classification": run["classification"],
                "field_only": (
                    {k: v for k, v in run["field_only"].items() if k != "feature_names"}
                    if run["field_only"] else None
                ),
                "summary_csv": run["summary_csv"],
                "auxiliary_csvs": run["auxiliary_csvs"],
            } for run in runs]
            for arm, runs in evaluated.items()
        },
        "gates": gates,
        "overall_gate_pass": overall,
        "failure_reason": failure_reason,
        "failure_detail": failure_detail,
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
        "test_auc": auc_by_arm,
        "auc_mean_se": evidence["auc_mean_se"],
        "test_accuracy": acc_by_arm,
        "delay_line_test_auc": delay["test_auc"],
        "failure_reason": failure_reason,
        "failure_detail": failure_detail,
        "driven_feature_count": len(driven_features),
    }, indent=2, allow_nan=True, default=json_default))
    print(f"evidence_file={destination}")
    print(f"markdown_file={args.markdown}")


if __name__ == "__main__":
    main()
