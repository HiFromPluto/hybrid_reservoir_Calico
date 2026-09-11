#!/usr/bin/env python3
"""Evaluate waveform-classification gates on the frozen Stage 6 dish.

Stage 6 NARMA-10 remains PASS. BenchA Mackey-Glass remains PASS.
BenchA Lorenz remains FAIL. Track B remains FAIL. Kinetics are not
retuned. Do not switch the gate to accuracy after seeing numbers.
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
WINDOWS_PER_BLOCK = 5
NUM_BLOCKS = 40
N_CLASSES = 3
CLASS_NAMES = ("sine", "square", "triangle")
CHANCE_ACC = 1.0 / 3.0
CHANCE_AUC = 0.5
EXPECTED_AUX_ROWS = NUM_WINDOWS * 16
LAST_SAMPLE = ("199", "15", "299.95")
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
SEEDS = (101, 202, 303)
TRIANGLE = (0.05, 0.25, 0.45, 0.25, 0.05)

CLASS_VECTOR = (
    "0,1,0,1,2,1,2,0,2,0,0,1,0,1,1,1,2,1,2,2,0,1,0,0,2,2,0,0,2,1,1,2,2,1,2,0,2,0,0,1"
)
CLASS_SHA = "c2155819f0fb512b770db53d44ef89e2577de3025f81141166278a8903446c29"
U_SHA = "c97a7dcb92c663892d16d50bfe342003b021b173211d9246b3c86c1e9734e0f8"

BIOLOGY_MEAN_PREFIXES = ("Receiver_R_", "Lum_Mean_")
BIOLOGY_LAST_PREFIXES = ("Input_Driven_Death_",)
BROWNIAN_MEAN_PREFIXES = ("Den_",)
BROWNIAN_LAST_PREFIXES = ()
FIELD_MEAN_PREFIXES = ("AHL_uM_",)
BIOLOGY_FORBIDDEN = (
    "AHL_uM_", "pH_", "Den_", "Att_", "Rep_", "Fraction_q", "Lum_Sum_",
    "Birth_", "Clamp_Death_", "OOB_Death_",
)
EXCLUDED_FROM_RIDGE = (
    "voxel AHL_uM (biology ridge)",
    "voxel Den (biology ridge)",
    "voxel pH",
    "Att",
    "raw u",
    "occupancy Fraction_q_gt_0_5",
    "Births",
    "Total_Deaths",
    "Clamp_Deaths",
    "OOB",
    "Population",
    "Lum_Sum",
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
TRACKB_DIR = HERE.parent / "BSimReservoirPlanTrackB"
STAGE9_DIR = HERE.parent / "BSimReservoirPlanStage9"


def json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not json serializable: {type(obj)}")


def sha256_u(values):
    payload = ",".join(f"{value:.12f}" for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def clip_u(value):
    return min(0.5, max(0.0, value))


def template_u(cls, w):
    if cls == 0:
        return clip_u(0.25 + 0.25 * math.sin(2.0 * math.pi * w / 5.0))
    if cls == 1:
        return clip_u(0.45 if w < 3 else 0.05)
    if cls == 2:
        return clip_u(TRIANGLE[w])
    raise ValueError(f"unknown class {cls}")


def load_sequence(path):
    values = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(float(line))
    if len(values) != NUM_WINDOWS:
        raise ValueError(f"{path} has {len(values)} values; expected {NUM_WINDOWS}")
    return values


def load_labels(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    if len(rows) != NUM_WINDOWS:
        raise ValueError(f"labels have {len(rows)} rows; expected {NUM_WINDOWS}")
    y = np.array([int(row["y"]) for row in rows], dtype=int)
    block = np.array([int(row["block"]) for row in rows], dtype=int)
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    return y, block, u, rows


def class_balance(y, start, end, name):
    slice_y = y[start:end]
    counts = [int(np.sum(slice_y == k)) for k in range(N_CLASSES)]
    return {
        "split": name,
        "windows": f"{start}..{end - 1}",
        "n": int(len(slice_y)),
        "sine": counts[0],
        "square": counts[1],
        "triangle": counts[2],
    }


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
        "last_sample_pass": last_ok if Path(path).name in {"results.csv", "voxels.csv"} else True,
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


def macro_ovr_auc(y_true, scores):
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    per_class = []
    for k in range(N_CLASSES):
        per_class.append(roc_auc((y_true == k).astype(int), scores[:, k]))
    return float(np.mean(per_class)), per_class


def confusion_counts(y_true, pred):
    matrix = np.zeros((N_CLASSES, N_CLASSES), dtype=int)
    for true, predicted in zip(y_true, pred):
        matrix[int(true), int(predicted)] += 1
    return matrix.tolist()


def accuracy_argmax(y_true, scores):
    pred = np.argmax(np.asarray(scores), axis=1)
    y_true = np.asarray(y_true)
    acc = float(np.mean(pred == y_true))
    return acc, pred, confusion_counts(y_true, pred)


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


def ovr_scores(X, y, X_pred, lam):
    columns = []
    weights = []
    for k in range(N_CLASSES):
        yk = (np.asarray(y) == k).astype(float)
        wk = ridge_fit(X, yk, lam)
        weights.append(wk)
        columns.append(ridge_predict(X_pred, wk))
    return np.column_stack(columns), weights


def select_lambda(X, y):
    """Fit 40..127, pick lambda on 128..149 only. Test is never in the val slice."""
    inner_train_end = TRAIN - INNER_VAL
    X_inner, y_inner = X[:inner_train_end], y[:inner_train_end]
    X_val, y_val = X[inner_train_end:TRAIN], y[inner_train_end:TRAIN]
    (X_inner_z, X_val_z), _, _ = standardize(X_inner, X_val)
    scored = []
    grid_scores = []
    for lam in RIDGE_GRID:
        scores_val, _ = ovr_scores(X_inner_z, y_inner, X_val_z, lam)
        macro, per_class = macro_ovr_auc(y_val, scores_val)
        scored.append((-macro, -lam, lam, macro))
        grid_scores.append({
            "lambda": lam,
            "val_macro_ovr_auc": macro,
            "val_per_class": per_class,
        })
    scored.sort()
    return scored[0][2], grid_scores


def block_pooled_auc(y_test, scores_test):
    if len(y_test) != TEST:
        raise ValueError(f"test length {len(y_test)}; expected {TEST}")
    n_blocks = TEST // WINDOWS_PER_BLOCK
    pooled_scores = []
    pooled_y = []
    for block in range(n_blocks):
        sl = slice(block * WINDOWS_PER_BLOCK, (block + 1) * WINDOWS_PER_BLOCK)
        pooled_scores.append(np.mean(scores_test[sl], axis=0))
        labels = np.unique(y_test[sl])
        if len(labels) != 1:
            raise ValueError("test block is not a constant label")
        pooled_y.append(int(labels[0]))
    pooled_scores = np.vstack(pooled_scores)
    macro, per_class = macro_ovr_auc(np.array(pooled_y), pooled_scores)
    return macro, per_class


def fit_eval(X, y, lam):
    X_train, y_train = X[:TRAIN], y[:TRAIN]
    X_test, y_test = X[TRAIN:], y[TRAIN:]
    (X_train_z, X_test_z), _, _ = standardize(X_train, X_test)
    train_scores, _ = ovr_scores(X_train_z, y_train, X_train_z, lam)
    test_scores, _ = ovr_scores(X_train_z, y_train, X_test_z, lam)
    train_macro, train_per = macro_ovr_auc(y_train, train_scores)
    test_macro, test_per = macro_ovr_auc(y_test, test_scores)
    test_acc, test_pred, confusion = accuracy_argmax(y_test, test_scores)
    pooled_macro, pooled_per = block_pooled_auc(y_test, test_scores)
    return {
        "lambda": lam,
        "train_macro_ovr_auc": train_macro,
        "test_macro_ovr_auc": test_macro,
        "test_per_class_ovr_auc": {
            CLASS_NAMES[k]: test_per[k] for k in range(N_CLASSES)
        },
        "test_accuracy_argmax": test_acc,
        "chance_accuracy": CHANCE_ACC,
        "test_confusion": confusion,
        "test_pred": [int(v) for v in test_pred],
        "block_pooled_macro_ovr_auc": pooled_macro,
        "block_pooled_per_class_ovr_auc": {
            CLASS_NAMES[k]: pooled_per[k] for k in range(N_CLASSES)
        },
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


def overlaps_chance(mean, se, chance=CHANCE_AUC):
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


def overall_pass(path):
    if not Path(path).exists():
        return False
    text = Path(path).read_text(encoding="utf-8")
    return bool(re.search(r"^## Overall:\s*PASS\s*$", text, re.MULTILINE))


def overall_is(path, expected):
    if not Path(path).exists():
        return False
    text = Path(path).read_text(encoding="utf-8")
    return bool(re.search(rf"^## Overall:\s*{expected}\s*$", text, re.MULTILINE))


def heading_overall(path, heading, expected):
    if not Path(path).exists():
        return False
    text = Path(path).read_text(encoding="utf-8")
    pattern = rf"^## {re.escape(heading)}\s*\n(?:.*\n)*?^### Overall:\s*{expected}\s*$"
    return bool(re.search(pattern, text, re.MULTILINE))


def java_untouched():
    stage6 = (STAGE6_DIR / "BSimReservoirPlanStage6.java").read_text(encoding="utf-8")
    stage6_voxel = (STAGE6_DIR / "VoxelAnalyzer.java").read_text(encoding="utf-8")
    bencha = (BENCHA_DIR / "BSimReservoirPlanBenchA.java").read_text(encoding="utf-8")
    bencha2 = (BENCHA2_DIR / "BSimReservoirPlanBenchA2.java").read_text(encoding="utf-8")
    trackb = (TRACKB_DIR / "BSimReservoirPlanTrackB.java").read_text(encoding="utf-8")
    return (
        "package BSimReservoirPlanStage6;" in stage6
        and "NARMA-10" in stage6
        and "package BSimReservoirPlanStage6;" in stage6_voxel
        and "Track B is not started" in bencha
        and "Track B is not started" in bencha2
        and "SKIP=50" in bencha2
        and "package BSimReservoirPlanTrackB;" in trackb
        and "PROD_RATE = 1e6" in trackb
    )


def waveform_constraints():
    java = (HERE / "BSimReservoirPlanWaveform.java").read_text(encoding="utf-8")
    if any(token in java for token in FORBIDDEN_JAVA):
        return False
    required = (
        "new Vector3d(500, 250, 5)",
        "new Vector3d(300, 375, 5)",
        "AHL_SOURCE_RATE = 128000000",
        "ACID_SOURCE_RATE = 2.0e11",
        "RECEIVER_K_UM = 1.6",
        "RECEIVER_HILL_N = 2.0",
        "RECEIVER_TAU_S = 15.0",
        "FLOW_SPEED = 0.0",
        "GROWTH_RATE = 4.0 * Math.PI / 1800.0",
        "CARRYING_CAPACITY = 2000",
        "BOUND_X = 1000.0, BOUND_Y = 500.0",
        "input_ahl_waveform200.txt",
    )
    if any(token not in java for token in required):
        return False
    if "new Vector3d(150, 100, 5)" in java or "new Vector3d(850, 400, 5)" in java:
        return False
    return java.count("new BSimChemicalField(") == 4


def write_markdown(path, evidence):
    gates = evidence["gates"]
    overall = "PASS" if evidence["overall_gate_pass"] else "FAIL"
    lines = [
        "# Waveform classification gate evidence",
        "",
        "Waveform classification copies frozen Stage 6 / BenchA and tests",
        "sine / square / triangle AHL templates against Brownian density,",
        "a silent-source control, and a Stage 9-style field-only AHL",
        "baseline. Kinetics were not retuned. Stage 6 NARMA-10 remains",
        "**PASS**. BenchA Mackey-Glass remains **PASS**. BenchA Lorenz",
        "remains **FAIL**. Track B remains **FAIL**. Accuracy vs 1/3 is",
        "printed, not a gate.",
        "",
        f"## Overall: {overall}",
        "",
        evidence["failure_reason"] if not evidence["overall_gate_pass"]
        else "Driven biology beat Brownian and field-only AHL; silent sat at chance.",
        "",
        evidence["failure_detail"],
        "",
        "| Gate | Result |",
        "|---|---|",
        f"| 1. Driven macro OVR AUC > Brownian, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_brownian_pass'] else 'FAIL'}** |",
        f"| 2. Driven macro OVR AUC > field-only AHL, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_field_pass'] else 'FAIL'}** |",
        f"| 3. Silent not clearly above 0.5; silent not ≈ driven if driven has skill | **{'PASS' if gates['silent_chance_pass'] else 'FAIL'}** |",
        f"| 4. CSV 200/3200/3200; last sample 199;15;299.95 | **{'PASS' if gates['csv_validation_pass'] else 'FAIL'}** |",
        f"| 5. Stage 6 PASS; BenchA MG PASS; BenchA Lorenz FAIL; Track B FAIL; Java untouched | **{'PASS' if gates['provenance_pass'] else 'FAIL'}** |",
        "",
        "## Frozen task",
        "",
        "40 blocks × 5 windows. Class vector SHA-256 "
        f"`{evidence['class_vector_sha256']}`. Templates and the shuffle",
        "seed were not changed after seeing AUC.",
        "",
        f"Class vector: `{evidence['class_vector']}`",
        "",
        "| Split | Windows | n | sine | square | triangle |",
        "|---|---|---|---|---|---|",
    ]
    for row in evidence["class_balance"]:
        lines.append(
            f"| {row['split']} | {row['windows']} | {row['n']} | "
            f"{row['sine']} | {row['square']} | {row['triangle']} |"
        )
    lines.extend([
        "",
        f"Feature counts: biology `{evidence['driven_feature_count']}`, "
        f"Brownian `{evidence['brownian_feature_count']}`, "
        f"field-only `{evidence['field_feature_count']}`.",
        "",
        "## Test macro OVR AUC (primary)",
        "",
        "Chance macro OVR AUC = 0.5. Chance accuracy = 1/3. Do not use",
        "Track B's binary 0.5 as the accuracy floor.",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 | mean ± s.e. |",
        "|---|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "driven", "field"):
        values = evidence["test_macro_ovr_auc"][arm]
        mean = evidence["auc_mean_se"][arm]["mean"]
        se = evidence["auc_mean_se"][arm]["se"]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} | {mean:.4f} ± {se:.4f} |"
        )
    lines.extend([
        "",
        "## Per-class OVR AUC (driven / field-only)",
        "",
        "| Arm | seed | sine | square | triangle |",
        "|---|---|---|---|---|",
    ])
    for arm in ("driven", "field", "brownian", "silent"):
        for seed, per_class in zip(SEEDS, evidence["test_per_class_ovr_auc"][arm]):
            lines.append(
                f"| {arm} | {seed} | {per_class['sine']:.4f} | "
                f"{per_class['square']:.4f} | {per_class['triangle']:.4f} |"
            )
    lines.extend([
        "",
        "Accuracy (argmax of three scores) is printed, not the gate.",
        "Compared to 1/3, not to 0.5.",
        "",
        "| Arm | seed 101 acc | seed 202 acc | seed 303 acc |",
        "|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "driven", "field"):
        values = evidence["test_accuracy"][arm]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} |"
        )
    lines.extend([
        "",
        "## Confusion counts (test windows, argmax)",
        "",
        "Rows = true class (sine, square, triangle). Columns = predicted.",
        "",
    ])
    for arm in ("driven", "field", "brownian", "silent"):
        lines.append(f"### {arm}")
        lines.append("")
        for seed, matrix in zip(SEEDS, evidence["test_confusion"][arm]):
            lines.append(f"seed {seed}: `{matrix}`")
        lines.append("")
    lines.extend([
        "## Block-pooled macro OVR AUC (diagnostic, not a gate)",
        "",
        "Mean of 5 window scores per test block (10 scores). Do not replace",
        "window AUC with this after seeing numbers.",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 |",
        "|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "driven", "field"):
        values = evidence["block_pooled_auc"][arm]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} |"
        )
    lines.extend([
        "",
        "## Driven feature list",
        "",
        "- `Receiver_R_*` / Mean_q, 20×10",
        "- `Lum_Mean_*` / Mean_L, 20×10",
        "- `Input_Driven_Death_*`, 4×2",
        "",
        f"Printed count: {evidence['driven_feature_count']} features.",
        "",
        f"Field-only: `AHL_uM_*` ({evidence['field_feature_count']}).",
        "",
        "Excluded from the biology ridge: " + ", ".join(EXCLUDED_FROM_RIDGE) + ".",
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
        "See `PROTOCOL.md`. Frozen Stage 6 dish. Block-boundary splits.",
        "Lambda fit 40..127, pick on 128..149 only. If gate 2 fails,",
        "Overall is FAIL. ACs were not added. Templates were not changed.",
        "The gate was not switched to accuracy.",
        "",
        "CSV completeness is the on-disk files, last sample `199;15;299.95`.",
        "Silent CSVs are Stage 6 copies. Silent was not rerun.",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brownian", nargs=3, default=[
        str(HERE / "results" / f"waveform_brownian_seed{seed}") for seed in SEEDS
    ])
    parser.add_argument("--silent", nargs=3, default=[
        str(HERE / "results" / f"waveform_silent_seed{seed}") for seed in SEEDS
    ])
    parser.add_argument("--driven", nargs=3, default=[
        str(HERE / "results" / f"waveform_driven_seed{seed}") for seed in SEEDS
    ])
    parser.add_argument("--labels", default=str(HERE / "waveform_labels.csv"))
    parser.add_argument("--expected-windows", type=int, default=NUM_WINDOWS)
    parser.add_argument("--expected-aux-rows", type=int, default=EXPECTED_AUX_ROWS)
    parser.add_argument("--evidence", default=str(HERE / "results" / "waveform_gate_evidence.json"))
    parser.add_argument("--markdown", default=str(HERE / "results" / "GATE_EVIDENCE.md"))
    args = parser.parse_args()

    u_file = load_sequence(HERE / "input_ahl_waveform200.txt")
    digest = sha256_u(u_file)
    if digest != U_SHA:
        raise ValueError(f"AHL SHA-256 {digest} does not match frozen {U_SHA}")

    y, block, u_csv, label_rows = load_labels(args.labels)
    reconstructed = ",".join(str(int(y[p * WINDOWS_PER_BLOCK])) for p in range(NUM_BLOCKS))
    class_sha = hashlib.sha256(reconstructed.encode("ascii")).hexdigest()
    if reconstructed != CLASS_VECTOR or class_sha != CLASS_SHA:
        raise ValueError("class vector does not match the frozen PROTOCOL.md declaration")
    for n, row in enumerate(label_rows):
        if int(row["block"]) != n // WINDOWS_PER_BLOCK:
            raise ValueError(f"block index mismatch at window {n}")
        if int(row["y"]) != int(y[n]):
            raise ValueError(f"label mismatch at window {n}")
        expected_u = template_u(int(y[n]), n % WINDOWS_PER_BLOCK)
        if abs(float(row["u"]) - expected_u) > 1e-12:
            raise ValueError(f"template u mismatch at window {n}")
        if abs(u_file[n] - expected_u) > 1e-12:
            raise ValueError(f"AHL file u mismatch at window {n}")
        if abs(u_csv[n] - expected_u) > 1e-12:
            raise ValueError(f"labels.csv u mismatch at window {n}")

    balances = [
        class_balance(y, 0, WASHOUT, "washout"),
        class_balance(y, WASHOUT, WASHOUT + TRAIN, "train"),
        class_balance(y, WASHOUT + TRAIN, NUM_WINDOWS, "test"),
        class_balance(y, 0, NUM_WINDOWS, "all"),
    ]
    print("CLASS_VECTOR")
    print(reconstructed)
    print(f"class_vector_sha256={class_sha}")
    print(f"u_sha256={digest}")
    print("CLASS_BALANCE")
    for row in balances:
        print(
            f"{row['split']} windows={row['windows']} n={row['n']} "
            f"sine={row['sine']} square={row['square']} triangle={row['triangle']}"
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
        arm: [run["classification"]["test_macro_ovr_auc"] for run in runs]
        for arm, runs in evaluated.items()
    }
    acc_by_arm = {
        arm: [run["classification"]["test_accuracy_argmax"] for run in runs]
        for arm, runs in evaluated.items()
    }
    pooled_by_arm = {
        arm: [run["classification"]["block_pooled_macro_ovr_auc"] for run in runs]
        for arm, runs in evaluated.items()
    }
    per_class_by_arm = {
        arm: [run["classification"]["test_per_class_ovr_auc"] for run in runs]
        for arm, runs in evaluated.items()
    }
    confusion_by_arm = {
        arm: [run["classification"]["test_confusion"] for run in runs]
        for arm, runs in evaluated.items()
    }
    field_auc = [run["field_only"]["test_macro_ovr_auc"] for run in evaluated["driven"] if run["field_only"]]
    field_acc = [
        run["field_only"]["test_accuracy_argmax"] for run in evaluated["driven"] if run["field_only"]
    ]
    field_pooled = [
        run["field_only"]["block_pooled_macro_ovr_auc"]
        for run in evaluated["driven"] if run["field_only"]
    ]
    field_per_class = [
        run["field_only"]["test_per_class_ovr_auc"]
        for run in evaluated["driven"] if run["field_only"]
    ]
    field_confusion = [
        run["field_only"]["test_confusion"]
        for run in evaluated["driven"] if run["field_only"]
    ]
    auc_by_arm["field"] = field_auc
    acc_by_arm["field"] = field_acc
    pooled_by_arm["field"] = field_pooled
    per_class_by_arm["field"] = field_per_class
    confusion_by_arm["field"] = field_confusion
    mean_se_by_arm = {arm: mean_se(values) for arm, values in auc_by_arm.items()}

    driven_features = evaluated["driven"][0]["feature_names"]
    expected_prefixes = ("Receiver_R_", "Lum_Mean_", "Input_Driven_Death_")
    feature_ok = (
        all(name.startswith(expected_prefixes) for name in driven_features)
        and not any(name.startswith(BIOLOGY_FORBIDDEN) for name in driven_features)
        and len(driven_features) == 408
        and all(run["n_features"] == 408 for run in evaluated["driven"] + evaluated["silent"])
        and all(run["n_features"] == 200 for run in evaluated["brownian"])
        and all(run["field_only"]["n_features"] == 200 for run in evaluated["driven"])
        and sum(name.startswith("Receiver_R_") for name in driven_features) == 200
        and sum(name.startswith("Lum_Mean_") for name in driven_features) == 200
        and sum(name.startswith("Input_Driven_Death_") for name in driven_features) == 8
    )
    print(f"driven_feature_count={len(driven_features)}")
    print(f"brownian_feature_count={evaluated['brownian'][0]['n_features']}")
    print(f"silent_feature_count={evaluated['silent'][0]['n_features']}")
    print(f"field_feature_count={evaluated['driven'][0]['field_only']['n_features']}")
    print("DRIVEN_FEATURE_LIST")
    for name in driven_features:
        print(name)

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
    silent_not_above_chance = silent_mean - silent_se <= CHANCE_AUC
    driven_has_skill = driven_mean - driven_se > CHANCE_AUC
    silent_overlap_driven = intervals_overlap(silent_mean, silent_se, driven_mean, driven_se)
    silent_not_approx_driven = (not silent_overlap_driven) if driven_has_skill else True
    silent_chance_pass = silent_not_above_chance and silent_not_approx_driven

    provenance_pass = (
        overall_pass(STAGE6_DIR / "results" / "GATE_EVIDENCE.md")
        and heading_overall(BENCHA_DIR / "results" / "GATE_EVIDENCE.md", "Mackey-Glass", "PASS")
        and heading_overall(BENCHA_DIR / "results" / "GATE_EVIDENCE.md", "Lorenz'63", "FAIL")
        and overall_is(TRACKB_DIR / "results" / "GATE_EVIDENCE.md", "FAIL")
        and java_untouched()
        and waveform_constraints()
        and feature_ok
        and (STAGE9_DIR / "PROTOCOL.md").exists()
    )

    gates = {
        "replicate_count_pass": all(len(values) == 3 for values in auc_by_arm.values()),
        "driven_beats_brownian_pass": driven_beats_brownian,
        "driven_beats_field_pass": driven_beats_field,
        "silent_chance_pass": silent_chance_pass,
        "csv_validation_pass": csv_pass,
        "provenance_pass": provenance_pass,
        "frozen_feature_list_pass": feature_ok,
    }
    overall = all(
        gates[key] for key in (
            "driven_beats_brownian_pass",
            "driven_beats_field_pass",
            "silent_chance_pass",
            "csv_validation_pass",
            "provenance_pass",
        )
    )
    if not gates["driven_beats_field_pass"]:
        failure_reason = (
            "field_wins: driven macro OVR AUC did not beat the field-only "
            "AHL baseline in every seed with non-overlapping mean±s.e. "
            "The plume is the waveform. Do not add ACs. Do not mix in "
            "Att/Rep. Do not raise source rates. Do not make the three "
            "templates more similar. Do not switch to accuracy vs 0.334."
        )
    elif not gates["driven_beats_brownian_pass"]:
        failure_reason = (
            "brownian_wins: driven macro OVR AUC did not beat Brownian "
            "density in every seed with non-overlapping mean±s.e."
        )
    elif not gates["silent_chance_pass"]:
        failure_reason = (
            "silent_skill: silent macro OVR AUC is clearly above 0.5 or "
            "overlaps driven while driven has skill."
        )
    elif not gates["csv_validation_pass"]:
        failure_reason = (
            "csv_incomplete: a summary/sample/voxel file is missing rows, is "
            "ragged, or last sample is not 199;15;299.95."
        )
    elif not gates["provenance_pass"]:
        failure_reason = (
            "provenance: Stage 6/BenchA/Track B labels were reopened, frozen "
            "Java was touched, or glucose/Danino/vesicles/extra ACs leaked in."
        )
    else:
        failure_reason = ""

    if not driven_beats_field:
        failure_detail = (
            "Gate 2 fails. Driven test macro OVR AUC mean ± s.e. is "
            f"{driven_mean:.4f} ± {driven_se:.4f}; field-only is "
            f"{field_mean:.4f} ± {field_se:.4f}. The AHL plume already "
            "ranks the waveform. Kinetics, source rates, AC sites, and "
            "templates were not changed after seeing this. Accuracy vs "
            "1/3 was not used as a gate."
        )
    elif not driven_beats_brownian:
        failure_detail = (
            "Gate 1 fails. Driven did not beat Brownian density with "
            "non-overlapping mean±s.e. in every seed."
        )
    elif not silent_chance_pass:
        failure_detail = (
            "Gate 3 fails. Silent macro OVR AUC is clearly above 0.5 or "
            "overlaps driven while driven has skill."
        )
    elif overall:
        failure_detail = (
            "No failure. Driven beat Brownian and field-only AHL; silent "
            "was not clearly above chance. Frozen 408-feature contract used. "
            "Accuracy vs 1/3 is diagnostic only."
        )
    else:
        failure_detail = failure_reason

    evidence = {
        "schema": "BSimReservoirPlanWaveform-gates-v1",
        "stage6_status": "PASS",
        "bencha_lorenz_status": "FAIL",
        "bencha_mg_status": "PASS",
        "trackb_status": "FAIL",
        "class_vector": reconstructed,
        "class_vector_sha256": class_sha,
        "u_sha256": digest,
        "protocol": {
            "num_windows": NUM_WINDOWS,
            "washout": WASHOUT,
            "train": TRAIN,
            "test": TEST,
            "inner_val": INNER_VAL,
            "lambda_validation_slice": "X[inner_train_end:TRAIN] after washout strip",
            "ridge_grid": list(RIDGE_GRID),
            "primary_metric": "window test macro one-vs-rest AUC",
            "chance_accuracy": CHANCE_ACC,
            "chance_macro_ovr_auc": CHANCE_AUC,
            "classes": list(CLASS_NAMES),
        },
        "class_balance": balances,
        "analysis_channels_driven": [
            "Receiver_R_* / Mean_q (20x10)",
            "Lum_Mean_* / Mean_L (20x10)",
            "Input_Driven_Death_* (4x2)",
        ],
        "driven_feature_count": len(driven_features),
        "brownian_feature_count": evaluated["brownian"][0]["n_features"],
        "field_feature_count": evaluated["driven"][0]["field_only"]["n_features"],
        "excluded_from_ridge": list(EXCLUDED_FROM_RIDGE),
        "test_macro_ovr_auc": auc_by_arm,
        "test_per_class_ovr_auc": per_class_by_arm,
        "test_accuracy": acc_by_arm,
        "test_confusion": confusion_by_arm,
        "block_pooled_auc": pooled_by_arm,
        "auc_mean_se": {
            arm: {"mean": mean, "se": se} for arm, (mean, se) in mean_se_by_arm.items()
        },
        "lambdas": {
            "brownian": [run["classification"]["lambda"] for run in evaluated["brownian"]],
            "silent": [run["classification"]["lambda"] for run in evaluated["silent"]],
            "driven": [run["classification"]["lambda"] for run in evaluated["driven"]],
            "field": [run["field_only"]["lambda"] for run in evaluated["driven"]],
        },
        "files": {
            arm: [{
                "path": run["path"],
                "n_features": run["n_features"],
                "classification": {
                    k: v for k, v in run["classification"].items() if k != "test_pred"
                },
                "field_only": (
                    {k: v for k, v in run["field_only"].items()
                     if k not in {"feature_names", "test_pred"}}
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
        "test_macro_ovr_auc": auc_by_arm,
        "auc_mean_se": evidence["auc_mean_se"],
        "test_accuracy": acc_by_arm,
        "test_per_class_ovr_auc": per_class_by_arm,
        "test_confusion": confusion_by_arm,
        "block_pooled_auc": pooled_by_arm,
        "failure_reason": failure_reason,
        "driven_feature_count": len(driven_features),
        "brownian_feature_count": evidence["brownian_feature_count"],
        "field_feature_count": evidence["field_feature_count"],
    }, indent=2, allow_nan=True, default=json_default))
    print(f"evidence_file={destination}")
    print(f"markdown_file={args.markdown}")
    print(f"Overall: {'PASS' if overall else 'FAIL'}")
    print("MACRO_OVR_AUC_TABLE")
    for arm in ("brownian", "silent", "driven", "field"):
        values = auc_by_arm[arm]
        mean, se = mean_se_by_arm[arm]
        print(
            f"{arm} {values[0]:.4f} {values[1]:.4f} {values[2]:.4f} "
            f"mean {mean:.4f} se {se:.4f}"
        )


if __name__ == "__main__":
    main()
