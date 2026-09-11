#!/usr/bin/env python3
"""Evaluate Waveform2b distinct-orbit temporal-order classification.

Does not rewrite Waveform GATE_EVIDENCE.md or Waveform2. Kinetics are
not retuned. Primary numbers are block-level. u-only must PASS before BSim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np

WASHOUT_BLOCKS = 8
TRAIN_BLOCKS = 26
TEST_BLOCKS = 16
INNER_VAL_BLOCKS = 6
WINDOWS_PER_BLOCK = 8
NUM_BLOCKS = 50
WASHOUT = WASHOUT_BLOCKS * WINDOWS_PER_BLOCK
TRAIN = TRAIN_BLOCKS * WINDOWS_PER_BLOCK
TEST = TEST_BLOCKS * WINDOWS_PER_BLOCK
INNER_VAL = INNER_VAL_BLOCKS * WINDOWS_PER_BLOCK
NUM_WINDOWS = WASHOUT + TRAIN + TEST
N_CLASSES = 3
CLASS_NAMES = ("ORDER_UNI", "ORDER_DOWN", "ORDER_UP")
CHANCE_ACC = 1.0 / 3.0
CHANCE_AUC = 0.5
EXPECTED_AUX_ROWS = NUM_WINDOWS * 16
LAST_SAMPLE = ("399", "15", "299.95")
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
SEEDS = (101, 202, 303)
SCOUT_SEED = 101

S2 = 0.25 * math.sqrt(2.0) / 2.0
LO = round(0.25 - S2, 12)
MID = round(0.25, 12)
HI = round(0.25 + S2, 12)
MINV = round(0.0, 12)
MAXV = round(0.5, 12)
SHARED_MULTISET = (MINV, LO, LO, MID, MID, HI, HI, MAXV)
TEMPLATES = (
    (MID, HI, MAXV, HI, MID, LO, MINV, LO),
    (MAXV, HI, HI, MID, MID, LO, LO, MINV),
    (MINV, LO, LO, MID, MID, HI, HI, MAXV),
)

CLASS_VECTOR = (
    "0,2,2,0,1,1,0,1,0,1,1,0,1,2,2,0,1,0,1,0,2,0,1,0,2,0,1,0,2,1,2,1,2,2,"
    "0,1,0,0,1,2,0,1,2,2,2,0,0,1,1,2"
)
PHASE_VECTOR = (
    "2,5,3,0,2,4,6,3,2,4,2,0,7,7,0,5,0,7,6,3,1,1,1,7,4,4,5,6,3,3,6,6,2,5,"
    "6,6,7,3,4,3,4,1,7,4,6,5,2,2,3,2"
)
CLASS_SHA = "26b44c1c2ee3932becc7fdebfc72179a2872bb87baca975624355c0f315772e6"
PHASE_SHA = "01f972847d1dad3c3b6591ff7fab7537b8f4b020da40f6f35376dfa02e97b7bc"
U_SHA = "339e8eb8c6a89510bfe1ef0f38938cc4407c6470e526e5fcf7c2ec734b9d5cdf"

WASHOUT_COUNTS = (3, 3, 2)
TRAIN_COUNTS = (9, 9, 8)
TEST_COUNTS = (6, 5, 5)

K_HILL = 1.6
TAU_R = 15.0
TAU_L = 1500.0
OCC_DEAD_R = 0.05
OCC_DEAD_CORR = 0.5

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
WAVEFORM_DIR = HERE.parent / "BSimReservoirPlanWaveform"
NARMA_DIR = HERE.parent / "BSimReservoirPlanNarma10b"


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


def sha256_ascii_ints(values):
    payload = ",".join(str(int(v)) for v in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def cyclic_shifts(seq):
    seq = tuple(seq)
    return {tuple(round(v, 12) for v in seq[p:] + seq[:p]) for p in range(len(seq))}


def clip_u(value):
    return min(MAXV, max(MINV, value))


def template_u(cls, phase, w):
    return clip_u(TEMPLATES[cls][(w + phase) % WINDOWS_PER_BLOCK])


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
    y = np.array([int(row["class"]) for row in rows], dtype=int)
    block = np.array([int(row["block"]) for row in rows], dtype=int)
    phase = np.array([int(row["phase"]) for row in rows], dtype=int)
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    return y, block, phase, u, rows


def verify_frozen_inputs():
    u_file = load_sequence(HERE / "input_ahl_waveform2b_400.txt")
    digest = sha256_u(u_file)
    if digest != U_SHA:
        raise ValueError(f"AHL SHA-256 {digest} does not match frozen {U_SHA}")
    y, block, phase, u_csv, label_rows = load_labels(HERE / "waveform2b_labels.csv")
    reconstructed_y = ",".join(str(int(y[p * WINDOWS_PER_BLOCK])) for p in range(NUM_BLOCKS))
    reconstructed_p = ",".join(str(int(phase[p * WINDOWS_PER_BLOCK])) for p in range(NUM_BLOCKS))
    class_sha = hashlib.sha256(reconstructed_y.encode("ascii")).hexdigest()
    phase_sha = hashlib.sha256(reconstructed_p.encode("ascii")).hexdigest()
    if reconstructed_y != CLASS_VECTOR or class_sha != CLASS_SHA:
        raise ValueError("class vector does not match frozen PROTOCOL.md")
    if reconstructed_p != PHASE_VECTOR or phase_sha != PHASE_SHA:
        raise ValueError("phase vector does not match frozen PROTOCOL.md")
    for n, row in enumerate(label_rows):
        if int(row["block"]) != n // WINDOWS_PER_BLOCK:
            raise ValueError(f"block index mismatch at window {n}")
        expected_u = round(template_u(int(y[n]), int(phase[n]), n % WINDOWS_PER_BLOCK), 12)
        if abs(float(row["u"]) - expected_u) > 1e-12:
            raise ValueError(f"template u mismatch at window {n}")
        if abs(u_file[n] - expected_u) > 1e-12:
            raise ValueError(f"AHL file u mismatch at window {n}")
        if abs(u_csv[n] - expected_u) > 1e-12:
            raise ValueError(f"labels.csv u mismatch at window {n}")
    acid = load_sequence(HERE / "input_acid_held05_400.txt")
    if any(abs(v - 0.5) > 1e-12 for v in acid):
        raise ValueError("acid file is not held at 0.5")
    orbits = [cyclic_shifts(t) for t in TEMPLATES]
    for i in range(N_CLASSES):
        for j in range(i + 1, N_CLASSES):
            if orbits[i] & orbits[j]:
                raise ValueError(
                    f"orbit gate failed: {CLASS_NAMES[i]} intersects {CLASS_NAMES[j]}"
                )
    return {
        "u": np.asarray(u_file, dtype=float),
        "y_window": y,
        "block": block,
        "phase": phase,
        "class_vector": reconstructed_y,
        "phase_vector": reconstructed_p,
        "class_sha": class_sha,
        "phase_sha": phase_sha,
        "u_sha": digest,
    }


def block_labels(y_window):
    labels = []
    for b in range(NUM_BLOCKS):
        sl = y_window[b * WINDOWS_PER_BLOCK:(b + 1) * WINDOWS_PER_BLOCK]
        uniq = np.unique(sl)
        if len(uniq) != 1:
            raise ValueError(f"block {b} is not a constant class")
        labels.append(int(uniq[0]))
    return np.array(labels, dtype=int)


def u_blocks(u):
    return np.asarray(u, dtype=float).reshape(NUM_BLOCKS, WINDOWS_PER_BLOCK)


def block_moments_from_u(block_u):
    mean = block_u.mean(axis=1)
    var = block_u.var(axis=1, ddof=0)
    power = np.mean(block_u * block_u, axis=1)
    mn = block_u.min(axis=1)
    mx = block_u.max(axis=1)
    return {
        "MEAN_ONLY": mean.reshape(-1, 1),
        "POWER_ONLY": power.reshape(-1, 1),
        "VARIANCE_ONLY": var.reshape(-1, 1),
        "MOMENTS": np.column_stack([mean, var, power, mn, mx, mx - mn]),
        "RAW_U_8": block_u.copy(),
        "START_U": block_u[:, :1].copy(),
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
    per_class = [roc_auc((y_true == k).astype(int), scores[:, k]) for k in range(N_CLASSES)]
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


def select_lambda_windows(X, y):
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


def select_lambda_blocks(X, y):
    inner_train_end = TRAIN_BLOCKS - INNER_VAL_BLOCKS
    X_inner, y_inner = X[:inner_train_end], y[:inner_train_end]
    X_val, y_val = X[inner_train_end:TRAIN_BLOCKS], y[inner_train_end:TRAIN_BLOCKS]
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


def pool_block_scores(y_windows, scores_windows, n_windows):
    if len(y_windows) != n_windows:
        raise ValueError(f"window count {len(y_windows)}; expected {n_windows}")
    n_blocks = n_windows // WINDOWS_PER_BLOCK
    pooled_scores = []
    pooled_y = []
    for block in range(n_blocks):
        sl = slice(block * WINDOWS_PER_BLOCK, (block + 1) * WINDOWS_PER_BLOCK)
        pooled_scores.append(np.mean(scores_windows[sl], axis=0))
        labels = np.unique(y_windows[sl])
        if len(labels) != 1:
            raise ValueError("block is not a constant label")
        pooled_y.append(int(labels[0]))
    pooled_scores = np.vstack(pooled_scores)
    pooled_y = np.array(pooled_y, dtype=int)
    macro, per_class = macro_ovr_auc(pooled_y, pooled_scores)
    acc, pred, confusion = accuracy_argmax(pooled_y, pooled_scores)
    return {
        "macro_ovr_auc": macro,
        "per_class_ovr_auc": {CLASS_NAMES[k]: per_class[k] for k in range(N_CLASSES)},
        "accuracy": acc,
        "confusion": confusion,
        "pred": [int(v) for v in pred],
        "n_blocks": n_blocks,
    }


def fit_eval_windows(X, y, lam):
    X_train, y_train = X[:TRAIN], y[:TRAIN]
    X_test, y_test = X[TRAIN:], y[TRAIN:]
    (X_train_z, X_test_z), _, _ = standardize(X_train, X_test)
    train_scores, _ = ovr_scores(X_train_z, y_train, X_train_z, lam)
    test_scores, _ = ovr_scores(X_train_z, y_train, X_test_z, lam)
    train_macro, train_per = macro_ovr_auc(y_train, train_scores)
    test_macro, test_per = macro_ovr_auc(y_test, test_scores)
    test_acc, test_pred, confusion = accuracy_argmax(y_test, test_scores)
    block = pool_block_scores(y_test, test_scores, TEST)
    return {
        "lambda": lam,
        "train_window_macro_ovr_auc": train_macro,
        "test_window_macro_ovr_auc": test_macro,
        "test_window_per_class_ovr_auc": {
            CLASS_NAMES[k]: test_per[k] for k in range(N_CLASSES)
        },
        "test_window_accuracy": test_acc,
        "test_window_confusion": confusion,
        "block_macro_ovr_auc": block["macro_ovr_auc"],
        "block_per_class_ovr_auc": block["per_class_ovr_auc"],
        "block_accuracy": block["accuracy"],
        "block_confusion": block["confusion"],
        "chance_accuracy": CHANCE_ACC,
    }


def fit_eval_blocks(X, y, lam):
    X_train, y_train = X[:TRAIN_BLOCKS], y[:TRAIN_BLOCKS]
    X_test, y_test = X[TRAIN_BLOCKS:], y[TRAIN_BLOCKS:]
    (X_train_z, X_test_z), _, _ = standardize(X_train, X_test)
    test_scores, _ = ovr_scores(X_train_z, y_train, X_test_z, lam)
    test_macro, test_per = macro_ovr_auc(y_test, test_scores)
    test_acc, test_pred, confusion = accuracy_argmax(y_test, test_scores)
    return {
        "lambda": lam,
        "block_macro_ovr_auc": test_macro,
        "block_per_class_ovr_auc": {
            CLASS_NAMES[k]: test_per[k] for k in range(N_CLASSES)
        },
        "block_accuracy": test_acc,
        "block_confusion": confusion,
        "n_features": int(X.shape[1]),
    }


def evaluate_window_readout(X_all, y_all, windows):
    if windows != list(range(NUM_WINDOWS)):
        raise ValueError(
            f"expected windows 0..{NUM_WINDOWS - 1}, got {windows[:3]}..{windows[-1:]}"
        )
    X = X_all[WASHOUT:]
    y = y_all[WASHOUT:]
    lam, grid_scores = select_lambda_windows(X, y)
    metrics = fit_eval_windows(X, y, lam)
    metrics["lambda_grid"] = grid_scores
    metrics["n_features"] = int(X.shape[1])
    return metrics


def evaluate_block_readout(X_all_blocks, y_all_blocks):
    X = X_all_blocks[WASHOUT_BLOCKS:]
    y = y_all_blocks[WASHOUT_BLOCKS:]
    lam, grid_scores = select_lambda_blocks(X, y)
    metrics = fit_eval_blocks(X, y, lam)
    metrics["lambda_grid"] = grid_scores
    return metrics


def mean_se(values):
    array = np.array(values, dtype=float)
    mean = float(np.mean(array))
    se = float(np.std(array, ddof=1) / math.sqrt(len(array))) if len(array) > 1 else float("nan")
    return mean, se


def columns_with_prefix(header, prefixes):
    names = []
    for prefix in prefixes:
        names.extend(name for name in header if name.startswith(prefix))
    return names


def validate_csv(path, expected_rows, check_last_sample=False):
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
        "last_sample_pass": last_ok if check_last_sample else True,
    }


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


def block_mean_features(X_windows):
    n_windows = X_windows.shape[0]
    if n_windows != NUM_WINDOWS:
        raise ValueError(f"expected {NUM_WINDOWS} windows, got {n_windows}")
    return X_windows.reshape(NUM_BLOCKS, WINDOWS_PER_BLOCK, X_windows.shape[1]).mean(axis=1)


def pearson(a, b):
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return 0.0 if den < 1e-15 else float(np.sum(a * b) / den)


def hill(c):
    c = np.asarray(c, dtype=float)
    return c * c / (K_HILL * K_HILL + c * c)


def load_voxel_arrays(path):
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    ahl_i = [header.index(f"AHL_uM_{i}") for i in range(200)]
    den_i = [header.index(f"Den_{i}") for i in range(200)]
    r_i = [header.index(f"Receiver_R_{i}") for i in range(200)]
    l_i = [header.index(f"Lum_Mean_{i}") for i in range(200)]
    win = header.index("Window")
    samp = header.index("Sample")
    tcol = header.index("TimeInWindow_s")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return {
        "window": data[:, win].astype(int),
        "sample": data[:, samp].astype(int),
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
        er = math.exp(-delta / TAU_R)
        el = math.exp(-delta / TAU_L)
        old_r, old_l = R, L
        R = h + (old_r - h) * er
        L = h + (old_l - h) * el + (old_r - h) * TAU_R / (TAU_R - TAU_L) * (er - el)
        out_r[index], out_l[index] = R, L
    return out_r, out_l


def occupancy(voxels, u):
    r_win = window_mean_matrix(voxels["receiver"], voxels["window"])
    a_win = window_mean_matrix(voxels["ahl"], voxels["window"])
    mean_r_w = np.mean(r_win, axis=1)
    mean_a_w = np.mean(a_win, axis=1)
    mean_r = float(np.mean(mean_r_w))
    r_u = pearson(mean_r_w, u)
    if abs(r_u) < OCC_DEAD_CORR or mean_r < OCC_DEAD_R:
        label = "DEAD"
    else:
        label = "ALIVE"
    return {
        "occupancy": label,
        "mean_R": mean_r,
        "r_meanR_u": r_u,
        "mean_AHL": float(np.mean(mean_a_w)),
    }


def summary_stats(run_dir):
    path = Path(run_dir) / "window_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    pop = np.array([float(row["Population"]) for row in rows], dtype=float)
    return {"mean_pop": float(np.mean(pop))}


def ahl_finite_nonneg(voxels_path):
    voxels = load_voxel_arrays(voxels_path)
    ahl = voxels["ahl"]
    finite = bool(np.isfinite(ahl).all())
    nonnegative = bool(np.min(ahl) >= -1e-15)
    return {
        "finite": finite,
        "nonnegative": nonnegative,
        "min": float(np.min(ahl)),
        "max": float(np.max(ahl)),
        "mean": float(np.mean(ahl)),
        "pass": finite and nonnegative,
    }


def u_only_baselines(frozen):
    y_block = block_labels(frozen["y_window"])
    features = block_moments_from_u(u_blocks(frozen["u"]))
    rows = []
    pass_scalar = True
    raw_ok = True
    leak = False
    for name, X in features.items():
        metrics = evaluate_block_readout(X, y_block)
        auc = metrics["block_macro_ovr_auc"]
        scalar = name != "RAW_U_8"
        if scalar:
            within = abs(auc - 0.5) <= 0.15
            if auc >= 0.70:
                leak = True
                within = False
            pass_scalar = pass_scalar and within
            criterion = "within_0.15_of_0.5_and_lt_0.70"
            passed = within
        else:
            passed = auc >= 0.90
            raw_ok = passed
            criterion = "ge_0.90"
        rows.append({
            "baseline": name,
            "n_features": metrics["n_features"],
            "lambda": metrics["lambda"],
            "block_macro_ovr_auc": auc,
            "block_accuracy": metrics["block_accuracy"],
            "ORDER_UNI_auc": metrics["block_per_class_ovr_auc"]["ORDER_UNI"],
            "ORDER_DOWN_auc": metrics["block_per_class_ovr_auc"]["ORDER_DOWN"],
            "ORDER_UP_auc": metrics["block_per_class_ovr_auc"]["ORDER_UP"],
            "confusion": metrics["block_confusion"],
            "criterion": criterion,
            "pass": passed,
        })
    overall = pass_scalar and raw_ok and not leak
    return {
        "rows": rows,
        "pass_scalar": pass_scalar,
        "raw_ok": raw_ok,
        "leak": leak,
        "overall_pass": overall,
        "moments_auc": next(r["block_macro_ovr_auc"] for r in rows if r["baseline"] == "MOMENTS"),
        "mean_auc": next(r["block_macro_ovr_auc"] for r in rows if r["baseline"] == "MEAN_ONLY"),
        "power_auc": next(r["block_macro_ovr_auc"] for r in rows if r["baseline"] == "POWER_ONLY"),
        "variance_auc": next(r["block_macro_ovr_auc"] for r in rows if r["baseline"] == "VARIANCE_ONLY"),
        "raw_auc": next(r["block_macro_ovr_auc"] for r in rows if r["baseline"] == "RAW_U_8"),
        "start_auc": next(r["block_macro_ovr_auc"] for r in rows if r["baseline"] == "START_U"),
    }


def write_u_only(audit):
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    csv_path = results / "u_only_baselines.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "baseline", "n_features", "lambda", "block_macro_ovr_auc",
            "block_accuracy", "ORDER_UNI_auc", "ORDER_DOWN_auc", "ORDER_UP_auc",
            "confusion", "criterion", "pass",
        ])
        writer.writeheader()
        for row in audit["rows"]:
            out = dict(row)
            out["confusion"] = json.dumps(row["confusion"])
            writer.writerow(out)
    status = "PASS" if audit["overall_pass"] else "FAIL"
    lines = [
        "# Waveform2b u-only baselines",
        "",
        "Block-level one-vs-rest ridge on the frozen distinct-orbit",
        "equal-histogram templates. Hashes were recorded in PROTOCOL.md",
        "before this ridge. Waveform2 remains TASK_VOID. This is not",
        "reservoir skill.",
        "",
        f"## Result: {status}",
        "",
        f"Scalar-cue family pass: {audit['pass_scalar']}.",
        f"RAW_U_8 pass: {audit['raw_ok']}.",
        f"Scalar leak (>=0.70): {audit['leak']}.",
        "",
        "| Baseline | Features | Lambda | Block macro AUC | Accuracy | ORDER_UNI | ORDER_DOWN | ORDER_UP | Pass |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in audit["rows"]:
        lines.append(
            f"| {row['baseline']} | {row['n_features']} | {row['lambda']:g} | "
            f"{row['block_macro_ovr_auc']:.4f} | {row['block_accuracy']:.4f} | "
            f"{row['ORDER_UNI_auc']:.4f} | {row['ORDER_DOWN_auc']:.4f} | "
            f"{row['ORDER_UP_auc']:.4f} | {row['pass']} |"
        )
    lines.extend([
        "",
        "Chance macro AUC = 0.5. MEAN/POWER/VARIANCE/MOMENTS/START_U must",
        "lie within 0.15 of 0.5. RAW_U_8 must be >= 0.90.",
        "",
        "Orbit intersections were empty before hashes. If RAW_U_8 is",
        "below 0.90, the three orders are not distinct and the task is void.",
        "Do not then permute templates after seeing reservoir scores.",
        "",
    ])
    (results / "u_only_baselines.md").write_text("\n".join(lines), encoding="utf-8")
    return csv_path


def patch_protocol_u_only(audit):
    path = HERE / "PROTOCOL.md"
    text = path.read_text(encoding="utf-8")
    status = "PASS" if audit["overall_pass"] else "FAIL"
    replacement = (
        "### u-only audit result\n"
        "\n"
        f"**{status}.** Block-level test macro OVR AUC:\n"
        "\n"
        f"- MEAN_ONLY {audit['mean_auc']:.4f}\n"
        f"- POWER_ONLY {audit['power_auc']:.4f}\n"
        f"- VARIANCE_ONLY {audit['variance_auc']:.4f}\n"
        f"- MOMENTS {audit['moments_auc']:.4f}\n"
        f"- RAW_U_8 {audit['raw_auc']:.4f}\n"
        f"- START_U {audit['start_auc']:.4f}\n"
        "\n"
        "See `results/u_only_baselines.md`.\n"
    )
    if not audit["overall_pass"]:
        if not audit["raw_ok"]:
            replacement += (
                "\nSTOP. RAW_U_8 is below 0.90. The three orders are not\n"
                "distinct. The task is void. Do not run BSim. Do not rewrite\n"
                "templates after that score. Do not retune K, n, tau, source,\n"
                "clamp, mortality, flow, or layout.\n"
            )
        elif audit["leak"] or not audit["pass_scalar"]:
            replacement += (
                "\nSTOP. A scalar baseline leaked class identity. Do not run BSim.\n"
            )
    else:
        replacement += "\nBSim is authorized to the smoke, then seed-101 scout.\n"
    pattern = r"### u-only audit result\n\nPending\. Hashes above were recorded before this ridge\.\n"
    if re.search(pattern, text) is None:
        text = re.sub(
            r"### u-only audit result\n.*?(?=\n## Ridge)",
            replacement + "\n",
            text,
            count=1,
            flags=re.S,
        )
    else:
        text = re.sub(pattern, replacement + "\n", text, count=1)
    path.write_text(text, encoding="utf-8")


def waveform_gate_untouched():
    path = WAVEFORM_DIR / "results" / "GATE_EVIDENCE.md"
    if not path.exists():
        return True
    text = path.read_text(encoding="utf-8")
    return "Waveform2" not in text


def java_constraints():
    java_path = HERE / "BSimReservoirPlanWaveform2b.java"
    if not java_path.exists():
        return False
    java = java_path.read_text(encoding="utf-8")
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
        "input_ahl_waveform2b_400.txt",
        "NUM_WINDOWS = 400",
    )
    if any(token not in java for token in required):
        return False
    if "new Vector3d(150, 100, 5)" in java or "new Vector3d(850, 400, 5)" in java:
        return False
    return java.count("new BSimChemicalField(") == 4


def read_run(path, arm):
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
        "summary_csv": validate_csv(summary_path, NUM_WINDOWS),
        "auxiliary_csvs": {
            "results.csv": validate_csv(results_path, EXPECTED_AUX_ROWS, True),
            "voxels.csv": validate_csv(voxels_path, EXPECTED_AUX_ROWS, True),
        },
    }


def evaluate_arm_run(run, y_window, y_block):
    matrix = run["matrix"]
    cls = evaluate_window_readout(matrix["X"], y_window, matrix["windows"])
    block_feat = evaluate_block_readout(block_mean_features(matrix["X"]), y_block)
    field = None
    if run["field"] is not None:
        field = evaluate_window_readout(run["field"]["X"], y_window, run["field"]["windows"])
        field["feature_names"] = run["field"]["feature_names"]
        field["n_features"] = len(run["field"]["feature_names"])
        field["block_mean_features"] = evaluate_block_readout(
            block_mean_features(run["field"]["X"]), y_block
        )
    return {
        "path": run["path"],
        "arm": run["arm"],
        "n_features": len(matrix["feature_names"]),
        "feature_names": matrix["feature_names"],
        "classification": cls,
        "block_mean_feature_ridge": block_feat,
        "field_only": field,
        "summary_csv": run["summary_csv"],
        "auxiliary_csvs": run["auxiliary_csvs"],
    }


def surrogate_metrics(driven_dir, y_window, y_block):
    voxels = load_voxel_arrays(Path(driven_dir) / "voxels.csv")
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
    return {
        "unmasked": evaluate_window_readout(unmasked, y_window, list(range(NUM_WINDOWS))),
        "masked": evaluate_window_readout(masked, y_window, list(range(NUM_WINDOWS))),
        "unmasked_block_features": evaluate_block_readout(block_mean_features(unmasked), y_block),
        "masked_block_features": evaluate_block_readout(block_mean_features(masked), y_block),
    }


def available_run(arm, seed):
    path = HERE / "results" / f"waveform2b_{arm}_seed{seed}"
    if (path / "voxels.csv").exists() and (path / "window_summary.csv").exists():
        return path
    return None


def write_scores_csv(path, rows):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def smoke_check(smoke_dir):
    directory = Path(smoke_dir)
    voxels = directory / "voxels.csv"
    summary = directory / "window_summary.csv"
    if not voxels.exists():
        return {"pass": False, "reason": f"missing {voxels}"}
    ahl = ahl_finite_nonneg(voxels)
    java_ok = java_constraints()
    return {
        "pass": ahl["pass"] and java_ok,
        "ahl": ahl,
        "java_ok": java_ok,
        "summary_rows": validate_csv(summary, None and 0) if False else (
            validate_csv(summary, 2) if summary.exists() else {}
        ),
        "dish": "1000x500x10",
        "source": "(500,250,5)",
        "FLOW_SPEED": 0.0,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--u-only", action="store_true")
    parser.add_argument("--smoke", default=str(HERE / "results" / "waveform2b_smoke"))
    parser.add_argument("--seeds", nargs="+", type=int, default=None)
    args = parser.parse_args()

    frozen = verify_frozen_inputs()
    y_block = block_labels(frozen["y_window"])
    print("FROZEN_INPUTS")
    print(f"class_vector={frozen['class_vector']}")
    print(f"phase_vector={frozen['phase_vector']}")
    print(f"class_vector_sha256={frozen['class_sha']}")
    print(f"phase_vector_sha256={frozen['phase_sha']}")
    print(f"u_sha256={frozen['u_sha']}")
    print(f"waveform_GATE_EVIDENCE_untouched={waveform_gate_untouched()}")

    audit = u_only_baselines(frozen)
    write_u_only(audit)
    patch_protocol_u_only(audit)
    print("U_ONLY")
    for row in audit["rows"]:
        print(
            f"{row['baseline']} auc={row['block_macro_ovr_auc']:.4f} "
            f"acc={row['block_accuracy']:.4f} pass={row['pass']}"
        )
    print(f"u_only_overall={audit['overall_pass']}")
    if args.u_only:
        print(f"U_ONLY_{'PASS' if audit['overall_pass'] else 'FAIL'}")
        if not audit["overall_pass"]:
            raise SystemExit(1)
        return

    if not audit["overall_pass"]:
        print("STOP: u-only audit FAIL. Do not interpret BSim as Waveform2b skill.")
        raise SystemExit(1)

    seeds = tuple(args.seeds) if args.seeds else SEEDS
    rows = []
    evaluated = {}
    occupancy_rows = {}
    for arm in ("brownian", "silent", "driven"):
        evaluated[arm] = []
        for seed in seeds:
            path = available_run(arm, seed)
            if path is None:
                continue
            run = read_run(path, arm)
            scored = evaluate_arm_run(run, frozen["y_window"], y_block)
            evaluated[arm].append((seed, scored))
            occ = None
            stats = summary_stats(path)
            if arm == "driven":
                voxels = load_voxel_arrays(path / "voxels.csv")
                occ = occupancy(voxels, frozen["u"])
                occupancy_rows[seed] = {**occ, **stats}
                sur = surrogate_metrics(path, frozen["y_window"], y_block)
                scored["surrogates"] = sur
            cls = scored["classification"]
            field = scored["field_only"]
            rows.append({
                "arm": arm,
                "seed": seed,
                "n_features": scored["n_features"],
                "lambda": cls["lambda"],
                "block_macro_ovr_auc": cls["block_macro_ovr_auc"],
                "block_accuracy": cls["block_accuracy"],
                "window_macro_ovr_auc": cls["test_window_macro_ovr_auc"],
                "window_accuracy": cls["test_window_accuracy"],
                "block_mean_feature_auc": scored["block_mean_feature_ridge"]["block_macro_ovr_auc"],
                "ORDER_UNI_auc": cls["block_per_class_ovr_auc"]["ORDER_UNI"],
                "ORDER_DOWN_auc": cls["block_per_class_ovr_auc"]["ORDER_DOWN"],
                "ORDER_UP_auc": cls["block_per_class_ovr_auc"]["ORDER_UP"],
                "occupancy": occ["occupancy"] if occ else "",
                "mean_R": occ["mean_R"] if occ else "",
                "r_meanR_u": occ["r_meanR_u"] if occ else "",
                "mean_AHL": occ["mean_AHL"] if occ else "",
                "mean_pop": stats["mean_pop"],
                "field_block_auc": field["block_macro_ovr_auc"] if field else "",
                "unmasked_block_auc": (
                    scored.get("surrogates", {}).get("unmasked", {}).get("block_macro_ovr_auc", "")
                ),
                "masked_block_auc": (
                    scored.get("surrogates", {}).get("masked", {}).get("block_macro_ovr_auc", "")
                ),
            })
            print(
                f"{arm} seed={seed} block_auc={cls['block_macro_ovr_auc']:.4f} "
                f"window_auc={cls['test_window_macro_ovr_auc']:.4f}"
            )

    write_scores_csv(HERE / "results" / "waveform2b_scores.csv", rows)
    (HERE / "results" / "waveform2b_gate_evidence.json").write_text(
        json.dumps({
            "u_only": audit,
            "occupancy": occupancy_rows,
            "rows": rows,
            "java_ok": java_constraints(),
            "waveform_gate_untouched": waveform_gate_untouched(),
        }, indent=2, allow_nan=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    print(f"scores={HERE / 'results' / 'waveform2b_scores.csv'}")


if __name__ == "__main__":
    main()
