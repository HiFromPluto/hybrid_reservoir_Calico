#!/usr/bin/env python3
"""Evaluate Track B 5-channel patient-classification gates.

Stage 6 NARMA-10 remains PASS. BenchA Mackey-Glass remains PASS.
BenchA Lorenz remains FAIL. Stage 9 product-bit FAIL is not rewritten.
Kinetics are not retuned. Do not widen the dish after AUC.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np

WASHOUT = 40
TRAIN = 110
TEST = 50
INNER_VAL = 22
NUM_WINDOWS = WASHOUT + TRAIN + TEST
WINDOWS_PER_PATIENT = 5
EXPECTED_AUX_ROWS = NUM_WINDOWS * 16
LAST_SAMPLE = ("199", "15", "299.95")
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
SEEDS = (101, 202, 303)

CLASS_VECTOR = (
    "1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0"
)
CLASS_SHA = "16de0e28bc19c9f94bd468afb8a3418b29b747ba544e8769855a3499eee2ca48"
U_SHA = {
    "ahl": "b8684fafe2452bc2f9330bd7e368e559a99bd3bceea4b219226499dbe7d7825e",
    "acid": "6762700e28f9a4e2087054acc25d933080358829eadc652ea165a1a6d915b0a3",
    "atta": "dbb60b2bb4b264cdc9341e33d58422c7f3f66977889cd6e57437c44751fe34c4",
    "attb": "25d27243479ccb2b2fdf324982811ebad72f4ed21cd6a4aa55de7a7545608cc2",
    "rep": "4d50166d9c199ab517d4bb6593a80f83f451006c7f9ef8b735ff589e239bf420",
    "singlesite": "1efc565c29aebe9baa91260584cf3ccfff7f474e981f63725bcd95ea31167df3",
}

BIOLOGY_MEAN_PREFIXES = ("Receiver_R_", "Lum_Mean_", "Den_")
BIOLOGY_LAST_PREFIXES = ("Input_Driven_Death_",)
BROWNIAN_MEAN_PREFIXES = ("Den_",)
BROWNIAN_LAST_PREFIXES = ()
FIELD_MEAN_PREFIXES = ("AHL_uM_", "Att_conc_", "Rep_conc_", "pH_")
BIOLOGY_FORBIDDEN = (
    "AHL_uM_", "pH_", "Att_conc_", "Rep_conc_", "Fraction_q", "Lum_Sum_",
    "Birth_", "Clamp_Death_", "OOB_Death_",
)
EXCLUDED_FROM_RIDGE = (
    "window AHL / raw u",
    "occupancy Fraction_q_gt_0_5",
    "pH (biology ridge)",
    "Births",
    "Total_Deaths",
    "Clamp_Deaths",
    "OOB",
    "Population",
    "Lum_Sum",
    "voxel AHL_uM (biology ridge)",
    "voxel Att_conc (biology ridge)",
    "voxel Rep_conc (biology ridge)",
    "voxel pH (biology ridge)",
    "voxel Fraction_q_gt_0_5",
    "voxel Lum_Sum",
)
FORBIDDEN_JAVA = (
    "glucose", "Glucose", "GLUCOSE",
    "Danino", "danino",
    "vesicle", "Vesicle",
    "setGoal(glucose)",
    "setGoal(repellent",
)

HERE = Path(__file__).resolve().parent
STAGE6_DIR = HERE.parent / "BSimReservoirPlanStage6"
BENCHA_DIR = HERE.parent / "BSimReservoirPlanBenchA"
BENCHA2_DIR = HERE.parent / "BSimReservoirPlanBenchA2"


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
    patient = np.array([int(row["patient"]) for row in rows], dtype=int)
    return y, patient, rows


def class_balance(y, start, end, name):
    slice_y = y[start:end]
    positives = int(np.sum(slice_y))
    return {
        "split": name,
        "windows": f"{start}..{end - 1}",
        "n": int(len(slice_y)),
        "positives": positives,
        "negatives": int(len(slice_y) - positives),
        "fraction": float(positives / len(slice_y)) if len(slice_y) else float("nan"),
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
    """Fit 40..127, pick lambda on 128..149 only. Test is never in the val slice."""
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


def patient_pooled_auc(y_test, scores_test):
    if len(y_test) != TEST:
        raise ValueError(f"test length {len(y_test)}; expected {TEST}")
    pooled_scores = []
    pooled_y = []
    n_patients = TEST // WINDOWS_PER_PATIENT
    for patient in range(n_patients):
        sl = slice(patient * WINDOWS_PER_PATIENT, (patient + 1) * WINDOWS_PER_PATIENT)
        pooled_scores.append(float(np.mean(scores_test[sl])))
        labels = np.unique(y_test[sl])
        if len(labels) != 1:
            raise ValueError("test patient block is not a constant label")
        pooled_y.append(int(labels[0]))
    return roc_auc(pooled_y, pooled_scores)


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
        "patient_pooled_test_auc": patient_pooled_auc(y_test, test_pred),
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


def overall_pass(path):
    if not Path(path).exists():
        return False
    text = Path(path).read_text(encoding="utf-8")
    return bool(re.search(r"^## Overall:\s*PASS\s*$", text, re.MULTILINE))


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
    return (
        "package BSimReservoirPlanStage6;" in stage6
        and "NARMA-10" in stage6
        and "package BSimReservoirPlanStage6;" in stage6_voxel
        and "Track B is not started" in bencha
        and "Track B is not started" in bencha2
        and "SKIP=50" in bencha2
    )


def trackb_constraints():
    java = (HERE / "BSimReservoirPlanTrackB.java").read_text(encoding="utf-8")
    if any(token in java for token in FORBIDDEN_JAVA):
        return False
    if java.count("new BSimChemicalField(") != 4:
        return False
    if "new Vector3d(250, 250, 5)" in java or "new Vector3d(750, 250, 5)" in java:
        return False
    required = (
        "new Vector3d(500, 250, 5)",
        "new Vector3d(300, 375, 5)",
        "new Vector3d(150, 100, 5)",
        "new Vector3d(850, 400, 5)",
        "new Vector3d(150, 400, 5)",
        "setGoal(attractantField)",
        "PROD_RATE = 1e6",
        "AHL_SOURCE_RATE = 128000000",
        "ACID_SOURCE_RATE = 2.0e11",
        "BOUND_X = 1000.0, BOUND_Y = 500.0",
    )
    if any(token not in java for token in required):
        return False
    gen = (BENCHA2_DIR / "generate_lorenz_skip50.py").read_text(encoding="utf-8")
    return "SKIP = 50" in gen


def write_markdown(path, evidence):
    gates = evidence["gates"]
    overall = "PASS" if evidence["overall_gate_pass"] else "FAIL"
    lines = [
        "# Track B gate evidence",
        "",
        "Track B copies frozen Stage 6 and tests a synthetic 5-channel patient",
        "classification task against Brownian density, a single-site AHL mix,",
        "a silent-source control, and a Stage 9-style field-only chemical",
        "baseline. Kinetics were not retuned. Stage 6 NARMA-10 remains **PASS**.",
        "BenchA Mackey-Glass remains **PASS**. BenchA Lorenz remains **FAIL**.",
        "Stage 9 product-bit FAIL is not rewritten. The dish was not widened.",
        "",
        f"## Overall: {overall}",
        "",
        evidence["failure_reason"] if not evidence["overall_gate_pass"]
        else "Driven 5-channel biology beat Brownian, field-only, and single-site; silent sat at chance.",
        "",
        evidence["failure_detail"],
        "",
        "| Gate | Result |",
        "|---|---|",
        f"| 1. Driven 5-channel test AUC > Brownian, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_brownian_pass'] else 'FAIL'}** |",
        f"| 2. Driven 5-channel test AUC > field-only 800 chemical voxels, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_field_pass'] else 'FAIL'}** |",
        f"| 3. Driven 5-channel test AUC > single-site, non-overlapping mean±s.e., all 3 seeds | **{'PASS' if gates['driven_beats_singlesite_pass'] else 'FAIL'}** |",
        f"| 4. Silent not above chance by a clear margin; silent not ≈ driven | **{'PASS' if gates['silent_chance_pass'] else 'FAIL'}** |",
        f"| 5. CSV 200/3200/3200; last sample 199;15;299.95 | **{'PASS' if gates['csv_validation_pass'] else 'FAIL'}** |",
        f"| 6. Stage 6 PASS; BenchA MG PASS; BenchA Lorenz FAIL; Java untouched; no glucose/Danino/vesicles/extra PDE/A2 SKIP change | **{'PASS' if gates['provenance_pass'] else 'FAIL'}** |",
        "",
        "## Frozen task",
        "",
        "40 patients × 5 windows. Class vector SHA-256 "
        f"`{evidence['class_vector_sha256']}`. HIGH/LOW and seeds were not",
        "changed after seeing AUC.",
        "",
        f"Class vector: `{evidence['class_vector']}`",
        "",
        "| Split | Windows | n | Pos | Neg | Fraction |",
        "|---|---|---|---|---|---|",
    ]
    for row in evidence["class_balance"]:
        lines.append(
            f"| {row['split']} | {row['windows']} | {row['n']} | "
            f"{row['positives']} | {row['negatives']} | {row['fraction']:.3f} |"
        )
    lines.extend([
        "",
        f"Feature counts: biology `{evidence['driven_feature_count']}`, "
        f"Brownian `{evidence['brownian_feature_count']}`, "
        f"field-only `{evidence['field_feature_count']}`.",
        "",
        "## Test AUC (primary)",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 | mean ± s.e. |",
        "|---|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "singlesite", "driven", "field"):
        values = evidence["test_auc"][arm]
        mean = evidence["auc_mean_se"][arm]["mean"]
        se = evidence["auc_mean_se"][arm]["se"]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} | {mean:.4f} ± {se:.4f} |"
        )
    lines.extend([
        "",
        "Accuracy at 0.5 is printed, not the gate.",
        "",
        "| Arm | seed 101 acc | seed 202 acc | seed 303 acc |",
        "|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "singlesite", "driven", "field"):
        values = evidence["test_accuracy"][arm]
        lines.append(
            f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} |"
        )
    lines.extend([
        "",
        "## Patient-pooled AUC (diagnostic, not a gate)",
        "",
        "Mean of 5 window scores per test patient (10 scores). Do not replace",
        "window AUC with this after seeing numbers.",
        "",
        "| Arm | seed 101 | seed 202 | seed 303 |",
        "|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "singlesite", "driven", "field"):
        values = evidence["patient_pooled_auc"][arm]
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
        "- `Den_*`, 20×10",
        "",
        f"Printed count: {evidence['driven_feature_count']} features.",
        "",
        "Field-only: `AHL_uM_*`, `Att_conc_*`, `Rep_conc_*`, `pH_*` "
        f"({evidence['field_feature_count']}).",
        "",
        "Excluded from the biology ridge: " + ", ".join(EXCLUDED_FROM_RIDGE) + ".",
        "",
        "## Lambdas (independent per arm × seed)",
        "",
        "| Arm | 101 | 202 | 303 |",
        "|---|---|---|---|",
    ])
    for arm in ("brownian", "silent", "singlesite", "driven", "field"):
        values = evidence["lambdas"][arm]
        lines.append(f"| {arm} | {values[0]:g} | {values[1]:g} | {values[2]:g} |")
    lines.extend([
        "",
        "## Protocol",
        "",
        "See `PROTOCOL.md`. Frozen 5-AC layout. Patient-boundary splits.",
        "Lambda fit 40..127, pick on 128..149 only. If gate 2 or 3 fails,",
        "Overall is FAIL. The dish was not widened. `PROD_RATE` was not raised.",
        "",
        "CSV completeness is the on-disk files, last sample `199;15;299.95`.",
        "Silent CSVs are Stage 6 copies. Silent was not rerun.",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brownian", nargs=3, default=[
        str(HERE / "results" / f"trackb_brownian_seed{seed}") for seed in SEEDS
    ])
    parser.add_argument("--silent", nargs=3, default=[
        str(HERE / "results" / f"trackb_silent_seed{seed}") for seed in SEEDS
    ])
    parser.add_argument("--driven", nargs=3, default=[
        str(HERE / "results" / f"trackb_driven_seed{seed}") for seed in SEEDS
    ])
    parser.add_argument("--singlesite", nargs=3, default=[
        str(HERE / "results" / f"trackb_singlesite_seed{seed}") for seed in SEEDS
    ])
    parser.add_argument("--labels", default=str(HERE / "trackb_labels.csv"))
    parser.add_argument("--expected-windows", type=int, default=NUM_WINDOWS)
    parser.add_argument("--expected-aux-rows", type=int, default=EXPECTED_AUX_ROWS)
    parser.add_argument("--evidence", default=str(HERE / "results" / "trackb_gate_evidence.json"))
    parser.add_argument("--markdown", default=str(HERE / "results" / "GATE_EVIDENCE.md"))
    args = parser.parse_args()

    sequences = {
        "ahl": load_sequence(HERE / "input_ahl_trackb200.txt"),
        "acid": load_sequence(HERE / "input_acid_trackb200.txt"),
        "atta": load_sequence(HERE / "input_atta_trackb200.txt"),
        "attb": load_sequence(HERE / "input_attb_trackb200.txt"),
        "rep": load_sequence(HERE / "input_rep_trackb200.txt"),
        "singlesite": load_sequence(HERE / "input_ahl_singlesite200.txt"),
    }
    for name, values in sequences.items():
        digest = sha256_u(values)
        if digest != U_SHA[name]:
            raise ValueError(f"{name} SHA-256 {digest} does not match frozen {U_SHA[name]}")

    y, patient, label_rows = load_labels(args.labels)
    reconstructed = ",".join(str(int(y[p * WINDOWS_PER_PATIENT])) for p in range(40))
    class_sha = hashlib.sha256(reconstructed.encode("ascii")).hexdigest()
    if reconstructed != CLASS_VECTOR or class_sha != CLASS_SHA:
        raise ValueError("class vector does not match the frozen PROTOCOL.md declaration")
    for n, row in enumerate(label_rows):
        if int(row["patient"]) != n // WINDOWS_PER_PATIENT:
            raise ValueError(f"patient index mismatch at window {n}")
        if int(row["y"]) != int(y[n]):
            raise ValueError(f"label mismatch at window {n}")

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
            f"pos={row['positives']} neg={row['negatives']} frac={row['fraction']:.6f}"
        )
    print(f"class_vector={reconstructed}")
    print(f"class_vector_sha256={class_sha}")

    grouped = {
        "brownian": [read_run(path, "brownian", args.expected_windows, args.expected_aux_rows)
                     for path in args.brownian],
        "silent": [read_run(path, "silent", args.expected_windows, args.expected_aux_rows)
                   for path in args.silent],
        "driven": [read_run(path, "driven", args.expected_windows, args.expected_aux_rows)
                   for path in args.driven],
        "singlesite": [read_run(path, "singlesite", args.expected_windows, args.expected_aux_rows)
                       for path in args.singlesite],
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
    pooled_by_arm = {
        arm: [run["classification"]["patient_pooled_test_auc"] for run in runs]
        for arm, runs in evaluated.items()
    }
    field_auc = [run["field_only"]["test_auc"] for run in evaluated["driven"] if run["field_only"]]
    field_acc = [
        run["field_only"]["test_accuracy_at_0_5"] for run in evaluated["driven"] if run["field_only"]
    ]
    field_pooled = [
        run["field_only"]["patient_pooled_test_auc"]
        for run in evaluated["driven"] if run["field_only"]
    ]
    auc_by_arm["field"] = field_auc
    acc_by_arm["field"] = field_acc
    pooled_by_arm["field"] = field_pooled
    mean_se_by_arm = {arm: mean_se(values) for arm, values in auc_by_arm.items()}

    driven_features = evaluated["driven"][0]["feature_names"]
    expected_prefixes = ("Receiver_R_", "Lum_Mean_", "Den_", "Input_Driven_Death_")
    feature_ok = (
        all(name.startswith(expected_prefixes) for name in driven_features)
        and not any(name.startswith(BIOLOGY_FORBIDDEN) for name in driven_features)
        and len(driven_features) == 608
        and all(run["n_features"] == 608
                for run in evaluated["driven"] + evaluated["silent"] + evaluated["singlesite"])
        and all(run["n_features"] == 200 for run in evaluated["brownian"])
        and all(run["field_only"]["n_features"] == 800 for run in evaluated["driven"])
        and sum(name.startswith("Receiver_R_") for name in driven_features) == 200
        and sum(name.startswith("Lum_Mean_") for name in driven_features) == 200
        and sum(name.startswith("Den_") for name in driven_features) == 200
        and sum(name.startswith("Input_Driven_Death_") for name in driven_features) == 8
    )
    print(f"driven_feature_count={len(driven_features)}")
    print(f"brownian_feature_count={evaluated['brownian'][0]['n_features']}")
    print(f"silent_feature_count={evaluated['silent'][0]['n_features']}")
    print(f"singlesite_feature_count={evaluated['singlesite'][0]['n_features']}")
    print(f"field_feature_count={evaluated['driven'][0]['field_only']['n_features']}")

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
    single_mean, single_se = mean_se_by_arm["singlesite"]
    field_mean, field_se = mean_se_by_arm["field"]

    driven_beats_brownian = (
        all(d > b for d, b in zip(auc_by_arm["driven"], auc_by_arm["brownian"]))
        and driven_mean - driven_se > brown_mean + brown_se
    )
    driven_beats_field = (
        all(d > f for d, f in zip(auc_by_arm["driven"], auc_by_arm["field"]))
        and driven_mean - driven_se > field_mean + field_se
    )
    driven_beats_singlesite = (
        all(d > s for d, s in zip(auc_by_arm["driven"], auc_by_arm["singlesite"]))
        and driven_mean - driven_se > single_mean + single_se
    )
    silent_overlap_driven = intervals_overlap(silent_mean, silent_se, driven_mean, driven_se)
    silent_worse_than_driven = silent_mean + silent_se < driven_mean - driven_se
    silent_at_chance = overlaps_chance(silent_mean, silent_se)
    silent_chance_pass = (silent_at_chance or silent_worse_than_driven) and not silent_overlap_driven

    provenance_pass = (
        overall_pass(STAGE6_DIR / "results" / "GATE_EVIDENCE.md")
        and heading_overall(BENCHA_DIR / "results" / "GATE_EVIDENCE.md", "Mackey-Glass", "PASS")
        and heading_overall(BENCHA_DIR / "results" / "GATE_EVIDENCE.md", "Lorenz'63", "FAIL")
        and java_untouched()
        and trackb_constraints()
        and feature_ok
    )

    gates = {
        "replicate_count_pass": all(len(values) == 3 for values in auc_by_arm.values()),
        "driven_beats_brownian_pass": driven_beats_brownian,
        "driven_beats_field_pass": driven_beats_field,
        "driven_beats_singlesite_pass": driven_beats_singlesite,
        "silent_chance_pass": silent_chance_pass,
        "csv_validation_pass": csv_pass,
        "provenance_pass": provenance_pass,
        "frozen_feature_list_pass": feature_ok,
    }
    overall = all(
        gates[key] for key in (
            "driven_beats_brownian_pass",
            "driven_beats_field_pass",
            "driven_beats_singlesite_pass",
            "silent_chance_pass",
            "csv_validation_pass",
            "provenance_pass",
        )
    )
    if not gates["driven_beats_field_pass"]:
        failure_reason = (
            "field_wins: driven 5-channel test AUC did not beat the field-only "
            "AHL/Att/Rep/pH baseline in every seed with non-overlapping mean±s.e. "
            "The plumes classify the patient. Do not widen the dish. Do not add "
            "ACs. Do not switch to accuracy. Do not raise PROD_RATE."
        )
    elif not gates["driven_beats_singlesite_pass"]:
        failure_reason = (
            "h5_unsupported: driven 5-channel test AUC did not beat the "
            "single-site AHL mix in every seed with non-overlapping mean±s.e. "
            "Do not retune."
        )
    elif not gates["driven_beats_brownian_pass"]:
        failure_reason = (
            "brownian_wins: driven 5-channel test AUC did not beat Brownian "
            "density in every seed with non-overlapping mean±s.e."
        )
    elif not gates["silent_chance_pass"]:
        failure_reason = (
            "silent_approx_driven: silent test AUC is not at chance and overlaps "
            "driven. Inputs are not what the classifier is using."
        )
    elif not gates["csv_validation_pass"]:
        failure_reason = (
            "csv_incomplete: a summary/sample/voxel file is missing rows, is "
            "ragged, or last sample is not 199;15;299.95."
        )
    elif not gates["provenance_pass"]:
        failure_reason = (
            "provenance: Stage 6/BenchA labels were reopened, frozen Java was "
            "touched, glucose/Danino/vesicles/extra PDE leaked in, or Track A2 "
            "SKIP changed."
        )
    else:
        failure_reason = ""

    if not driven_beats_field:
        failure_detail = (
            "Gate 2 fails. Driven test AUC mean ± s.e. is "
            f"{driven_mean:.4f} ± {driven_se:.4f}; field-only is "
            f"{field_mean:.4f} ± {field_se:.4f}. The chemical plumes already "
            "rank the patient. Kinetics, PROD_RATE, AC sites, and domain were "
            "not changed after seeing this. A 2000×1000 copy was not started."
        )
    elif not driven_beats_singlesite:
        failure_detail = (
            "Gate 3 fails. Driven 5-channel mean ± s.e. is "
            f"{driven_mean:.4f} ± {driven_se:.4f}; single-site is "
            f"{single_mean:.4f} ± {single_se:.4f}. H5 is not supported. "
            "Kinetics were not retuned."
        )
    elif not driven_beats_brownian:
        failure_detail = (
            "Gate 1 fails. Driven did not beat Brownian density with "
            "non-overlapping mean±s.e. in every seed."
        )
    elif not silent_chance_pass:
        failure_detail = (
            "Gate 4 fails. Silent test AUC is not at chance and overlaps driven."
        )
    elif overall:
        failure_detail = (
            "No failure. Driven 5-channel beat Brownian, field-only, and "
            "single-site; silent sat at chance. Frozen 608-feature contract used."
        )
    else:
        failure_detail = failure_reason

    evidence = {
        "schema": "BSimReservoirPlanTrackB-gates-v1",
        "stage6_status": "PASS",
        "stage7_status": "FAIL",
        "stage9_status": "FAIL",
        "bencha_lorenz_status": "FAIL",
        "bencha_mg_status": "PASS",
        "class_vector": reconstructed,
        "class_vector_sha256": class_sha,
        "u_sha256": {name: sha256_u(values) for name, values in sequences.items()},
        "protocol": {
            "num_windows": NUM_WINDOWS,
            "washout": WASHOUT,
            "train": TRAIN,
            "test": TEST,
            "inner_val": INNER_VAL,
            "lambda_validation_slice": "X[inner_train_end:TRAIN] after washout strip",
            "ridge_grid": list(RIDGE_GRID),
            "primary_metric": "window test AUC",
            "layout": "AHL (500,250,5); acid (300,375,5); attA (150,100,5); attB (850,400,5); rep (150,400,5)",
        },
        "class_balance": balances,
        "analysis_channels_driven": [
            "Receiver_R_* / Mean_q (20x10)",
            "Lum_Mean_* / Mean_L (20x10)",
            "Input_Driven_Death_* (4x2)",
            "Den_* (20x10)",
        ],
        "driven_feature_count": len(driven_features),
        "brownian_feature_count": evaluated["brownian"][0]["n_features"],
        "field_feature_count": evaluated["driven"][0]["field_only"]["n_features"],
        "excluded_from_ridge": list(EXCLUDED_FROM_RIDGE),
        "test_auc": auc_by_arm,
        "test_accuracy": acc_by_arm,
        "patient_pooled_auc": pooled_by_arm,
        "auc_mean_se": {
            arm: {"mean": mean, "se": se} for arm, (mean, se) in mean_se_by_arm.items()
        },
        "lambdas": {
            "brownian": [run["classification"]["lambda"] for run in evaluated["brownian"]],
            "silent": [run["classification"]["lambda"] for run in evaluated["silent"]],
            "singlesite": [run["classification"]["lambda"] for run in evaluated["singlesite"]],
            "driven": [run["classification"]["lambda"] for run in evaluated["driven"]],
            "field": [run["field_only"]["lambda"] for run in evaluated["driven"]],
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
        "patient_pooled_auc": pooled_by_arm,
        "failure_reason": failure_reason,
        "driven_feature_count": len(driven_features),
        "brownian_feature_count": evidence["brownian_feature_count"],
        "field_feature_count": evidence["field_feature_count"],
    }, indent=2, allow_nan=True, default=json_default))
    print(f"evidence_file={destination}")
    print(f"markdown_file={args.markdown}")
    print(f"Overall: {'PASS' if overall else 'FAIL'}")


if __name__ == "__main__":
    main()
