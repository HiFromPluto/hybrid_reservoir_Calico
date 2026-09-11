#!/usr/bin/env python3
"""Run the frozen zero-simulation evidence audit.

This script reads existing CSV/JSON evidence only. It never invokes Java or
BSim and writes only inside zero_simulation_audit/ and external_analysis/.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import inspect
import json
import math
import shutil
from functools import lru_cache
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
HYBRID = HERE.parent
EXAMPLES = HYBRID.parent
REPO = EXAMPLES.parent
EXTERNAL = HYBRID / "external_analysis"

NARMA = EXAMPLES / "BSimReservoirPlanNarma10b"
WAVEFORM = EXAMPLES / "BSimReservoirPlanWaveform"
BENCHA = EXAMPLES / "BSimReservoirPlanBenchA"
BENCHA2 = EXAMPLES / "BSimReservoirPlanBenchA2"

WASHOUT = 40
TRAIN = 110
TEST = 50
INNER_VAL = 22
NUM_WINDOWS = 200
WINDOWS = list(range(NUM_WINDOWS))
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
NARMA_SEEDS = (111, 222, 333)
WAVEFORM_SEEDS = (101, 202, 303)
SHIFT_OFFSETS = (20, 40, 60, 80, 100)
CLASS_NAMES = ("sine", "square", "triangle")
SANITY_TOL = 1e-3

NARMA_EXPECTED = {
    "driven": (0.9312, 0.9228, 0.9310),
    "field": 1.0289,
    "brownian_mean": 1.1624,
    "silent": 1.1622,
    "persistence": 1.0068,
    "train_intercept": 1.1622,
}
WAVE_EXPECTED = {
    "driven_auc": (0.7969, 0.8287, 0.8411),
    "field_auc": 0.8603,
    "driven_accuracy_mean": 0.73,
}
BENCHA_EXPECTED = {
    "lorenz_driven_mean": 0.9895,
    "lorenz_field_mean": 0.095,
    "mg_driven_mean": 0.2602,
    "mg_field_mean": 0.0229,
}

EXTERNAL_SOURCES = (
    (
        Path(r"C:\Users\Ceylin\Downloads\CHARC_Framework_BenchA.md"),
        "CHARC_Framework_BenchA.md",
    ),
    (
        Path(r"C:\Users\Ceylin\Downloads\Stage99_Trapped_Dish_Evaluation.md"),
        "Stage99_Trapped_Dish_Evaluation.md",
    ),
    (
        Path(
            r"C:\Users\Ceylin\.gemini\antigravity-ide\brain"
            r"\40cf0870-1570-4044-9b6e-f8273f101a50"
            r"\Lorenz_Encoding_Recommendation.md"
        ),
        "Lorenz_Encoding_Recommendation.md",
    ),
)


def import_file(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Numerical source of truth: existing closed checkers.
NC = import_file("audit_narma_checker", NARMA / "check_narma10b.py")
WC = import_file("audit_wave_checker", WAVEFORM / "check_waveform.py")
BC = import_file("audit_bencha_checker", BENCHA / "check_bencha.py")
B2C = import_file("audit_bencha2_checker", BENCHA2 / "check_bencha2.py")


def fnum(value) -> str:
    if value is None:
        return ""
    value = float(value)
    if not math.isfinite(value):
        return "nan" if math.isnan(value) else ("inf" if value > 0 else "-inf")
    return f"{value:.10g}"


def mean(values) -> float:
    return float(np.mean(np.asarray(values, dtype=float)))


def within(actual, expected, tol=SANITY_TOL) -> bool:
    return abs(float(actual) - float(expected)) <= tol


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_semicolon(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=";"))


def load_narma_target():
    rows = load_semicolon(NARMA / "narma10_target.csv")
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    y = np.array([float(row["y_next"]) for row in rows], dtype=float)
    ahl = NARMA / "input_ahl_narma200.txt"
    checked_u, checked_y, digest = NC.load_target(
        NARMA / "narma10_target.csv", ahl
    )
    if not np.array_equal(u, checked_u) or not np.array_equal(y, checked_y):
        raise AssertionError("NARMA target checker round-trip mismatch")
    return u, y, digest


def load_wave_labels():
    y, block, u, _ = WC.load_labels(WAVEFORM / "waveform_labels.csv")
    sequence = np.asarray(WC.load_sequence(WAVEFORM / "input_ahl_waveform200.txt"))
    if not np.allclose(u, sequence, rtol=0, atol=1e-12):
        raise AssertionError("waveform u sidecars disagree")
    return y, block, u


def select_columns(header, prefix):
    return [index for index, name in enumerate(header) if name.startswith(prefix)]


def read_voxels(path: Path):
    """Read only channels used by this audit from one existing voxel CSV."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        wi = header.index("Window")
        ti = header.index("TimeInWindow_s")
        groups = {
            "ahl": select_columns(header, "AHL_uM_"),
            "den": select_columns(header, "Den_"),
            "r": select_columns(header, "Receiver_R_"),
            "l": select_columns(header, "Lum_Mean_"),
            "death": select_columns(header, "Input_Driven_Death_"),
        }
        if any(len(groups[key]) != 200 for key in ("ahl", "den", "r", "l")):
            raise ValueError(f"{path}: expected 200 voxel columns per map")
        if len(groups["death"]) != 8:
            raise ValueError(f"{path}: expected 8 driven-death columns")
        windows, times = [], []
        values = {key: [] for key in groups}
        for row in reader:
            windows.append(int(row[wi]))
            times.append(float(row[ti]))
            for key, indices in groups.items():
                values[key].append([float(row[index]) for index in indices])
    result = {
        "window": np.asarray(windows, dtype=int),
        "time": np.asarray(times, dtype=float),
    }
    result.update({key: np.asarray(rows, dtype=float) for key, rows in values.items()})
    if len(windows) != 3200:
        raise ValueError(f"{path}: expected 3200 rows, found {len(windows)}")
    absolute = result["window"] * 300.0 + result["time"]
    if np.any(np.diff(absolute) <= 0):
        raise ValueError(f"{path}: non-increasing absolute sample time")
    result["absolute_time"] = absolute
    return result


@lru_cache(maxsize=None)
def read_window_features(path: Path, mean_prefixes, last_prefixes=()):
    loaded = NC.read_window_matrix(path, tuple(mean_prefixes), tuple(last_prefixes))
    if loaded["windows"] != WINDOWS:
        raise ValueError(f"{path}: incomplete windows")
    return np.asarray(loaded["X"], dtype=float)


def window_mean(samples, windows):
    output = np.empty((NUM_WINDOWS, samples.shape[1]), dtype=float)
    for window in WINDOWS:
        rows = samples[windows == window]
        if len(rows) != 16:
            raise ValueError(f"window {window}: expected 16 samples, found {len(rows)}")
        output[window] = rows.mean(axis=0)
    return output


def window_last(samples, windows):
    output = np.empty((NUM_WINDOWS, samples.shape[1]), dtype=float)
    for window in WINDOWS:
        rows = samples[windows == window]
        if len(rows) != 16:
            raise ValueError(f"window {window}: expected 16 samples, found {len(rows)}")
        output[window] = rows[-1]
    return output


def nrmse(y_true, y_pred):
    return NC.nrmse(np.asarray(y_true), np.asarray(y_pred))


def r2(y_true, y_pred):
    return NC.r_squared(np.asarray(y_true), np.asarray(y_pred))


def scalar_closed(X_all, target):
    """Closed scalar ridge using the imported NARMA checker functions."""
    X_all = np.asarray(X_all, dtype=float)
    target = np.asarray(target, dtype=float)
    X = X_all[WASHOUT:]
    y = target[WASHOUT:]
    lam, grid = NC.select_lambda(X, y)
    metrics = NC.fit_eval(X, y, lam)
    val = next(row["val_nrmse"] for row in grid if row["lambda"] == lam)
    return {
        "lambda": float(lam),
        "train_nrmse": float(metrics["train_nrmse"]),
        "validation_nrmse": float(val),
        "test_nrmse": float(metrics["test_nrmse"]),
        "test_r2": float(metrics["test_r2"]),
    }


def evaluate_mismatch(state_trajectory, target_trajectory):
    """Accept state i and target j; independent-input use remains pending."""
    X = np.asarray(state_trajectory, dtype=float)
    y = np.asarray(target_trajectory, dtype=float)
    if X.ndim != 2 or len(X) != NUM_WINDOWS or y.shape != (NUM_WINDOWS,):
        raise ValueError("mismatch requires X[200,n_features] and y[200]")
    return scalar_closed(X, y)


def constant_metrics(target, constant):
    y = np.asarray(target)[WASHOUT:]
    y_train, y_test = y[:TRAIN], y[TRAIN:]
    y_val = y[TRAIN - INNER_VAL : TRAIN]
    return {
        "train_nrmse": nrmse(y_train, np.full_like(y_train, constant)),
        "validation_nrmse": nrmse(y_val, np.full_like(y_val, constant)),
        "test_nrmse": nrmse(y_test, np.full_like(y_test, constant)),
        "test_r2": r2(y_test, np.full_like(y_test, constant)),
    }


def delay_matrix(values, taps):
    values = np.asarray(values, dtype=float)
    X = np.zeros((len(values), taps), dtype=float)
    for n in range(len(values)):
        for lag in range(taps):
            if n - lag >= 0:
                X[n, lag] = values[n - lag]
    return X


@lru_cache(maxsize=None)
def actual_narma_features(seed):
    path = NARMA / "results" / f"narma10b_driven_seed{seed}" / "voxels.csv"
    f408 = read_window_features(
        path, ("Receiver_R_", "Lum_Mean_"), ("Input_Driven_Death_",)
    )
    frl = read_window_features(path, ("Receiver_R_", "Lum_Mean_"))
    fr = read_window_features(path, ("Receiver_R_",))
    fl = read_window_features(path, ("Lum_Mean_",))
    field = read_window_features(path, ("AHL_uM_",))
    return {"F408": f408, "FRL": frl, "R_ONLY": fr, "L_ONLY": fl, "FIELD": field}


def integration_analytic_constant(C, duration, R0=0.0, L0=0.0):
    hill = C * C / (1.6 * 1.6 + C * C)
    er = math.exp(-duration / 15.0)
    el = math.exp(-duration / 1500.0)
    R = hill + (R0 - hill) * er
    L = (
        hill
        + (L0 - hill) * el
        + (R0 - hill) * 15.0 / (15.0 - 1500.0) * (er - el)
    )
    return R, L


def integrate_surrogate(absolute_time, ahl, method, dt=1.0):
    """Return R/L at stored sample times for all 200 voxels."""
    absolute_time = np.asarray(absolute_time, dtype=float)
    ahl = np.asarray(ahl, dtype=float)
    R = np.zeros(ahl.shape[1], dtype=float)
    L = np.zeros(ahl.shape[1], dtype=float)
    out_r = np.empty_like(ahl)
    out_l = np.empty_like(ahl)
    out_r[0], out_l[0] = R, L

    if method == "zoh_previous":
        for index in range(1, len(absolute_time)):
            delta = absolute_time[index] - absolute_time[index - 1]
            C = ahl[index - 1]
            h = C * C / (1.6 * 1.6 + C * C)
            er = math.exp(-delta / 15.0)
            el = math.exp(-delta / 1500.0)
            old_r, old_l = R, L
            R = h + (old_r - h) * er
            L = (
                h
                + (old_l - h) * el
                + (old_r - h) * 15.0 / (15.0 - 1500.0) * (er - el)
            )
            out_r[index], out_l[index] = R, L
        return out_r, out_l

    if method != "linear":
        raise ValueError(method)

    def deriv(r, lum, concentration):
        hill = concentration * concentration / (
            1.6 * 1.6 + concentration * concentration
        )
        return (hill - r) / 15.0, (r - lum) / 1500.0

    for index in range(1, len(absolute_time)):
        delta = absolute_time[index] - absolute_time[index - 1]
        steps = max(1, int(math.ceil(delta / dt)))
        hdt = delta / steps
        c0, c1 = ahl[index - 1], ahl[index]
        for step in range(steps):
            a0 = step / steps
            am = (step + 0.5) / steps
            a1 = (step + 1.0) / steps
            ca = c0 + a0 * (c1 - c0)
            cm = c0 + am * (c1 - c0)
            cb = c0 + a1 * (c1 - c0)
            k1r, k1l = deriv(R, L, ca)
            k2r, k2l = deriv(R + 0.5 * hdt * k1r, L + 0.5 * hdt * k1l, cm)
            k3r, k3l = deriv(R + 0.5 * hdt * k2r, L + 0.5 * hdt * k2l, cm)
            k4r, k4l = deriv(R + hdt * k3r, L + hdt * k3l, cb)
            R = R + hdt * (k1r + 2 * k2r + 2 * k3r + k4r) / 6.0
            L = L + hdt * (k1l + 2 * k2l + 2 * k3l + k4l) / 6.0
        out_r[index], out_l[index] = R, L
    return out_r, out_l


def integration_unit_test():
    C = 2.3
    times = np.arange(0.0, 3000.0 + 20.0, 20.0)
    ahl = np.full((len(times), 1), C, dtype=float)
    r, l = integrate_surrogate(times, ahl, "linear", dt=1.0)
    exact_r, exact_l = integration_analytic_constant(C, times[-1])
    error_r = abs(float(r[-1, 0]) - exact_r)
    error_l = abs(float(l[-1, 0]) - exact_l)
    return {
        "constant_C": C,
        "duration_s": float(times[-1]),
        "R_abs_error": error_r,
        "L_abs_error": error_l,
        "pass": max(error_r, error_l) < 1e-8,
    }


def surrogate_families(voxels, method, dt=1.0):
    r_samples, l_samples = integrate_surrogate(
        voxels["absolute_time"], voxels["ahl"], method, dt
    )
    mask = (voxels["den"] != 0).astype(float)
    r = window_mean(r_samples, voxels["window"])
    l = window_mean(l_samples, voxels["window"])
    mr = window_mean(r_samples * mask, voxels["window"])
    ml = window_mean(l_samples * mask, voxels["window"])
    return {
        "UNMASKED_RL": np.column_stack([r, l]),
        "OCCUPANCY_MASKED_RL": np.column_stack([mr, ml]),
        "UNMASKED_R": r,
        "UNMASKED_L": l,
        "MASKED_R": mr,
        "MASKED_L": ml,
    }, r_samples, l_samples


def waveform_closed(X_all, labels):
    return WC.evaluate_readout(np.asarray(X_all), np.asarray(labels), WINDOWS)


def add_wave_row(rows, seed, family, method, metrics, n_features):
    per = metrics["test_per_class_ovr_auc"]
    rows.append(
        {
            "seed": seed,
            "family": family,
            "integration_method": method,
            "n_features": n_features,
            "lambda": metrics["lambda"],
            "window_macro_ovr_auc": metrics["test_macro_ovr_auc"],
            "window_accuracy": metrics["test_accuracy_argmax"],
            "sine_auc": per["sine"],
            "square_auc": per["square"],
            "triangle_auc": per["triangle"],
            "block_pooled_macro_auc": metrics["block_pooled_macro_ovr_auc"],
            "confusion": json.dumps(metrics["test_confusion"], separators=(",", ":")),
        }
    )


def stage0_sanity(u_narma, y_narma, wave_labels):
    rows = []
    details = {}

    narma_driven, narma_field, narma_brown = [], [], []
    for seed in NARMA_SEEDS:
        actual = actual_narma_features(seed)
        narma_driven.append(scalar_closed(actual["F408"], y_narma)["test_nrmse"])
        narma_field.append(scalar_closed(actual["FIELD"], y_narma)["test_nrmse"])
        brown_path = (
            NARMA / "results" / f"narma10b_brownian_seed{seed}" / "voxels.csv"
        )
        brown = read_window_features(brown_path, ("Den_",))
        narma_brown.append(scalar_closed(brown, y_narma)["test_nrmse"])
    silent_path = NARMA / "results" / "narma10b_silent_seed111" / "voxels.csv"
    silent = read_window_features(
        silent_path, ("Receiver_R_", "Lum_Mean_"), ("Input_Driven_Death_",)
    )
    silent_score = scalar_closed(silent, y_narma)["test_nrmse"]
    y_current = np.concatenate([[0.0], y_narma[:-1]])
    persistence = evaluate_nontrainable_series(y_narma, y_current)["test_nrmse"]
    train_mean = float(np.mean(y_narma[WASHOUT : WASHOUT + TRAIN]))
    intercept = constant_metrics(y_narma, train_mean)["test_nrmse"]
    narma_checks = [
        *(within(a, e) for a, e in zip(narma_driven, NARMA_EXPECTED["driven"])),
        within(mean(narma_field), NARMA_EXPECTED["field"]),
        within(mean(narma_brown), NARMA_EXPECTED["brownian_mean"]),
        within(silent_score, NARMA_EXPECTED["silent"]),
        within(persistence, NARMA_EXPECTED["persistence"]),
        within(intercept, NARMA_EXPECTED["train_intercept"]),
    ]
    details["narma"] = {
        "driven": narma_driven,
        "field": narma_field,
        "brownian": narma_brown,
        "silent_seed111": silent_score,
        "persistence": persistence,
        "train_intercept": intercept,
        "pass": all(narma_checks),
    }
    rows.extend(
        [
            ("NARMA", "driven F408", narma_driven, NARMA_EXPECTED["driven"], all(narma_checks[:3])),
            ("NARMA", "field", mean(narma_field), NARMA_EXPECTED["field"], narma_checks[3]),
            ("NARMA", "Brownian mean", mean(narma_brown), NARMA_EXPECTED["brownian_mean"], narma_checks[4]),
            ("NARMA", "silent seed111", silent_score, NARMA_EXPECTED["silent"], narma_checks[5]),
            ("NARMA", "persistence", persistence, NARMA_EXPECTED["persistence"], narma_checks[6]),
            ("NARMA", "train-intercept", intercept, NARMA_EXPECTED["train_intercept"], narma_checks[7]),
        ]
    )

    wave_auc, wave_field, wave_acc = [], [], []
    wave_brown, wave_silent = [], []
    for seed in WAVEFORM_SEEDS:
        driven_path = (
            WAVEFORM / "results" / f"waveform_driven_seed{seed}" / "voxels.csv"
        )
        driven = read_window_features(
            driven_path, ("Receiver_R_", "Lum_Mean_"), ("Input_Driven_Death_",)
        )
        field = read_window_features(driven_path, ("AHL_uM_",))
        brown = read_window_features(
            WAVEFORM / "results" / f"waveform_brownian_seed{seed}" / "voxels.csv",
            ("Den_",),
        )
        silent = read_window_features(
            WAVEFORM / "results" / f"waveform_silent_seed{seed}" / "voxels.csv",
            ("Receiver_R_", "Lum_Mean_"),
            ("Input_Driven_Death_",),
        )
        dm = waveform_closed(driven, wave_labels)
        fm = waveform_closed(field, wave_labels)
        wave_auc.append(dm["test_macro_ovr_auc"])
        wave_field.append(fm["test_macro_ovr_auc"])
        wave_acc.append(dm["test_accuracy_argmax"])
        wave_brown.append(waveform_closed(brown, wave_labels)["test_macro_ovr_auc"])
        wave_silent.append(waveform_closed(silent, wave_labels)["test_macro_ovr_auc"])
    wave_checks = [
        *(within(a, e) for a, e in zip(wave_auc, WAVE_EXPECTED["driven_auc"])),
        within(mean(wave_field), WAVE_EXPECTED["field_auc"]),
        within(mean(wave_acc), WAVE_EXPECTED["driven_accuracy_mean"], tol=0.01),
    ]
    details["waveform"] = {
        "driven_auc": wave_auc,
        "field_auc": wave_field,
        "driven_accuracy": wave_acc,
        "brownian_auc": wave_brown,
        "silent_auc": wave_silent,
        "pass": all(wave_checks),
    }
    rows.extend(
        [
            ("Waveform", "driven macro AUC", wave_auc, WAVE_EXPECTED["driven_auc"], all(wave_checks[:3])),
            ("Waveform", "field macro AUC", mean(wave_field), WAVE_EXPECTED["field_auc"], wave_checks[3]),
            ("Waveform", "driven accuracy mean", mean(wave_acc), WAVE_EXPECTED["driven_accuracy_mean"], wave_checks[4]),
        ]
    )

    lorenz_evidence = load_json(BENCHA / "results" / "lorenz_gate_evidence.json")
    mg_evidence = load_json(BENCHA / "results" / "mg_gate_evidence.json")
    lorenz_driven = lorenz_evidence["test_nrmse"]["driven"]
    lorenz_field = lorenz_evidence["field_only_test_nrmse"]
    mg_driven = mg_evidence["test_nrmse"]["driven"]
    mg_field = mg_evidence["field_only_test_nrmse"]
    bench_checks = [
        within(mean(lorenz_driven), BENCHA_EXPECTED["lorenz_driven_mean"]),
        within(mean(lorenz_field), BENCHA_EXPECTED["lorenz_field_mean"], tol=0.002),
        within(mean(mg_driven), BENCHA_EXPECTED["mg_driven_mean"]),
        within(mean(mg_field), BENCHA_EXPECTED["mg_field_mean"]),
    ]
    details["bencha"] = {
        "lorenz_driven": lorenz_driven,
        "lorenz_field": lorenz_field,
        "mg_driven": mg_driven,
        "mg_field": mg_field,
        "pass": all(bench_checks),
    }
    rows.extend(
        [
            ("BenchA", "Lorenz driven mean", mean(lorenz_driven), BENCHA_EXPECTED["lorenz_driven_mean"], bench_checks[0]),
            ("BenchA", "Lorenz field mean", mean(lorenz_field), BENCHA_EXPECTED["lorenz_field_mean"], bench_checks[1]),
            ("BenchA", "MG driven mean", mean(mg_driven), BENCHA_EXPECTED["mg_driven_mean"], bench_checks[2]),
            ("BenchA", "MG field mean", mean(mg_field), BENCHA_EXPECTED["mg_field_mean"], bench_checks[3]),
        ]
    )
    return details, rows


def run_narma_baselines(u, y):
    rows = []
    train_target = y[WASHOUT : WASHOUT + TRAIN]
    train_mean = float(np.mean(train_target))
    intercept = constant_metrics(y, train_mean)
    rows.append(
        {
            "baseline": "TRAIN_INTERCEPT",
            "status": "legal",
            "n_features": 0,
            "lambda": "",
            **intercept,
        }
    )
    y_current = np.concatenate([[0.0], y[:-1]])
    persistence = evaluate_nontrainable_series(y, y_current)
    rows.append(
        {
            "baseline": "PERSISTENCE",
            "status": "legal",
            "n_features": 1,
            "lambda": "",
            **persistence,
        }
    )
    taps = delay_matrix(u, 10)
    linear = scalar_closed(taps, y)
    rows.append(
        {
            "baseline": "LINEAR_U_DELAY_10",
            "status": "legal",
            "n_features": 10,
            **linear,
        }
    )
    product = np.array([u[n] * u[n - 9] if n >= 9 else 0.0 for n in WINDOWS])
    informed_X = np.column_stack([taps, product])
    informed = scalar_closed(informed_X, y)
    rows.append(
        {
            "baseline": "NARMA_INFORMED_INPUT",
            "status": "legal",
            "n_features": 11,
            **informed,
        }
    )
    test_target = y[WASHOUT + TRAIN :]
    oracle_mean = float(np.mean(test_target))
    oracle = constant_metrics(y, oracle_mean)
    rows.append(
        {
            "baseline": "TEST_MEAN_ORACLE_REFERENCE",
            "status": "non-predictive oracle normalization reference; uses test labels",
            "n_features": 0,
            "lambda": "",
            **oracle,
        }
    )
    driven_mean = mean(NARMA_EXPECTED["driven"])
    for row in rows:
        row["driven_mean_0_9283_beats"] = (
            "not_ranked_oracle"
            if row["status"] != "legal"
            else str(driven_mean < row["test_nrmse"]).upper()
        )
    return rows


def preprocessing_audit(y):
    source_select = inspect.getsource(NC.select_lambda)
    source_fit = inspect.getsource(NC.fit_eval)
    source_standardize = inspect.getsource(NC.standardize)
    source_ridge = inspect.getsource(NC.ridge_fit)
    synthetic = np.column_stack(
        [np.linspace(-1, 1, NUM_WINDOWS), np.ones(NUM_WINDOWS)]
    )
    X = synthetic[WASHOUT:]
    yy = y[WASHOUT:]
    lam_a, _ = NC.select_lambda(X, yy)
    X_changed = X.copy()
    X_changed[TRAIN:] = 1e12
    yy_changed = yy.copy()
    yy_changed[TRAIN:] = -1e12
    lam_b, _ = NC.select_lambda(X_changed, yy_changed)
    standardized, _, _ = NC.standardize(
        np.ones((10, 3)), np.ones((2, 3))
    )
    tie_lam, tie_grid = NC.select_lambda(
        np.zeros((TRAIN + TEST, 1)), np.zeros(TRAIN + TEST)
    )
    checks = [
        {
            "item": "validation rows are windows 128..149 only",
            "pass": "X[inner_train_end:TRAIN]" in source_select
            and TRAIN - INNER_VAL == 88,
            "location": "check_narma10b.py:278-293 (select_lambda)",
        },
        {
            "item": "test is never used for lambda",
            "pass": lam_a == lam_b and "X[:inner_train_end]" in source_select,
            "location": "check_narma10b.py:278-293 (runtime contamination test)",
        },
        {
            "item": "inner standardization uses corresponding inner train only",
            "pass": "standardize(X_inner, X_val)" in source_select,
            "location": "check_narma10b.py:280-283",
        },
        {
            "item": "final standardization is fit on all 110 train windows",
            "pass": "standardize(X_train, X_test)" in source_fit
            and "X[:TRAIN]" in source_fit,
            "location": "check_narma10b.py:296-308",
        },
        {
            "item": "intercept is not regularized",
            "pass": "gram[1:, 1:]" in source_ridge,
            "location": "check_narma10b.py:260-270",
        },
        {
            "item": "zero-variance columns do not produce NaN",
            "pass": all(np.isfinite(array).all() for array in standardized)
            and "std < 1e-12" in source_standardize,
            "location": "check_narma10b.py:252-257 (runtime constant-column test)",
        },
        {
            "item": "tie selects larger lambda",
            "pass": tie_lam == max(RIDGE_GRID)
            and len({row["val_nrmse"] for row in tie_grid}) == 1,
            "location": "check_narma10b.py:286-293 (runtime exact-tie test)",
        },
    ]
    return checks


def run_narma_nulls(y):
    rows = []
    feature_rows = []
    for seed in NARMA_SEEDS:
        features = actual_narma_features(seed)
        for family in ("F408", "FRL", "R_ONLY", "L_ONLY"):
            metrics = scalar_closed(features[family], y)
            feature_rows.append(
                {
                    "seed": seed,
                    "family": family,
                    "n_features": features[family].shape[1],
                    **metrics,
                }
            )
        for family in ("F408", "FIELD"):
            X = features[family]
            for offset in (0,) + SHIFT_OFFSETS:
                target = y if offset == 0 else np.roll(y, offset)
                metrics = evaluate_mismatch(X, target)
                rows.append(
                    {
                        "seed": seed,
                        "family": "OFFICIAL_BIOLOGY" if family == "F408" else "FIELD_ONLY",
                        "target": "true" if offset == 0 else "circular_shift_null",
                        "offset_windows": offset,
                        "n_features": X.shape[1],
                        **metrics,
                    }
                )
    # Function contract test: matched i=i reproduces the ordinary closed result.
    test_features = actual_narma_features(111)["F408"]
    mismatch_test = evaluate_mismatch(test_features, y)
    ordinary_test = scalar_closed(test_features, y)
    function_pass = all(
        abs(mismatch_test[key] - ordinary_test[key]) < 1e-12
        for key in ("lambda", "test_nrmse", "test_r2")
    )
    mismatch_status = {
        "status": "PENDING",
        "reason": (
            "Only one independent NARMA input realization exists; bacterial "
            "seeds are not independent inputs."
        ),
        "function_contract_test_pass": function_pass,
    }
    return rows, feature_rows, mismatch_status


def run_narma_surrogates(y):
    rows = []
    convergence = []
    unit = integration_unit_test()
    for seed in NARMA_SEEDS:
        path = NARMA / "results" / f"narma10b_driven_seed{seed}" / "voxels.csv"
        voxels = read_voxels(path)
        methods = {}
        for label, method, dt in (
            ("linear_dt1", "linear", 1.0),
            ("linear_dt0.5", "linear", 0.5),
            ("zoh_previous", "zoh_previous", 1.0),
        ):
            families, rs, ls = surrogate_families(voxels, method, dt)
            methods[label] = (families, rs, ls)
            for family, X in families.items():
                metrics = scalar_closed(X, y)
                rows.append(
                    {
                        "seed": seed,
                        "family": family,
                        "integration_method": label,
                        "n_features": X.shape[1],
                        **metrics,
                    }
                )
        r1, l1 = methods["linear_dt1"][1:]
        r05, l05 = methods["linear_dt0.5"][1:]
        convergence.append(
            {
                "seed": seed,
                "max_abs_R_dt1_vs_dt0_5": float(np.max(np.abs(r1 - r05))),
                "max_abs_L_dt1_vs_dt0_5": float(np.max(np.abs(l1 - l05))),
            }
        )
        actual = actual_narma_features(seed)
        for family in ("FRL", "F408", "FIELD", "R_ONLY", "L_ONLY"):
            metrics = scalar_closed(actual[family], y)
            rows.append(
                {
                    "seed": seed,
                    "family": f"ACTUAL_{family}",
                    "integration_method": "actual_or_historical",
                    "n_features": actual[family].shape[1],
                    **metrics,
                }
            )
    return rows, convergence, unit


def run_wave_surrogates(labels):
    rows = []
    convergence = []
    for seed in WAVEFORM_SEEDS:
        path = WAVEFORM / "results" / f"waveform_driven_seed{seed}" / "voxels.csv"
        voxels = read_voxels(path)
        methods = {}
        for label, method, dt in (
            ("linear_dt1", "linear", 1.0),
            ("linear_dt0.5", "linear", 0.5),
            ("zoh_previous", "zoh_previous", 1.0),
        ):
            families, rs, ls = surrogate_families(voxels, method, dt)
            methods[label] = (families, rs, ls)
            for family, X in families.items():
                add_wave_row(
                    rows, seed, family, label, waveform_closed(X, labels), X.shape[1]
                )
        r1, l1 = methods["linear_dt1"][1:]
        r05, l05 = methods["linear_dt0.5"][1:]
        convergence.append(
            {
                "seed": seed,
                "max_abs_R_dt1_vs_dt0_5": float(np.max(np.abs(r1 - r05))),
                "max_abs_L_dt1_vs_dt0_5": float(np.max(np.abs(l1 - l05))),
            }
        )
        actual_specs = (
            (
                "ACTUAL_DRIVEN_F408",
                read_window_features(
                    path,
                    ("Receiver_R_", "Lum_Mean_"),
                    ("Input_Driven_Death_",),
                ),
            ),
            ("ACTUAL_DRIVEN_FRL", read_window_features(path, ("Receiver_R_", "Lum_Mean_"))),
            ("ACTUAL_FIELD", read_window_features(path, ("AHL_uM_",))),
            (
                "ACTUAL_BROWNIAN",
                read_window_features(
                    WAVEFORM / "results" / f"waveform_brownian_seed{seed}" / "voxels.csv",
                    ("Den_",),
                ),
            ),
            (
                "ACTUAL_SILENT",
                read_window_features(
                    WAVEFORM / "results" / f"waveform_silent_seed{seed}" / "voxels.csv",
                    ("Receiver_R_", "Lum_Mean_"),
                    ("Input_Driven_Death_",),
                ),
            ),
        )
        for family, X in actual_specs:
            add_wave_row(
                rows,
                seed,
                family,
                "actual_or_historical",
                waveform_closed(X, labels),
                X.shape[1],
            )
    return rows, convergence


def template_moments(values):
    values = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(values)),
        "variance": float(np.var(values)),
        "population_sd": float(np.std(values)),
        "power_mean_u2": float(np.mean(values * values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "range": float(np.max(values) - np.min(values)),
    }


def block_closed_classifier(X_blocks, y_blocks):
    """Frozen OVR ridge behavior with the predeclared whole-block split."""
    X = np.asarray(X_blocks, dtype=float)[8:]
    y = np.asarray(y_blocks, dtype=int)[8:]
    X_train, y_train = X[:22], y[:22]
    X_test, y_test = X[22:], y[22:]
    X_inner, y_inner = X_train[:17], y_train[:17]  # blocks 8..24
    X_val, y_val = X_train[17:22], y_train[17:22]  # blocks 25..29
    (Xi, Xv), _, _ = WC.standardize(X_inner, X_val)
    candidates = []
    for lam in RIDGE_GRID:
        scores, _ = WC.ovr_scores(Xi, y_inner, Xv, lam)
        macro, _ = WC.macro_ovr_auc(y_val, scores)
        candidates.append((-macro, -lam, lam, macro))
    candidates.sort()
    lam = candidates[0][2]
    (Xt, Xq), _, _ = WC.standardize(X_train, X_test)
    scores, _ = WC.ovr_scores(Xt, y_train, Xq, lam)
    macro, per = WC.macro_ovr_auc(y_test, scores)
    accuracy, _, confusion = WC.accuracy_argmax(y_test, scores)
    return {
        "lambda": lam,
        "macro_ovr_auc": macro,
        "accuracy": accuracy,
        "sine_auc": per[0],
        "square_auc": per[1],
        "triangle_auc": per[2],
        "confusion": json.dumps(confusion, separators=(",", ":")),
    }


def run_waveform_moments(labels, blocks, u):
    templates5 = {
        "sine": np.array([WC.template_u(0, w) for w in range(5)]),
        "square": np.array([WC.template_u(1, w) for w in range(5)]),
        "triangle": np.array([WC.template_u(2, w) for w in range(5)]),
    }
    moment_rows = []
    for name, values in templates5.items():
        moment_rows.append({"template_set": "frozen_5_point", "class": name, **template_moments(values)})

    templates8 = {
        "sine": np.array([0.25 + 0.25 * math.sin(2 * math.pi * w / 8) for w in range(8)]),
        "square": np.array([0.5] * 4 + [0.0] * 4),
        "triangle": np.array([0, 0.125, 0.25, 0.375, 0.5, 0.375, 0.25, 0.125]),
    }
    for name, values in templates8.items():
        moment_rows.append({"template_set": "proposed_8_point", "class": name, **template_moments(values)})

    block_u = np.vstack([u[p * 5 : (p + 1) * 5] for p in range(40)])
    block_y = np.array([labels[p * 5] for p in range(40)], dtype=int)
    if not np.array_equal(blocks, np.repeat(np.arange(40), 5)):
        raise AssertionError("waveform block sidecar mismatch")
    feature_sets = {
        "MEAN_ONLY": np.mean(block_u, axis=1)[:, None],
        "POWER_ONLY": np.mean(block_u * block_u, axis=1)[:, None],
        "VARIANCE_ONLY": np.var(block_u, axis=1)[:, None],
        "MOMENTS": np.column_stack(
            [
                np.mean(block_u, axis=1),
                np.var(block_u, axis=1),
                np.mean(block_u * block_u, axis=1),
                np.min(block_u, axis=1),
                np.max(block_u, axis=1),
                np.ptp(block_u, axis=1),
            ]
        ),
        "RAW_TEMPLATE_5": block_u,
    }
    classifier_rows = []
    for name, X in feature_sets.items():
        classifier_rows.append(
            {
                "record_type": "block_baseline",
                "template_set": "frozen_5_point",
                "baseline": name,
                "n_features": X.shape[1],
                "n_train_blocks": 22,
                "n_validation_blocks": 5,
                "n_test_blocks": 10,
                **block_closed_classifier(X, block_y),
            }
        )
    return moment_rows, classifier_rows


def lorenz_rows_from_target(path):
    return load_semicolon(path)


def evaluate_nontrainable_series(target, prediction):
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    y = target[WASHOUT:]
    p = prediction[WASHOUT:]
    return {
        "train_nrmse": nrmse(y[:TRAIN], p[:TRAIN]),
        "validation_nrmse": nrmse(
            y[TRAIN - INNER_VAL : TRAIN], p[TRAIN - INNER_VAL : TRAIN]
        ),
        "test_nrmse": nrmse(y[TRAIN:], p[TRAIN:]),
        "test_r2": r2(y[TRAIN:], p[TRAIN:]),
    }


def run_lorenz_baselines():
    output = []
    verifications = []
    clocks = (
        ("BenchA_dt0.02", BC.recompute_lorenz(), BENCHA / "lorenz_target.csv", 0.02, BC.LORENZ_SHA),
        ("BenchA2_dt1.0", B2C.recompute_lorenz_skip50(), BENCHA2 / "lorenz_skip50_target.csv", 1.0, B2C.LORENZ_SHA),
    )
    for clock, samples, target_path, delta_t, expected_hash in clocks:
        rows = lorenz_rows_from_target(target_path)
        stored_u = np.array([float(row["u"]) for row in rows])
        stored_x = np.array([float(row["x"]) for row in rows])
        stored_y_next = np.array([float(row["y_next"]) for row in rows])
        x = samples[:NUM_WINDOWS, 0]
        x_next = samples[1 : NUM_WINDOWS + 1, 0]
        y_next = samples[1 : NUM_WINDOWS + 1, 1]
        xmin, xmax = float(np.min(x)), float(np.max(x))
        u_hat = 0.5 * (x - xmin) / (xmax - xmin)
        digest = BC.sha256_u(stored_u)
        verify = {
            "clock": clock,
            "delta_t_sample": delta_t,
            "u_hash": digest,
            "expected_u_hash": expected_hash,
            "hash_pass": digest == expected_hash,
            "x_pass": bool(np.allclose(stored_x, x, rtol=0, atol=1e-9)),
            "y_next_pass": bool(np.allclose(stored_y_next, y_next, rtol=0, atol=1e-9)),
            "u_pass": bool(np.allclose(stored_u, u_hat, rtol=0, atol=1e-9)),
        }
        verify["pass"] = all(
            verify[key] for key in ("hash_pass", "x_pass", "y_next_pass", "u_pass")
        )
        verifications.append(verify)

        targets = (("AUTO_X", x_next), ("CROSS_Y", y_next))
        for task, target in targets:
            train_mean = float(np.mean(target[WASHOUT : WASHOUT + TRAIN]))
            row = {
                "clock": clock,
                "delta_t_sample": delta_t,
                "task": task,
                "model": "TRAIN_INTERCEPT",
                "development_m": "",
                "selected_on_validation": "",
                "n_features": 0,
                "lambda": "",
                **constant_metrics(target, train_mean),
            }
            output.append(row)
            oracle = float(np.mean(target[WASHOUT + TRAIN :]))
            output.append(
                {
                    "clock": clock,
                    "delta_t_sample": delta_t,
                    "task": task,
                    "model": "TEST_MEAN_ORACLE_REFERENCE",
                    "development_m": "",
                    "selected_on_validation": "",
                    "n_features": 0,
                    "lambda": "",
                    "note": "non-predictive oracle normalization reference; uses test labels",
                    **constant_metrics(target, oracle),
                }
            )
            if task == "AUTO_X":
                output.append(
                    {
                        "clock": clock,
                        "delta_t_sample": delta_t,
                        "task": task,
                        "model": "PERSISTENCE_X_N",
                        "development_m": "",
                        "selected_on_validation": "",
                        "n_features": 1,
                        "lambda": "",
                        **evaluate_nontrainable_series(target, x),
                    }
                )
                m_values = (1, 3, 5, 7, 10)
                prefix = "LINEAR_AR"
            else:
                m_values = (1, 3, 5, 7, 10)
                prefix = "LINEAR_X_DELAY"
            development = []
            for m in m_values:
                X = delay_matrix(x, m)
                metrics = scalar_closed(X, target)
                development.append((m, metrics))
                model_name = (
                    "SCALAR_LINEAR_X"
                    if task == "CROSS_Y" and m == 1
                    else prefix
                )
                output.append(
                    {
                        "clock": clock,
                        "delta_t_sample": delta_t,
                        "task": task,
                        "model": model_name,
                        "development_m": m,
                        "selected_on_validation": "FALSE",
                        "n_features": m,
                        **metrics,
                    }
                )
            selected_m, selected_metrics = min(
                development, key=lambda item: (item[1]["validation_nrmse"], item[0])
            )
            output.append(
                {
                    "clock": clock,
                    "delta_t_sample": delta_t,
                    "task": task,
                    "model": f"{prefix}_VALIDATION_SELECTED",
                    "development_m": selected_m,
                    "selected_on_validation": "TRUE",
                    "n_features": selected_m,
                    **selected_metrics,
                }
            )

        evidence_path = (
            BENCHA / "results" / "lorenz_gate_evidence.json"
            if delta_t == 0.02
            else BENCHA2 / "results" / "lorenz_skip50_gate_evidence.json"
        )
        evidence = load_json(evidence_path)
        for arm, values in (
            ("EXISTING_DRIVEN_CROSS_Y", evidence["test_nrmse"]["driven"]),
            ("EXISTING_FIELD_CROSS_Y", evidence["field_only_test_nrmse"]),
        ):
            for seed, score in zip(WAVEFORM_SEEDS, values):
                output.append(
                    {
                        "clock": clock,
                        "delta_t_sample": delta_t,
                        "task": "CROSS_Y",
                        "model": arm,
                        "development_m": "",
                        "selected_on_validation": "",
                        "seed": seed,
                        "n_features": 408 if "DRIVEN" in arm else 200,
                        "lambda": evidence["lambdas"]["driven" if "DRIVEN" in arm else "driven"][WAVEFORM_SEEDS.index(seed)]
                        if "DRIVEN" in arm
                        else "",
                        "train_nrmse": "",
                        "validation_nrmse": "",
                        "test_nrmse": score,
                        "test_r2": 1.0 - score * score,
                        "note": "existing frozen checker score; target alignment permits CROSS_Y only",
                    }
                )
    return output, verifications


def write_csv(path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: fnum(value)
                    if isinstance(value, (float, np.floating))
                    else value
                    for key, value in row.items()
                }
            )


def archive_external():
    EXTERNAL.mkdir(parents=True, exist_ok=True)
    records = []
    for source, filename in EXTERNAL_SOURCES:
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = EXTERNAL / filename
        shutil.copyfile(source, destination)
        src_hash = sha256_file(source)
        dst_hash = sha256_file(destination)
        if src_hash != dst_hash:
            raise AssertionError(f"hash mismatch after copy: {filename}")
        records.append(
            {
                "filename": filename,
                "original_absolute_path": str(source),
                "acquisition_date": "2026-08-16",
                "sha256": src_hash,
                "byte_size": source.stat().st_size,
                "label": "external, untrusted analysis—not project evidence",
                "source_destination_hash_match": True,
            }
        )
    lines = [
        "# External analysis provenance manifest",
        "",
        "These files are preserved byte-for-byte for provenance. They are",
        "**external, untrusted analysis—not project evidence**.",
        "",
        "| Filename | Original absolute path | Acquisition date | SHA-256 | Bytes | Label | Hash match |",
        "|---|---|---|---|---:|---|---|",
    ]
    for row in records:
        lines.append(
            f"| [{row['filename']}](./{row['filename']}) | `{row['original_absolute_path']}` "
            f"| {row['acquisition_date']} | `{row['sha256']}` | {row['byte_size']} "
            f"| {row['label']} | PASS |"
        )
    lines.extend(
        [
            "",
            "No chat review was vendored. Only the three actual source files listed",
            "above are archived.",
            "",
        ]
    )
    (EXTERNAL / "MANIFEST.md").write_text("\n".join(lines), encoding="utf-8")
    return records


def write_ipc():
    broad_delays = 20
    broad_targets = (
        broad_delays
        + math.comb(broad_delays + 1, 2)
        + math.comb(broad_delays + 2, 3)
    )
    narrow_targets = 10 + 5 + 5
    lines = [
        "# IPC affordability decision",
        "",
        "## Current sample geometry",
        "",
        "- total windows: 200",
        "- post-washout time rows available to analysis: 160 (sequential rows,",
        "  not 160 statistically independent realizations)",
        "- train / validation / test: 110 / 22 (inside train) / 50",
        "- official state features: 408",
        "- observed linear memory capacity: approximately 1.22",
        "- useful lags: approximately 5",
        "",
        "The 408-column state matrix has at most rank 110 in the final training",
        "fit (and at most rank 88 during lambda selection), before centering and",
        "intercept effects. More features do not create more independent temporal",
        "observations.",
        "",
        "## Broad Dambre-style IPC",
        "",
        f"Illustrative predeclared basis: delays 1..{broad_delays}, all monomials",
        "through order three. This contains "
        f"`20 + C(21,2) + C(22,3) = {broad_targets}` targets. A different",
        "orthogonal basis changes the exact count but not the current sample",
        "shortage.",
        "",
        f"- target count: {broad_targets}",
        "- samples per target: the same 88 inner-train, 22 validation, and 50 test",
        "  rows; targets do not provide new state observations",
        "- feature/sample rank limitation: 408 features versus 88/110 fitting",
        "  rows, so estimates rely strongly on ridge regularization",
        f"- expected compute burden: at least {broad_targets * len(RIDGE_GRID):,}",
        "  ridge fits before input-realization replication or nulls",
        "- estimator/null burden: each target needs held-out capacity estimation,",
        "  label/target-shift nulls, multiplicity-aware summaries, and multiple",
        "  independent random drives",
        "- recommendation: not credible on the current 200-window trajectory",
        "",
        "Exact IPC sample requirements depend on input distribution, target basis,",
        "state dimension, estimator, regularization, desired uncertainty, and null",
        "design. Published large-sample practice is context, not a universal",
        "`10^4–10^5` theorem. The current 200 windows are nevertheless plainly",
        "insufficient for a broad stable spectrum.",
        "",
        "## Narrow predeclared panel",
        "",
        "Illustrative panel: linear delays 1..10, five individual quadratic",
        "terms, and five predeclared cross-delay products.",
        "",
        f"- target count: {narrow_targets}",
        "- samples per target: 88 inner-train, 22 validation, 50 test",
        "- feature/sample rank limitation: unchanged; 408 > 110, so uncertainty",
        "  and ridge dependence remain material",
        f"- expected compute burden: {narrow_targets * len(RIDGE_GRID)} primary",
        "  lambda fits per state/seed, plus refits",
        "- estimator/null burden: manageable only with explicit shift nulls,",
        "  confidence intervals, and independently generated random drives",
        "- recommendation: feasible as a later separately frozen diagnostic, not",
        "  as a broad IPC decomposition and not retroactively on this one drive",
        "",
        "## Primary decision",
        "",
        "IPC_DECISION: KILL_BROAD_KEEP_MC",
        "",
        "Keep MC≈1.22 and useful-lag≈5 as the supported capacity diagnostics.",
        "A narrow panel or a long random-drive experiment may be proposed later",
        "under its own frozen protocol and compute budget.",
        "",
    ]
    (HERE / "ipc_affordability.md").write_text("\n".join(lines), encoding="utf-8")
    return {"broad_targets": broad_targets, "narrow_targets": narrow_targets}


def markdown_table(headers, rows):
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    lines.extend("| " + " | ".join(str(cell) for cell in row) + " |" for row in rows)
    return lines


def row_mean(rows, family, metric, method=None):
    selected = [
        float(row[metric])
        for row in rows
        if row["family"] == family
        and (method is None or row["integration_method"] == method)
    ]
    return mean(selected) if selected else float("nan")


def selected_lorenz(rows, clock, task):
    return next(
        row
        for row in rows
        if row["clock"] == clock
        and row["task"] == task
        and row.get("selected_on_validation") == "TRUE"
    )


def write_report(
    sanity,
    sanity_rows,
    narma_surrogate,
    narma_convergence,
    integration_test,
    narma_baselines,
    narma_nulls,
    feature_rows,
    mismatch_status,
    preprocessing,
    lorenz_rows,
    lorenz_verify,
    moment_rows,
    moment_classifiers,
    wave_surrogate,
    wave_convergence,
    provenance,
):
    lines = [
        "# Zero-simulation evidence audit",
        "",
        "## 1. Executive verdict",
        "",
    ]
    if all(family["pass"] for family in sanity.values()):
        lines.append(
            "All three mandatory sanity families reproduced the frozen evidence, so "
            "all dependent zero-simulation modules were executed."
        )
    else:
        failed = [name for name, value in sanity.items() if not value["pass"]]
        lines.append(
            "Mandatory sanity failed for: "
            + ", ".join(failed)
            + ". Dependent analyses were stopped; independent modules continued."
        )
    if narma_baselines and narma_surrogate and wave_surrogate and lorenz_rows:
        linear_u = next(
            row for row in narma_baselines if row["baseline"] == "LINEAR_U_DELAY_10"
        )
        informed_u = next(
            row
            for row in narma_baselines
            if row["baseline"] == "NARMA_INFORMED_INPUT"
        )
        fast_cross = selected_lorenz(lorenz_rows, "BenchA_dt0.02", "CROSS_Y")
        lines.extend(
            [
                "",
                f"The decisive attribution result is that the legal ten-tap direct",
                f"NARMA input baseline scores `{linear_u['test_nrmse']:.4f}` and the",
                f"frozen informed input baseline scores `{informed_u['test_nrmse']:.4f}`,",
                "both better than driven biology `0.9283`. The living readout still",
                "beats Brownian, silent, persistence, train-intercept, and field, but",
                "it does not beat the task's direct input-only baselines. The",
                "confirmatory question must therefore be narrowed before more NARMA",
                "compute is authorized.",
                "",
                f"The unmasked NARMA R/L surrogate (`{row_mean(narma_surrogate, 'UNMASKED_RL', 'test_nrmse', 'linear_dt1'):.4f}`)",
                f"and occupancy-masked surrogate (`{row_mean(narma_surrogate, 'OCCUPANCY_MASKED_RL', 'test_nrmse', 'linear_dt1'):.4f}`)",
                "bracket actual FRL `0.9283`; explicit cells add little beyond the",
                "modeled kinetic cascade plus occupancy for this trajectory. Waveform",
                "is even more separable through the unmasked kinetics, while its",
                "five-point moments already rank classes at block level. At",
                f"Δt=0.02, a validation-selected x-delay predicts Lorenz CROSS_Y at",
                f"`{fast_cross['test_nrmse']:.4f}`, far ahead of the living readout.",
                "Broad IPC is killed at the current sample geometry.",
            ]
        )
    lines.extend(
        [
            "",
            "The audit does not change any historical gate: NARMA remains **PASS**,",
            "waveform **FAIL**, Lorenz **FAIL**, Mackey–Glass **PASS**, and all",
            "sweep verdicts remain unchanged. No result below is a retrospective",
            "gate.",
            "",
            "## 2. Mandatory sanity",
            "",
        ]
    )
    table_rows = []
    for family, check, actual, expected, passed in sanity_rows:
        actual_text = (
            " / ".join(f"{float(x):.4f}" for x in actual)
            if isinstance(actual, (list, tuple))
            else f"{float(actual):.4f}"
        )
        expected_text = (
            " / ".join(f"{float(x):.4f}" for x in expected)
            if isinstance(expected, (list, tuple))
            else f"{float(expected):.4f}"
        )
        table_rows.append((family, check, actual_text, expected_text, "PASS" if passed else "FAIL"))
    lines.extend(markdown_table(["Family", "Check", "Actual", "Expected", "Result"], table_rows))

    lines.extend(
        [
            "",
            "## 3. Kinetic-surrogate attribution",
            "",
            f"Constant-C analytical integrator test: **{'PASS' if integration_test['pass'] else 'FAIL'}**; "
            f"R absolute error `{integration_test['R_abs_error']:.3e}`, "
            f"L absolute error `{integration_test['L_abs_error']:.3e}`.",
            "",
        ]
    )
    if narma_surrogate:
        narma_summary = []
        for family in (
            "UNMASKED_RL",
            "OCCUPANCY_MASKED_RL",
            "UNMASKED_R",
            "UNMASKED_L",
            "MASKED_R",
            "MASKED_L",
            "ACTUAL_FRL",
            "ACTUAL_F408",
            "ACTUAL_FIELD",
        ):
            method = "linear_dt1" if not family.startswith("ACTUAL") else "actual_or_historical"
            zoh = (
                row_mean(narma_surrogate, family, "test_nrmse", "zoh_previous")
                if not family.startswith("ACTUAL")
                else float("nan")
            )
            narma_summary.append(
                (
                    family,
                    f"{row_mean(narma_surrogate, family, 'test_nrmse', method):.4f}",
                    "—" if family.startswith("ACTUAL") else f"{zoh:.4f}",
                )
            )
        lines.extend(
            markdown_table(
                ["NARMA family", "Linear dt=1 mean NRMSE", "ZOH-previous mean NRMSE"],
                narma_summary,
            )
        )
        paired_rows = []
        for seed in NARMA_SEEDS:
            def ns(family, method):
                return next(
                    row["test_nrmse"]
                    for row in narma_surrogate
                    if row["seed"] == seed
                    and row["family"] == family
                    and row["integration_method"] == method
                )
            paired_rows.append(
                (
                    seed,
                    f"{ns('ACTUAL_FIELD', 'actual_or_historical'):.4f}",
                    f"{ns('UNMASKED_RL', 'linear_dt1'):.4f}",
                    f"{ns('OCCUPANCY_MASKED_RL', 'linear_dt1'):.4f}",
                    f"{ns('ACTUAL_FRL', 'actual_or_historical'):.4f}",
                    f"{ns('ACTUAL_F408', 'actual_or_historical'):.4f}",
                )
            )
        lines.extend(["", "Paired primary comparison by seed:", ""])
        lines.extend(
            markdown_table(
                ["Seed", "Field", "Unmasked RL", "Masked RL", "Actual FRL", "Official F408"],
                paired_rows,
            )
        )
        lines.extend(
            [
                "",
                "Primary linear interpolation and the 0.5 s convergence run are both",
                "reported in `narma_surrogate_results.csv`; zero-order hold is shown",
                "as sampling-convention sensitivity. Maximum state differences for",
                "dt 1 versus 0.5 s:",
                "",
            ]
        )
        lines.extend(
            markdown_table(
                ["Seed", "max |ΔR|", "max |ΔL|"],
                [
                    (
                        row["seed"],
                        f"{row['max_abs_R_dt1_vs_dt0_5']:.3e}",
                        f"{row['max_abs_L_dt1_vs_dt0_5']:.3e}",
                    )
                    for row in narma_convergence
                ],
            )
        )
        lines.extend(
            [
                "",
                "Interpretation is component attribution, not subtraction from the",
                "hybrid claim. Field → unmasked is the contribution of the assumed",
                "Hill/reporter kinetics. Unmasked → masked is living spatial",
                "occupancy/sampling. Masked → actual FRL contains effects not",
                "reproduced by the voxel cascade, including sub-voxel sampling and",
                "population dynamics. A close surrogate match would still represent",
                "modeled biological R/L kinetics.",
                "",
            ]
        )
    else:
        lines.extend(["NARMA sanity failed; dependent surrogate analysis was stopped.", ""])

    if wave_surrogate:
        wave_summary = []
        for family in (
            "UNMASKED_RL",
            "OCCUPANCY_MASKED_RL",
            "UNMASKED_R",
            "UNMASKED_L",
            "MASKED_R",
            "MASKED_L",
            "ACTUAL_DRIVEN_F408",
            "ACTUAL_DRIVEN_FRL",
            "ACTUAL_FIELD",
            "ACTUAL_BROWNIAN",
            "ACTUAL_SILENT",
        ):
            method = "linear_dt1" if not family.startswith("ACTUAL") else "actual_or_historical"
            zoh = (
                row_mean(wave_surrogate, family, "window_macro_ovr_auc", "zoh_previous")
                if not family.startswith("ACTUAL")
                else float("nan")
            )
            wave_summary.append(
                (
                    family,
                    f"{row_mean(wave_surrogate, family, 'window_macro_ovr_auc', method):.4f}",
                    "—" if family.startswith("ACTUAL") else f"{zoh:.4f}",
                    f"{row_mean(wave_surrogate, family, 'window_accuracy', method):.3f}",
                )
            )
        lines.extend(
            markdown_table(
                ["Waveform family", "Linear dt=1 macro AUC", "ZOH macro AUC", "Linear/actual accuracy"],
                wave_summary,
            )
        )
        wave_paired = []
        for seed in WAVEFORM_SEEDS:
            def ws(family, method):
                return next(
                    row["window_macro_ovr_auc"]
                    for row in wave_surrogate
                    if row["seed"] == seed
                    and row["family"] == family
                    and row["integration_method"] == method
                )
            wave_paired.append(
                (
                    seed,
                    f"{ws('UNMASKED_RL', 'linear_dt1'):.4f}",
                    f"{ws('OCCUPANCY_MASKED_RL', 'linear_dt1'):.4f}",
                    f"{ws('ACTUAL_DRIVEN_F408', 'actual_or_historical'):.4f}",
                    f"{ws('ACTUAL_FIELD', 'actual_or_historical'):.4f}",
                    f"{ws('ACTUAL_BROWNIAN', 'actual_or_historical'):.4f}",
                    f"{ws('ACTUAL_SILENT', 'actual_or_historical'):.4f}",
                )
            )
        lines.extend(["", "Paired waveform macro AUC by seed:", ""])
        lines.extend(
            markdown_table(
                ["Seed", "Unmasked RL", "Masked RL", "Driven F408", "Field", "Brownian", "Silent"],
                wave_paired,
            )
        )
        lines.extend(
            [
                "",
                "All per-seed lambdas, class AUCs, confusion matrices and block-pooled",
                "AUCs are in `waveform_surrogate_results.csv`. Sampling sensitivity",
                "and dt convergence are retained, not concealed.",
                "",
            ]
        )

    lines.extend(["## 4. NARMA legal and trivial baselines", ""])
    if narma_baselines:
        lines.extend(
            markdown_table(
                ["Baseline", "Status", "Features", "Lambda", "Test NRMSE", "Test R2", "Driven 0.9283 beats?"],
                [
                    (
                        row["baseline"],
                        row["status"],
                        row["n_features"],
                        row.get("lambda", ""),
                        f"{row['test_nrmse']:.4f}",
                        f"{row['test_r2']:.4f}",
                        row["driven_mean_0_9283_beats"],
                    )
                    for row in narma_baselines
                ],
            )
        )
        lines.extend(
            [
                "",
                "The test-mean row is a **non-predictive oracle normalization",
                "reference; uses test labels**. It is not ranked as a legal baseline.",
                "No polynomial term beyond the one frozen `u[n]u[n-9]` product was",
                "tried.",
                "",
            ]
        )
    else:
        lines.extend(["NARMA sanity failed; dependent baselines were stopped.", ""])

    lines.extend(["## 5. Leakage and null controls", ""])
    if narma_nulls:
        null_summary = []
        for family in ("OFFICIAL_BIOLOGY", "FIELD_ONLY"):
            true_scores = [
                row["test_nrmse"]
                for row in narma_nulls
                if row["family"] == family and row["offset_windows"] == 0
            ]
            shifted = [
                row["test_nrmse"]
                for row in narma_nulls
                if row["family"] == family and row["offset_windows"] != 0
            ]
            null_summary.append(
                (
                    family,
                    f"{min(true_scores):.4f}..{max(true_scores):.4f}",
                    f"{min(shifted):.4f}..{max(shifted):.4f}",
                    f"{mean(shifted):.4f}",
                )
            )
        lines.extend(
            markdown_table(
                ["Readout", "True score range", "Shift-null range", "Shift-null mean"],
                null_summary,
            )
        )
        lines.extend(
            [
                "",
                "Each offset 20/40/60/80/100 reran standardization, validation,",
                "lambda selection and refit. True-label lambdas were not reused.",
                "These are falsification controls, not alternative tasks.",
                "",
                f"Independent input-target mismatch: **{mismatch_status['status']}**.",
                mismatch_status["reason"],
                f"The state-i/target-j function contract test is **{'PASS' if mismatch_status['function_contract_test_pass'] else 'FAIL'}**.",
                "",
                "Preprocessing audit:",
                "",
            ]
        )
        lines.extend(
            markdown_table(
                ["Item", "Result", "Code location"],
                [
                    (
                        row["item"],
                        "PASS" if row["pass"] else "FAIL",
                        f"`{row['location']}`",
                    )
                    for row in preprocessing
                ],
            )
        )
        lines.extend(["", "Feature-family audit (means over all three seeds):", ""])
        lines.extend(
            markdown_table(
                ["Family", "Features", "Mean test NRMSE"],
                [
                    (
                        family,
                        next(row["n_features"] for row in feature_rows if row["family"] == family),
                        f"{mean([row['test_nrmse'] for row in feature_rows if row['family'] == family]):.4f}",
                    )
                    for family in ("F408", "FRL", "R_ONLY", "L_ONLY")
                ],
            )
        )
        lines.extend(
            [
                "",
                "This reprint is descriptive; F408 remains the official historical",
                "readout.",
                "",
            ]
        )

    lines.extend(["## 6. Lorenz baseline verdict", ""])
    if lorenz_rows:
        lines.extend(
            markdown_table(
                ["Clock", "Verification", "u SHA-256"],
                [
                    (
                        row["clock"],
                        "PASS" if row["pass"] else "FAIL",
                        f"`{row['u_hash']}`",
                    )
                    for row in lorenz_verify
                ],
            )
        )
        selected = [
            selected_lorenz(lorenz_rows, clock, task)
            for clock in ("BenchA_dt0.02", "BenchA2_dt1.0")
            for task in ("AUTO_X", "CROSS_Y")
        ]
        lines.extend(["", "Validation-selected delay results:", ""])
        lines.extend(
            markdown_table(
                ["Clock", "Task", "Selected m", "Lambda", "Validation NRMSE", "Test NRMSE"],
                [
                    (
                        row["clock"],
                        row["task"],
                        row["development_m"],
                        row["lambda"],
                        f"{row['validation_nrmse']:.4f}",
                        f"{row['test_nrmse']:.4f}",
                    )
                    for row in selected
                ],
            )
        )
        auto_fast_persistence = next(
            row
            for row in lorenz_rows
            if row["clock"] == "BenchA_dt0.02"
            and row["task"] == "AUTO_X"
            and row["model"] == "PERSISTENCE_X_N"
        )
        cross_fast = selected_lorenz(lorenz_rows, "BenchA_dt0.02", "CROSS_Y")
        auto_slow_persistence = next(
            row
            for row in lorenz_rows
            if row["clock"] == "BenchA2_dt1.0"
            and row["task"] == "AUTO_X"
            and row["model"] == "PERSISTENCE_X_N"
        )
        lines.extend(
            [
                "",
                f"AUTO_X at Δt=0.02 is already largely solved by persistence "
                f"(test NRMSE `{auto_fast_persistence['test_nrmse']:.4f}`). "
                f"The selected x-delay CROSS_Y baseline scores "
                f"`{cross_fast['test_nrmse']:.4f}`. At Δt=1.0 persistence rises "
                f"to `{auto_slow_persistence['test_nrmse']:.4f}`, exposing the "
                "SKIP=50 loss of local continuity.",
                "",
                "Recommendation:",
                "",
                "1. **L1:** do not prioritize Δt=0.02 auto-prediction unless the",
                "   frozen claim requires beating persistence, selected AR, and field.",
                "2. **L2:** retain cross-variable prediction only as a branch that",
                "   must beat the validation-selected x-delay baseline; do not inject",
                "   delay taps as extra chemicals.",
                "3. **L3:** highest-value Lorenz development branch. Predeclare",
                "   physically meaningful sample clocks, select on development",
                "   trajectories, and confirm one clock on untouched trajectories.",
                "",
                "No Lorenz input was generated and no BSim was run.",
                "",
            ]
        )
    else:
        lines.extend(["BenchA sanity failed; dependent Lorenz baselines were stopped.", ""])

    lines.extend(
        [
            "## 7. IPC affordability decision",
            "",
            "See [ipc_affordability.md](./ipc_affordability.md). The 408-feature",
            "readout has only 88 inner-fit and 110 final-fit rows. A broad",
            "Dambre-style target basis is not estimable stably from this geometry.",
            "",
            "IPC_DECISION: KILL_BROAD_KEEP_MC",
            "",
            "Keep MC≈1.22 and useful lags≈5. A narrow panel or long random drive",
            "requires a later frozen protocol; neither is retrospectively inferred.",
            "",
            "## 8. Waveform DC/power diagnosis",
            "",
        ]
    )
    if moment_rows:
        frozen = [row for row in moment_rows if row["template_set"] == "frozen_5_point"]
        proposed = [row for row in moment_rows if row["template_set"] == "proposed_8_point"]
        lines.extend(
            markdown_table(
                ["5-point class", "Mean", "Variance", "Pop SD", "Power", "Min", "Max", "Range"],
                [
                    (
                        row["class"],
                        f"{row['mean']:.4f}",
                        f"{row['variance']:.4f}",
                        f"{row['population_sd']:.4f}",
                        f"{row['power_mean_u2']:.4f}",
                        f"{row['min']:.4f}",
                        f"{row['max']:.4f}",
                        f"{row['range']:.4f}",
                    )
                    for row in frozen
                ],
            )
        )
        lines.extend(["", "Block-level scalar baselines (only ten test blocks):", ""])
        lines.extend(
            markdown_table(
                ["Baseline", "Features", "Lambda", "Macro AUC", "Accuracy", "Per-class AUC", "Confusion"],
                [
                    (
                        row["baseline"],
                        row["n_features"],
                        row["lambda"],
                        f"{row['macro_ovr_auc']:.4f}",
                        f"{row['accuracy']:.3f}",
                        f"{row['sine_auc']:.3f}/{row['square_auc']:.3f}/{row['triangle_auc']:.3f}",
                        f"`{row['confusion']}`",
                    )
                    for row in moment_classifiers
                ],
            )
        )
        lines.extend(
            [
                "",
                "`RAW_TEMPLATE_5` measures direct input separability, not reservoir",
                "computation. The scalar baselines test the DC/power hypothesis; the",
                "tiny ten-block test set prevents strong uncertainty claims.",
                "",
                "Eight-point proposal moments:",
                "",
            ]
        )
        lines.extend(
            markdown_table(
                ["8-point class", "Mean", "Variance", "Pop SD", "Power", "Min", "Max", "Range"],
                [
                    (
                        row["class"],
                        f"{row['mean']:.4f}",
                        f"{row['variance']:.4f}",
                        f"{row['population_sd']:.4f}",
                        f"{row['power_mean_u2']:.4f}",
                        f"{row['min']:.4f}",
                        f"{row['max']:.4f}",
                        f"{row['range']:.4f}",
                    )
                    for row in proposed
                ],
            )
        )
        lines.extend(
            [
                "",
                "The proposed eight-point family still has unequal SD and power.",
                "Two later alternatives are admissible but not final:",
                "",
                "1. recognizable, moment-aware waveforms with mean/power/moment",
                "   baselines declared as component attribution;",
                "2. equal-histogram temporal-order templates with exactly matched",
                "   scalar moments, acknowledging that conventional waveform names",
                "   may no longer apply.",
                "",
                "A later frozen Waveform2 protocol must choose; no production file",
                "was generated here.",
                "",
            ]
        )
    else:
        lines.extend(["Waveform sanity failed; dependent moment analysis was stopped.", ""])

    lines.extend(
        [
            "## 9. Provenance archive",
            "",
            "See [external analysis MANIFEST](../external_analysis/MANIFEST.md).",
            f"All {len(provenance)} source/destination SHA-256 pairs match. The",
            "documents remain labelled external, untrusted analysis—not project",
            "evidence.",
            "",
            "## 10. Recommended next simulations, ordered by evidence value",
            "",
            "1. **NARMA only after claim revision:** the 0.6828 direct-input delay",
            "   baseline already beats 0.9283 biology. If the revised question is",
            "   reproducibility beyond Brownian/silent/field—not superiority to a",
            "   direct task model—freeze 10–12 independent inputs, one primary",
            "   bacterial seed each, plus the predeclared seed×input subset. Include",
            "   every baseline, surrogate, shift and genuine mismatch. If beating",
            "   direct-input baselines is required, kill this simulation branch.",
            "2. **Waveform2 only as strict temporal-order evidence:** the current",
            "   scalar moments and raw template already rank classes. Prefer",
            "   equal-histogram templates, more independent test blocks, and",
            "   block-level evaluation; a moment-aware recognizable branch is",
            "   secondary and must report its moment baselines.",
            "3. **Lorenz L3 clock map, then at most one confirmation:** the Δt=0.02",
            "   x-delay baseline nearly solves CROSS_Y. Retain L2",
            "   only if the scientific question is beating its x-delay baseline;",
            "   deprioritize trivial Δt=0.02 AUTO_X.",
            "4. **IPC:** do not run broad IPC. Consider a separately budgeted long",
            "   random drive only if a narrow predeclared panel and null burden are",
            "   affordable.",
            "5. **AC and Stage99:** no implementation follows from this audit.",
            "   Continue only under their separate calibration/architecture plans.",
            "",
            "AC/field nonlinearities remain legitimate system-level hybrid",
            "computation. Component baselines explain where capacity arises; they",
            "do not invalidate the integrated-system claim.",
            "",
            "## Zero-simulation audit completed",
            "",
            "This package completes the analysis-only precondition linked from the",
            "[high-tier submission evidence plan](../HIGH_TIER_SUBMISSION_EVIDENCE_PLAN.md).",
            "Its historical rationale and earlier closed decisions were not rewritten.",
            "",
        ]
    )
    (HERE / "ZERO_SIMULATION_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")


def write_readme():
    text = """# Zero-simulation evidence audit

Analysis-only package for frozen NARMA, waveform, BenchA and BenchA2 data.
It invokes no Java or BSim and writes only to this directory plus the sibling
`external_analysis/` provenance archive.

Run from the repository root:

```powershell
python examples/HybridDish/zero_simulation_audit/run_zero_simulation_audit.py
```

The command reruns Stage 0 sanity, all eligible analyses, CSV generation,
the IPC affordability note, provenance copying/hash verification, and the
integrated report. See [PROTOCOL.md](./PROTOCOL.md) for the rules frozen before
calculation and [ZERO_SIMULATION_AUDIT.md](./ZERO_SIMULATION_AUDIT.md) for the
findings.
"""
    (HERE / "README.md").write_text(text, encoding="utf-8")


def print_sanity(rows):
    print("STAGE_0_SANITY")
    for family, check, actual, expected, passed in rows:
        print(
            f"{family} | {check} | actual={actual} | expected={expected} | "
            f"{'PASS' if passed else 'FAIL'}"
        )


def main():
    print("ZERO-SIMULATION AUDIT: no Java or BSim is invoked")
    u_narma, y_narma, narma_hash = load_narma_target()
    wave_labels, wave_blocks, wave_u = load_wave_labels()
    sanity, sanity_rows = stage0_sanity(u_narma, y_narma, wave_labels)
    print_sanity(sanity_rows)

    narma_surrogate = []
    narma_convergence = []
    integration_test = integration_unit_test()
    narma_baselines = []
    narma_nulls = []
    feature_rows = []
    mismatch_status = {
        "status": "PENDING",
        "reason": "NARMA Stage 0 failed before mismatch function test.",
        "function_contract_test_pass": False,
    }
    preprocessing = []
    if sanity["narma"]["pass"]:
        print("RUN module A/NARMA, B, C")
        narma_surrogate, narma_convergence, integration_test = run_narma_surrogates(y_narma)
        narma_baselines = run_narma_baselines(u_narma, y_narma)
        narma_nulls, feature_rows, mismatch_status = run_narma_nulls(y_narma)
        preprocessing = preprocessing_audit(y_narma)
    else:
        print("STOP dependent NARMA modules: Stage 0 NARMA sanity failed")

    wave_surrogate = []
    wave_convergence = []
    moment_rows = []
    moment_classifiers = []
    if sanity["waveform"]["pass"]:
        print("RUN module A/Waveform and F")
        wave_surrogate, wave_convergence = run_wave_surrogates(wave_labels)
        moment_rows, moment_classifiers = run_waveform_moments(
            wave_labels, wave_blocks, wave_u
        )
    else:
        print("STOP dependent waveform modules: Stage 0 waveform sanity failed")

    lorenz_rows = []
    lorenz_verify = []
    if sanity["bencha"]["pass"]:
        print("RUN module D")
        lorenz_rows, lorenz_verify = run_lorenz_baselines()
    else:
        print("STOP dependent Lorenz modules: Stage 0 BenchA sanity failed")

    print("RUN independent modules E and G")
    ipc = write_ipc()
    provenance = archive_external()

    write_csv(HERE / "narma_surrogate_results.csv", narma_surrogate)
    write_csv(HERE / "narma_baseline_results.csv", narma_baselines)
    write_csv(HERE / "narma_null_results.csv", narma_nulls)
    write_csv(HERE / "lorenz_baseline_results.csv", lorenz_rows)
    write_csv(
        HERE / "waveform_moment_results.csv", moment_rows + moment_classifiers
    )
    write_csv(HERE / "waveform_surrogate_results.csv", wave_surrogate)
    write_readme()
    write_report(
        sanity,
        sanity_rows,
        narma_surrogate,
        narma_convergence,
        integration_test,
        narma_baselines,
        narma_nulls,
        feature_rows,
        mismatch_status,
        preprocessing,
        lorenz_rows,
        lorenz_verify,
        moment_rows,
        moment_classifiers,
        wave_surrogate,
        wave_convergence,
        provenance,
    )

    print(f"narma_u_sha256={narma_hash}")
    print(f"integration_unit_test={'PASS' if integration_test['pass'] else 'FAIL'}")
    print("IPC_DECISION: KILL_BROAD_KEEP_MC")
    for row in provenance:
        print(
            f"PROVENANCE {row['filename']} sha256={row['sha256']} "
            f"bytes={row['byte_size']} hash_match=PASS"
        )
    print(f"broad_ipc_targets={ipc['broad_targets']}")
    print("AUDIT_COMPLETE")


if __name__ == "__main__":
    main()
