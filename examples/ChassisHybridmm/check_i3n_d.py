#!/usr/bin/env python3
"""Run the frozen, post-hoc ChassisHybridmm I3n-D diagnostics."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
NARMA_DIR = HERE.parent / "BSimReservoirPlanNarma10b"
PROTOCOL = HERE / "I3N_D_PROTOCOL.md"
EVIDENCE_PATH = RESULTS / "i3n_d_evidence.json"
STANDING_PATH = HERE / "I3N_D_STANDING.md"
SEEDS = (111, 222, 333)
NUM_WINDOWS = 200
NUM_SAMPLES = 16
NUM_BINS = 200
NUM_DEATHS = 8
WASHOUT = 40
TRAIN = 110
INNER_TRAIN = 88
RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)


def load_narma_module():
    path = NARMA_DIR / "check_narma10b.py"
    spec = importlib.util.spec_from_file_location("closed_narma10b_checker", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NARMA = load_narma_module()


def run_dir(seed: int) -> Path:
    return RESULTS / f"i3n_driven_seed{seed}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mean_se(values: list[float]) -> tuple[float, float]:
    data = np.asarray(values, dtype=float)
    return float(np.mean(data)), float(np.std(data, ddof=1) / np.sqrt(len(data)))


def require_protocol() -> None:
    text = PROTOCOL.read_text(encoding="utf-8")
    required = (
        "Frozen 2026-08-25 before computing any I3n-D NRMSE",
        "windows `128..149`",
        "missing standardized coordinate to `0`",
        "cell-mean \\(R/L\\) ridge",
        "not a NARMA Overall",
    )
    missing = [phrase for phrase in required if phrase not in text]
    if missing:
        raise ValueError(f"I3n-D protocol is absent or incomplete: {missing}")


def read_summary(seed: int) -> dict:
    path = run_dir(seed) / "window_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    windows = [int(row["Window"]) for row in rows]
    if windows != list(range(NUM_WINDOWS)):
        raise ValueError(f"{path}: expected windows 0..199")
    mean_r = np.asarray([float(row["Mean_R"]) for row in rows], dtype=float)
    mean_l = np.asarray([float(row["Mean_L"]) for row in rows], dtype=float)
    input_u = np.asarray([float(row["AHL_Input"]) for row in rows], dtype=float)
    population = np.asarray([int(row["Population"]) for row in rows], dtype=int)
    return {
        "path": str(path),
        "X": np.column_stack([mean_r, mean_l]),
        "windows": windows,
        "mean_r": mean_r,
        "input_u": input_u,
        "population": population,
    }


def indexed_columns(header: list[str], prefix: str, count: int) -> list[int]:
    names = [f"{prefix}{index}" for index in range(count)]
    missing = [name for name in names if name not in header]
    if missing:
        raise ValueError(f"missing columns with prefix {prefix}: {missing[:3]}")
    return [header.index(name) for name in names]


def read_masked_voxels(seed: int) -> dict:
    path = run_dir(seed) / "voxels.csv"
    sums_r = np.zeros((NUM_WINDOWS, NUM_BINS), dtype=float)
    sums_l = np.zeros((NUM_WINDOWS, NUM_BINS), dtype=float)
    occupied_samples = np.zeros((NUM_WINDOWS, NUM_BINS), dtype=int)
    deaths = np.zeros((NUM_WINDOWS, NUM_DEATHS), dtype=float)
    samples_per_window = np.zeros(NUM_WINDOWS, dtype=int)
    sample_occupied_counts: list[int] = []

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        window_idx = header.index("Window")
        den_idx = indexed_columns(header, "Den_", NUM_BINS)
        r_idx = indexed_columns(header, "Receiver_R_", NUM_BINS)
        l_idx = indexed_columns(header, "Lum_Mean_", NUM_BINS)
        death_idx = indexed_columns(header, "Input_Driven_Death_", NUM_DEATHS)
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"{path}:{row_number}: ragged row")
            window = int(row[window_idx])
            if window < 0 or window >= NUM_WINDOWS:
                raise ValueError(f"{path}:{row_number}: invalid window {window}")
            density = np.fromiter((float(row[index]) for index in den_idx), dtype=float)
            occupied = density > 0.0
            r_values = np.fromiter((float(row[index]) for index in r_idx), dtype=float)
            l_values = np.fromiter((float(row[index]) for index in l_idx), dtype=float)
            sums_r[window, occupied] += r_values[occupied]
            sums_l[window, occupied] += l_values[occupied]
            occupied_samples[window, occupied] += 1
            deaths[window] = np.fromiter(
                (float(row[index]) for index in death_idx), dtype=float
            )
            samples_per_window[window] += 1
            sample_occupied_counts.append(int(np.sum(occupied)))

    if not np.all(samples_per_window == NUM_SAMPLES):
        raise ValueError(
            f"{path}: samples/window range "
            f"{samples_per_window.min()}..{samples_per_window.max()}, expected 16"
        )

    observed_bins = occupied_samples > 0
    masked_r = np.full((NUM_WINDOWS, NUM_BINS), np.nan, dtype=float)
    masked_l = np.full((NUM_WINDOWS, NUM_BINS), np.nan, dtype=float)
    masked_r[observed_bins] = sums_r[observed_bins] / occupied_samples[observed_bins]
    masked_l[observed_bins] = sums_l[observed_bins] / occupied_samples[observed_bins]
    X = np.column_stack([masked_r, masked_l, deaths])
    observed = np.column_stack(
        [
            observed_bins,
            observed_bins,
            np.ones((NUM_WINDOWS, NUM_DEATHS), dtype=bool),
        ]
    )
    return {
        "path": str(path),
        "X": X,
        "observed": observed,
        "windows": list(range(NUM_WINDOWS)),
        "sample_occupied_counts": sample_occupied_counts,
        "window_occupied_counts": np.sum(observed_bins, axis=1).astype(int),
    }


def masked_standardize(
    X_fit: np.ndarray,
    observed_fit: np.ndarray,
    X_other: np.ndarray,
    observed_other: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    spatial_count = 2 * NUM_BINS
    means = np.zeros(X_fit.shape[1], dtype=float)
    stds = np.ones(X_fit.shape[1], dtype=float)
    keep = np.ones(X_fit.shape[1], dtype=bool)

    for column in range(spatial_count):
        values = X_fit[observed_fit[:, column], column]
        if len(values) == 0:
            keep[column] = False
            continue
        means[column] = float(np.mean(values))
        std = float(np.std(values))
        if std < 1e-12:
            keep[column] = False
            continue
        stds[column] = std

    for column in range(spatial_count, X_fit.shape[1]):
        values = X_fit[:, column]
        means[column] = float(np.mean(values))
        std = float(np.std(values))
        stds[column] = std if std >= 1e-12 else 1.0

    def convert(X: np.ndarray, observed: np.ndarray) -> np.ndarray:
        selected = np.flatnonzero(keep)
        converted = np.zeros((len(X), len(selected)), dtype=float)
        for destination, source in enumerate(selected):
            present = observed[:, source]
            converted[present, destination] = (
                X[present, source] - means[source]
            ) / stds[source]
        return converted

    return convert(X_fit, observed_fit), convert(X_other, observed_other), keep


def evaluate_masked(X_all: np.ndarray, observed_all: np.ndarray, target: np.ndarray) -> dict:
    X = X_all[WASHOUT:]
    observed = observed_all[WASHOUT:]
    y = target[WASHOUT:]
    X_inner = X[:INNER_TRAIN]
    observed_inner = observed[:INNER_TRAIN]
    X_val = X[INNER_TRAIN:TRAIN]
    observed_val = observed[INNER_TRAIN:TRAIN]
    y_inner = y[:INNER_TRAIN]
    y_val = y[INNER_TRAIN:TRAIN]
    inner_z, val_z, inner_keep = masked_standardize(
        X_inner, observed_inner, X_val, observed_val
    )

    candidates = []
    grid_scores = []
    for lam in RIDGE_GRID:
        weights = NARMA.ridge_fit(inner_z, y_inner, lam)
        prediction = NARMA.ridge_predict(val_z, weights)
        score = float(NARMA.nrmse(y_val, prediction))
        candidates.append((score, -lam, lam))
        grid_scores.append({"lambda": lam, "val_nrmse": score})
    candidates.sort()
    selected_lambda = float(candidates[0][2])

    X_train = X[:TRAIN]
    observed_train = observed[:TRAIN]
    X_test = X[TRAIN:]
    observed_test = observed[TRAIN:]
    train_z, test_z, final_keep = masked_standardize(
        X_train, observed_train, X_test, observed_test
    )
    weights = NARMA.ridge_fit(train_z, y[:TRAIN], selected_lambda)
    train_prediction = NARMA.ridge_predict(train_z, weights)
    test_prediction = NARMA.ridge_predict(test_z, weights)
    return {
        "lambda": selected_lambda,
        "lambda_grid": grid_scores,
        "train_nrmse": float(NARMA.nrmse(y[:TRAIN], train_prediction)),
        "test_nrmse": float(NARMA.nrmse(y[TRAIN:], test_prediction)),
        "test_r2": float(NARMA.r_squared(y[TRAIN:], test_prediction)),
        "inner_retained_features": int(np.sum(inner_keep)),
        "final_retained_features": int(np.sum(final_keep)),
        "final_retained_spatial_features": int(np.sum(final_keep[: 2 * NUM_BINS])),
        "death_features": NUM_DEATHS,
        "missing_spatial_fraction_all_windows": float(
            1.0 - np.mean(observed_all[:, : 2 * NUM_BINS])
        ),
        "missing_spatial_fraction_train_windows": float(
            1.0 - np.mean(observed_all[WASHOUT : WASHOUT + TRAIN, : 2 * NUM_BINS])
        ),
        "missing_spatial_fraction_test_windows": float(
            1.0 - np.mean(observed_all[WASHOUT + TRAIN :, : 2 * NUM_BINS])
        ),
    }


def validate_frozen_i3n_evidence(evidence: dict) -> None:
    if evidence.get("overall") != "FAIL":
        raise ValueError("frozen I3n evidence no longer says FAIL")
    expected = {
        "driven": 1.1980451014454812,
        "brownian": 1.162790759205085,
        "silent": 1.1621883174398704,
    }
    for arm, value in expected.items():
        actual = float(evidence["mean_nrmse"][arm])
        if not math.isclose(actual, value, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"frozen I3n {arm} mean changed: {actual} != {value}")


def write_standing(evidence: dict) -> None:
    masked = evidence["diagnostic_scores"]["masked_spatial"]
    scalar = evidence["diagnostic_scores"]["cell_mean_rl"]
    nulls = evidence["frozen_i3n_context"]["mean_nrmse"]
    lines = [
        "# ChassisHybridmm I3n-D standing",
        "",
        "**I3n-D: POST-HOC DIAGNOSTIC — NO NARMA OVERALL**",
        "",
        "I3n remains **OVERALL FAIL**. Its frozen unmasked driven mean is",
        f"`{nulls['driven']:.4f}` versus Brownian `{nulls['brownian']:.4f}`",
        f"and silent `{nulls['silent']:.4f}`. Field AHL remains",
        f"`{evidence['frozen_i3n_context']['field_mean_nrmse']:.4f}`.",
        "",
        "**I3n-D answer:** the FAIL comparison does not survive either",
        "declared diagnostic readout. This establishes that A1 was load-bearing;",
        "it does not promote either post-hoc score to an I3n System result.",
        "",
        "Every I3n-D number below is post-hoc. No simulation was rerun, no",
        "seed was dropped, and no frozen I3n CSV was edited.",
        "",
        "## Diagnostic scores",
        "",
        "| Post-hoc readout | seed 111 | seed 222 | seed 333 | mean ± s.e. | Beats both frozen null means? |",
        "|---|---:|---:|---:|---:|---|",
        (
            f"| Occupancy-masked spatial R/L + 8 deaths | "
            f"{masked['test_nrmse'][0]:.4f} | {masked['test_nrmse'][1]:.4f} | "
            f"{masked['test_nrmse'][2]:.4f} | {masked['mean_nrmse']:.4f} ± "
            f"{masked['se_nrmse']:.4f} | {'YES' if masked['beats_both_nulls'] else 'NO'} |"
        ),
        (
            f"| Cell-mean R/L (2 features) | "
            f"{scalar['test_nrmse'][0]:.4f} | {scalar['test_nrmse'][1]:.4f} | "
            f"{scalar['test_nrmse'][2]:.4f} | {scalar['mean_nrmse']:.4f} ± "
            f"{scalar['se_nrmse']:.4f} | {'YES' if scalar['beats_both_nulls'] else 'NO'} |"
        ),
        f"| Frozen Brownian Den null | — | — | — | {nulls['brownian']:.4f} | context |",
        f"| Frozen silent null | — | — | — | {nulls['silent']:.4f} | context |",
        "",
        "Lambda was selected only on windows `128..149`. Empty bin samples",
        "were missing, not raw `R=L=0`; missing standardized coordinates were",
        "neutral training-mean coordinates. The masked score is not the frozen",
        "I3n 408-D score.",
        "",
        "Masked final retained feature counts (including eight death columns):",
        "`"
        + " / ".join(str(value) for value in masked["final_retained_features"])
        + "` for seeds 111/222/333. Spatial missing fractions over all windows:",
        "`"
        + " / ".join(f"{value:.3f}" for value in masked["missing_spatial_fraction_all_windows"])
        + "`.",
        "",
        "## Population and map occupancy",
        "",
        "| Seed | sample occupied bins, median [range] | window occupied bins, median [range] | N window 0 → 150 → 199 | corr(Mean_R, u) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for seed in SEEDS:
        row = evidence["descriptive"][str(seed)]
        lines.append(
            f"| {seed} | {row['sample_occupied_bins']['median']:.1f} "
            f"[{row['sample_occupied_bins']['minimum']}–{row['sample_occupied_bins']['maximum']}] | "
            f"{row['window_occupied_bins']['median']:.1f} "
            f"[{row['window_occupied_bins']['minimum']}–{row['window_occupied_bins']['maximum']}] | "
            f"{row['population']['window_0']} → {row['population']['window_150']} → "
            f"{row['population']['window_199']} | {row['corr_mean_r_u']:.3f} |"
        )
    lines.extend(
        [
            "",
            "The full 200-value per-window occupied-bin and population traces are",
            "stored in `results/i3n_d_evidence.json`.",
            "",
            "## Frozen interpretation applied",
            "",
            evidence["interpretation"],
            "",
            "This result does not rewrite I3n, does not resolve the clock mismatch",
            "(A2), and does not establish delayed full-dish map occupancy (A3).",
            "The next named job must freeze cell occupancy, map occupancy, and empty",
            "encoding before code.",
        ]
    )
    STANDING_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    require_protocol()
    _, target, input_digest = NARMA.load_target(
        HERE / "narma10_target.csv", HERE / "input_ahl_narma200.txt"
    )
    frozen_path = RESULTS / "i3n_evidence.json"
    frozen = json.loads(frozen_path.read_text(encoding="utf-8"))
    validate_frozen_i3n_evidence(frozen)

    input_paths = [
        path
        for seed in SEEDS
        for path in (
            run_dir(seed) / "voxels.csv",
            run_dir(seed) / "window_summary.csv",
        )
    ]
    hashes_before = {str(path.relative_to(HERE)): sha256(path) for path in input_paths}

    masked_scores = []
    scalar_scores = []
    masked_details = []
    scalar_details = []
    descriptive = {}
    for seed in SEEDS:
        summary = read_summary(seed)
        voxels = read_masked_voxels(seed)
        masked_result = evaluate_masked(voxels["X"], voxels["observed"], target)
        scalar_result = NARMA.evaluate_readout(
            summary["X"], target, summary["windows"]
        )
        masked_scores.append(float(masked_result["test_nrmse"]))
        scalar_scores.append(float(scalar_result["test_nrmse"]))
        masked_details.append(masked_result)
        scalar_details.append(scalar_result)

        sample_counts = np.asarray(voxels["sample_occupied_counts"], dtype=int)
        window_counts = np.asarray(voxels["window_occupied_counts"], dtype=int)
        population = summary["population"]
        descriptive[str(seed)] = {
            "sample_occupied_bins": {
                "median": float(np.median(sample_counts)),
                "minimum": int(np.min(sample_counts)),
                "maximum": int(np.max(sample_counts)),
            },
            "window_occupied_bins": {
                "values": window_counts.tolist(),
                "median": float(np.median(window_counts)),
                "minimum": int(np.min(window_counts)),
                "maximum": int(np.max(window_counts)),
            },
            "population": {
                "values": population.tolist(),
                "window_0": int(population[0]),
                "window_150": int(population[150]),
                "window_199": int(population[199]),
                "minimum": int(np.min(population)),
                "maximum": int(np.max(population)),
            },
            "corr_mean_r_u": float(
                np.corrcoef(summary["mean_r"], summary["input_u"])[0, 1]
            ),
        }

    hashes_after = {str(path.relative_to(HERE)): sha256(path) for path in input_paths}
    if hashes_before != hashes_after:
        raise RuntimeError("a frozen I3n CSV changed while I3n-D was running")

    masked_mean, masked_se = mean_se(masked_scores)
    scalar_mean, scalar_se = mean_se(scalar_scores)
    brownian_mean = float(frozen["mean_nrmse"]["brownian"])
    silent_mean = float(frozen["mean_nrmse"]["silent"])
    masked_beats = masked_mean < brownian_mean and masked_mean < silent_mean
    scalar_beats = scalar_mean < brownian_mean and scalar_mean < silent_mean
    if not masked_beats and not scalar_beats:
        interpretation = (
            "Neither post-hoc readout beats both frozen null means. Empty-zero "
            "encoding is not sufficient to explain I3n FAIL; A2 and A3 remain."
        )
    elif masked_beats and scalar_beats:
        interpretation = (
            "Both post-hoc readouts beat both frozen null means. Empty-zero "
            "encoding contributed to I3n FAIL, but I3n remains FAIL and A2/A3 remain."
        )
    else:
        winner = "masked spatial" if masked_beats else "cell-mean R/L"
        loser = "cell-mean R/L" if masked_beats else "masked spatial"
        interpretation = (
            f"The diagnostics split: {winner} beats both frozen null means, while "
            f"{loser} does not. Zero encoding contributed for one declared readout; "
            "the favorable readout is not promoted to an I3n PASS, and A2/A3 remain."
        )

    evidence = {
        "schema": "ChassisHybridmm-I3n-D-post-hoc-v1",
        "protocol": str(PROTOCOL),
        "status": "POST_HOC_DIAGNOSTIC_NO_OVERALL",
        "seeds": list(SEEDS),
        "input_u_sha256": input_digest,
        "ridge_hygiene": {
            "washout_windows": "0..39",
            "train_windows": "40..149",
            "lambda_validation_windows": "128..149 only",
            "test_windows": "150..199",
            "lambda_grid": list(RIDGE_GRID),
            "intercept_regularized": False,
            "test_used_for_selection_or_standardization": False,
        },
        "frozen_i3n_context": {
            "overall": frozen["overall"],
            "test_nrmse": frozen["test_nrmse"],
            "mean_nrmse": frozen["mean_nrmse"],
            "field_mean_nrmse": frozen["field_mean_nrmse"],
        },
        "diagnostic_scores": {
            "masked_spatial": {
                "label": "post-hoc occupancy-masked spatial R/L plus 8 deaths",
                "test_nrmse": masked_scores,
                "mean_nrmse": masked_mean,
                "se_nrmse": masked_se,
                "lambdas": [float(item["lambda"]) for item in masked_details],
                "final_retained_features": [
                    item["final_retained_features"] for item in masked_details
                ],
                "final_retained_spatial_features": [
                    item["final_retained_spatial_features"] for item in masked_details
                ],
                "missing_spatial_fraction_all_windows": [
                    item["missing_spatial_fraction_all_windows"]
                    for item in masked_details
                ],
                "details": masked_details,
                "beats_both_nulls": masked_beats,
            },
            "cell_mean_rl": {
                "label": "post-hoc cell-mean R/L two-feature ridge",
                "test_nrmse": scalar_scores,
                "mean_nrmse": scalar_mean,
                "se_nrmse": scalar_se,
                "lambdas": [float(item["lambda"]) for item in scalar_details],
                "details": scalar_details,
                "beats_both_nulls": scalar_beats,
            },
        },
        "descriptive": descriptive,
        "interpretation": interpretation,
        "frozen_csv_sha256": hashes_before,
        "frozen_csv_unchanged": hashes_before == hashes_after,
        "i3n_overall_rewritten": False,
    }
    EVIDENCE_PATH.write_text(
        json.dumps(evidence, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    write_standing(evidence)

    print("I3n-D POST-HOC DIAGNOSTIC — NO OVERALL")
    print(
        "masked spatial: "
        + " / ".join(f"{value:.4f}" for value in masked_scores)
        + f" mean={masked_mean:.4f}"
    )
    print(
        "cell-mean R/L: "
        + " / ".join(f"{value:.4f}" for value in scalar_scores)
        + f" mean={scalar_mean:.4f}"
    )
    print(interpretation)
    print("I3n OVERALL remains FAIL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
