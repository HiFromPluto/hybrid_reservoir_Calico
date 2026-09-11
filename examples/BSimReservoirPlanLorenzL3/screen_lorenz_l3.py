#!/usr/bin/env python3
"""Module 1: Lorenz L3 trivial-baseline screen. No Java, no BSim.

Hashes in results/lorenz_l3_input_hashes.csv must exist before this
script ranks anything. Frozen labels:

  TRIVIAL:   best legal test NRMSE <= 0.30
  DESTROYED: best legal test NRMSE >= 0.95
  SURVIVE:   0.30 < best legal test NRMSE < 0.95

k=1 must be TRIVIAL and k=50 DESTROYED if zero-sim selected-delay
test NRMSE reproduces within 1e-3. If that sanity fails, stop.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EXAMPLES / "BSimReservoirPlanBenchA"))
import check_bencha as BC  # noqa: E402
from generate_lorenz_l3 import integrate as integrate_lorenz  # noqa: E402

WASHOUT = 40
TRAIN = 110
TEST = 50
INNER_VAL = 22
NUM_WINDOWS = 200
M_GRID = (1, 3, 5, 7, 10)
TRIVIAL_MAX = 0.30
DESTROYED_MIN = 0.95
SANITY_TOL = 1e-3

CLOCKS = (
    (1, "anchor", "BenchA_k1"),
    (5, "development", "L3_k5"),
    (10, "development", "L3_k10"),
    (20, "development", "L3_k20"),
    (40, "development", "L3_k40"),
    (50, "anchor", "BenchA2_k50"),
)

ZERO_SIM_SELECTED_TEST = {
    (1, "AUTO_X"): 0.0014,
    (1, "CROSS_Y"): 0.0118,
    (50, "AUTO_X"): 1.0047,
    (50, "CROSS_Y"): 1.0019,
}

BENCHA = EXAMPLES / "BSimReservoirPlanBenchA"
BENCHA2 = EXAMPLES / "BSimReservoirPlanBenchA2"


def sha256_u(u: np.ndarray) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def load_hashes() -> dict:
    path = HERE / "results" / "lorenz_l3_input_hashes.csv"
    if not path.exists():
        raise SystemExit(
            "STOP: hashes file missing. Run generate_lorenz_l3.py before ranking."
        )
    rows = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows[int(row["k"])] = row
    for skip, _, _ in CLOCKS:
        if skip not in rows:
            raise SystemExit(f"STOP: hash missing for k={skip}")
    return rows


def load_stored_u(path: Path) -> np.ndarray:
    values = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values.append(float(stripped))
    return np.array(values, dtype=float)


def delay_matrix(values: np.ndarray, taps: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    X = np.zeros((len(values), taps), dtype=float)
    for n in range(len(values)):
        for lag in range(taps):
            if n - lag >= 0:
                X[n, lag] = values[n - lag]
    return X


def scalar_closed(X_all: np.ndarray, target: np.ndarray) -> dict:
    X_all = np.asarray(X_all, dtype=float)
    target = np.asarray(target, dtype=float)
    X = X_all[WASHOUT:]
    y = target[WASHOUT:]
    lam, grid = BC.select_lambda(X, y)
    metrics = BC.fit_eval(X, y, lam)
    val = next(row["val_nrmse"] for row in grid if row["lambda"] == lam)
    return {
        "lambda": float(lam),
        "train_nrmse": float(metrics["train_nrmse"]),
        "validation_nrmse": float(val),
        "test_nrmse": float(metrics["test_nrmse"]),
        "test_r2": float(metrics["test_r2"]),
    }


def constant_metrics(target: np.ndarray, constant: float) -> dict:
    y = np.asarray(target)[WASHOUT:]
    y_train, y_test = y[:TRAIN], y[TRAIN:]
    y_val = y[TRAIN - INNER_VAL:TRAIN]
    pred_train = np.full_like(y_train, constant)
    pred_val = np.full_like(y_val, constant)
    pred_test = np.full_like(y_test, constant)
    return {
        "lambda": "",
        "train_nrmse": BC.nrmse(y_train, pred_train),
        "validation_nrmse": BC.nrmse(y_val, pred_val),
        "test_nrmse": BC.nrmse(y_test, pred_test),
        "test_r2": BC.r_squared(y_test, pred_test),
    }


def evaluate_nontrainable(target: np.ndarray, prediction: np.ndarray) -> dict:
    y = np.asarray(target, dtype=float)[WASHOUT:]
    p = np.asarray(prediction, dtype=float)[WASHOUT:]
    return {
        "lambda": "",
        "train_nrmse": BC.nrmse(y[:TRAIN], p[:TRAIN]),
        "validation_nrmse": BC.nrmse(
            y[TRAIN - INNER_VAL:TRAIN], p[TRAIN - INNER_VAL:TRAIN]
        ),
        "test_nrmse": BC.nrmse(y[TRAIN:], p[TRAIN:]),
        "test_r2": BC.r_squared(y[TRAIN:], p[TRAIN:]),
    }


def integrate(skip: int) -> np.ndarray:
    return integrate_lorenz(skip)


def clock_paths(skip: int) -> tuple[Path, Path]:
    if skip == 1:
        return BENCHA / "input_ahl_lorenz200.txt", BENCHA / "lorenz_target.csv"
    if skip == 50:
        return (
            BENCHA2 / "input_ahl_lorenz_skip50.txt",
            BENCHA2 / "lorenz_skip50_target.csv",
        )
    return HERE / f"input_ahl_lorenz_k{skip}.txt", HERE / f"lorenz_target_k{skip}.csv"


def label_from_best_test(score: float) -> str:
    if score <= TRIVIAL_MAX:
        return "TRIVIAL"
    if score >= DESTROYED_MIN:
        return "DESTROYED"
    return "SURVIVE"


def fmt(value) -> str:
    if value == "" or value is None:
        return ""
    if isinstance(value, str):
        return value
    return f"{float(value):.10g}"


def write_markdown(path: Path, hash_rows: dict, model_rows: list, pair_rows: list,
                   sanity: dict, disposition: str, selected: dict | None) -> None:
    lines = [
        "# Lorenz L3 — sample-clock screen (Module 1)",
        "",
        "Track E3. No Java, no BSim in this module. Hashes were recorded",
        "by `generate_lorenz_l3.py` before any baseline ranking.",
        "BenchA Lorenz FAIL and BenchA2 SKIP=50 FAIL are not rewritten.",
        "",
        "Frozen labels, declared before the new-k numbers:",
        "",
        "- TRIVIAL: best legal test NRMSE ≤ 0.30",
        "- DESTROYED: best legal test NRMSE ≥ 0.95",
        "- SURVIVE: 0.30 < best legal test NRMSE < 0.95",
        "",
        "Legal AUTO_X: train-intercept, persistence x[n], validation-selected",
        "linear AR on [x[n],…,x[n−m+1]] for m∈{1,3,5,7,10}.",
        "Legal CROSS_Y: train-intercept, persistence y[n] (teacher-forced),",
        "validation-selected linear x-delay on the same m grid.",
        "Test-mean is an oracle reference, not ranked.",
        "m is selected on validation only; ties take the smaller m.",
        "Lambda ties take the larger lambda. Closed NARMA ridge.",
        "",
        f"## Disposition: {disposition}",
        "",
    ]
    if selected:
        lines.extend([
            "Selected living scout pair (Module-1 numbers only; hardest",
            "remaining linear SURVIVE clock = largest best-legal validation",
            "NRMSE):",
            "",
            f"- k = {selected['k']}",
            f"- Δt = {selected['dt_sample']}",
            f"- target = {selected['task']}",
            f"- label = {selected['label']}",
            f"- best legal validation NRMSE = {selected['best_legal_val_nrmse']:.4f}",
            f"- best legal test NRMSE = {selected['best_legal_test_nrmse']:.4f}",
            f"- selected delay m = {selected['selected_m']}",
            "",
        ])
    else:
        lines.extend([
            "No new (k, target) SURVIVEd. No living BSim.",
            "",
        ])

    lines.extend([
        "## Input hashes (before ranking)",
        "",
        "| k | Δt | Role | u SHA-256 | xmin | xmax | Regenerated |",
        "|---|---|---|---|---|---|---|",
    ])
    for skip, role, _ in CLOCKS:
        row = hash_rows[skip]
        lines.append(
            f"| {skip} | {row['dt_sample']} | {role} | `{row['u_sha256']}` | "
            f"{row['xmin']} | {row['xmax']} | {row['regenerated']} |"
        )

    lines.extend([
        "",
        "## Anchor sanity versus zero-simulation selected delay",
        "",
        f"Zero-sim reproduction within {SANITY_TOL:g}: "
        f"**{'PASS' if sanity['pass'] else 'FAIL'}**.",
        "",
        "| Clock | Task | This screen | Zero-sim | |Δ| | Match | Label |",
        "|---|---|---|---|---|---|---|",
    ])
    for row in sanity["rows"]:
        lines.append(
            f"| k={row['k']} | {row['task']} | {row['got']:.4f} | "
            f"{row['expected']:.4f} | {row['delta']:.4g} | {row['match']} | "
            f"{row['label']} |"
        )
    lines.extend([
        "",
        f"k=1 TRIVIAL: {sanity['k1_trivial']}. k=50 DESTROYED: {sanity['k50_destroyed']}.",
        "",
        "## Clock × target labels",
        "",
        "| k | Δt | Target | Best legal model | m | λ | Val NRMSE | Test NRMSE | Label |",
        "|---|---|---|---|---|---|---|---|---|",
    ])
    for row in pair_rows:
        lines.append(
            f"| {row['k']} | {row['dt_sample']} | {row['task']} | "
            f"{row['best_legal_model']} | {row['selected_m']} | {fmt(row['best_legal_lambda'])} | "
            f"{row['best_legal_val_nrmse']:.4f} | {row['best_legal_test_nrmse']:.4f} | "
            f"{row['label']} |"
        )
    lines.extend([
        "",
        "## Full legal and development table",
        "",
        "Every m is reported. Only the validation-selected delay enters the",
        "legal set with intercept and persistence. Oracle is not ranked.",
        "",
        "| k | Δt | Task | Model | m | Selected | λ | Train | Val | Test | R² | Legal |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ])
    for row in model_rows:
        lines.append(
            f"| {row['k']} | {row['dt_sample']} | {row['task']} | {row['model']} | "
            f"{fmt(row['development_m'])} | {row['selected_on_validation']} | "
            f"{fmt(row['lambda'])} | {row['train_nrmse']:.4f} | "
            f"{row['validation_nrmse']:.4f} | {row['test_nrmse']:.4f} | "
            f"{row['test_r2']:.4f} | {row['legal']} |"
        )
    lines.extend([
        "",
        "## Selection rule (frozen before new-k scores)",
        "",
        "Consider SURVIVE pairs only. Choose the pair whose best legal",
        "*validation* NRMSE is largest. At most one living BSim pair.",
        "Do not select on a dish score. Tie-break, frozen: larger k, then",
        "CROSS_Y over AUTO_X.",
        "",
        "Module 2, if authorized, uses the claim dish (CENTER / FLOW=0),",
        "75 s pulse / 300 s window, one AHL channel, acid held 0.5.",
        "Do not retune K, n, tau, rate, clamp, mortality, flow, or layout.",
        "Do not inject y or z. Do not use acid as a Lorenz channel.",
        "k=1 and k=50 are not rerun in BSim.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    hash_rows = load_hashes()
    print("HASHES_PRESENT_BEFORE_RANKING")
    for skip, _, name in CLOCKS:
        print(f"  k={skip} {name} sha256={hash_rows[skip]['u_sha256']}")

    model_rows = []
    pair_rows = []
    for skip, role, name in CLOCKS:
        ahl_path, target_path = clock_paths(skip)
        samples = integrate(skip)
        pack_u = samples[:NUM_WINDOWS, 0]
        xmin = float(np.min(pack_u))
        xmax = float(np.max(pack_u))
        u = 0.5 * (pack_u - xmin) / (xmax - xmin)
        stored_u = load_stored_u(ahl_path)
        digest = sha256_u(stored_u)
        expected_hash = hash_rows[skip]["u_sha256"]
        if digest != expected_hash:
            raise SystemExit(
                f"STOP: k={skip} file hash {digest} != recorded {expected_hash}"
            )
        if not np.allclose(stored_u, u, rtol=0, atol=1e-9):
            raise SystemExit(f"STOP: k={skip} stored u does not match affine RK4")
        x = samples[:NUM_WINDOWS, 0]
        y_now = samples[:NUM_WINDOWS, 1]
        x_next = samples[1:NUM_WINDOWS + 1, 0]
        y_next = samples[1:NUM_WINDOWS + 1, 1]
        dt_sample = skip * 0.02
        tasks = (
            ("AUTO_X", x_next, x, "LINEAR_AR"),
            ("CROSS_Y", y_next, y_now, "LINEAR_X_DELAY"),
        )
        for task, target, persist_series, delay_prefix in tasks:
            train_mean = float(np.mean(target[WASHOUT:WASHOUT + TRAIN]))
            intercept = constant_metrics(target, train_mean)
            intercept_row = {
                "k": skip,
                "dt_sample": f"{dt_sample:.2f}",
                "clock": name,
                "role": role,
                "task": task,
                "model": "TRAIN_INTERCEPT",
                "development_m": "",
                "selected_on_validation": "",
                "n_features": 0,
                "legal": "TRUE",
                **intercept,
            }
            model_rows.append(intercept_row)
            oracle = float(np.mean(target[WASHOUT + TRAIN:]))
            model_rows.append({
                "k": skip,
                "dt_sample": f"{dt_sample:.2f}",
                "clock": name,
                "role": role,
                "task": task,
                "model": "TEST_MEAN_ORACLE_REFERENCE",
                "development_m": "",
                "selected_on_validation": "",
                "n_features": 0,
                "legal": "FALSE",
                "note": "non-predictive oracle; uses test labels; not ranked",
                **constant_metrics(target, oracle),
            })
            persist_name = "PERSISTENCE_X_N" if task == "AUTO_X" else "PERSISTENCE_Y_N"
            persist = evaluate_nontrainable(target, persist_series)
            persist_row = {
                "k": skip,
                "dt_sample": f"{dt_sample:.2f}",
                "clock": name,
                "role": role,
                "task": task,
                "model": persist_name,
                "development_m": "",
                "selected_on_validation": "",
                "n_features": 1,
                "legal": "TRUE",
                **persist,
            }
            model_rows.append(persist_row)

            development = []
            for m in M_GRID:
                metrics = scalar_closed(delay_matrix(x, m), target)
                development.append((m, metrics))
                model_name = delay_prefix
                if task == "CROSS_Y" and m == 1:
                    model_name = "SCALAR_LINEAR_X"
                model_rows.append({
                    "k": skip,
                    "dt_sample": f"{dt_sample:.2f}",
                    "clock": name,
                    "role": role,
                    "task": task,
                    "model": model_name,
                    "development_m": m,
                    "selected_on_validation": "FALSE",
                    "n_features": m,
                    "legal": "FALSE",
                    **metrics,
                })
            selected_m, selected_metrics = min(
                development,
                key=lambda item: (item[1]["validation_nrmse"], item[0]),
            )
            selected_row = {
                "k": skip,
                "dt_sample": f"{dt_sample:.2f}",
                "clock": name,
                "role": role,
                "task": task,
                "model": f"{delay_prefix}_VALIDATION_SELECTED",
                "development_m": selected_m,
                "selected_on_validation": "TRUE",
                "n_features": selected_m,
                "legal": "TRUE",
                **selected_metrics,
            }
            model_rows.append(selected_row)

            legal = [
                ("TRAIN_INTERCEPT", intercept_row),
                (persist_name, persist_row),
                (selected_row["model"], selected_row),
            ]
            best_name, best_row = min(legal, key=lambda item: item[1]["test_nrmse"])
            best_val_name, best_val_row = min(
                legal, key=lambda item: item[1]["validation_nrmse"]
            )
            label = label_from_best_test(best_row["test_nrmse"])
            pair_rows.append({
                "k": skip,
                "dt_sample": f"{dt_sample:.2f}",
                "clock": name,
                "role": role,
                "task": task,
                "selected_m": selected_m,
                "selected_lambda": selected_metrics["lambda"],
                "selected_val_nrmse": selected_metrics["validation_nrmse"],
                "selected_test_nrmse": selected_metrics["test_nrmse"],
                "best_legal_model": best_name,
                "best_legal_lambda": best_row["lambda"],
                "best_legal_val_nrmse": best_val_row["validation_nrmse"],
                "best_legal_val_model": best_val_name,
                "best_legal_test_nrmse": best_row["test_nrmse"],
                "label": label,
            })

    sanity_rows = []
    sanity_ok = True
    for skip, task in ((1, "AUTO_X"), (1, "CROSS_Y"), (50, "AUTO_X"), (50, "CROSS_Y")):
        pair = next(row for row in pair_rows if row["k"] == skip and row["task"] == task)
        expected = ZERO_SIM_SELECTED_TEST[(skip, task)]
        got = pair["selected_test_nrmse"]
        delta = abs(got - expected)
        match = delta <= SANITY_TOL
        sanity_ok = sanity_ok and match
        sanity_rows.append({
            "k": skip,
            "task": task,
            "got": got,
            "expected": expected,
            "delta": delta,
            "match": "PASS" if match else "FAIL",
            "label": pair["label"],
        })
        print(
            f"SANITY k={skip} {task} selected_test={got:.6f} "
            f"zero_sim={expected:.4f} d={delta:.4g} {'PASS' if match else 'FAIL'}"
        )

    k1_labels = {row["task"]: row["label"] for row in pair_rows if row["k"] == 1}
    k50_labels = {row["task"]: row["label"] for row in pair_rows if row["k"] == 50}
    k1_trivial = all(label == "TRIVIAL" for label in k1_labels.values())
    k50_destroyed = all(label == "DESTROYED" for label in k50_labels.values())
    sanity = {
        "pass": sanity_ok,
        "rows": sanity_rows,
        "k1_trivial": k1_trivial,
        "k50_destroyed": k50_destroyed,
    }
    if not sanity_ok:
        results = HERE / "results"
        results.mkdir(parents=True, exist_ok=True)
        write_outputs(hash_rows, model_rows, pair_rows, sanity, "SANITY_FAIL_STOP", None)
        raise SystemExit(
            "STOP: k=1/k=50 selected-delay test NRMSE did not reproduce "
            "zero-sim within 1e-3. Thresholds were not changed."
        )
    if not k1_trivial or not k50_destroyed:
        write_outputs(hash_rows, model_rows, pair_rows, sanity, "SANITY_FAIL_STOP", None)
        raise SystemExit(
            "STOP: k=1 was not TRIVIAL or k=50 was not DESTROYED after "
            "zero-sim reproduction. Thresholds were not changed."
        )

    survivors = [
        row for row in pair_rows
        if row["role"] == "development" and row["label"] == "SURVIVE"
    ]
    selected = None
    if not survivors:
        disposition = "NO_INTERMEDIATE_CLOCK"
        print("L3_DISPOSITION: NO_INTERMEDIATE_CLOCK")
    else:
        selected = max(
            survivors,
            key=lambda row: (
                row["best_legal_val_nrmse"],
                row["k"],
                1 if row["task"] == "CROSS_Y" else 0,
            ),
        )
        disposition = (
            f"SCOUT_{selected['task']}_K{selected['k']}"
        )
        print(
            f"L3_DISPOSITION: {disposition} "
            f"best_legal_val={selected['best_legal_val_nrmse']:.4f} "
            f"best_legal_test={selected['best_legal_test_nrmse']:.4f}"
        )

    for row in pair_rows:
        print(
            f"LABEL k={row['k']} {row['task']} {row['label']} "
            f"best_test={row['best_legal_test_nrmse']:.4f} "
            f"best_val={row['best_legal_val_nrmse']:.4f} "
            f"model={row['best_legal_model']} m={row['selected_m']}"
        )

    write_outputs(hash_rows, model_rows, pair_rows, sanity, disposition, selected)


def write_outputs(hash_rows, model_rows, pair_rows, sanity, disposition, selected):
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    csv_path = results / "lorenz_l3_clock_screen.csv"
    fieldnames = [
        "k", "dt_sample", "clock", "role", "task", "model", "development_m",
        "selected_on_validation", "n_features", "lambda", "train_nrmse",
        "validation_nrmse", "test_nrmse", "test_r2", "legal", "note",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in model_rows:
            out = dict(row)
            out.setdefault("note", "")
            writer.writerow(out)
    pair_csv = results / "lorenz_l3_clock_labels.csv"
    with pair_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0].keys()))
        writer.writeheader()
        writer.writerows(pair_rows)
    md_path = results / "LORENZ_L3_CLOCK_SCREEN.md"
    write_markdown(md_path, hash_rows, model_rows, pair_rows, sanity, disposition, selected)
    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
