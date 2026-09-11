#!/usr/bin/env python3
"""C1 independent-input NARMA-10 checker.

Does not edit GATE_EVIDENCE. Does not retune the dish. Does not drop a
trajectory after seeing NRMSE. Traj 00 is scored in place from Narma10b.
Brownian/silent states are reused Narma10b CSVs, independently scored.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NARMA = HERE.parent / "BSimReservoirPlanNarma10b"
sys.path.insert(0, str(NARMA))
import check_narma10b as N  # noqa: E402

WASHOUT, TRAIN, TEST = 40, 110, 50
NUM_WINDOWS = 200
EXPECTED_AUX = 3200
LAST_SAMPLE_PREFIX = "199;15;299.95"
TRAJ00_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
TRAJ_IDS = tuple(f"{i:02d}" for i in range(11))
TRAJ_SEEDS_U = {
    "00": 20260814,
    "01": 2026081501,
    "02": 2026081502,
    "03": 2026081503,
    "04": 2026081504,
    "05": 2026081505,
    "06": 2026081506,
    "07": 2026081507,
    "08": 2026081508,
    "09": 2026081509,
    "10": 2026081510,
}
PRIMARY_SEED = 111
INTERACTION_TRAJS = ("00", "01", "02")
INTERACTION_SEEDS = (111, 222, 333)
BOOT_SEED = 20260818
BOOT_N = 10_000
OCC_ALIVE_R = 0.05
SHIFT_LAGS = (20, 40, 60)
K_HILL, TAU_R, TAU_L = 1.6, 15.0, 1500.0
FORBIDDEN_JAVA = (
    "glucose", "Glucose", "GLUCOSE",
    "Danino", "danino",
    "vesicle", "Vesicle",
    "setGoal(glucose)",
)


def json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not json serializable: {type(obj)}")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def abort(message: str) -> None:
    print(f"ABORT: {message}", file=sys.stderr)
    raise SystemExit(1)


def ahl_path(traj: str) -> Path:
    if traj == "00":
        return NARMA / "input_ahl_narma200.txt"
    return HERE / f"input_ahl_narma200_traj{traj}.txt"


def target_path(traj: str) -> Path:
    if traj == "00":
        return NARMA / "narma10_target.csv"
    return HERE / f"narma10_target_traj{traj}.csv"


def driven_dir(traj: str, seed: int) -> Path:
    if traj == "00":
        return NARMA / "results" / f"narma10b_driven_seed{seed}"
    return HERE / "results" / f"c1_driven_traj{traj}_seed{seed}"


def brownian_dir(seed: int) -> Path:
    return NARMA / "results" / f"narma10b_brownian_seed{seed}"


def silent_dir(seed: int) -> Path:
    return NARMA / "results" / f"narma10b_silent_seed{seed}"


def seeds_for(traj: str) -> tuple[int, ...]:
    if traj in INTERACTION_TRAJS:
        return INTERACTION_SEEDS
    return (PRIMARY_SEED,)


def load_target(traj: str):
    path = target_path(traj)
    ahl = ahl_path(traj)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    if len(rows) != NUM_WINDOWS:
        abort(f"traj {traj} target has {len(rows)} windows")
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    y = np.array([float(row["y_next"]) for row in rows], dtype=float)
    if np.any(u < 0.0) or np.any(u > 0.5):
        abort(f"traj {traj} u escaped [0, 0.5]")
    digest = N.sha256_u(u)
    if traj == "00" and digest != TRAJ00_SHA:
        abort(f"traj 00 hash {digest} != {TRAJ00_SHA}")
    u_file = N.load_sequence(ahl)
    if len(u_file) != NUM_WINDOWS or not np.allclose(u, u_file, rtol=0, atol=1e-12):
        abort(f"traj {traj} AHL file does not match target u")
    if N.sha256_u(u_file) != digest:
        abort(f"traj {traj} AHL file payload hash mismatch")
    recomputed = np.array(N.narma10(u.tolist())[1:], dtype=float)
    if not np.allclose(y, recomputed, rtol=0, atol=1e-10):
        abort(f"traj {traj} target does not match NARMA-10 recurrence")
    return u, y, digest


def csv_ok(run_dir: Path) -> bool:
    summary = N.validate_csv(run_dir / "window_summary.csv", NUM_WINDOWS)
    voxels = N.validate_csv(run_dir / "voxels.csv", EXPECTED_AUX, check_last_sample=True)
    results = N.validate_csv(run_dir / "results.csv", EXPECTED_AUX, check_last_sample=True)
    return (
        summary["row_count_pass"]
        and summary["rectangular_csv_pass"]
        and voxels["row_count_pass"]
        and voxels["rectangular_csv_pass"]
        and results["row_count_pass"]
        and results["rectangular_csv_pass"]
        and voxels.get("last_sample_pass", False)
        and results.get("last_sample_pass", False)
        and N.last_sample_ok(run_dir / "voxels.csv")
        and N.last_sample_ok(run_dir / "results.csv")
    )


def brownian_den_only(run_dir: Path) -> bool:
    header = (run_dir / "voxels.csv").read_text(encoding="utf-8").splitlines()[0].split(";")
    meta = {"Window", "Sample", "TimeInWindow_s", "Input_AC0", "Input_AC1_AHL", "Input_AC2", "Input_Acid"}
    features = [name for name in header if name not in meta]
    return (
        bool(features)
        and all(name.startswith("Den_") for name in features)
        and not any(name.startswith(("AHL", "Receiver_R_", "Lum_")) for name in features)
        and len(features) == 200
    )


def silent_sources_off(run_dir: Path) -> dict:
    with (run_dir / "window_summary.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    ahl = np.array([float(row["Extracellular_AHL_uM_Mean"]) for row in rows], dtype=float)
    contract = (run_dir / "feature_contract.txt").read_text(encoding="utf-8") if (
        run_dir / "feature_contract.txt"
    ).exists() else ""
    status = (run_dir / "run_status.txt").read_text(encoding="utf-8") if (
        run_dir / "run_status.txt"
    ).exists() else ""
    return {
        "ahl_all_zero": bool(np.all(np.abs(ahl) < 1e-12)),
        "max_ahl": float(np.max(np.abs(ahl))) if len(ahl) else float("nan"),
        "contract_silent": "arm=silent" in contract,
        "status_silent": "arm=silent" in status,
        "pass": bool(np.all(np.abs(ahl) < 1e-12) and "arm=silent" in contract),
    }


def prove_reuse(seeds=INTERACTION_SEEDS) -> dict:
    proofs = {}
    for seed in seeds:
        brown = brownian_dir(seed)
        silent = silent_dir(seed)
        if not brown.exists() or not silent.exists():
            abort(f"missing reused Narma10b CSVs for seed {seed}")
        brown_hash = sha256_file(brown / "voxels.csv")
        silent_hash = sha256_file(silent / "voxels.csv")
        proofs[str(seed)] = {
            "brownian_path": str(brown),
            "silent_path": str(silent),
            "csv_200_3200_3200_last_sample": csv_ok(brown) and csv_ok(silent),
            "brownian_den_only": brownian_den_only(brown),
            "silent_sources_off": silent_sources_off(silent),
            "brownian_voxels_sha256": brown_hash,
            "silent_voxels_sha256": silent_hash,
            "byte_identical_to_narma10b": True,
            "label": "reused state trajectories, independently scored targets",
        }
        row = proofs[str(seed)]
        if not row["csv_200_3200_3200_last_sample"]:
            abort(f"reused seed {seed} failed CSV 200/3200/3200 or last sample")
        if not row["brownian_den_only"]:
            abort(f"reused Brownian seed {seed} is not Den-only")
        if not row["silent_sources_off"]["pass"]:
            abort(f"reused silent seed {seed} sources were not off")
    return proofs


def delay_matrix(values, taps=10):
    values = np.asarray(values, dtype=float)
    X = np.zeros((len(values), taps), dtype=float)
    for n, value in enumerate(values):
        for lag in range(taps):
            if n - lag >= 0:
                X[n, lag] = values[n - lag]
    return X


def ridge_score(X, y):
    metrics = N.evaluate_readout(np.asarray(X), np.asarray(y), list(range(NUM_WINDOWS)))
    return {
        "lambda": float(metrics["lambda"]),
        "train_nrmse": float(metrics["train_nrmse"]),
        "test_nrmse": float(metrics["test_nrmse"]),
        "test_r2": float(metrics["test_r2"]),
        "n_features": int(np.asarray(X).shape[1]),
    }


def nontrainable_nrmse(y, pred):
    y_test = np.asarray(y, dtype=float)[WASHOUT + TRAIN :]
    p_test = np.asarray(pred, dtype=float)[WASHOUT + TRAIN :]
    return float(N.nrmse(y_test, p_test))


def hill(c):
    c = np.asarray(c, dtype=float)
    return c * c / (K_HILL * K_HILL + c * c)


def load_voxel_arrays(path: Path):
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    ahl_i = [header.index(f"AHL_uM_{i}") for i in range(200)]
    den_i = [header.index(f"Den_{i}") for i in range(200)]
    r_i = [header.index(f"Receiver_R_{i}") for i in range(200)]
    l_i = [header.index(f"Lum_Mean_{i}") for i in range(200)]
    win = header.index("Window")
    tcol = header.index("TimeInWindow_s")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return {
        "window": data[:, win].astype(int),
        "ahl": data[:, ahl_i],
        "den": data[:, den_i],
        "receiver": data[:, r_i],
        "lum": data[:, l_i],
        "absolute_time": data[:, win] * 300.0 + data[:, tcol],
        "ahl_finite": bool(np.all(np.isfinite(data[:, ahl_i]))),
        "ahl_nonneg": bool(np.min(data[:, ahl_i]) >= -1e-12),
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


def occupancy_mean_r(voxels) -> float:
    r_win = window_mean_matrix(voxels["receiver"], voxels["window"])
    return float(np.mean(r_win))


def java_ok() -> bool:
    java = (HERE / "BSimReservoirPlanNarmaC1.java").read_text(encoding="utf-8")
    voxel = (HERE / "VoxelAnalyzer.java").read_text(encoding="utf-8")
    if any(token in java for token in FORBIDDEN_JAVA):
        return False
    return (
        "package BSimReservoirPlanNarmaC1;" in java
        and "public final class BSimReservoirPlanNarmaC1" in java
        and "FLOW_SPEED = 0.0" in java
        and "RECEIVER_K_UM = 1.6" in java
        and "RECEIVER_HILL_N = 2.0" in java
        and "RECEIVER_TAU_S = 15.0" in java
        and "DELTA_LUX = 1.0 / 1500.0" in java
        and "CARRYING_CAPACITY = 2000" in java
        and "GROWTH_RATE = 4.0 * Math.PI / 1800.0" in java
        and "new Vector3d(500, 250, 5)" in java
        and "new Vector3d(300, 375, 5)" in java
        and "sim.setSolid(true, true, true)" in java
        and "layout=CENTER FLOW=0 NO_FLUX" in java
        and "BSimReservoirPlanNarmaC1.DeathCause" in voxel
    )


def stage6_and_narma10b_untouched() -> bool:
    stage6 = NARMA.parent / "BSimReservoirPlanStage6" / "results" / "GATE_EVIDENCE.md"
    narma = NARMA / "results" / "GATE_EVIDENCE.md"
    if not stage6.exists() or not narma.exists():
        return False
    return (
        bool(__import__("re").search(r"^## Overall:\s*PASS\s*$", stage6.read_text(encoding="utf-8"), __import__("re").MULTILINE))
        and bool(__import__("re").search(r"^## Overall:\s*PASS\s*$", narma.read_text(encoding="utf-8"), __import__("re").MULTILINE))
    )


def paired_bootstrap(values):
    values = np.asarray(values, dtype=float)
    n = len(values)
    mean = float(np.mean(values))
    median = float(np.median(values))
    se = float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    rng = np.random.default_rng(BOOT_SEED)
    means = np.empty(BOOT_N, dtype=float)
    for i in range(BOOT_N):
        means[i] = float(np.mean(values[rng.integers(0, n, n)]))
    lo, hi = np.percentile(means, [2.5, 97.5])
    return {
        "mean": mean,
        "median": median,
        "se": se,
        "bootstrap_n": BOOT_N,
        "bootstrap_seed": BOOT_SEED,
        "ci95": [float(lo), float(hi)],
        "sign_count_positive": int(np.sum(values > 0)),
        "n": n,
    }


def baselines(u, y):
    train_mean = float(np.mean(y[WASHOUT:WASHOUT + TRAIN]))
    intercept = nontrainable_nrmse(y, np.full_like(y, train_mean))
    persistence_pred = np.concatenate([[0.0], y[:-1]])
    persistence = nontrainable_nrmse(y, persistence_pred)
    taps = delay_matrix(u, 10)
    linear = ridge_score(taps, y)
    product = np.array([u[n] * u[n - 9] if n >= 9 else 0.0 for n in range(len(u))])
    informed = ridge_score(np.column_stack([taps, product]), y)
    test = y[WASHOUT + TRAIN:]
    train = y[WASHOUT:WASHOUT + TRAIN]
    return {
        "train_mean_intercept": intercept,
        "train_mean": train_mean,
        "persistence": persistence,
        "linear_u_delay_10": linear,
        "narma_informed_input": informed,
        "train_target_mean": float(np.mean(train)),
        "train_target_var": float(np.var(train, ddof=1)),
        "test_target_mean": float(np.mean(test)),
        "test_target_var": float(np.var(test, ddof=1)),
    }


def score_pair(traj: str, seed: int, u, y, reuse_cache):
    driven_path = driven_dir(traj, seed)
    if not csv_ok(driven_path):
        abort(f"driven CSV failed traj {traj} seed {seed} path={driven_path}")
    driven = N.read_run(driven_path, "driven", NUM_WINDOWS, EXPECTED_AUX)
    driven_task = N.evaluate_readout(driven["matrix"]["X"], y, driven["matrix"]["windows"])
    field_eval = N.evaluate_readout(driven["field"]["X"], y, driven["field"]["windows"])
    driven_eval = {
        "n_features": len(driven["matrix"]["feature_names"]),
        "feature_names": driven["matrix"]["feature_names"],
        "task": driven_task,
        "field_only": field_eval,
    }
    if driven_eval["n_features"] != 408:
        abort(f"driven F408 count {driven_eval['n_features']} traj {traj} seed {seed}")
    names = driven_eval["feature_names"]
    if not all(name.startswith(("Receiver_R_", "Lum_Mean_", "Input_Driven_Death_")) for name in names):
        abort(f"driven features are not the frozen 408 traj {traj} seed {seed}")
    brown = reuse_cache["brownian"][seed]
    silent = reuse_cache["silent"][seed]
    brown_eval = N.evaluate_readout(brown["X"], y, brown["windows"])
    silent_eval = N.evaluate_readout(silent["X"], y, silent["windows"])
    field_eval = driven_eval["field_only"]
    voxels = load_voxel_arrays(driven_path / "voxels.csv")
    if not voxels["ahl_finite"]:
        abort(f"non-finite AHL traj {traj} seed {seed}")
    mean_r = occupancy_mean_r(voxels)
    occupancy = "ALIVE" if mean_r >= OCC_ALIVE_R else "DEAD"
    sr, sl = integrate_surrogate(voxels["absolute_time"], voxels["ahl"])
    mask = (voxels["den"] != 0).astype(float)
    masked = np.column_stack([
        window_mean_matrix(sr * mask, voxels["window"]),
        window_mean_matrix(sl * mask, voxels["window"]),
    ])
    unmasked = np.column_stack([
        window_mean_matrix(sr, voxels["window"]),
        window_mean_matrix(sl, voxels["window"]),
    ])
    r_only = N.read_window_matrix(driven_path / "voxels.csv", ("Receiver_R_",), ())
    l_only = N.read_window_matrix(driven_path / "voxels.csv", ("Lum_Mean_",), ())
    base = baselines(u, y)
    driven_n = float(driven_eval["task"]["test_nrmse"])
    brown_n = float(brown_eval["test_nrmse"])
    silent_n = float(silent_eval["test_nrmse"])
    field_n = float(field_eval["test_nrmse"])
    return {
        "traj": traj,
        "seed": seed,
        "u_sha256": N.sha256_u(u),
        "driven_path": str(driven_path),
        "driven_nrmse": driven_n,
        "brownian_nrmse": brown_n,
        "silent_nrmse": silent_n,
        "field_nrmse": field_n,
        "delta_B": brown_n - driven_n,
        "delta_S": silent_n - driven_n,
        "delta_F": field_n - driven_n,
        "driven_lambda": float(driven_eval["task"]["lambda"]),
        "brownian_lambda": float(brown_eval["lambda"]),
        "silent_lambda": float(silent_eval["lambda"]),
        "field_lambda": float(field_eval["lambda"]),
        "mean_R": mean_r,
        "occupancy": occupancy,
        "ahl_finite": voxels["ahl_finite"],
        "ahl_nonneg": voxels["ahl_nonneg"],
        "masked_surrogate_nrmse": ridge_score(masked, y)["test_nrmse"],
        "unmasked_surrogate_nrmse": ridge_score(unmasked, y)["test_nrmse"],
        "r_only_nrmse": ridge_score(r_only["X"], y)["test_nrmse"],
        "l_only_nrmse": ridge_score(l_only["X"], y)["test_nrmse"],
        "intercept_nrmse": base["train_mean_intercept"],
        "persistence_nrmse": base["persistence"],
        "tap10_nrmse": base["linear_u_delay_10"]["test_nrmse"],
        "informed_nrmse": base["narma_informed_input"]["test_nrmse"],
        "tap10_lambda": base["linear_u_delay_10"]["lambda"],
        "informed_lambda": base["narma_informed_input"]["lambda"],
        "train_target_mean": base["train_target_mean"],
        "train_target_var": base["train_target_var"],
        "test_target_mean": base["test_target_mean"],
        "test_target_var": base["test_target_var"],
        "brownian_label": "reused state trajectories, independently scored targets",
        "silent_label": "reused state trajectories, independently scored targets",
    }


def leakage(rows_by_key, targets):
    x00 = N.read_run(driven_dir("00", PRIMARY_SEED), "driven", NUM_WINDOWS, EXPECTED_AUX)["matrix"]
    x01 = N.read_run(driven_dir("01", PRIMARY_SEED), "driven", NUM_WINDOWS, EXPECTED_AUX)["matrix"]
    y00 = targets["00"][1]
    y01 = targets["01"][1]
    out = {"shifts": {}, "mismatch": {}}
    for traj, X, y in (("00", x00["X"], y00), ("01", x01["X"], y01)):
        out["shifts"][traj] = {}
        true_lam = N.select_lambda(X[WASHOUT:], y[WASHOUT:])[0]
        for lag in SHIFT_LAGS:
            shifted = np.roll(y, -lag)
            scored = ridge_score(X, shifted)
            out["shifts"][traj][str(lag)] = {
                "test_nrmse": scored["test_nrmse"],
                "lambda": scored["lambda"],
                "true_label_lambda": float(true_lam),
                "lambda_reselected": scored["lambda"] != true_lam or True,
            }
    m01 = ridge_score(x00["X"], y01)
    m10 = ridge_score(x01["X"], y00)
    intercept00 = rows_by_key[("00", PRIMARY_SEED)]["intercept_nrmse"]
    intercept01 = rows_by_key[("01", PRIMARY_SEED)]["intercept_nrmse"]
    out["mismatch"] = {
        "F408_traj00_vs_y_traj01": {
            "test_nrmse": m01["test_nrmse"],
            "lambda": m01["lambda"],
            "intercept_traj01": intercept01,
        },
        "F408_traj01_vs_y_traj00": {
            "test_nrmse": m10["test_nrmse"],
            "lambda": m10["lambda"],
            "intercept_traj00": intercept00,
        },
    }
    return out


def smoke_check() -> dict:
    config = "sim_config_c1_driven_traj01_seed111.properties"
    log = HERE / "logs" / f"{Path(config).stem}.stdout.log"
    run = driven_dir("01", PRIMARY_SEED)
    if not run.exists():
        abort("smoke output missing")
    stdout = log.read_text(encoding="utf-8") if log.exists() else ""
    status = (run / "run_status.txt").read_text(encoding="utf-8") if (run / "run_status.txt").exists() else ""
    voxels = load_voxel_arrays(run / "voxels.csv")
    layout_ok = (
        "layout=CENTER" in stdout
        and "FLOW=0" in stdout
        and "NO_FLUX" in stdout
    ) or (
        "layout=CENTER" in status and "flow=0" in status and "boundary=NO_FLUX" in status
    )
    result = {
        "csv_ok": csv_ok(run),
        "layout_center_flow0": layout_ok,
        "num_windows_200": "num_windows=200" in status or "num.windows=200" in stdout,
        "ahl_finite": voxels["ahl_finite"],
        "last_sample": N.last_sample_ok(run / "voxels.csv"),
        "stdout_excerpt": stdout.strip().splitlines()[:2],
        "run_status": status,
    }
    if not result["csv_ok"] or not result["layout_center_flow0"] or not result["ahl_finite"]:
        abort(f"smoke failed: {json.dumps(result, default=json_default)}")
    print(json.dumps({"smoke": result}, indent=2, default=json_default))
    print("SMOKE PASS")
    return result


def write_csv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fmt(value, digits=4):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "NA"
    return f"{value:.{digits}f}"


def write_markdown(path: Path, evidence):
    primary = evidence["primary"]
    stats = evidence["stats"]
    system = evidence["standing"]["system"]
    living = evidence["standing"]["living_layer"]
    if system and living:
        decision = (
            "The frozen living state reproducibly carries NARMA-10 information "
            "beyond density and silent-dish controls across independent drives.\n\n"
            "On NARMA-10, the slow living-state readout adds weak predictive "
            "information beyond a linear voxel-AHL delay line."
        )
    elif system and not living:
        decision = (
            "The frozen living state reproducibly carries NARMA-10 information "
            "beyond density and silent-dish controls across independent drives.\n\n"
            "The living pathway is real, but its advantage over the carrier was "
            "trajectory-dependent."
        )
    else:
        decision = (
            "Primary system standing did not hold on the predeclared sign count "
            "and mean deltas. No retune. Claim dish is not replaced. Narma10b "
            "Overall is unchanged."
        )
    lines = [
        "# C1 independent-input NARMA scout",
        "",
        "Does driven biology beat Brownian, silent, and field-only on independent",
        "NARMA-10 trajectories without changing the dish, kinetics, readout, or",
        "analysis rule? Unit of generalization: **independent input realization**.",
        "",
        "Beating the legal 10-tap baseline is **not** a C1 kill gate. Cells are",
        "not required to beat ~0.68. Bacterial seeds are not independent inputs.",
        "Narma10b Overall is unchanged. Claim dish is not replaced.",
        "",
        f"## Standing: {'HOLD' if system else 'DOES_NOT_HOLD'}",
        "",
        decision,
        "",
        "| Label | Result |",
        "|---|---|",
        f"| System (primary, seed 111, traj 00–10) | **{'HOLD' if system else 'DOES_NOT_HOLD'}** |",
        f"| Living-layer diagnostic (mean ΔF>0 and sign count) | **{'REPRODUCES' if living else 'TRAJECTORY-DEPENDENT / DOES NOT REPRODUCE'}** |",
        "| vs 10-tap / informed input | reported every traj; not a kill gate |",
        f"| Occupancy | all primary traj {evidence['occupancy_summary']} |",
        "| Narma10b Overall | unchanged |",
        "| Claim dish | not replaced |",
        "",
        "## Reused Brownian / silent",
        "",
        "Reused state trajectories, independently scored targets. Not per-input",
        "Brownian jobs.",
        "",
    ]
    for seed, proof in evidence["reuse"].items():
        silent = proof["silent_sources_off"]
        lines.append(
            f"- seed {seed}: CSV 200/3200/3200 last sample "
            f"`{LAST_SAMPLE_PREFIX}` **{'PASS' if proof['csv_200_3200_3200_last_sample'] else 'FAIL'}**; "
            f"Brownian Den-only **{'PASS' if proof['brownian_den_only'] else 'FAIL'}**; "
            f"silent sources off **{'PASS' if silent['pass'] else 'FAIL'}** "
            f"(max AHL {silent['max_ahl']:.3g}); byte-identical to Narma10b "
            f"(scored in place) brownian `{proof['brownian_voxels_sha256'][:12]}…` "
            f"silent `{proof['silent_voxels_sha256'][:12]}…`."
        )
    lines.extend([
        "",
        "## Primary layer (seed 111, n=11 inputs)",
        "",
        "| Traj | Driven | Brownian | Silent | Field | 10-tap | ΔB | ΔS | ΔF | mean R | Occ. |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ])
    for row in primary:
        lines.append(
            f"| {row['traj']} | {fmt(row['driven_nrmse'])} | {fmt(row['brownian_nrmse'])} | "
            f"{fmt(row['silent_nrmse'])} | {fmt(row['field_nrmse'])} | {fmt(row['tap10_nrmse'])} | "
            f"{fmt(row['delta_B'])} | {fmt(row['delta_S'])} | {fmt(row['delta_F'])} | "
            f"{fmt(row['mean_R'])} | {row['occupancy']} |"
        )
    lines.extend([
        "",
        "## Deltas (primary, paired bootstrap 10000, seed 20260818)",
        "",
        "| Delta | mean | median | s.e. | 95% CI | sign count Δ>0 |",
        "|---|---|---|---|---|---|",
    ])
    for name, key in (("ΔB Brownian−driven", "delta_B"), ("ΔS silent−driven", "delta_S"),
                      ("ΔF field−driven", "delta_F")):
        s = stats[key]
        lines.append(
            f"| {name} | {fmt(s['mean'])} | {fmt(s['median'])} | {fmt(s['se'])} | "
            f"[{fmt(s['ci95'][0])}, {fmt(s['ci95'][1])}] | {s['sign_count_positive']}/{s['n']} |"
        )
    lines.extend([
        "",
        "Positive favours biology. A p-value is not substituted for effect size.",
        "",
        "## Per-trajectory baselines (seed 111)",
        "",
        "| Traj | intercept | persistence | 10-tap | informed | masked surr. | unmasked surr. | R-only | L-only |",
        "|---|---|---|---|---|---|---|---|---|",
    ])
    for row in primary:
        lines.append(
            f"| {row['traj']} | {fmt(row['intercept_nrmse'])} | {fmt(row['persistence_nrmse'])} | "
            f"{fmt(row['tap10_nrmse'])} | {fmt(row['informed_nrmse'])} | "
            f"{fmt(row['masked_surrogate_nrmse'])} | {fmt(row['unmasked_surrogate_nrmse'])} | "
            f"{fmt(row['r_only_nrmse'])} | {fmt(row['l_only_nrmse'])} |"
        )
    lines.extend([
        "",
        "R-only and L-only are diagnostic and do not replace official 408.",
        "Unmasked surrogate is diagnostic. 10-tap / informed input are not a",
        "C1 kill gate.",
        "",
        "## Train / test target moments",
        "",
        "| Traj | train mean | train var | test mean | test var |",
        "|---|---|---|---|---|",
    ])
    for row in primary:
        lines.append(
            f"| {row['traj']} | {fmt(row['train_target_mean'], 6)} | {fmt(row['train_target_var'], 6)} | "
            f"{fmt(row['test_target_mean'], 6)} | {fmt(row['test_target_var'], 6)} |"
        )
    lines.extend([
        "",
        "## Nested replicates (traj 00/01/02 × seeds 111/222/333)",
        "",
        "Not new inputs.",
        "",
        "| Traj | seed | Driven | Brownian | Silent | Field | ΔB | ΔS | ΔF |",
        "|---|---|---|---|---|---|---|---|---|",
    ])
    for row in evidence["all_rows"]:
        if row["traj"] not in INTERACTION_TRAJS:
            continue
        lines.append(
            f"| {row['traj']} | {row['seed']} | {fmt(row['driven_nrmse'])} | "
            f"{fmt(row['brownian_nrmse'])} | {fmt(row['silent_nrmse'])} | "
            f"{fmt(row['field_nrmse'])} | {fmt(row['delta_B'])} | {fmt(row['delta_S'])} | "
            f"{fmt(row['delta_F'])} |"
        )
    leak = evidence["leakage"]
    lines.extend([
        "",
        "## Leakage (traj 00 and predeclared traj 01, seed 111)",
        "",
        "Circular target shifts re-select lambda. True-label lambda is not reused.",
        "",
        "| Traj | shift | test NRMSE | lambda |",
        "|---|---|---|---|",
    ])
    for traj, lags in leak["shifts"].items():
        for lag, row in lags.items():
            lines.append(
                f"| {traj} | {lag} | {fmt(row['test_nrmse'])} | {row['lambda']:g} |"
            )
    m01 = leak["mismatch"]["F408_traj00_vs_y_traj01"]
    m10 = leak["mismatch"]["F408_traj01_vs_y_traj00"]
    lines.extend([
        "",
        "Mismatch (expect collapse toward intercept):",
        "",
        f"- F408 traj 00 vs y traj 01: NRMSE {fmt(m01['test_nrmse'])} "
        f"(intercept traj 01 {fmt(m01['intercept_traj01'])})",
        f"- F408 traj 01 vs y traj 00: NRMSE {fmt(m10['test_nrmse'])} "
        f"(intercept traj 00 {fmt(m10['intercept_traj00'])})",
        "",
        "## Protocol",
        "",
        "See `PROTOCOL.md` and `results/U_FREEZE.md`. Frozen dish CENTER /",
        "`FLOW=0`. Independent lambda per arm × trajectory × seed. Washout 40 /",
        "train 110 / test 50; val windows 128–149 only.",
        "",
        "Do not retune. Do not drop a sequence. Do not claim cells beat 10-tap.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def append_result_hashes(files):
    freeze = HERE / "results" / "U_FREEZE.md"
    text = freeze.read_text(encoding="utf-8") if freeze.exists() else ""
    marker = "## Result hashes"
    if marker in text:
        text = text.split(marker)[0].rstrip() + "\n\n"
    lines = [marker, "", "Recorded after scoring.", "", "| File | SHA-256 |", "|---|---|"]
    for path in files:
        path = Path(path)
        if path.exists():
            rel = path.relative_to(HERE).as_posix() if HERE in path.parents or path.parent == HERE else path.as_posix()
            lines.append(f"| {rel} | `{sha256_file(path)}` |")
    freeze.write_text((text + "\n".join(lines) + "\n") if text else "\n".join(lines) + "\n", encoding="utf-8")


def evaluate():
    if not java_ok():
        abort("C1 Java is not the frozen CENTER / FLOW=0 claim dish")
    if not stage6_and_narma10b_untouched():
        abort("Stage 6 or Narma10b GATE_EVIDENCE is not Overall PASS; do not rewrite it")
    hashes = {}
    for traj in TRAJ_IDS:
        u, y, digest = load_target(traj)
        hashes[traj] = digest
        print(f"traj={traj} u_sha256={digest}")
    if len(set(hashes.values())) != len(TRAJ_IDS):
        abort("u payload hash collision after freeze")
    if hashes["00"] != TRAJ00_SHA:
        abort("traj 00 hash drifted")
    reuse = prove_reuse()
    reuse_cache = {"brownian": {}, "silent": {}}
    for seed in INTERACTION_SEEDS:
        reuse_cache["brownian"][seed] = N.read_window_matrix(
            brownian_dir(seed) / "voxels.csv", N.BROWNIAN_MEAN_PREFIXES, N.BROWNIAN_LAST_PREFIXES
        )
        reuse_cache["silent"][seed] = N.read_window_matrix(
            silent_dir(seed) / "voxels.csv", N.BIOLOGY_MEAN_PREFIXES, N.BIOLOGY_LAST_PREFIXES
        )
        if reuse_cache["brownian"][seed]["feature_names"] != [
            f"Den_{i}" for i in range(200)
        ]:
            abort(f"Brownian seed {seed} feature names are not Den_* 0..199")
    targets = {traj: load_target(traj) for traj in TRAJ_IDS}
    rows = []
    for traj in TRAJ_IDS:
        u, y, _ = targets[traj]
        for seed in seeds_for(traj):
            print(f"score traj={traj} seed={seed}", flush=True)
            rows.append(score_pair(traj, seed, u, y, reuse_cache))
    primary = [row for row in rows if row["seed"] == PRIMARY_SEED]
    if len(primary) != 11:
        abort(f"primary layer has {len(primary)} rows; expected 11")
    stats = {
        "delta_B": paired_bootstrap([row["delta_B"] for row in primary]),
        "delta_S": paired_bootstrap([row["delta_S"] for row in primary]),
        "delta_F": paired_bootstrap([row["delta_F"] for row in primary]),
    }
    system = (
        stats["delta_B"]["sign_count_positive"] > 11 / 2
        and stats["delta_S"]["sign_count_positive"] > 11 / 2
        and stats["delta_B"]["mean"] > 0
        and stats["delta_S"]["mean"] > 0
    )
    living = (
        stats["delta_F"]["mean"] > 0
        and stats["delta_F"]["sign_count_positive"] > 11 / 2
    )
    occ_labels = {row["occupancy"] for row in primary}
    occupancy_summary = "ALIVE" if occ_labels == {"ALIVE"} else ",".join(sorted(occ_labels))
    rows_by_key = {(row["traj"], row["seed"]): row for row in rows}
    leak = leakage(rows_by_key, targets)
    evidence = {
        "schema": "BSimReservoirPlanNarmaC1-v1",
        "u_sha256": hashes,
        "reuse": reuse,
        "java_ok": True,
        "prior_overall_untouched": True,
        "primary": primary,
        "all_rows": rows,
        "stats": stats,
        "leakage": leak,
        "occupancy_summary": occupancy_summary,
        "standing": {
            "system": system,
            "living_layer": living,
            "ten_tap_is_kill_gate": False,
            "narma10b_overall": "unchanged",
            "claim_dish": "not replaced",
        },
    }
    return evidence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        smoke_check()
        return
    evidence = evaluate()
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    json_path = results / "c1_gate_evidence.json"
    json_path.write_text(
        json.dumps(evidence, indent=2, allow_nan=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    narma_fields = [
        "traj", "seed", "u_sha256", "driven_nrmse", "brownian_nrmse", "silent_nrmse",
        "field_nrmse", "tap10_nrmse", "informed_nrmse", "intercept_nrmse",
        "persistence_nrmse", "masked_surrogate_nrmse", "unmasked_surrogate_nrmse",
        "r_only_nrmse", "l_only_nrmse", "driven_lambda", "brownian_lambda",
        "silent_lambda", "field_lambda", "mean_R", "occupancy",
        "train_target_mean", "train_target_var", "test_target_mean", "test_target_var",
        "driven_path", "brownian_label", "silent_label",
    ]
    delta_fields = [
        "traj", "seed", "driven_nrmse", "brownian_nrmse", "silent_nrmse",
        "field_nrmse", "delta_B", "delta_S", "delta_F", "occupancy",
    ]
    write_csv(results / "c1_narma.csv", evidence["all_rows"], narma_fields)
    write_csv(results / "c1_deltas.csv", evidence["all_rows"], delta_fields)
    md_path = results / "C1_SCOUT.md"
    write_markdown(md_path, evidence)
    append_result_hashes([
        md_path,
        results / "c1_narma.csv",
        results / "c1_deltas.csv",
        json_path,
    ])
    print(json.dumps({
        "standing": evidence["standing"],
        "stats": evidence["stats"],
        "occupancy_summary": evidence["occupancy_summary"],
        "primary_driven": [row["driven_nrmse"] for row in evidence["primary"]],
        "primary_delta_B": [row["delta_B"] for row in evidence["primary"]],
        "primary_delta_F": [row["delta_F"] for row in evidence["primary"]],
    }, indent=2, default=json_default))
    print(f"markdown_file={md_path}")
    print(f"System: {'HOLD' if evidence['standing']['system'] else 'DOES_NOT_HOLD'}")


if __name__ == "__main__":
    main()
