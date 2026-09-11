#!/usr/bin/env python3
"""Lorenz L3 living scout checker. Selected pair is frozen: k=10 AUTO_X.

Does not rewrite BenchA / BenchA2 GATE_EVIDENCE. Does not retune K, n,
tau, rate, clamp, mortality, flow, or layout. Does not switch k or
target after NRMSE. Seeds 202/303 are authorized only after seed 101
beats the best legal baseline, Brownian, and silent, and is not
occupancy-DEAD.
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

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
STAGE6 = EXAMPLES / "BSimReservoirPlanStage6"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EXAMPLES / "BSimReservoirPlanBenchA"))
import check_bencha as BC  # noqa: E402
from generate_lorenz_l3 import integrate as integrate_lorenz  # noqa: E402
from screen_lorenz_l3 import (  # noqa: E402
    delay_matrix,
    scalar_closed,
    constant_metrics,
    evaluate_nontrainable,
    M_GRID,
)

WASHOUT, TRAIN, TEST, INNER_VAL = 40, 110, 50, 22
NUM_WINDOWS = 200
EXPECTED_AUX = 3200
LAST_SAMPLE_PREFIX = "199;15;299.95"
SELECTED_K = 10
SELECTED_TASK = "AUTO_X"
SELECTED_SHA = "69b861f256fb39527c037b994885d119ccc8ab2873a3d2bca59f7a8802c44aff"
SCOUT_SEED = 101
CONFIRM_SEEDS = (202, 303)
K_HILL = 1.6
TAU_R = 15.0
TAU_L = 1500.0
OCC_DEAD_R = 0.05
OCC_DEAD_CORR = 0.5
WEAK_BAND = 0.90
WEAK_MARGIN = 0.03
RX, RY = 20, 10

SILENT_DIRS = {
    101: STAGE6 / "results" / "stage6_silent_seed101",
    202: STAGE6 / "results" / "stage6_silent_seed202",
    303: STAGE6 / "results" / "stage6_silent_seed303",
}

FORBIDDEN_JAVA = (
    "glucose", "Glucose", "GLUCOSE",
    "Danino", "danino",
    "vesicle", "Vesicle",
    "setGoal(glucose)",
)


def sha256_u(u):
    return hashlib.sha256(",".join(f"{v:.12f}" for v in u).encode("ascii")).hexdigest()


def load_selected_target():
    path = HERE / "lorenz_target_k10.csv"
    ahl = HERE / "input_ahl_lorenz_k10.txt"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    if len(rows) != NUM_WINDOWS:
        raise SystemExit(f"target has {len(rows)} rows")
    u = np.array([float(row["u"]) for row in rows], dtype=float)
    x = np.array([float(row["x"]) for row in rows], dtype=float)
    x_next = np.array([float(row["x_next"]) for row in rows], dtype=float)
    y_next = np.array([float(row["y_next"]) for row in rows], dtype=float)
    digest = sha256_u(u)
    if digest != SELECTED_SHA:
        raise SystemExit(f"STOP u hash {digest} != frozen {SELECTED_SHA}")
    stored = []
    for line in ahl.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            stored.append(float(s))
    if not np.allclose(u, np.array(stored), rtol=0, atol=1e-9):
        raise SystemExit("STOP AHL file does not match target u")
    samples = integrate_lorenz(SELECTED_K)
    if not np.allclose(x, samples[:NUM_WINDOWS, 0], rtol=0, atol=1e-9):
        raise SystemExit("STOP x does not match RK4 k=10")
    if not np.allclose(x_next, samples[1:NUM_WINDOWS + 1, 0], rtol=0, atol=1e-9):
        raise SystemExit("STOP x_next does not match RK4 k=10")
    if not np.allclose(y_next, samples[1:NUM_WINDOWS + 1, 1], rtol=0, atol=1e-9):
        raise SystemExit("STOP y_next drifted; AUTO_X does not use it")
    y_now = samples[:NUM_WINDOWS, 1]
    return u, x, x_next, y_now, y_next


def legal_baselines(x, x_next, y_now):
    target = x_next
    train_mean = float(np.mean(target[WASHOUT:WASHOUT + TRAIN]))
    intercept = constant_metrics(target, train_mean)
    persist = evaluate_nontrainable(target, x)
    development = []
    for m in M_GRID:
        metrics = scalar_closed(delay_matrix(x, m), target)
        development.append((m, metrics))
    selected_m, selected = min(
        development, key=lambda item: (item[1]["validation_nrmse"], item[0])
    )
    legal = {
        "TRAIN_INTERCEPT": intercept,
        "PERSISTENCE_X_N": persist,
        "LINEAR_AR_VALIDATION_SELECTED": selected,
    }
    best_name = min(legal, key=lambda name: legal[name]["test_nrmse"])
    return {
        "legal": legal,
        "selected_m": selected_m,
        "best_name": best_name,
        "best_test": legal[best_name]["test_nrmse"],
        "oracle": constant_metrics(target, float(np.mean(target[WASHOUT + TRAIN:]))),
    }


def ridge_score(X, y):
    metrics = BC.evaluate_readout(np.asarray(X), np.asarray(y), list(range(NUM_WINDOWS)))
    val = next(row["val_nrmse"] for row in metrics["lambda_grid"] if row["lambda"] == metrics["lambda"])
    return {
        "lambda": metrics["lambda"],
        "train_nrmse": metrics["train_nrmse"],
        "validation_nrmse": val,
        "test_nrmse": metrics["test_nrmse"],
        "test_r2": metrics["test_r2"],
        "n_features": np.asarray(X).shape[1],
    }


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
    tcol = header.index("TimeInWindow_s")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return {
        "window": data[:, win].astype(int),
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


def occupancy(voxels, u):
    r_win = window_mean_matrix(voxels["receiver"], voxels["window"])
    mean_r_w = np.mean(r_win, axis=1)
    mean_r = float(np.mean(mean_r_w))
    r_u = pearson(mean_r_w, u)
    if abs(r_u) < OCC_DEAD_CORR or mean_r < OCC_DEAD_R:
        label = "DEAD"
    else:
        label = "ALIVE"
    ahl_finite = bool(np.all(np.isfinite(voxels["ahl"])))
    return {
        "occupancy": label,
        "mean_R": mean_r,
        "r_meanR_u": r_u,
        "mean_AHL": float(np.mean(voxels["ahl"])),
        "mean_L": float(np.mean(voxels["lum"])),
        "ahl_finite": ahl_finite,
    }


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


def csv_ok(run_dir):
    d = Path(run_dir)
    s = BC.validate_csv(d / "window_summary.csv", NUM_WINDOWS)
    r = BC.validate_csv(d / "results.csv", EXPECTED_AUX, check_last_sample=True)
    v = BC.validate_csv(d / "voxels.csv", EXPECTED_AUX, check_last_sample=True)
    return (
        s["row_count_pass"] and s["rectangular_csv_pass"]
        and r["row_count_pass"] and r["rectangular_csv_pass"] and r.get("last_sample_pass")
        and v["row_count_pass"] and v["rectangular_csv_pass"] and v.get("last_sample_pass")
    )


def living_beats(driven, other):
    if not math.isfinite(driven) or not math.isfinite(other):
        return False
    if driven >= WEAK_BAND and other >= WEAK_BAND:
        return (other - driven) > WEAK_MARGIN
    return driven < other


def java_ok():
    java = (HERE / "BSimReservoirPlanLorenzL3.java").read_text(encoding="utf-8")
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
        "input_ahl_lorenz_k10.txt",
        "NUM_WINDOWS = 200",
        "PULSE_DURATION = 75",
    )
    return all(token in java for token in required) and java.count("new BSimChemicalField(") == 4


def gate_evidence_untouched():
    texts = []
    for path in (
        EXAMPLES / "BSimReservoirPlanBenchA" / "results" / "GATE_EVIDENCE.md",
        EXAMPLES / "BSimReservoirPlanBenchA2" / "results" / "GATE_EVIDENCE.md",
    ):
        if path.exists():
            texts.append(path.read_text(encoding="utf-8"))
    joined = "\n".join(texts)
    return "Lorenz L3" not in joined and "AUTO_X_K10" not in joined


def score_seed(seed, u, target, baselines):
    driven_dir = HERE / "results" / f"l3_k10_driven_seed{seed}"
    brown_dir = HERE / "results" / f"l3_k10_brownian_seed{seed}"
    silent_dir = SILENT_DIRS[seed]
    if not (driven_dir / "voxels.csv").exists():
        return None
    if not csv_ok(driven_dir):
        raise SystemExit(f"STOP CSV completeness failed for {driven_dir}")
    if not csv_ok(brown_dir):
        raise SystemExit(f"STOP CSV completeness failed for {brown_dir}")
    if not csv_ok(silent_dir):
        raise SystemExit(f"STOP CSV completeness failed for reused silent {silent_dir}")
    voxels = load_voxel_arrays(driven_dir / "voxels.csv")
    occ = occupancy(voxels, u)
    driven = BC.evaluate_arm_run(
        BC.read_run(driven_dir, "driven", NUM_WINDOWS, EXPECTED_AUX), u, target
    )
    brown = BC.evaluate_arm_run(
        BC.read_run(brown_dir, "brownian", NUM_WINDOWS, EXPECTED_AUX), u, target
    )
    silent = BC.evaluate_arm_run(
        BC.read_run(silent_dir, "silent", NUM_WINDOWS, EXPECTED_AUX), u, target
    )
    frl = BC.read_window_matrix(driven_dir / "voxels.csv", ("Receiver_R_", "Lum_Mean_"), ())
    r_only = BC.read_window_matrix(driven_dir / "voxels.csv", ("Receiver_R_",), ())
    l_only = BC.read_window_matrix(driven_dir / "voxels.csv", ("Lum_Mean_",), ())
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
    f408 = driven["task"]["test_nrmse"]
    field = driven["field_only"]["test_nrmse"]
    brown_n = brown["task"]["test_nrmse"]
    silent_n = silent["task"]["test_nrmse"]
    best_legal = baselines["best_test"]
    row = {
        "seed": seed,
        "silent_reused": str(silent_dir),
        **occ,
        "csv_ok": True,
        "F408_nrmse": f408,
        "F408_r2": driven["task"]["test_r2"],
        "F408_lambda": driven["task"]["lambda"],
        "field_nrmse": field,
        "FRL_nrmse": ridge_score(frl["X"], target)["test_nrmse"],
        "R_only_nrmse": ridge_score(r_only["X"], target)["test_nrmse"],
        "L_only_nrmse": ridge_score(l_only["X"], target)["test_nrmse"],
        "unmasked_RL_nrmse": ridge_score(unmasked, target)["test_nrmse"],
        "masked_RL_nrmse": ridge_score(masked, target)["test_nrmse"],
        "brownian_nrmse": brown_n,
        "silent_nrmse": silent_n,
        "best_legal_nrmse": best_legal,
        "best_legal_model": baselines["best_name"],
        "beats_legal": living_beats(f408, best_legal),
        "beats_brownian": living_beats(f408, brown_n),
        "beats_silent": living_beats(f408, silent_n),
        "field_beats_driven": field < f408,
    }
    row["living_layer_pass"] = (
        occ["occupancy"] != "DEAD"
        and row["beats_legal"]
        and row["beats_brownian"]
        and row["beats_silent"]
    )
    return row


def write_scout(baselines, rows, living_pass, authorize_more):
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    csv_path = results / "lorenz_l3_scout.csv"
    if rows:
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    legal = baselines["legal"]
    lines = [
        "# Lorenz L3 living scout",
        "",
        "Selected pair was frozen from Module 1 before Java:",
        "**k=10, Δt=0.20, AUTO_X** (`x[n] → x[n+1]`).",
        "u SHA-256 "
        f"`{SELECTED_SHA}`.",
        "Claim dish CENTER / FLOW=0. 75 s pulse / 300 s window.",
        "One AHL channel. Acid held 0.5. Silent reused from Stage 6.",
        "BenchA / BenchA2 GATE_EVIDENCE not edited. k=1 and k=50 were",
        "not rerun in BSim. Kinetics were not retuned.",
        "",
        f"## Living-layer (seed 101): {'PASS' if living_pass else 'FAIL'}",
        "",
        "Driven F408 must beat the Module-1 best legal baseline,",
        "Brownian, and silent. If both scores are ≥ 0.90, margin must",
        "exceed 0.03. Occupancy DEAD stops 202/303. Field is diagnostic.",
        "",
        "## Module-1 legal baselines on this target",
        "",
        f"Selected delay m={baselines['selected_m']}. "
        f"Best legal = {baselines['best_name']} "
        f"test NRMSE {baselines['best_test']:.4f}.",
        "",
        "| Model | λ | Train | Val | Test | R² |",
        "|---|---|---|---|---|---|",
    ]
    for name, metrics in legal.items():
        lam = metrics["lambda"]
        lam_s = f"{lam:g}" if lam != "" else ""
        lines.append(
            f"| {name} | {lam_s} | {metrics['train_nrmse']:.4f} | "
            f"{metrics['validation_nrmse']:.4f} | {metrics['test_nrmse']:.4f} | "
            f"{metrics['test_r2']:.4f} |"
        )
    lines.extend([
        "",
        "| Model | Test NRMSE | Ranked |",
        "|---|---|---|",
        f"| TEST_MEAN_ORACLE_REFERENCE | {baselines['oracle']['test_nrmse']:.4f} | no |",
        "",
        "## Scout rows",
        "",
        "| Seed | Occupancy | mean_R | r(R,u) | F408 | Field | FRL | R | L | Unmasked | Masked | Brownian | Silent | Beats legal | Living |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ])
    for row in rows:
        lines.append(
            f"| {row['seed']} | {row['occupancy']} | {row['mean_R']:.4f} | "
            f"{row['r_meanR_u']:.3f} | {row['F408_nrmse']:.4f} | "
            f"{row['field_nrmse']:.4f} | {row['FRL_nrmse']:.4f} | "
            f"{row['R_only_nrmse']:.4f} | {row['L_only_nrmse']:.4f} | "
            f"{row['unmasked_RL_nrmse']:.4f} | {row['masked_RL_nrmse']:.4f} | "
            f"{row['brownian_nrmse']:.4f} | {row['silent_nrmse']:.4f} | "
            f"{row['beats_legal']} | {row['living_layer_pass']} |"
        )
    if not rows:
        lines.append("| (no driven voxels yet) | | | | | | | | | | | | | | |")
    seed101 = next((row for row in rows if row["seed"] == 101), None)
    lines.extend(["", "## Decision", ""])
    if seed101 is None:
        lines.append("Seed 101 voxels are not present. Smoke/production still required.")
    elif seed101["occupancy"] == "DEAD":
        lines.append(
            "Seed 101 occupancy **DEAD**. Keep the row. Do not raise rate or K. "
            "Do not run 202/303. Do not retune. Do not switch target or k."
        )
    elif not seed101["beats_legal"]:
        lines.append(
            "Seed 101 does not beat the best legal baseline. Stop. "
            "Do not run 202/303. Do not retune. Do not switch target or k."
        )
    elif not living_pass:
        lines.append(
            "Seed 101 does not beat Brownian and silent under the frozen "
            "margin rule. Stop. Do not run 202/303."
        )
    elif authorize_more:
        lines.append(
            "Seed 101 living-layer PASS. Seeds 202 and 303 are authorized "
            "(driven + Brownian each; silent reused from Stage 6). "
            "No best-seed selection. Do not freeze a new claim dish."
        )
    else:
        lines.append("Seed 101 living-layer PASS recorded. Confirmation seeds scored below.")
    lines.extend([
        "",
        f"GATE_EVIDENCE untouched: {gate_evidence_untouched()}.",
        f"Java claim-dish constraints: {java_ok()}.",
        "",
    ])
    (results / "LORENZ_L3_SCOUT.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {results / 'LORENZ_L3_SCOUT.md'}")


def check_smoke():
    log = HERE / "logs" / "l3_smoke.stdout.log"
    text = ""
    if log.exists():
        raw = log.read_bytes()
        for encoding in ("utf-8-sig", "utf-16", "utf-16-le", "cp1252", "utf-8"):
            try:
                text = raw.decode(encoding)
                break
            except UnicodeDecodeError:
                continue
        text = text.replace("\x00", "")
    smoke_dir = HERE / "results" / "l3_smoke"
    voxels = smoke_dir / "voxels.csv"
    ahl_finite = False
    if voxels.exists():
        data = np.loadtxt(voxels, delimiter=";", skiprows=1)
        # AHL columns start after Window;Sample;TimeInWindow_s;...; look at header
        with voxels.open(encoding="utf-8") as handle:
            header = handle.readline().strip().split(";")
        ahl_i = [i for i, name in enumerate(header) if name.startswith("AHL_uM_")]
        if ahl_i:
            ahl_finite = bool(np.all(np.isfinite(data[:, ahl_i])))
    ok = (
        "dish=1000x500x10" in text
        and "FLOW_SPEED=0.0" in text
        and "AHL_source=(500,250,5)" in text
        and ("AHL_finite=true" in text or ahl_finite)
    )
    print("SMOKE dish/source/FLOW=0/finite AHL:", "PASS" if ok else "FAIL")
    print(text.strip() or "(no smoke log)")
    print("AHL_finite_from_voxels:", ahl_finite)
    if not ok:
        raise SystemExit("STOP smoke print missing dish/source/FLOW=0/finite AHL")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        check_smoke()
        return
    if not java_ok():
        raise SystemExit("STOP Java claim-dish constraints failed")
    if not gate_evidence_untouched():
        raise SystemExit("STOP GATE_EVIDENCE was edited")
    u, x, x_next, y_now, y_next = load_selected_target()
    del y_next
    baselines = legal_baselines(x, x_next, y_now)
    print(
        f"SELECTED k={SELECTED_K} {SELECTED_TASK} "
        f"best_legal={baselines['best_name']} "
        f"test={baselines['best_test']:.4f} m={baselines['selected_m']}"
    )
    rows = []
    seed101 = score_seed(SCOUT_SEED, u, x_next, baselines)
    if seed101:
        rows.append(seed101)
        print(
            f"seed 101 occupancy={seed101['occupancy']} "
            f"F408={seed101['F408_nrmse']:.4f} "
            f"legal={seed101['best_legal_nrmse']:.4f} "
            f"brown={seed101['brownian_nrmse']:.4f} "
            f"silent={seed101['silent_nrmse']:.4f} "
            f"living={seed101['living_layer_pass']}"
        )
        living_pass = seed101["living_layer_pass"]
        authorize = living_pass
        if living_pass:
            for seed in CONFIRM_SEEDS:
                scored = score_seed(seed, u, x_next, baselines)
                if scored:
                    rows.append(scored)
                    print(
                        f"seed {seed} occupancy={scored['occupancy']} "
                        f"F408={scored['F408_nrmse']:.4f} "
                        f"living={scored['living_layer_pass']}"
                    )
                else:
                    print(f"seed {seed} not present (authorized, not yet run)")
        else:
            print("STOP 202/303: seed 101 did not pass living-layer / occupancy")
            authorize = False
        write_scout(baselines, rows, living_pass, authorize and len(rows) == 1)
    else:
        print("seed 101 voxels not present")
        write_scout(baselines, [], False, False)


if __name__ == "__main__":
    main()
