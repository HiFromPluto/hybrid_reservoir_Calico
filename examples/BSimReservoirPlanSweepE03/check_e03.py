#!/usr/bin/env python3
"""E0.3 layout-by-flow scout: claim sanity, then occupancy/ridge/surrogates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NARMA = HERE.parent / "BSimReservoirPlanNarma10b"
SWEEPS5 = HERE.parent / "BSimReservoirPlanSweepS5"
sys.path.insert(0, str(NARMA))
import check_narma10b as N  # noqa: E402

WASHOUT, TRAIN, TEST = 40, 110, 50
NUM_WINDOWS = 200
NARMA_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
RX, RY = 20, 10
K_HILL = 1.6


def hill(c):
    c = np.asarray(c, dtype=float)
    return c * c / (K_HILL * K_HILL + c * c)


def pearson(a, b):
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return 0.0 if den < 1e-15 else float(np.sum(a * b) / den)


def constant_nrmse(y, constant):
    y_test = y[WASHOUT + TRAIN :]
    pred = np.full_like(y_test, constant)
    return N.nrmse(y_test, pred)


def delay_matrix(u, taps=10):
    X = np.zeros((len(u), taps), dtype=float)
    for n, value in enumerate(u):
        for lag in range(taps):
            if n - lag >= 0:
                X[n, lag] = u[n - lag]
    return X


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


def occupancy(voxels, u):
    r_win = window_mean_matrix(voxels["receiver"], voxels["window"])
    a_win = window_mean_matrix(voxels["ahl"], voxels["window"])
    h_win = hill(a_win)
    mean_r_w = np.mean(r_win, axis=1)
    mean_a_w = np.mean(a_win, axis=1)
    frac = np.mean(r_win > 0.5, axis=1)
    weak = float(np.mean(h_win < 0.05))
    sat = float(np.mean(h_win > 0.90))
    mean_r = float(np.mean(mean_r_w))
    r_u = pearson(mean_r_w, u)
    if abs(r_u) < 0.5 or mean_r < 0.05:
        label = "DEAD"
    elif mean_r > 0.8 or float(np.mean(frac)) > 0.8:
        label = "SATURATED"
    else:
        label = "ALIVE"
    wall = np.zeros((RX, RY), dtype=bool)
    wall[[0, -1], :] = True
    wall[:, [0, -1]] = True
    ahl_map = a_win.reshape(NUM_WINDOWS, RX, RY)
    wall_mean = float(np.mean(ahl_map[:, wall]))
    interior_mean = float(np.mean(ahl_map[:, ~wall]))
    ratio = wall_mean / interior_mean if interior_mean > 0 else float("inf")
    coverage = float(np.mean(h_win >= 0.10))
    return {
        "occupancy": label,
        "r_meanR_u": r_u,
        "mean_R": mean_r,
        "mean_AHL": float(np.mean(mean_a_w)),
        "mean_L": float(np.mean(voxels["lum"])),
        "frac_R_gt_0_5": float(np.mean(frac)),
        "weak_H_lt_0p05": weak,
        "sat_H_gt_0p90": sat,
        "mean_H": float(np.mean(h_win)),
        "wall_interior": ratio,
        "plume_coverage": coverage,
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
        er = math.exp(-delta / 15.0)
        el = math.exp(-delta / 1500.0)
        old_r, old_l = R, L
        R = h + (old_r - h) * er
        L = h + (old_l - h) * el + (old_r - h) * 15.0 / (15.0 - 1500.0) * (er - el)
        out_r[index], out_l[index] = R, L
    return out_r, out_l


def ridge_score(X, y):
    metrics = N.evaluate_readout(np.asarray(X), np.asarray(y), list(range(NUM_WINDOWS)))
    return {
        "lambda": metrics["lambda"],
        "test_nrmse": metrics["test_nrmse"],
        "test_r2": metrics["test_r2"],
        "n_features": np.asarray(X).shape[1],
    }


def summary_stats(run_dir):
    path = Path(run_dir) / "window_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    def col(name):
        return np.array([float(row[name]) for row in rows], dtype=float)
    return {
        "mean_pop": float(np.mean(col("Population"))),
        "births": float(np.sum(col("Births"))),
        "clamp_deaths": float(np.sum(col("Clamp_Deaths"))),
        "acid_deaths": float(np.sum(col("Input_Driven_Deaths"))),
        "oob_deaths": float(np.sum(col("OOB_Deaths"))),
    }


def csv_ok(run_dir, expected_windows=200, expected_aux=3200):
    d = Path(run_dir)
    s = N.validate_csv(d / "window_summary.csv", expected_windows)
    r = N.validate_csv(d / "results.csv", expected_aux, check_last_sample=True)
    v = N.validate_csv(d / "voxels.csv", expected_aux, check_last_sample=True)
    return s["row_count_pass"] and r["row_count_pass"] and v["row_count_pass"] and r.get("last_sample_pass") and v.get("last_sample_pass")


def sanity():
    u, y, digest = N.load_target(HERE / "narma10_target.csv", HERE / "input_ahl_narma200.txt")
    if digest != NARMA_SHA:
        raise SystemExit(f"ABORT u hash {digest}")
    driven = N.read_run(NARMA / "results" / "narma10b_driven_seed111", "driven", NUM_WINDOWS, 3200)
    brown = N.read_run(NARMA / "results" / "narma10b_brownian_seed111", "brownian", NUM_WINDOWS, 3200)
    silent = N.read_run(NARMA / "results" / "narma10b_silent_seed111", "silent", NUM_WINDOWS, 3200)
    d = N.evaluate_arm_run(driven, u, y)
    b = N.evaluate_arm_run(brown, u, y)
    s = N.evaluate_arm_run(silent, u, y)
    intercept = constant_nrmse(y, float(np.mean(y[WASHOUT:WASHOUT + TRAIN])))
    persistence = N.nrmse(y[WASHOUT + TRAIN :], y[WASHOUT + TRAIN - 1 : WASHOUT + TRAIN + TEST - 1])
    checks = {
        "driven_F408": (d["task"]["test_nrmse"], 0.9312),
        "field": (d["field_only"]["test_nrmse"], 1.0289),
        "brownian": (b["task"]["test_nrmse"], 1.1625),
        "intercept": (intercept, 1.1622),
        "persistence": (persistence, 1.0068),
    }
    print("CLAIM SANITY seed 111")
    failed = False
    for name, (got, expected) in checks.items():
        ok = abs(got - expected) <= 1e-3
        print(f"  {name:12} {got:.4f} expected {expected:.4f} {'PASS' if ok else 'FAIL'}")
        failed = failed or not ok
    if failed:
        raise SystemExit("STOP claim sanity failed; do not interpret new conditions")
    print("SANITY PASS")
    return u, y


CONDITIONS = [
    ("PRI_CENTER_F0p0", 500.0, 0.0, NARMA / "results" / "narma10b_driven_seed111", True),
    ("PRI_CENTER_F0p25", 500.0, 0.25, HERE / "results" / "e03_center_f0p25_driven_seed111", False),
    ("PRI_CENTER_F0p5", 500.0, 0.5, HERE / "results" / "e03_center_f0p5_driven_seed111", False),
    ("PRI_CENTER_F0p666667", 500.0, 0.666667, HERE / "results" / "e03_center_f0p666667_driven_seed111", False),
    ("PRI_CENTER_F1p0", 500.0, 1.0, HERE / "results" / "e03_center_f1p0_driven_seed111", False),
    ("PRI_UPSTREAM_CENTER_F0p0", 200.0, 0.0, HERE / "results" / "e03_upstream_f0p0_driven_seed111", False),
    ("PRI_UPSTREAM_CENTER_F0p666667", 200.0, 0.666667, HERE / "results" / "e03_upstream_f0p666667_driven_seed111", False),
    ("PRI_DOWNSTREAM_CENTER_F0p0", 800.0, 0.0, HERE / "results" / "e03_downstream_f0p0_driven_seed111", False),
    ("PRI_DOWNSTREAM_CENTER_F0p666667", 800.0, 0.666667, HERE / "results" / "e03_downstream_f0p666667_driven_seed111", False),
]


def null_dir(flow, arm):
    if flow == 0.0:
        return NARMA / "results" / f"narma10b_{arm}_seed111"
    token = str(flow).replace(".", "p")
    return HERE / "results" / f"e03_f{token}_{arm}_seed111"


def analyze(u, y):
    taps = delay_matrix(u, 10)
    product = np.array([u[n] * u[n - 9] if n >= 9 else 0.0 for n in range(len(u))])
    informed = np.column_stack([taps, product])
    baselines = {
        "LINEAR_U_DELAY_10": ridge_score(taps, y),
        "NARMA_INFORMED_INPUT": ridge_score(informed, y),
        "TRAIN_INTERCEPT": {
            "lambda": "",
            "test_nrmse": constant_nrmse(y, float(np.mean(y[WASHOUT:WASHOUT + TRAIN]))),
            "test_r2": float("nan"),
            "n_features": 0,
        },
        "PERSISTENCE": {
            "lambda": "",
            "test_nrmse": N.nrmse(y[WASHOUT + TRAIN :], y[WASHOUT + TRAIN - 1 : WASHOUT + TRAIN + TEST - 1]),
            "test_r2": float("nan"),
            "n_features": 1,
        },
    }
    rows = []
    for cid, x, flow, driven_dir, reused in CONDITIONS:
        if not (Path(driven_dir) / "voxels.csv").exists():
            print(f"missing {driven_dir}")
            continue
        if not csv_ok(driven_dir):
            raise SystemExit(f"STOP CSV completeness failed for {driven_dir}")
        voxels = load_voxel_arrays(Path(driven_dir) / "voxels.csv")
        occ = occupancy(voxels, u)
        stats = summary_stats(driven_dir)
        driven = N.read_run(driven_dir, "driven", NUM_WINDOWS, 3200)
        scored = N.evaluate_arm_run(driven, u, y)
        frl = N.read_window_matrix(Path(driven_dir) / "voxels.csv", ("Receiver_R_", "Lum_Mean_"), ())
        r_only = N.read_window_matrix(Path(driven_dir) / "voxels.csv", ("Receiver_R_",), ())
        l_only = N.read_window_matrix(Path(driven_dir) / "voxels.csv", ("Lum_Mean_",), ())
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
        silent_dir = null_dir(flow, "silent")
        brown_dir = null_dir(flow, "brownian")
        silent = N.evaluate_arm_run(N.read_run(silent_dir, "silent", NUM_WINDOWS, 3200), u, y)
        brown = N.evaluate_arm_run(N.read_run(brown_dir, "brownian", NUM_WINDOWS, 3200), u, y)
        transit = float("inf") if flow == 0 else 1000.0 / flow
        courant = flow * 0.05 / 20.0
        row = {
            "ConditionID": cid,
            "ReusedClaimRun": reused,
            "AHL_x": x,
            "Flow_um_s": flow,
            "Boundary": "NO_FLUX" if flow == 0 else "OUTFLOW",
            "Transit_s": transit,
            "Courant": courant,
            **occ,
            **stats,
            "driven_F408_nrmse": scored["task"]["test_nrmse"],
            "driven_F408_r2": scored["task"]["test_r2"],
            "field_nrmse": scored["field_only"]["test_nrmse"],
            "FRL_nrmse": ridge_score(frl["X"], y)["test_nrmse"],
            "R_only_nrmse": ridge_score(r_only["X"], y)["test_nrmse"],
            "L_only_nrmse": ridge_score(l_only["X"], y)["test_nrmse"],
            "unmasked_RL_nrmse": ridge_score(unmasked, y)["test_nrmse"],
            "masked_RL_nrmse": ridge_score(masked, y)["test_nrmse"],
            "silent_nrmse": silent["task"]["test_nrmse"],
            "brownian_nrmse": brown["task"]["test_nrmse"],
            "direct10_nrmse": baselines["LINEAR_U_DELAY_10"]["test_nrmse"],
            "informed_nrmse": baselines["NARMA_INFORMED_INPUT"]["test_nrmse"],
            "intercept_nrmse": baselines["TRAIN_INTERCEPT"]["test_nrmse"],
            "persistence_nrmse": baselines["PERSISTENCE"]["test_nrmse"],
        }
        exclusions = []
        if occ["occupancy"] in ("DEAD", "SATURATED"):
            exclusions.append(occ["occupancy"])
        if occ["mean_R"] < 0.02 or occ["weak_H_lt_0p05"] > 0.90:
            exclusions.append("WEAK_RECEIVER")
        if occ["mean_H"] > 0.90 or occ["sat_H_gt_0p90"] > 0.90:
            exclusions.append("SATURATED_H")
        if stats["mean_pop"] < 400:
            exclusions.append("POPULATION_COLLAPSE")
        row["Exclusions"] = "|".join(exclusions)
        rows.append(row)
        print(cid, occ["occupancy"], f"R={occ['mean_R']:.4f}", f"F408={row['driven_F408_nrmse']:.4f}")

    # historical 8 um/s: cite SweepS5 voxels; do not rerun BSim
    s5 = SWEEPS5 / "results" / "s5_flow_driven_seed111"
    s5_brown = SWEEPS5 / "results" / "s5_flow_brownian_seed111"
    if (s5 / "voxels.csv").exists():
        voxels = load_voxel_arrays(s5 / "voxels.csv")
        occ = occupancy(voxels, u)
        driven = N.read_run(s5, "driven", NUM_WINDOWS, 3200)
        scored = N.evaluate_arm_run(driven, u, y)
        frl = N.read_window_matrix(s5 / "voxels.csv", ("Receiver_R_", "Lum_Mean_"), ())
        r_only = N.read_window_matrix(s5 / "voxels.csv", ("Receiver_R_",), ())
        l_only = N.read_window_matrix(s5 / "voxels.csv", ("Lum_Mean_",), ())
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
        if (s5_brown / "voxels.csv").exists():
            brown_nrmse = N.evaluate_arm_run(
                N.read_run(s5_brown, "brownian", NUM_WINDOWS, 3200), u, y
            )["task"]["test_nrmse"]
        else:
            brown_nrmse = 1.1623
        rows.append({
            "ConditionID": "HISTORICAL_CENTER_F8_ANCHOR",
            "ReusedClaimRun": True,
            "AHL_x": 500.0,
            "Flow_um_s": 8.0,
            "Boundary": "OUTFLOW",
            "Transit_s": 125.0,
            "Courant": 0.02,
            **occ,
            **summary_stats(s5),
            "driven_F408_nrmse": scored["task"]["test_nrmse"],
            "driven_F408_r2": scored["task"]["test_r2"],
            "field_nrmse": scored["field_only"]["test_nrmse"],
            "FRL_nrmse": ridge_score(frl["X"], y)["test_nrmse"],
            "R_only_nrmse": ridge_score(r_only["X"], y)["test_nrmse"],
            "L_only_nrmse": ridge_score(l_only["X"], y)["test_nrmse"],
            "unmasked_RL_nrmse": ridge_score(unmasked, y)["test_nrmse"],
            "masked_RL_nrmse": ridge_score(masked, y)["test_nrmse"],
            "silent_nrmse": float("nan"),
            "brownian_nrmse": brown_nrmse,
            "direct10_nrmse": baselines["LINEAR_U_DELAY_10"]["test_nrmse"],
            "informed_nrmse": baselines["NARMA_INFORMED_INPUT"]["test_nrmse"],
            "intercept_nrmse": baselines["TRAIN_INTERCEPT"]["test_nrmse"],
            "persistence_nrmse": baselines["PERSISTENCE"]["test_nrmse"],
            "Exclusions": "DEAD|HISTORICAL_DEAD_FLOW_ANCHOR_NO_NEW_RUN",
        })
        print("HISTORICAL_CENTER_F8_ANCHOR", occ["occupancy"], f"R={occ['mean_R']:.4f} (SweepS5, not rerun)")

    out = HERE / "results"
    out.mkdir(exist_ok=True)
    csv_path = out / "e0_3_layout_flow_scout.csv"
    if rows:
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print("wrote", csv_path)
    return rows, baselines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sanity-only", action="store_true")
    args = parser.parse_args()
    u, y = sanity()
    if args.sanity_only:
        return
    analyze(u, y)


if __name__ == "__main__":
    main()
