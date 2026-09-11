#!/usr/bin/env python3
"""E5.1 isolated one-AHL CRP scout checker.

Prints the four claims. Refuses to retune K, n, tau, Jmax, clamp,
mortality, flow, or layout. Does not classify Y. SYNTHETIC_COHORT stays.
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
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MAT = REPO / "examples" / "HybridDish" / "external_analysis" / "synthetic_biomarker_BSim_inputs.mat"

DEV = [1, 4, 6, 7, 9, 13, 16, 18, 20, 22, 24, 25]
CONF = [8, 12, 14, 15, 17, 19, 21, 23]
DEV_SHA = "83acc7102e532112c961bb45a1327504290ffddea1c06e659bfb5e68be032915"
CONF_SHA = "c90cca8f5627c737b38024b4863415e46a384d8715e6d06a5823e8cb93ba3787"
LAMBDA_GRID = np.array([1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6])
NUM_WINDOWS = 48
EXPECTED_AUX = 768
LAST_SAMPLE_PREFIX = "47;15;299.95"
OCC_DEAD_R = 0.05
OCC_DEAD_CORR = 0.5
WEAK_MARGIN = 0.03
T0, T_ALIGN, T_LAST = 0, 9, 46
K_HILL = 1.6
FORBIDDEN_JAVA = (
    "glucose", "Glucose", "GLUCOSE",
    "Danino", "danino",
    "vesicle", "Vesicle",
    "setGoal(glucose)",
)

BIOLOGY_MEAN = ("Receiver_R_", "Lum_Mean_")
BIOLOGY_LAST = ("Input_Driven_Death_",)
FIELD_MEAN = ("AHL_uM_",)
BROWN_MEAN = ("Den_",)


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def abort(msg):
    print(f"ABORT: {msg}", file=sys.stderr)
    raise SystemExit(1)


def refuse_retune(reason):
    print("REFUSE_RETUNE: do not raise Jmax, K, n, tau, clamp, mortality,")
    print("flow, or layout. Do not redraw the scout. Do not classify Y.")
    print(reason)


def nrmse_pop(ytrue, ypred):
    ytrue = np.asarray(ytrue, dtype=float).ravel()
    ypred = np.asarray(ypred, dtype=float).ravel()
    rmse = np.sqrt(np.mean((ytrue - ypred) ** 2))
    denom = ytrue.std(ddof=0)
    return float(rmse / denom) if denom > 0 else float("inf")


def pearson(a, b):
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return 0.0 if den < 1e-15 else float(np.sum(a * b) / den)


def standardize_fit(X):
    mu = X.mean(axis=0)
    sg = X.std(axis=0, ddof=0)
    sg = np.where(sg < 1e-12, 1.0, sg)
    return (X - mu) / sg, mu, sg


def ridge_fit(X, y, lam):
    n = X.shape[1]
    Xb = np.column_stack([np.ones(len(X)), X])
    G = Xb.T @ Xb
    G[1:, 1:] = G[1:, 1:] + lam * np.eye(n)
    try:
        return np.linalg.solve(G, Xb.T @ y)
    except np.linalg.LinAlgError:
        w, *_ = np.linalg.lstsq(G, Xb.T @ y, rcond=None)
        return w


def ridge_predict(X, w):
    return np.column_stack([np.ones(len(X)), X]) @ w


def load_u(p):
    path = HERE / "input" / f"input_ahl_p{p:04d}.txt"
    vals = [float(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    return np.asarray(vals, float)


def load_target(p):
    path = HERE / "input" / f"target_crp_p{p:04d}.txt"
    vals = [float(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    return np.asarray(vals, float)


def last_sample_ok(path):
    lines = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines()[1:] if ln.strip()]
    return bool(lines) and lines[-1].startswith(LAST_SAMPLE_PREFIX)


def validate_csv(path, expected_rows, check_last=False):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, [])
        widths = [len(row) for row in reader]
    out = {
        "path": str(path),
        "row_count": len(widths),
        "row_count_pass": len(widths) == expected_rows,
        "rectangular": bool(header) and all(w == len(header) for w in widths),
    }
    if check_last:
        out["last_sample_pass"] = last_sample_ok(path)
    return out


def columns_with_prefix(header, prefixes):
    names = []
    for prefix in prefixes:
        names.extend(n for n in header if n.startswith(prefix))
    return names


def read_window_matrix(voxels_path, mean_prefixes, last_prefixes=()):
    path = Path(voxels_path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        mean_cols = columns_with_prefix(header, mean_prefixes)
        last_cols = columns_with_prefix(header, last_prefixes)
        feature_names = mean_cols + last_cols
        window_idx = header.index("Window")
        mean_idx = [header.index(n) for n in mean_cols]
        last_idx = [header.index(n) for n in last_cols]
        sums, counts, lasts = {}, {}, {}
        for row in reader:
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
    rows = []
    for window in windows:
        parts = []
        if mean_idx:
            parts.append(sums[window] / counts[window])
        if last_idx:
            parts.append(lasts[window])
        rows.append(np.concatenate(parts) if parts else np.zeros(0))
    return {"windows": windows, "X": np.vstack(rows), "feature_names": feature_names}


def load_voxel_arrays(path):
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    ahl_i = [header.index(f"AHL_uM_{i}") for i in range(200)]
    den_i = [header.index(f"Den_{i}") for i in range(200)]
    r_i = [header.index(f"Receiver_R_{i}") for i in range(200)]
    l_i = [header.index(f"Lum_Mean_{i}") for i in range(200)]
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    win = header.index("Window")
    tcol = header.index("TimeInWindow_s")
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


def hill(c):
    c = np.asarray(c, dtype=float)
    return c * c / (K_HILL * K_HILL + c * c)


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


def occupancy_from_voxels(voxels, u):
    r_win = window_mean_matrix(voxels["receiver"], voxels["window"])
    mean_r_w = np.mean(r_win, axis=1)
    mean_r = float(np.mean(mean_r_w))
    r_u = pearson(mean_r_w, u)
    dead = abs(r_u) < OCC_DEAD_CORR or mean_r < OCC_DEAD_R
    return {
        "occupancy": "DEAD" if dead else "ALIVE",
        "mean_R": mean_r,
        "r_meanR_u": r_u,
        "mean_AHL": float(np.mean(voxels["ahl"])),
        "ahl_finite": voxels["ahl_finite"],
        "ahl_nonneg": voxels["ahl_nonneg"],
    }


def stack_pairs(X_by_patient, patients, t_lo, t_hi):
    feats, targets, pids = [], [], []
    for p in patients:
        X = X_by_patient[p]
        y = load_target(p)
        for t in range(t_lo, t_hi + 1):
            feats.append(X[t])
            targets.append(y[t])
            pids.append(p)
    return np.asarray(feats, float), np.asarray(targets, float), np.asarray(pids, int)


def loo_lambda(X, y, pids, patients):
    best, best_lam = np.inf, LAMBDA_GRID[-1]
    grid = []
    for lam in LAMBDA_GRID:
        scores = []
        for held in patients:
            tr = pids != held
            va = pids == held
            Xz, mu, sg = standardize_fit(X[tr])
            w = ridge_fit(Xz, y[tr], lam)
            pred = ridge_predict((X[va] - mu) / sg, w)
            scores.append(nrmse_pop(y[va], pred))
        pooled = float(np.mean(scores))
        grid.append({"lambda": float(lam), "loo_pooled": pooled})
        if pooled < best - 1e-15 or (abs(pooled - best) <= 1e-15 and lam > best_lam):
            best, best_lam = pooled, lam
    return float(best_lam), float(best), grid


def fit_confirm(Xtr, ytr, ptr, Xte, yte, pte, lam):
    Xz, mu, sg = standardize_fit(Xtr)
    w = ridge_fit(Xz, ytr, lam)
    pred_tr = ridge_predict(Xz, w)
    pred_te = ridge_predict((Xte - mu) / sg, w)
    tr_scores = [nrmse_pop(ytr[ptr == p], pred_tr[ptr == p]) for p in np.unique(ptr)]
    te_scores = [nrmse_pop(yte[pte == p], pred_te[pte == p]) for p in np.unique(pte)]
    return {
        "lambda": float(lam),
        "development_pooled": float(np.mean(tr_scores)),
        "confirmation_pooled": float(np.mean(te_scores)),
        "confirmation_concat": nrmse_pop(yte, pred_te),
        "confirmation_per_patient": {
            str(int(p)): float(nrmse_pop(yte[pte == p], pred_te[pte == p]))
            for p in np.unique(pte)
        },
        "n_features": int(Xtr.shape[1]),
    }


def score_readout(X_by_patient, t_lo=T_ALIGN, t_hi=T_LAST):
    Xtr, ytr, ptr = stack_pairs(X_by_patient, DEV, t_lo, t_hi)
    Xte, yte, pte = stack_pairs(X_by_patient, CONF, t_lo, t_hi)
    lam, loo, grid = loo_lambda(Xtr, ytr, ptr, DEV)
    metrics = fit_confirm(Xtr, ytr, ptr, Xte, yte, pte, lam)
    metrics["loo_pooled_development"] = loo
    metrics["loo_grid"] = grid
    return metrics


def run_dir(kind, p=None):
    if kind == "driven":
        return HERE / "results" / f"e51_driven_p{p:04d}_seed111"
    if kind == "silent":
        return HERE / "results" / "e51_silent_seed111"
    if kind == "brownian":
        return HERE / "results" / "e51_brownian_seed111"
    raise ValueError(kind)


def java_ok():
    java = (HERE / "BSimReservoirPlanE51CRP.java").read_text(encoding="utf-8")
    voxel = (HERE / "VoxelAnalyzer.java").read_text(encoding="utf-8")
    if any(tok in java for tok in FORBIDDEN_JAVA):
        return False
    if "FLOW_SPEED = 0.0" not in java:
        return False
    if "RECEIVER_K_UM = 1.6" not in java:
        return False
    if "new Vector3d(500, 250, 5)" not in java:
        return False
    if "new Vector3d(300, 375, 5)" not in java:
        return False
    if "NUM_WINDOWS = 48" not in java:
        return False
    if "package BSimReservoirPlanE51CRP;" not in voxel:
        return False
    return True


def verify_hashes():
    if sha256_text(",".join(str(i) for i in DEV)) != DEV_SHA:
        abort("development scout hash mismatch")
    if sha256_text(",".join(str(i) for i in CONF)) != CONF_SHA:
        abort("confirmation scout hash mismatch")
    recorded = json.loads((HERE / "results" / "e51_input_hashes.json").read_text(encoding="utf-8"))
    for p in DEV + CONF:
        path = HERE / "input" / f"input_ahl_p{p:04d}.txt"
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        exp = recorded["files"][path.name]["file_sha256"]
        if got != exp:
            abort(f"input hash drift {path.name}")
    return recorded


def load_ceilings():
    return json.loads((HERE / "results" / "e51_input_ceilings.json").read_text(encoding="utf-8"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke-only", action="store_true")
    args = parser.parse_args()
    verify_hashes()
    ceilings = load_ceilings()
    if not java_ok():
        abort("Java failed frozen-dish / forbidden-token checks")

    if args.smoke_only:
        smoke = run_dir("driven", DEV[0])
        voxels_path = smoke / "voxels.csv"
        if not voxels_path.exists():
            abort(f"smoke voxels missing: {voxels_path}")
        csvs = validate_csv(voxels_path, EXPECTED_AUX, check_last=True)
        summary = validate_csv(smoke / "window_summary.csv", NUM_WINDOWS)
        voxels = load_voxel_arrays(voxels_path)
        u = load_u(DEV[0])
        occ = occupancy_from_voxels(voxels, u)
        print("SMOKE patient", DEV[0])
        print("u[0]", float(u[0]))
        print("last_sample", LAST_SAMPLE_PREFIX, csvs.get("last_sample_pass"))
        print("summary_rows", summary["row_count"], "voxel_rows", csvs["row_count"])
        print("AHL finite", occ["ahl_finite"], "nonneg", occ["ahl_nonneg"], "mean_AHL", occ["mean_AHL"])
        print("occupancy", occ["occupancy"], "mean_R", occ["mean_R"], "r(R,u)", occ["r_meanR_u"])
        print("COHORT_STATUS SYNTHETIC_COHORT")
        print("CENTER FLOW=0 one-AHL")
        return

    missing = []
    for p in DEV:
        if not (run_dir("driven", p) / "voxels.csv").exists():
            missing.append(f"driven p{p:04d}")
    for name in ("silent", "brownian"):
        if not (run_dir(name) / "voxels.csv").exists():
            missing.append(name)
    if missing:
        print("Module 0 hashes OK. BSim voxels missing:", ", ".join(missing))
        print("COHORT_STATUS SYNTHETIC_COHORT")
        print("Run isolated BSim before Module 2 scores.")
        raise SystemExit(2)

    occ_rows = []
    X408, Xfield, Xmask = {}, {}, {}
    dead_dev = 0

    def load_driven(p, split):
        d = run_dir("driven", p)
        vox_ok = validate_csv(d / "voxels.csv", EXPECTED_AUX, check_last=True)
        sum_ok = validate_csv(d / "window_summary.csv", NUM_WINDOWS)
        if not (vox_ok["row_count_pass"] and vox_ok["rectangular"] and vox_ok["last_sample_pass"] and sum_ok["row_count_pass"]):
            abort(f"CSV incomplete for driven p{p:04d}")
        voxels = load_voxel_arrays(d / "voxels.csv")
        u = load_u(p)
        occ = occupancy_from_voxels(voxels, u)
        occ["patient"] = p
        occ["split"] = split
        bio = read_window_matrix(d / "voxels.csv", BIOLOGY_MEAN, BIOLOGY_LAST)
        field = read_window_matrix(d / "voxels.csv", FIELD_MEAN, ())
        if len(bio["feature_names"]) != 408:
            abort(f"expected 408 features, got {len(bio['feature_names'])} for p{p:04d}")
        sr, sl = integrate_surrogate(voxels["absolute_time"], voxels["ahl"])
        mask = (voxels["den"] != 0).astype(float)
        X408[p] = bio["X"]
        Xfield[p] = field["X"]
        Xmask[p] = np.column_stack([
            window_mean_matrix(sr * mask, voxels["window"]),
            window_mean_matrix(sl * mask, voxels["window"]),
        ])
        return occ

    for p in DEV:
        occ = load_driven(p, "development")
        occ_rows.append(occ)
        if occ["occupancy"] == "DEAD":
            dead_dev += 1

    skip_conf = dead_dev >= 3
    if skip_conf:
        print("Development occupancy DEAD on", dead_dev, "of 12. Confirmation not scored.")
        refuse_retune("Occupancy DEAD. Do not raise Jmax.")
        ceilings = load_ceilings()
        persist = ceilings["persistence_aligned_t9_46"]["confirmation_pooled"]
        persist_all = ceilings["persistence_t0_46"]["confirmation_pooled"]
        crp10 = ceilings["crp_delay10"]["confirmation_pooled"]
        ch5 = ceilings["five_channel_delay10"]["confirmation_pooled"]
        payload = {
            "cohort_status": "SYNTHETIC_COHORT",
            "clinical": False,
            "skip_confirmation": True,
            "dead_development": dead_dev,
            "occupancy": occ_rows,
            "F408": "NOT_SCORED_CONFIRMATION_SKIPPED",
            "field": "NOT_SCORED_CONFIRMATION_SKIPPED",
            "ceilings": {
                "persistence_t0_46_conf": persist_all,
                "persistence_aligned_conf": persist,
                "crp_delay10_conf": crp10,
                "five_channel_delay10_conf": ch5,
            },
            "claims": {
                "system": "NOT_SCORED",
                "living_layer": "NOT_SCORED",
                "ceiling": "NOT_REQUIRED",
                "fair_input": "NOT_SCORED",
                "clinical": "FORBIDDEN",
                "no_story_move": True,
                "cohort_status": "SYNTHETIC_COHORT",
            },
        }
        (HERE / "results" / "e51_scout.json").write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        with (HERE / "results" / "e51_scout.csv").open("w", newline="", encoding="utf-8") as handle:
            w = csv.DictWriter(handle, fieldnames=["readout", "split", "nrmse_patient_pooled", "note"])
            w.writeheader()
            w.writerow({"readout": "F408", "split": "confirmation", "nrmse_patient_pooled": "", "note": "NOT_SCORED occupancy DEAD"})
            w.writerow({"readout": "field", "split": "confirmation", "nrmse_patient_pooled": "", "note": "NOT_SCORED occupancy DEAD"})
            w.writerow({"readout": "persistence_aligned", "split": "confirmation", "nrmse_patient_pooled": f"{persist:.10g}", "note": "Module 0 ceiling"})
            w.writerow({"readout": "crp_delay10", "split": "confirmation", "nrmse_patient_pooled": f"{crp10:.10g}", "note": "Module 0 ceiling"})
            w.writerow({"readout": "five_channel_delay10", "split": "confirmation", "nrmse_patient_pooled": f"{ch5:.10g}", "note": "Module 0 ceiling"})
        lines = [
            "# E5.1 isolated one-AHL CRP scout",
            "",
            "Seed 111 only. Task B only. `COHORT_STATUS: SYNTHETIC_COHORT`.",
            "Not a clinical dataset. Track B FAIL is not rewritten.",
            "Confirmation driven BSim was **not** run: development occupancy",
            f"was DEAD on {dead_dev} of 12 patients (`mean_R < 0.05` and/or",
            "`|r(mean_R, u)| < 0.5`). Keep DEAD rows. Do not raise `Jmax`.",
            "",
            "## Claims",
            "",
            "| Claim | Result |",
            "|---|---|",
            "| System (F408 vs Brownian and silent) | **NOT_SCORED** |",
            "| Living-layer (F408 vs field) | **NOT_SCORED** |",
            "| Ceiling (beat 5-ch delay-10) | **NOT_REQUIRED** |",
            "| Fair input | Module 0 ceilings only |",
            "| Clinical | **FORBIDDEN** |",
            "| NO_STORY_MOVE | yes (no living confirmation score) |",
            "",
            "## Development occupancy (seed 111)",
            "",
            "| Patient | Occupancy | mean_R | r(mean_R, u) | mean_AHL |",
            "|---|---|---|---|---|",
        ]
        for row in occ_rows:
            lines.append(
                f"| {row['patient']} | {row['occupancy']} | {row['mean_R']:.4f} | "
                f"{row['r_meanR_u']:.3f} | {row['mean_AHL']:.4f} |"
            )
        lines.extend([
            "",
            "## Module 0 confirmation ceilings (12+8, no BSim)",
            "",
            "| Baseline | Conf pooled NRMSE |",
            "|---|---|",
            f"| persistence `t=0..46` | {persist_all:.4f} |",
            f"| persistence aligned `t=9..46` | {persist:.4f} |",
            f"| CRP-only delay-10 | {crp10:.4f} |",
            f"| 5-channel delay-10 (ceiling) | {ch5:.4f} |",
            "",
            "F408 / field were not scored. Silent and Brownian nulls exist",
            "but were not used as a living-layer gate without confirmation",
            "driven CSVs. REFUSE_RETUNE.",
            "",
        ])
        (HERE / "results" / "E5_1_SCOUT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        print("COHORT_STATUS SYNTHETIC_COHORT")
        print("occupancy DEAD development", dead_dev, "/12")
        print("F408 NOT_SCORED")
        print("field NOT_SCORED")
        print("persistence_aligned", persist)
        print("crp_delay10", crp10)
        print("five_channel_delay10", ch5)
        raise SystemExit(0)

    for p in CONF:
        occ_rows.append(load_driven(p, "confirmation"))

    silent_bio = read_window_matrix(run_dir("silent") / "voxels.csv", BIOLOGY_MEAN, BIOLOGY_LAST)
    brown_den = read_window_matrix(run_dir("brownian") / "voxels.csv", BROWN_MEAN, ())
    Xsilent = {p: silent_bio["X"] for p in DEV + CONF}
    Xbrown = {p: brown_den["X"] for p in DEV + CONF}

    f408 = score_readout(X408)
    field = score_readout(Xfield)
    masked = score_readout(Xmask)
    silent = score_readout(Xsilent)
    brown = score_readout(Xbrown)
    f408_diag = score_readout(X408, T0, T_LAST)

    persist = ceilings["persistence_aligned_t9_46"]["confirmation_pooled"]
    persist_all = ceilings["persistence_t0_46"]["confirmation_pooled"]
    crp10 = ceilings["crp_delay10"]["confirmation_pooled"]
    ch5 = ceilings["five_channel_delay10"]["confirmation_pooled"]
    f = f408["confirmation_pooled"]
    fld = field["confirmation_pooled"]
    system_pass = f < brown["confirmation_pooled"] and f < silent["confirmation_pooled"]
    living_pass = f < fld
    delta = abs(f - fld)
    no_story = delta < WEAK_MARGIN and not living_pass
    system = "PASS" if system_pass else "FAIL"
    living = "PASS" if living_pass else "FAIL"
    ceiling = "NOT_REQUIRED"
    fair = {
        "vs_persistence_aligned": f - persist,
        "vs_crp_delay10": f - crp10,
    }
    clinical = "FORBIDDEN"

    claims = {
        "system": system,
        "living_layer": living,
        "ceiling": ceiling,
        "fair_input": fair,
        "clinical": clinical,
        "no_story_move": no_story,
        "cohort_status": "SYNTHETIC_COHORT",
    }

    payload = {
        "cohort_status": "SYNTHETIC_COHORT",
        "clinical": False,
        "last_sample_prefix": LAST_SAMPLE_PREFIX,
        "occupancy": occ_rows,
        "dead_development": dead_dev,
        "F408": f408,
        "field": field,
        "masked_surrogate": masked,
        "silent": silent,
        "brownian": brown,
        "F408_t0_46_diagnostic": f408_diag,
        "ceilings": {
            "persistence_t0_46_conf": persist_all,
            "persistence_aligned_conf": persist,
            "crp_delay10_conf": crp10,
            "five_channel_delay10_conf": ch5,
        },
        "claims": claims,
        "delta_F408_field": f - fld,
    }
    (HERE / "results" / "e51_scout.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )

    with (HERE / "results" / "e51_scout.csv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=["readout", "split", "nrmse_patient_pooled", "nrmse_concat", "lambda", "n_features"])
        w.writeheader()
        for name, row in (
            ("F408", f408),
            ("field", field),
            ("masked_surrogate", masked),
            ("silent", silent),
            ("brownian", brown),
        ):
            w.writerow({
                "readout": name,
                "split": "confirmation",
                "nrmse_patient_pooled": f"{row['confirmation_pooled']:.10g}",
                "nrmse_concat": f"{row['confirmation_concat']:.10g}",
                "lambda": f"{row['lambda']:g}",
                "n_features": row["n_features"],
            })
        w.writerow({
            "readout": "persistence_aligned",
            "split": "confirmation",
            "nrmse_patient_pooled": f"{persist:.10g}",
            "nrmse_concat": f"{ceilings['persistence_aligned_t9_46']['confirmation_concat_diagnostic']:.10g}",
            "lambda": "",
            "n_features": 1,
        })
        w.writerow({
            "readout": "crp_delay10",
            "split": "confirmation",
            "nrmse_patient_pooled": f"{crp10:.10g}",
            "nrmse_concat": f"{ceilings['crp_delay10']['confirmation_concat']:.10g}",
            "lambda": f"{ceilings['crp_delay10']['lambda']:g}",
            "n_features": 10,
        })
        w.writerow({
            "readout": "five_channel_delay10",
            "split": "confirmation",
            "nrmse_patient_pooled": f"{ch5:.10g}",
            "nrmse_concat": f"{ceilings['five_channel_delay10']['confirmation_concat']:.10g}",
            "lambda": f"{ceilings['five_channel_delay10']['lambda']:g}",
            "n_features": 50,
        })

    occ_alive = sum(1 for r in occ_rows if r["occupancy"] == "ALIVE")
    md = []
    md.append("# E5.1 isolated one-AHL CRP scout")
    md.append("")
    md.append("Seed 111 only. Task B only. `COHORT_STATUS: SYNTHETIC_COHORT`.")
    md.append("Not a clinical dataset. Track B FAIL is not rewritten.")
    md.append("K, n, tau, Jmax, clamp, mortality, flow, and layout were not retuned.")
    md.append("")
    md.append("## Claims")
    md.append("")
    md.append("| Claim | Result |")
    md.append("|---|---|")
    md.append(f"| System (F408 vs Brownian and silent) | **{system}** |")
    md.append(f"| Living-layer (F408 vs field) | **{living}** |")
    md.append("| Ceiling (beat 5-ch delay-10) | **NOT_REQUIRED** |")
    md.append("| Fair input (vs persistence and CRP delay-10) | report only |")
    md.append("| Clinical | **FORBIDDEN** |")
    md.append(f"| NO_STORY_MOVE | {'yes' if no_story else 'no'} |")
    md.append("")
    md.append("## Confirmation patient-pooled NRMSE (`t=9..46`)")
    md.append("")
    md.append("| Readout | Pooled NRMSE | Concat (diagnostic) |")
    md.append("|---|---|---|")
    md.append(f"| F408 | {f:.4f} | {f408['confirmation_concat']:.4f} |")
    md.append(f"| field AHL 20×10 | {fld:.4f} | {field['confirmation_concat']:.4f} |")
    md.append(f"| occupancy-masked RL surrogate | {masked['confirmation_pooled']:.4f} | {masked['confirmation_concat']:.4f} |")
    md.append(f"| reused Brownian Den | {brown['confirmation_pooled']:.4f} | {brown['confirmation_concat']:.4f} |")
    md.append(f"| reused silent F408 | {silent['confirmation_pooled']:.4f} | {silent['confirmation_concat']:.4f} |")
    md.append(f"| persistence aligned | {persist:.4f} | {ceilings['persistence_aligned_t9_46']['confirmation_concat_diagnostic']:.4f} |")
    md.append(f"| CRP-only delay-10 | {crp10:.4f} | {ceilings['crp_delay10']['confirmation_concat']:.4f} |")
    md.append(f"| 5-channel delay-10 (ceiling) | {ch5:.4f} | {ceilings['five_channel_delay10']['confirmation_concat']:.4f} |")
    md.append("")
    md.append(f"Occupancy ALIVE {occ_alive}/{len(occ_rows)}. Development DEAD {dead_dev}/12.")
    md.append(f"|F408 − field| = {delta:.4f}. F408 diagnostic `t=0..46` pooled {f408_diag['confirmation_pooled']:.4f}.")
    md.append("")
    md.append("Do not promote a new claim dish. Seeds 222/333 were not run.")
    md.append("")
    if not system_pass or not living_pass:
        md.append("REFUSE_RETUNE: do not raise Jmax or kinetics after these scores.")
        md.append("")
    (HERE / "results" / "E5_1_SCOUT.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print("COHORT_STATUS SYNTHETIC_COHORT")
    print("occupancy ALIVE", occ_alive, "/", len(occ_rows), "dev_DEAD", dead_dev)
    print("F408", f)
    print("field", fld)
    print("persistence_aligned", persist)
    print("crp_delay10", crp10)
    print("five_channel_delay10", ch5)
    print("system", system, "living_layer", living, "NO_STORY_MOVE", no_story)
    print("clinical FORBIDDEN")
    if not system_pass or not living_pass or no_story:
        refuse_retune("Report the scores. Do not retune.")


if __name__ == "__main__":
    main()
