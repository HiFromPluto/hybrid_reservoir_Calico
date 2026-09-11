#!/usr/bin/env python3
"""E5.2 occupancy-first CRP scout checker.

Prints occupancy, system, living-layer, and fair-input readouts.
Refuses to retune K, n, tau, Jmax, clamp, mortality, flow, or layout.
Does not classify Y. SYNTHETIC_COHORT stays. E5.1 NOT_SCORED stays.
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
E51 = REPO / "examples" / "BSimReservoirPlanE51CRP"

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
LIVING_ALIVE_MIN = 10
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


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def abort(msg):
    print(f"ABORT: {msg}", file=sys.stderr)
    raise SystemExit(1)


def refuse_retune(reason):
    print("REFUSE_RETUNE: do not raise Jmax, K, n, tau, clamp, mortality,")
    print("flow, or layout. Do not add a map. Do not classify Y.")
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
        "in_cohort": mean_r >= OCC_DEAD_R,
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
        "development_pooled": float(np.mean(tr_scores)) if tr_scores else float("nan"),
        "confirmation_pooled": float(np.mean(te_scores)) if te_scores else float("nan"),
        "confirmation_concat": nrmse_pop(yte, pred_te) if len(yte) else float("nan"),
        "confirmation_per_patient": {
            str(int(p)): float(nrmse_pop(yte[pte == p], pred_te[pte == p]))
            for p in np.unique(pte)
        },
        "n_features": int(Xtr.shape[1]),
        "n_train_patients": int(len(np.unique(ptr))),
        "n_test_patients": int(len(np.unique(pte))),
    }


def score_readout(X_by_patient, train_patients, test_patients, t_lo=T_ALIGN, t_hi=T_LAST):
    if not train_patients or not test_patients:
        return {"lambda": None, "development_pooled": None, "confirmation_pooled": None, "note": "empty split"}
    Xtr, ytr, ptr = stack_pairs(X_by_patient, train_patients, t_lo, t_hi)
    Xte, yte, pte = stack_pairs(X_by_patient, test_patients, t_lo, t_hi)
    lam, loo, grid = loo_lambda(Xtr, ytr, ptr, train_patients)
    metrics = fit_confirm(Xtr, ytr, ptr, Xte, yte, pte, lam)
    metrics["loo_pooled_development"] = loo
    metrics["loo_grid"] = grid
    return metrics


def driven_dir(p):
    return HERE / "results" / f"e52_driven_p{p:04d}_seed111"


def reused_null_dir(kind):
    e52 = HERE / "results" / f"e52_{kind}_seed111"
    if (e52 / "voxels.csv").exists():
        return e52
    return E51 / "results" / f"e51_{kind}_seed111"


def java_ok():
    java = (HERE / "BSimReservoirPlanE52CRP.java").read_text(encoding="utf-8")
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
    if "AHL_SOURCE_RATE = 128000000" not in java:
        return False
    if "package BSimReservoirPlanE52CRP;" not in voxel:
        return False
    return True


def load_freeze():
    path = HERE / "results" / "e52_encoding_freeze.json"
    if not path.exists():
        abort("ENCODING_FREEZE missing; run generate_e52_encoding.py first")
    return json.loads(path.read_text(encoding="utf-8"))


def verify_hashes(freeze):
    if sha256_text(",".join(str(i) for i in DEV)) != DEV_SHA:
        abort("development scout hash mismatch")
    if sha256_text(",".join(str(i) for i in CONF)) != CONF_SHA:
        abort("confirmation scout hash mismatch")
    recorded = json.loads((HERE / "results" / "e52_input_hashes.json").read_text(encoding="utf-8"))
    for spec in freeze["maps"]:
        map_id = spec["id"] if isinstance(spec, dict) and "id" in spec else spec.get("id")
        recs = recorded["maps"][map_id]
        for p in DEV:
            path = HERE / "input" / "maps" / map_id / f"input_ahl_p{p:04d}.txt"
            got = sha256_file(path)
            exp = recs[f"input_ahl_{map_id}_p{p:04d}.txt"]["file_sha256"]
            if got != exp:
                abort(f"map u hash drift {path}")
    if freeze["status"] == "ENCODING_FROZEN":
        for p in DEV + CONF:
            path = HERE / "input" / f"input_ahl_p{p:04d}.txt"
            got = sha256_file(path)
            exp = recorded["files"][path.name]["file_sha256"]
            if got != exp:
                abort(f"winner u hash drift {path.name}")
    return recorded


def file_manifest():
    items = {
        "BSimReservoirPlanE52CRP.java": HERE / "BSimReservoirPlanE52CRP.java",
        "VoxelAnalyzer.java": HERE / "VoxelAnalyzer.java",
        "check_e52.py": HERE / "check_e52.py",
        "generate_e52_encoding.py": HERE / "generate_e52_encoding.py",
        "PROTOCOL.md": HERE / "PROTOCOL.md",
        "ENCODING_FREEZE.md": HERE / "results" / "ENCODING_FREEZE.md",
        "e52_occupancy.csv": HERE / "results" / "e52_occupancy.csv",
    }
    out = {}
    for name, path in items.items():
        if path.exists():
            out[name] = sha256_file(path)
    cfg = HERE / "sim_config_e52.properties"
    if cfg.exists():
        out[cfg.name] = sha256_file(cfg)
    return out


def write_not_scored(reason, freeze, occ_rows, ceilings, extra_lines=()):
    persist = ceilings["persistence_aligned_t9_46"]["confirmation_pooled"]
    persist_all = ceilings["persistence_t0_46"]["confirmation_pooled"]
    crp10 = ceilings["crp_delay10"]["confirmation_pooled"]
    ch5 = ceilings["five_channel_delay10"]["confirmation_pooled"]
    payload = {
        "cohort_status": "SYNTHETIC_COHORT",
        "clinical": False,
        "e51_status": "NOT_SCORED_UNCHANGED",
        "skip_confirmation": True,
        "reason": reason,
        "encoding": freeze.get("status"),
        "winner": freeze.get("winner"),
        "occupancy": occ_rows,
        "F408": "NOT_SCORED",
        "field": "NOT_SCORED",
        "ceilings": {
            "persistence_t0_46_conf": persist_all,
            "persistence_aligned_conf": persist,
            "crp_delay10_conf": crp10,
            "five_channel_delay10_conf": ch5,
        },
        "claims": {
            "occupancy": reason,
            "system": "NOT_SCORED",
            "living_layer": "NOT_SCORED",
            "ceiling": "NOT_REQUIRED",
            "fair_input": "NOT_SCORED",
            "clinical": "FORBIDDEN",
            "task_a": "UNCHANGED",
            "track_b": "UNCHANGED",
            "e51": "NOT_SCORED_UNCHANGED",
            "claim_dish": "NOT_REPLACED",
            "cohort_status": "SYNTHETIC_COHORT",
        },
        "manifest_sha256": file_manifest(),
    }
    (HERE / "results" / "e52_scout.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    with (HERE / "results" / "e52_taskb.csv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=["readout", "split", "nrmse_patient_pooled", "note"])
        w.writeheader()
        w.writerow({"readout": "F408", "split": "confirmation", "nrmse_patient_pooled": "", "note": f"NOT_SCORED {reason}"})
        w.writerow({"readout": "field", "split": "confirmation", "nrmse_patient_pooled": "", "note": f"NOT_SCORED {reason}"})
        w.writerow({"readout": "persistence_aligned", "split": "confirmation", "nrmse_patient_pooled": f"{persist:.10g}", "note": "Module 0 ceiling"})
        w.writerow({"readout": "crp_delay10", "split": "confirmation", "nrmse_patient_pooled": f"{crp10:.10g}", "note": "Module 0 ceiling"})
        w.writerow({"readout": "five_channel_delay10", "split": "confirmation", "nrmse_patient_pooled": f"{ch5:.10g}", "note": "Module 0 ceiling"})
    lines = [
        "# E5.2 occupancy-first one-AHL CRP scout",
        "",
        "Seed 111 only. Task B only. `COHORT_STATUS: SYNTHETIC_COHORT`.",
        "Not a clinical dataset. Track B FAIL is not rewritten.",
        "E5.1 `NOT_SCORED` is unchanged. Claim dish is not replaced.",
        "",
        f"**{reason}**",
        "",
        "## Claims",
        "",
        "| Claim | Result |",
        "|---|---|",
        f"| Occupancy | **{reason}** |",
        "| System (F408 vs Brownian and silent) | **NOT_SCORED** |",
        "| Living-layer (F408 vs field) | **NOT_SCORED** |",
        "| Ceiling (beat 5-ch delay-10) | **NOT_REQUIRED** |",
        "| Fair input | Module 0 ceilings only |",
        "| Clinical | **FORBIDDEN** |",
        "| Task A / Track B / E5.1 | unchanged |",
        "| Claim dish | not replaced |",
        "",
        f"Winner / freeze: `{freeze.get('winner')}` / `{freeze.get('status')}`.",
        "",
        "## Module 0 confirmation ceilings (12+8, no living score)",
        "",
        "| Baseline | Conf pooled NRMSE |",
        "|---|---|",
        f"| persistence `t=0..46` | {persist_all:.4f} |",
        f"| persistence aligned `t=9..46` | {persist:.4f} |",
        f"| CRP-only delay-10 | {crp10:.4f} |",
        f"| 5-channel delay-10 (ceiling) | {ch5:.4f} |",
        "",
        "REFUSE_RETUNE: do not raise `Jmax` or invent Map 5.",
        "",
        *extra_lines,
    ]
    (HERE / "results" / "E5_2_SCOUT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("COHORT_STATUS SYNTHETIC_COHORT")
    print("E5.1 NOT_SCORED unchanged")
    print("F408 NOT_SCORED")
    print(reason)


def reuse_nulls_ok():
    notes = {}
    for kind in ("silent", "brownian"):
        d = reused_null_dir(kind)
        vox = d / "voxels.csv"
        summary = d / "window_summary.csv"
        if not vox.exists() or not summary.exists():
            abort(f"reused {kind} CSVs missing at {d}")
        csvs = validate_csv(vox, EXPECTED_AUX, check_last=True)
        sums = validate_csv(summary, NUM_WINDOWS)
        if not (csvs["row_count_pass"] and csvs["last_sample_pass"] and sums["row_count_pass"]):
            abort(f"reused {kind} last-sample / 48/768 mismatch: {d}")
        notes[kind] = {
            "path": str(vox),
            "sha256": sha256_file(vox),
            "rows": csvs["row_count"],
            "last_sample_pass": True,
            "reused_from_e51": "BSimReservoirPlanE51CRP" in str(vox),
        }
    return notes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke-only", action="store_true")
    args = parser.parse_args()
    freeze = load_freeze()
    verify_hashes(freeze)
    ceilings = json.loads((HERE / "results" / "e52_input_ceilings.json").read_text(encoding="utf-8"))
    if not java_ok():
        abort("Java failed frozen-dish / forbidden-token checks")
    loadmat(str(MAT))

    if freeze["status"] != "ENCODING_FROZEN":
        write_not_scored("ENCODING_FAIL", freeze, [], ceilings)
        refuse_retune("Encoding screen did not freeze a map. Do not raise Jmax.")
        raise SystemExit(0)

    if args.smoke_only:
        smoke = driven_dir(DEV[0])
        voxels_path = smoke / "voxels.csv"
        if not voxels_path.exists():
            abort(f"smoke voxels missing: {voxels_path}")
        csvs = validate_csv(voxels_path, EXPECTED_AUX, check_last=True)
        summary = validate_csv(smoke / "window_summary.csv", NUM_WINDOWS)
        voxels = load_voxel_arrays(voxels_path)
        u = load_u(DEV[0])
        occ = occupancy_from_voxels(voxels, u)
        print("SMOKE patient", DEV[0])
        print("map", freeze["winner"])
        print("u[0]", float(u[0]))
        print("CENTER FLOW=0")
        print("last_sample", LAST_SAMPLE_PREFIX, csvs.get("last_sample_pass"))
        print("summary_rows", summary["row_count"], "voxel_rows", csvs["row_count"])
        print("AHL finite", occ["ahl_finite"], "nonneg", occ["ahl_nonneg"], "mean_AHL", occ["mean_AHL"])
        print("occupancy", occ["occupancy"], "mean_R", occ["mean_R"], "r(R,u)", occ["r_meanR_u"])
        print("COHORT_STATUS SYNTHETIC_COHORT")
        return

    missing = []
    for p in DEV:
        if not (driven_dir(p) / "voxels.csv").exists():
            missing.append(f"driven p{p:04d}")
    if missing:
        print("Encoding freeze OK. BSim voxels missing:", ", ".join(missing))
        print("COHORT_STATUS SYNTHETIC_COHORT")
        print("Run isolated BSim before Module 2 scores.")
        raise SystemExit(2)

    reuse_nulls_ok()
    occ_rows = []
    X408, Xfield, Xmask = {}, {}, {}

    def load_driven(p, split):
        d = driven_dir(p)
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
        occ_rows.append(load_driven(p, "development"))
    alive_dev = sum(1 for r in occ_rows if r["occupancy"] == "ALIVE")
    if alive_dev < LIVING_ALIVE_MIN:
        extra = [
            "## Development occupancy (living, seed 111)",
            "",
            "| Patient | Occupancy | mean_R | r(mean_R, u) | mean_AHL | in_cohort |",
            "|---|---|---|---|---|---|",
        ]
        for row in occ_rows:
            extra.append(
                f"| {row['patient']} | {row['occupancy']} | {row['mean_R']:.4f} | "
                f"{row['r_meanR_u']:.3f} | {row['mean_AHL']:.4f} | {row['in_cohort']} |"
            )
        extra.extend(["", f"Living development `n_ALIVE={alive_dev}/12`. Confirmation not started.", ""])
        write_not_scored(
            f"NOT_SCORED living development n_ALIVE={alive_dev}/12",
            freeze,
            occ_rows,
            ceilings,
            extra,
        )
        with (HERE / "results" / "e52_occupancy_living.csv").open("w", newline="", encoding="utf-8") as handle:
            w = csv.DictWriter(
                handle,
                fieldnames=["map_id", "patient", "split", "occupancy", "mean_R", "r_meanR_u", "mean_AHL", "in_cohort", "source"],
            )
            w.writeheader()
            for row in occ_rows:
                w.writerow({
                    "map_id": freeze["winner"],
                    "patient": row["patient"],
                    "split": "development",
                    "occupancy": row["occupancy"],
                    "mean_R": f"{row['mean_R']:.10g}",
                    "r_meanR_u": f"{row['r_meanR_u']:.10g}",
                    "mean_AHL": f"{row['mean_AHL']:.10g}",
                    "in_cohort": row["in_cohort"],
                    "source": "living_BSim",
                })
        refuse_retune("Living occupancy missed n_ALIVE>=10. Do not switch maps.")
        raise SystemExit(0)

    missing_conf = [p for p in CONF if not (driven_dir(p) / "voxels.csv").exists()]
    if missing_conf:
        print("Living development n_ALIVE", alive_dev, "/12. Confirmation voxels missing:", missing_conf)
        print("COHORT_STATUS SYNTHETIC_COHORT")
        raise SystemExit(2)

    for p in CONF:
        occ_rows.append(load_driven(p, "confirmation"))

    silent_bio = read_window_matrix(reused_null_dir("silent") / "voxels.csv", BIOLOGY_MEAN, BIOLOGY_LAST)
    brown_den = read_window_matrix(reused_null_dir("brownian") / "voxels.csv", BROWN_MEAN, ())
    Xsilent = {p: silent_bio["X"] for p in DEV + CONF}
    Xbrown = {p: brown_den["X"] for p in DEV + CONF}

    f408 = score_readout(X408, DEV, CONF)
    field = score_readout(Xfield, DEV, CONF)
    masked = score_readout(Xmask, DEV, CONF)
    silent = score_readout(Xsilent, DEV, CONF)
    brown = score_readout(Xbrown, DEV, CONF)
    f408_diag = score_readout(X408, DEV, CONF, T0, T_LAST)

    occ_before_f408 = [{k: r[k] for k in ("patient", "split", "occupancy", "mean_R", "r_meanR_u", "in_cohort", "mean_AHL")} for r in occ_rows]
    dev_in = [r["patient"] for r in occ_rows if r["split"] == "development" and r["in_cohort"]]
    conf_in = [r["patient"] for r in occ_rows if r["split"] == "confirmation" and r["in_cohort"]]
    f408_occ = score_readout(X408, dev_in, conf_in) if dev_in and conf_in else {"note": "occupancy-conditioned empty"}

    persist = ceilings["persistence_aligned_t9_46"]["confirmation_pooled"]
    persist_all = ceilings["persistence_t0_46"]["confirmation_pooled"]
    crp10 = ceilings["crp_delay10"]["confirmation_pooled"]
    ch5 = ceilings["five_channel_delay10"]["confirmation_pooled"]
    f = f408["confirmation_pooled"]
    fld = field["confirmation_pooled"]
    system_pass = f < brown["confirmation_pooled"] and f < silent["confirmation_pooled"]
    living_delta = f - fld
    alive_conf = sum(1 for r in occ_rows if r["split"] == "confirmation" and r["occupancy"] == "ALIVE")
    system = "PASS" if system_pass else "FAIL"
    living = "REPORTED"
    claims = {
        "occupancy": f"ALIVE {alive_dev}/12 development, {alive_conf}/8 confirmation",
        "system": system,
        "living_layer": {"F408": f, "field": fld, "delta": living_delta},
        "ceiling": "NOT_REQUIRED",
        "fair_input": {"vs_persistence_aligned": f - persist, "vs_crp_delay10": f - crp10},
        "clinical": "FORBIDDEN",
        "task_a": "UNCHANGED",
        "track_b": "UNCHANGED",
        "e51": "NOT_SCORED_UNCHANGED",
        "claim_dish": "NOT_REPLACED",
        "cohort_status": "SYNTHETIC_COHORT",
    }
    payload = {
        "cohort_status": "SYNTHETIC_COHORT",
        "clinical": False,
        "last_sample_prefix": LAST_SAMPLE_PREFIX,
        "winner": freeze["winner"],
        "occupancy_before_F408": occ_before_f408,
        "F408": f408,
        "field": field,
        "masked_surrogate": masked,
        "silent": silent,
        "brownian": brown,
        "F408_t0_46_diagnostic": f408_diag,
        "F408_occupancy_conditioned": f408_occ,
        "occupancy_conditioned_patients": {"development": dev_in, "confirmation": conf_in},
        "ceilings": {
            "persistence_t0_46_conf": persist_all,
            "persistence_aligned_conf": persist,
            "crp_delay10_conf": crp10,
            "five_channel_delay10_conf": ch5,
        },
        "claims": claims,
        "manifest_sha256": file_manifest(),
        "null_reuse": reuse_nulls_ok(),
    }
    (HERE / "results" / "e52_scout.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    with (HERE / "results" / "e52_taskb.csv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=["readout", "split", "nrmse_patient_pooled", "nrmse_concat", "lambda", "n_features", "note"])
        w.writeheader()
        for name, row, note in (
            ("F408", f408, "primary all confirmation"),
            ("field", field, "primary"),
            ("masked_surrogate", masked, "primary"),
            ("silent", silent, "reused E5.1"),
            ("brownian", brown, "reused E5.1"),
        ):
            w.writerow({
                "readout": name,
                "split": "confirmation",
                "nrmse_patient_pooled": f"{row['confirmation_pooled']:.10g}",
                "nrmse_concat": f"{row['confirmation_concat']:.10g}",
                "lambda": f"{row['lambda']:g}",
                "n_features": row["n_features"],
                "note": note,
            })
        if isinstance(f408_occ, dict) and f408_occ.get("confirmation_pooled") is not None:
            w.writerow({
                "readout": "F408_occupancy_conditioned",
                "split": "confirmation_in_cohort",
                "nrmse_patient_pooled": f"{f408_occ['confirmation_pooled']:.10g}",
                "nrmse_concat": f"{f408_occ.get('confirmation_concat', float('nan')):.10g}",
                "lambda": f"{f408_occ['lambda']:g}",
                "n_features": f408_occ["n_features"],
                "note": "secondary; occupancy before F408; not the only number",
            })
        w.writerow({
            "readout": "persistence_aligned",
            "split": "confirmation",
            "nrmse_patient_pooled": f"{persist:.10g}",
            "nrmse_concat": f"{ceilings['persistence_aligned_t9_46']['confirmation_concat_diagnostic']:.10g}",
            "lambda": "",
            "n_features": 1,
            "note": "fair input",
        })
        w.writerow({
            "readout": "crp_delay10",
            "split": "confirmation",
            "nrmse_patient_pooled": f"{crp10:.10g}",
            "nrmse_concat": f"{ceilings['crp_delay10']['confirmation_concat']:.10g}",
            "lambda": f"{ceilings['crp_delay10']['lambda']:g}",
            "n_features": 10,
            "note": "fair input",
        })
        w.writerow({
            "readout": "five_channel_delay10",
            "split": "confirmation",
            "nrmse_patient_pooled": f"{ch5:.10g}",
            "nrmse_concat": f"{ceilings['five_channel_delay10']['confirmation_concat']:.10g}",
            "lambda": f"{ceilings['five_channel_delay10']['lambda']:g}",
            "n_features": 50,
            "note": "ceiling, not a living-layer bar",
        })

    living_occ_path = HERE / "results" / "e52_occupancy_living.csv"
    with living_occ_path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=["map_id", "patient", "split", "occupancy", "mean_R", "r_meanR_u", "mean_AHL", "in_cohort"])
        w.writeheader()
        for row in occ_rows:
            w.writerow({
                "map_id": freeze["winner"],
                "patient": row["patient"],
                "split": row["split"],
                "occupancy": row["occupancy"],
                "mean_R": f"{row['mean_R']:.10g}",
                "r_meanR_u": f"{row['r_meanR_u']:.10g}",
                "mean_AHL": f"{row['mean_AHL']:.10g}",
                "in_cohort": row["in_cohort"],
            })

    md = [
        "# E5.2 occupancy-first one-AHL CRP scout",
        "",
        "Seed 111 only. Task B only. `COHORT_STATUS: SYNTHETIC_COHORT`.",
        "Not a clinical dataset. Track B FAIL is not rewritten.",
        "E5.1 `NOT_SCORED` is unchanged. K, n, tau, Jmax were not retuned.",
        f"Frozen map: `{freeze['winner']}`.",
        "",
        "## Claims",
        "",
        "| Claim | Result |",
        "|---|---|",
        f"| Occupancy | **ALIVE {alive_dev}/12 development, {alive_conf}/8 confirmation** |",
        f"| System (F408 vs Brownian and silent) | **{system}** |",
        f"| Living-layer (F408 vs field) | F408 `{f:.4f}` vs field `{fld:.4f}` (Δ `{living_delta:.4f}`) |",
        "| Ceiling (beat 5-ch delay-10) | **NOT_REQUIRED** |",
        "| Fair input (vs persistence and CRP delay-10) | report only |",
        "| Clinical | **FORBIDDEN** |",
        "| Task A / Track B / E5.1 | unchanged |",
        "| Claim dish | not replaced |",
        "",
        "## Confirmation patient-pooled NRMSE (`t=9..46`, primary = all 8)",
        "",
        "| Readout | Pooled NRMSE | Concat (diagnostic) |",
        "|---|---|---|",
        f"| F408 | {f:.4f} | {f408['confirmation_concat']:.4f} |",
        f"| field AHL 20×10 | {fld:.4f} | {field['confirmation_concat']:.4f} |",
        f"| occupancy-masked RL surrogate | {masked['confirmation_pooled']:.4f} | {masked['confirmation_concat']:.4f} |",
        f"| reused Brownian Den | {brown['confirmation_pooled']:.4f} | {brown['confirmation_concat']:.4f} |",
        f"| reused silent F408 | {silent['confirmation_pooled']:.4f} | {silent['confirmation_concat']:.4f} |",
        f"| persistence aligned | {persist:.4f} | {ceilings['persistence_aligned_t9_46']['confirmation_concat_diagnostic']:.4f} |",
        f"| CRP-only delay-10 | {crp10:.4f} | {ceilings['crp_delay10']['confirmation_concat']:.4f} |",
        f"| 5-channel delay-10 (ceiling) | {ch5:.4f} | {ceilings['five_channel_delay10']['confirmation_concat']:.4f} |",
        "",
    ]
    if isinstance(f408_occ, dict) and f408_occ.get("confirmation_pooled") is not None:
        md.extend([
            "## Occupancy-conditioned (secondary, not the only number)",
            "",
            f"In-cohort = living `mean_R ≥ 0.05`, computed before F408. "
            f"Development {dev_in}; confirmation {conf_in}.",
            f"F408 occupancy-conditioned confirmation pooled `{f408_occ['confirmation_pooled']:.4f}`.",
            "",
        ])
    md.extend([
        f"F408 diagnostic `t=0..46` pooled {f408_diag['confirmation_pooled']:.4f}.",
        "",
        "## Living occupancy",
        "",
        "| Patient | Split | Occupancy | mean_R | r(mean_R, u) | mean_AHL | in_cohort |",
        "|---|---|---|---|---|---|---|",
    ])
    for row in occ_rows:
        md.append(
            f"| {row['patient']} | {row['split']} | {row['occupancy']} | {row['mean_R']:.4f} | "
            f"{row['r_meanR_u']:.3f} | {row['mean_AHL']:.4f} | {row['in_cohort']} |"
        )
    md.extend(["", "Silent / Brownian seed 111 reused from E5.1 after last-sample",
               "`47;15;299.95` and 48/768/768 matched. Driven CSVs were not reused.",
               "",
               "Occupancy-conditioned secondary equals primary here if every",
               "scored patient has living `mean_R ≥ 0.05`. It is labelled, not the",
               "only number.",
               "",
               "Do not promote a new claim dish. Seeds 222/333 were not run.", ""])
    if not system_pass:
        md.append("REFUSE_RETUNE: do not raise Jmax or kinetics after these scores.")
        md.append("")
    (HERE / "results" / "E5_2_SCOUT.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print("COHORT_STATUS SYNTHETIC_COHORT")
    print("E5.1 NOT_SCORED unchanged")
    print("map", freeze["winner"])
    print("occupancy ALIVE development", alive_dev, "/12 confirmation", alive_conf, "/8")
    print("F408", f)
    print("field", fld)
    print("persistence_aligned", persist)
    print("crp_delay10", crp10)
    print("five_channel_delay10", ch5)
    print("system", system, "clinical FORBIDDEN")
    if not system_pass:
        refuse_retune("Report the scores. Do not retune.")


if __name__ == "__main__":
    main()
