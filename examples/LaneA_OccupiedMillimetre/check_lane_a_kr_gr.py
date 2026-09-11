#!/usr/bin/env python3
"""Lane A KR/GR + Jaeger MC on new T=840 drives. Not a NARMA-map rewrite."""

from __future__ import annotations

import csv
import hashlib
import json
import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from check_lane_a_carrier import (  # noqa: E402
    N_BINS,
    apply_std,
    nrmse,
    ridge_fit,
    ridge_predict,
    standardize_fit,
)
from check_lane_a_narma10 import delay_u10, load_u, sha256_u  # noqa: E402

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_KR_GR.md"
PROTOCOL_JSON = HERE / "configs" / "kr_gr_protocol.json"
KR_U_FILE = HERE / "input_u_kr_iid840.txt"
GR_U_FILE = HERE / "input_u_gr_const840.txt"
LANEA_JAVA = HERE.parents[1] / "src" / "bsim" / "lanea"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
NARMA_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_NARMA10_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"
MC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_MEMORY_CAPACITY_STANDING.md"
METRIC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_METRIC_LOCK_STANDING.md"

NUM_WINDOWS = 840
SAMPLES_PER_WINDOW = 16
WASHOUT = 40
TRAIN = 560
INNER_VAL = 88
TEST = 240
OCC_LO, OCC_HI = 600, 839
RANK_DELTA = 1e-6
SIGMA_MAX_FLOOR = 1e-15
ZERO_VAR = 1e-12
K_MAX = 20
DELTA_MC = 1.0
ALIVE_MIN = 0.05
KR_SEED = 20260831
KR_SHA = "1134e3ac2031cf33f699341090b3d186f9a776a81cfe072b18e9d0c9553ba1b5"
GR_SHA = "b7503301c6a2fabb9c4324c44d7da57ebea8632442d8681e862886643c4a8f72"
NARMA_SHA = "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b"
NARMA10B_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
ARMS = ("RL", "FIELD", "SILENT", "DELAY_U_10")


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("sweep", "0.928", "narma10b", "charc_sweep")):
        raise SystemExit("LANE_A_KR_GR refuses sweep / Narma10b / 0.928")


def lagged_u(u, k):
    return [0.0 if (n - k) < 0 else u[n - k] for n in range(len(u))]


def r2_from_nrmse(score: float) -> float:
    return 1.0 - score * score


def window_means(rows, prefix, n_windows: int):
    by_w = {}
    for row in rows:
        w = int(row["window"])
        by_w.setdefault(w, []).append([float(row[f"{prefix}_{i}"]) for i in range(N_BINS)])
    out = []
    for w in range(n_windows):
        block = by_w.get(w)
        if not block or len(block) != SAMPLES_PER_WINDOW:
            raise SystemExit(f"window {w} has {0 if not block else len(block)} samples")
        p = len(block[0])
        mean = [sum(block[t][j] for t in range(SAMPLES_PER_WINDOW)) / SAMPLES_PER_WINDOW for j in range(p)]
        out.append(mean)
    return out


def concat(a, b):
    return [x + y for x, y in zip(a, b)]


def numerical_rank(X):
    X = np.asarray(X, dtype=float)
    n_rows, n_cols = X.shape
    std = X.std(axis=0)
    std = np.where(std < ZERO_VAR, 1.0, std)
    Z = (X - X.mean(axis=0)) / std
    sigma = np.linalg.svd(Z, compute_uv=False)
    sigma = np.clip(sigma, 0.0, None)
    sig_max = float(sigma[0]) if len(sigma) else 0.0
    if sig_max < SIGMA_MAX_FLOOR:
        rank = 0
    else:
        rank = int(np.sum(sigma > RANK_DELTA * sig_max))
    ceiling = int(min(n_rows, n_cols))
    return {
        "rank": rank,
        "S": int(n_rows),
        "N": int(n_cols),
        "ceiling": ceiling,
        "row_capped": bool(n_rows <= n_cols),
        "feature_ceiling": bool(rank == n_cols and n_rows > n_cols),
        "sigma_max": sig_max,
        "delta": RANK_DELTA,
    }


def select_lambda(X, y):
    inner_end = TRAIN - INNER_VAL
    X_inner, y_inner = X[:inner_end], y[:inner_end]
    X_val, y_val = X[inner_end:TRAIN], y[inner_end:TRAIN]
    best_lam = 1e6
    best_score = float("inf")
    for lam in (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6):
        mu, sg = standardize_fit(X_inner)
        w = ridge_fit(apply_std(X_inner, mu, sg), y_inner, lam)
        pred = ridge_predict(apply_std(X_val, mu, sg), w)
        score = nrmse(y_val, pred)
        if score < best_score - 1e-15 or (abs(score - best_score) <= 1e-15 and lam > best_lam):
            best_score = score
            best_lam = lam
    return best_lam


def fit_eval(X_all, y_all):
    X = X_all[WASHOUT:]
    y = y_all[WASHOUT:]
    if len(X) != TRAIN + TEST:
        raise SystemExit(f"post-washout rows {len(X)} want {TRAIN + TEST}")
    lam = select_lambda(X, y)
    Xtr, ytr = X[:TRAIN], y[:TRAIN]
    Xte, yte = X[TRAIN:], y[TRAIN:]
    mu, sg = standardize_fit(Xtr)
    w = ridge_fit(apply_std(Xtr, mu, sg), ytr, lam)
    train_n = nrmse(ytr, ridge_predict(apply_std(Xtr, mu, sg), w))
    test_n = nrmse(yte, ridge_predict(apply_std(Xte, mu, sg), w))
    return {
        "lambda": lam,
        "train_nrmse": train_n,
        "test_nrmse": test_n,
        "train_r2": r2_from_nrmse(train_n),
        "test_r2": r2_from_nrmse(test_n),
    }


def mc_fading(r2_by_k) -> float:
    return float(sum(max(0.0, r2_by_k[k]) for k in range(1, K_MAX + 1)))


def load_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if "LANE_A_KR_GR" not in names or "u" not in names:
            raise SystemExit(f"{path.name} missing LANE_A_KR_GR/u")
        return list(reader)


def band_mean(rows, key, lo, hi):
    vals = [float(r[key]) for r in rows if lo <= int(r["window"]) <= hi]
    return sum(vals) / len(vals) if vals else float("nan")


def require_parents() -> None:
    occ = OCC_STANDING.read_text(encoding="utf-8")
    nar = NARMA_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    mc = MC_STANDING.read_text(encoding="utf-8")
    met = METRIC_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: PASS**" not in nar or "LANE_A_NARMA10" not in nar:
        raise SystemExit("NARMA standing must remain system PASS")
    if "**Status: FAIL**" not in car or "LANE_A_CARRIER" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    if "**Status: PASS**" not in mc or "LANE_A_MEMORY_CAPACITY" not in mc:
        raise SystemExit("memory-capacity standing must remain PASS")
    if "**Status: PASS**" not in met or "LANE_A_METRIC_LOCK" not in met:
        raise SystemExit("metric-lock standing must remain PASS")


def require_protocol_frozen() -> dict:
    proto = PROTOCOL.read_text(encoding="utf-8")
    js_text = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js_text:
        raise SystemExit("KR_GR PROTOCOL not frozen")
    js = json.loads(js_text)
    if js.get("kr_u", {}).get("sha256") != KR_SHA or js.get("gr_u", {}).get("sha256") != GR_SHA:
        raise SystemExit("PROTOCOL json u hash mismatch")
    if js.get("T_windows") != NUM_WINDOWS or js.get("S_post_washout") != 800:
        raise SystemExit("T/S freeze drifted")
    if js.get("rank_delta") != RANK_DELTA:
        raise SystemExit("rank delta drifted")
    if js.get("reuse_narma_u_as_kr_stream") is not False:
        raise SystemExit("must not reuse NARMA u")
    if js.get("lane_a_charc_sweep_started") is not False or js.get("sweep_K_tau_J") is not False:
        raise SystemExit("must not start a parameter sweep")
    if js.get("J_max") != 128000000.0 or js.get("raise_J_max") is not False:
        raise SystemExit("J_max drifted")
    if js.get("retune_Hill") is not False:
        raise SystemExit("Hill must not be retuned")
    if js.get("motility") is not False or js.get("growth") is not False or js.get("death") is not False:
        raise SystemExit("motility/growth/death must stay OFF")
    if js.get("warmup_s") != 0.0:
        raise SystemExit("warmup must stay 0")
    if js.get("restage_0895_as_mc") is not False:
        raise SystemExit("must not restage 0.895 as MC")
    return js


def refuse_math_random() -> None:
    hits = [p.name for p in sorted(LANEA_JAVA.glob("*.java")) if "Math.random()" in p.read_text(encoding="utf-8")]
    if hits:
        raise SystemExit(f"Math.random() forbidden in src/bsim/lanea: {hits}")


def verify_u(path: Path, expected_sha: str, n: int, regen=None):
    u = load_u(path)
    if len(u) != n:
        raise SystemExit(f"{path.name} has {len(u)}")
    digest = sha256_u(u)
    if digest != expected_sha:
        raise SystemExit(f"{path.name} sha256 drifted: {digest}")
    if digest in (NARMA_SHA, NARMA10B_SHA):
        raise SystemExit(f"{path.name} reused a NARMA hash")
    if regen is not None and sha256_u(regen) != digest:
        raise SystemExit(f"{path.name} Random({KR_SEED}) does not reproduce frozen u")
    return u


def csv_pulse_u(rows):
    out = []
    for w in range(NUM_WINDOWS):
        block = [r for r in rows if int(r["window"]) == w]
        on = [float(r["u"]) for r in block if float(r["u"]) > 0]
        out.append(on[0] if on else 0.0)
    return out


def write_not_scored(test_r: float) -> None:
    summary = {
        "gate": "LaneA_KR_GR",
        "status_label": "LANE_A_KR_GR",
        "LANE_A_KR_GR": "NOT_SCORED",
        "occupancy": "DEAD",
        "test_mean_R": test_r,
        "note": "occupancy DEAD on the i.i.d. u. No rank. Do not raise J_max.",
    }
    (RESULTS / "lane_a_kr_gr_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("LANE_A_KR_GR=NOT_SCORED occupancy DEAD on the i.i.d. u. Do not raise J_max. No rank.")


def main() -> None:
    refuse_forbidden()
    require_parents()
    require_protocol_frozen()
    refuse_math_random()

    rng = random.Random(KR_SEED)
    kr_regen = [rng.uniform(0.0, 0.5) for _ in range(NUM_WINDOWS)]
    u_kr = verify_u(KR_U_FILE, KR_SHA, NUM_WINDOWS, kr_regen)
    u_gr = verify_u(GR_U_FILE, GR_SHA, NUM_WINDOWS)
    if any(abs(v - 0.25) > 1e-15 for v in u_gr):
        raise SystemExit("GR u is not constant 0.25")

    iid_path = RESULTS / "java_LANE_A_KR_GR_IID_DRIVEN.csv"
    const_path = RESULTS / "java_LANE_A_KR_GR_CONST_DRIVEN.csv"
    silent_path = RESULTS / "java_LANE_A_KR_GR_SILENT.csv"
    if not iid_path.is_file():
        raise SystemExit("IID driven CSV missing; run ant lane-a-kr-gr. Do not raise J.")
    iid = load_rows(iid_path)
    if len(iid) != NUM_WINDOWS * SAMPLES_PER_WINDOW:
        raise SystemExit(f"IID samples {len(iid)}")
    if any(abs(a - b) > 1e-12 for a, b in zip(csv_pulse_u(iid), u_kr)):
        raise SystemExit("IID CSV pulse amplitudes do not match frozen KR u")

    test_r = band_mean(iid, "mean_R", OCC_LO, OCC_HI)
    train_r = band_mean(iid, "mean_R", 40, 599)
    occupancy_alive = test_r >= ALIVE_MIN
    if not occupancy_alive:
        write_not_scored(test_r)
        raise SystemExit(0)

    if not const_path.is_file() or not silent_path.is_file():
        raise SystemExit("CONST or SILENT CSV missing after ALIVE occupancy; do not invent a rank")

    const = load_rows(const_path)
    silent = load_rows(silent_path)
    if len(const) != NUM_WINDOWS * SAMPLES_PER_WINDOW:
        raise SystemExit(f"CONST samples {len(const)}")
    if len(silent) != NUM_WINDOWS * SAMPLES_PER_WINDOW:
        raise SystemExit(f"SILENT samples {len(silent)}")
    if any(abs(a - b) > 1e-12 for a, b in zip(csv_pulse_u(const), u_gr)):
        raise SystemExit("CONST CSV pulse amplitudes do not match frozen GR u")

    silent_r = band_mean(silent, "mean_R", OCC_LO, OCC_HI)
    const_test_r = band_mean(const, "mean_R", OCC_LO, OCC_HI)
    if silent_r >= ALIVE_MIN:
        raise SystemExit(f"silent occupied like driven mean_R={silent_r}")

    iid_ahl = window_means(iid, "AHL", NUM_WINDOWS)
    iid_rmap = window_means(iid, "R", NUM_WINDOWS)
    iid_lmap = window_means(iid, "L", NUM_WINDOWS)
    const_ahl = window_means(const, "AHL", NUM_WINDOWS)
    const_rmap = window_means(const, "R", NUM_WINDOWS)
    const_lmap = window_means(const, "L", NUM_WINDOWS)
    sil_rmap = window_means(silent, "R", NUM_WINDOWS)
    sil_lmap = window_means(silent, "L", NUM_WINDOWS)

    kr_features = {
        "RL": concat(iid_rmap, iid_lmap),
        "FIELD": iid_ahl,
        "SILENT": concat(sil_rmap, sil_lmap),
        "DELAY_U_10": delay_u10(u_kr),
    }
    gr_features = {
        "RL": concat(const_rmap, const_lmap),
        "FIELD": const_ahl,
        "SILENT": concat(sil_rmap, sil_lmap),
        "DELAY_U_10": delay_u10(u_gr),
    }

    kr = {}
    gr = {}
    for arm in ARMS:
        kr[arm] = numerical_rank(kr_features[arm][WASHOUT:])
        gr[arm] = numerical_rank(gr_features[arm][WASHOUT:])
        if arm == "RL" and kr[arm]["row_capped"]:
            raise SystemExit("RL KR is row-capped; do not sell that integer as reservoir rank")

    curves = {}
    for arm in ARMS:
        X = kr_features[arm]
        by_k = {}
        for k in range(0, K_MAX + 1):
            by_k[k] = fit_eval(X, lagged_u(u_kr, k))
        test_r2 = {k: by_k[k]["test_r2"] for k in range(0, K_MAX + 1)}
        curves[arm] = {
            "test_r2": [test_r2[k] for k in range(0, K_MAX + 1)],
            "MC_fading": mc_fading(test_r2),
            "r2_k0": test_r2[0],
        }

    mc_rl = curves["RL"]["MC_fading"]
    mc_field = curves["FIELD"]["MC_fading"]
    mc_silent = curves["SILENT"]["MC_fading"]
    mc_delay = curves["DELAY_U_10"]["MC_fading"]

    summary = {
        "gate": "LaneA_KR_GR",
        "status_label": "LANE_A_KR_GR",
        "LANE_A_KR_GR": "SCORED",
        "occupancy": "ALIVE",
        "test_mean_R": test_r,
        "train_mean_R": train_r,
        "const_test_mean_R": const_test_r,
        "silent_mean_R": silent_r,
        "occupancy_parent": "PASS",
        "narma_parent": "PASS",
        "carrier_parent": "FAIL",
        "mc_parent": "PASS",
        "metric_lock_parent": "PASS",
        "metric_lock_narma_maps_kr_gr_still_forbidden": True,
        "narma_0895_is_not_this_mc": True,
        "parent_mc_2069_is_not_this_mc": True,
        "kr_u_sha256": KR_SHA,
        "gr_u_sha256": GR_SHA,
        "T": NUM_WINDOWS,
        "S": 800,
        "rank_delta": RANK_DELTA,
        "KR": {arm: kr[arm] for arm in ARMS},
        "GR": {arm: gr[arm] for arm in ARMS},
        "MC_fading": {
            "RL": mc_rl,
            "FIELD": mc_field,
            "SILENT": mc_silent,
            "DELAY_U_10": mc_delay,
        },
        "MC_test_r2": {arm: curves[arm]["test_r2"] for arm in ARMS},
        "field_wins_delayed_u": bool(mc_field > mc_rl),
        "delay_u_10_is_kill": False,
        "lane_a_charc_sweep_started": False,
        "dambre_ipc": False,
    }
    (RESULTS / "lane_a_kr_gr_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print("LANE_A_KR_GR=SCORED")
    print(f"occupancy ALIVE test_mean_R={test_r:.6g} const_R={const_test_r:.6g} silent_R={silent_r:.6g}")
    for arm in ARMS:
        print(
            f"  {arm} KR={kr[arm]['rank']} GR={gr[arm]['rank']} "
            f"S={kr[arm]['S']} N={kr[arm]['N']} ceiling={kr[arm]['ceiling']} "
            f"row_capped={kr[arm]['row_capped']} feature_ceiling={kr[arm]['feature_ceiling']} "
            f"MC_fading={curves[arm]['MC_fading']:.6g}"
        )
    print("Field winning delayed-u is allowed. DELAY_U_10 is a ceiling, not a kill.")
    print("NOT CHARC. METRIC_LOCK NARMA-map ban stands. No sweep.")


if __name__ == "__main__":
    main()
