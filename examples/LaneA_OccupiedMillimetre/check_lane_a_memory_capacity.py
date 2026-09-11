#!/usr/bin/env python3
"""Lane A Jaeger delay capacity. Reads existing NARMA CSVs only. Not CHARC."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

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
from check_lane_a_narma10 import (  # noqa: E402
    INNER_VAL,
    LAMBDA_GRID,
    NUM_WINDOWS,
    SAMPLES_PER_WINDOW,
    TRAIN,
    WASHOUT,
    delay_u10,
    load_u,
    sha256_u,
    window_means,
    concat,
)

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_MEMORY_CAPACITY.md"
PROTOCOL_JSON = HERE / "configs" / "memory_capacity_protocol.json"
U_FILE = HERE / "input_u_narma200.txt"
LANEA_JAVA = HERE.parents[1] / "src" / "bsim" / "lanea"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
NARMA_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_NARMA10_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"
AUDIT_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_AUDIT_STANDING.md"
SEED_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_SEED_REPLICATE_STANDING.md"

U_SHA = "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b"
NARMA10B_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
K_MAX = 20
DELTA_MC = 1.0
ALIVE_MIN = 0.05
K_HALF_THRESH = 0.5
ARMS = ("RL", "FIELD", "SILENT_RL", "DELAY_U_10")


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("charc", "ipc", "erickson", "lane-d", "laned")):
        raise SystemExit("Lane A memory capacity refuses CHARC/IPC/Lane D/Erickson")


def lagged_u(u, k):
    out = []
    for n in range(len(u)):
        j = n - k
        out.append(0.0 if j < 0 else u[j])
    return out


def r2_from_nrmse(score: float) -> float:
    return 1.0 - score * score


def select_lambda(X, y):
    inner_end = TRAIN - INNER_VAL
    X_inner = X[:inner_end]
    y_inner = y[:inner_end]
    X_val = X[inner_end:TRAIN]
    y_val = y[inner_end:TRAIN]
    best_lam = LAMBDA_GRID[-1]
    best_score = float("inf")
    for lam in LAMBDA_GRID:
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
    if len(X) != TRAIN + 50:
        raise SystemExit(f"post-washout rows {len(X)}")
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


def k_half(r2_by_k):
    for k in range(0, K_MAX + 1):
        if r2_by_k[k] < K_HALF_THRESH:
            return k
    return "never"


def load_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if "LANE_A_NARMA10" not in names or "u" not in names:
            raise SystemExit(f"{path.name} missing LANE_A_NARMA10/u")
        return list(reader)


def band_mean(rows, key, lo, hi):
    vals = [float(r[key]) for r in rows if lo <= int(r["window"]) <= hi]
    return sum(vals) / len(vals) if vals else float("nan")


def require_parents() -> None:
    occ = OCC_STANDING.read_text(encoding="utf-8")
    nar = NARMA_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    aud = AUDIT_STANDING.read_text(encoding="utf-8")
    seed = SEED_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: PASS**" not in nar or "LANE_A_NARMA10" not in nar:
        raise SystemExit("NARMA standing must remain system PASS")
    if "**Status: FAIL**" not in car or "LANE_A_CARRIER" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    if "**Status: SCOPE_NOTE**" not in aud:
        raise SystemExit("audit standing must remain SCOPE_NOTE")
    if "**Status: AXIS_HOLDS**" not in seed:
        raise SystemExit("seed-replicate standing must remain AXIS_HOLDS")


def require_protocol_frozen() -> dict:
    proto = PROTOCOL.read_text(encoding="utf-8")
    js_text = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js_text:
        raise SystemExit("MEMORY_CAPACITY PROTOCOL not frozen")
    js = json.loads(js_text)
    if js.get("u_sha256") != U_SHA:
        raise SystemExit("PROTOCOL json u hash mismatch")
    if js.get("delta_MC") != DELTA_MC:
        raise SystemExit("delta_MC drifted")
    if js.get("charc") is not False or js.get("dambre_ipc") is not False:
        raise SystemExit("must freeze charc/ipc false")
    if js.get("double_washout") is not False:
        raise SystemExit("Paper 1 double washout is forbidden")
    if js.get("hunt_seeds_202_303") is not False:
        raise SystemExit("must not hunt 202/303 after scores")
    if js.get("J_max") != 128000000.0:
        raise SystemExit("J_max drifted")
    if js.get("motility") is not False or js.get("growth") is not False or js.get("death") is not False:
        raise SystemExit("motility/growth/death must stay OFF")
    if js.get("warmup_s") != 0.0:
        raise SystemExit("warmup must stay 0")
    return js


def refuse_math_random() -> None:
    hits = []
    for path in sorted(LANEA_JAVA.glob("*.java")):
        text = path.read_text(encoding="utf-8")
        if "Math.random()" in text:
            hits.append(path.name)
    if hits:
        raise SystemExit(f"Math.random() forbidden in src/bsim/lanea: {hits}")


def main() -> None:
    refuse_forbidden()
    require_parents()
    require_protocol_frozen()
    refuse_math_random()

    u = load_u(U_FILE)
    if len(u) != NUM_WINDOWS:
        raise SystemExit(f"u has {len(u)}")
    digest = sha256_u(u)
    if digest != U_SHA:
        raise SystemExit(f"u sha256 drifted: {digest}")
    if digest == NARMA10B_SHA:
        raise SystemExit("copied Narma10b u hash")

    driven_path = RESULTS / "java_LANE_A_NARMA10_DRIVEN.csv"
    silent_path = RESULTS / "java_LANE_A_NARMA10_SILENT.csv"
    if not driven_path.is_file() or not silent_path.is_file():
        raise SystemExit(
            "required NARMA CSVs missing; rerun only the NARMA target, do not raise J"
        )

    driven = load_rows(driven_path)
    silent = load_rows(silent_path)
    if len(driven) != NUM_WINDOWS * SAMPLES_PER_WINDOW:
        raise SystemExit(f"driven samples {len(driven)}")
    if len(silent) != NUM_WINDOWS * SAMPLES_PER_WINDOW:
        raise SystemExit(f"silent samples {len(silent)}")

    csv_u = []
    for w in range(NUM_WINDOWS):
        block = [r for r in driven if int(r["window"]) == w]
        on = [float(r["u"]) for r in block if float(r["u"]) > 0]
        csv_u.append(on[0] if on else 0.0)
    if any(abs(a - b) > 1e-12 for a, b in zip(csv_u, u)):
        raise SystemExit("CSV pulse amplitudes do not match frozen u")

    test_r = band_mean(driven, "mean_R", 150, 199)
    train_r = band_mean(driven, "mean_R", 40, 149)
    silent_r = band_mean(silent, "mean_R", 150, 199)
    occupancy_alive = test_r >= ALIVE_MIN

    if not occupancy_alive:
        summary = {
            "gate": "LaneA_MEMORY_CAPACITY",
            "status_label": "LANE_A_MEMORY_CAPACITY",
            "LANE_A_MEMORY_CAPACITY": "NOT_SCORED",
            "system": "NOT_SCORED",
            "test_mean_R": test_r,
            "occupancy": "DEAD",
        }
        (RESULTS / "lane_a_memory_capacity_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )
        print("LANE_A_MEMORY_CAPACITY=NOT_SCORED occupancy DEAD on this u. Do not raise J_max.")
        raise SystemExit(0)

    if silent_r >= ALIVE_MIN:
        raise SystemExit(f"silent occupied like driven mean_R={silent_r}")

    ahl = window_means(driven, "AHL")
    rmap = window_means(driven, "R")
    lmap = window_means(driven, "L")
    s_r = window_means(silent, "R")
    s_l = window_means(silent, "L")
    features = {
        "RL": concat(rmap, lmap),
        "FIELD": ahl,
        "SILENT_RL": concat(s_r, s_l),
        "DELAY_U_10": delay_u10(u),
    }

    curves = {}
    for arm in ARMS:
        X = features[arm]
        by_k = {}
        for k in range(0, K_MAX + 1):
            y = lagged_u(u, k)
            by_k[k] = fit_eval(X, y)
        test_r2 = {k: by_k[k]["test_r2"] for k in range(0, K_MAX + 1)}
        curves[arm] = {
            "by_k": {str(k): by_k[k] for k in range(0, K_MAX + 1)},
            "test_r2": [test_r2[k] for k in range(0, K_MAX + 1)],
            "MC_fading": mc_fading(test_r2),
            "k_half": k_half(test_r2),
            "r2_k0": test_r2[0],
        }

    mc_rl = curves["RL"]["MC_fading"]
    mc_silent = curves["SILENT_RL"]["MC_fading"]
    mc_field = curves["FIELD"]["MC_fading"]
    mc_delay = curves["DELAY_U_10"]["MC_fading"]
    system_pass = mc_rl >= mc_silent + DELTA_MC
    verdict = "PASS" if system_pass else "FAIL"
    field_wins_sum = mc_field > mc_rl
    delay_wins_sum = mc_delay > mc_rl

    summary = {
        "gate": "LaneA_MEMORY_CAPACITY",
        "status_label": "LANE_A_MEMORY_CAPACITY",
        "LANE_A_MEMORY_CAPACITY": verdict,
        "system": verdict,
        "metric": "jaeger_delay_capacity_windowed_maps",
        "occupancy_parent": "PASS",
        "narma_parent": "PASS",
        "carrier_parent": "FAIL",
        "audit_parent": "SCOPE_NOTE",
        "seed_replicate_parent": "AXIS_HOLDS",
        "charc": False,
        "dambre_ipc": False,
        "kernel_rank": False,
        "paper1_mc_rewrite": False,
        "narma_nrmse_is_not_mc": True,
        "double_washout": False,
        "u_sha256": digest,
        "delta_MC": DELTA_MC,
        "k_min": 0,
        "k_max": K_MAX,
        "test_mean_R": test_r,
        "train_mean_R": train_r,
        "silent_mean_R": silent_r,
        "occupancy": "ALIVE",
        "motility": "OFF",
        "growth": "OFF",
        "death": "OFF",
        "J_max": 128000000.0,
        "warmup_s": 0.0,
        "bacterial_seed": 101,
        "frozen_before_traces": True,
        "MC_fading": {arm: curves[arm]["MC_fading"] for arm in ARMS},
        "k_half": {arm: curves[arm]["k_half"] for arm in ARMS},
        "r2_k0": {arm: curves[arm]["r2_k0"] for arm in ARMS},
        "test_r2": {arm: curves[arm]["test_r2"] for arm in ARMS},
        "lambda": {
            arm: [curves[arm]["by_k"][str(k)]["lambda"] for k in range(0, K_MAX + 1)]
            for arm in ARMS
        },
        "test_nrmse": {
            arm: [curves[arm]["by_k"][str(k)]["test_nrmse"] for k in range(0, K_MAX + 1)]
            for arm in ARMS
        },
        "field_wins_MC_sum": field_wins_sum,
        "delay_u_10_wins_MC_sum": delay_wins_sum,
        "delay_u_10_is_kill": False,
        "field_arm_is_gate": False,
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
        "lane_c_started": False,
        "ieee_axis_rewrite": False,
        "hunt_seeds_202_303": False,
    }
    (RESULTS / "lane_a_memory_capacity_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print("LANE_A_MEMORY_CAPACITY Jaeger delay capacity on frozen NARMA u.")
    print("Not CHARC. Not Dambre IPC. Not NARMA 0.895 restaged as MC. One washout.")
    print(f"u_sha256={digest} occupancy test_mean_R={test_r:.6g} train_mean_R={train_r:.6g} ALIVE")
    print(f"honesty motility=OFF growth=OFF death=OFF J_max=1.28e8 warmup=0 frozen_before_traces=true")
    print(f"silent_R={silent_r:.6g} delta_MC={DELTA_MC}")
    for arm in ARMS:
        c = curves[arm]
        print(
            f"  {arm} MC_fading={c['MC_fading']:.6g} k_1/2={c['k_half']} "
            f"R2(0)={c['r2_k0']:.6g} LANE_A_MEMORY_CAPACITY"
        )
        bits = " ".join(f"k{k}={c['test_r2'][k]:.4f}" for k in range(0, K_MAX + 1))
        print(f"    test_R2 {bits}")
    print(
        f"system={verdict} MC_RL={mc_rl:.6g} MC_SILENT={mc_silent:.6g} "
        f"MC_FIELD={mc_field:.6g} MC_DELAY10={mc_delay:.6g} "
        f"field_wins_sum={field_wins_sum} delay_wins_sum={delay_wins_sum} "
        f"LANE_A_MEMORY_CAPACITY={verdict}"
    )
    print(
        "Occupancy parent remains PASS. NARMA parent remains PASS. "
        "Carrier parent remains FAIL. Seed-replicate remains AXIS_HOLDS. "
        "DELAY_U_10 is a ceiling, not a kill. Gate 4 never kills."
    )
    if verdict == "PASS":
        print(
            "LANE_A_MEMORY_CAPACITY PASS system: RL fading-memory sum beats silent "
            "by frozen delta_MC. Not CHARC. Not IPC. Lead with the curve."
        )
    else:
        print(
            "LANE_A_MEMORY_CAPACITY FAIL system: RL did not beat silent by delta_MC. "
            "NARMA system PASS unchanged. Do not raise J_max."
        )


if __name__ == "__main__":
    main()
