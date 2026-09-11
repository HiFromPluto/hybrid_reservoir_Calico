#!/usr/bin/env python3
"""Lane A carrier maps audit. Not a living-layer reopen. Not NARMA."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from check_lane_a_carrier import (  # noqa: E402
    N_BINS,
    TRAIN,
    delay_features,
    fit_score,
    load_rows,
    mat_from,
    u_of,
    windows_of,
)

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_CARRIER_AUDIT.md"
PROTOCOL_JSON = HERE / "configs" / "carrier_audit_protocol.json"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"

OCC_WINDOWS = {4, 5, 6, 7}
OCCUPIED_EPS = 1e-12
SILENT_EPS = 1e-12
CENTRE_INDEX = 105
READOUT_X = 20
READOUT_Y = 10
SEED_X0, SEED_X1 = 300.0, 700.0
SEED_Y0, SEED_Y1 = 150.0, 350.0
VOXEL_X = 1000.0 / READOUT_X
VOXEL_Y = 500.0 / READOUT_Y

STANDING_3DEC = {
    "FIELD": 0.444,
    "R": 0.476,
    "L": 0.571,
    "SILENT_R": 1.000,
    "DELAY_U_10": 0.544,
}


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("narma", "charc", "lorenz", "mackey", "fig4b")):
        raise SystemExit("Lane A carrier audit refuses NARMA/CHARC/Lorenz/MG/Fig4b")


def round3(x: float) -> float:
    return float(f"{x:.3f}")


def pearson(a, b) -> float:
    n = len(a)
    if n == 0:
        return float("nan")
    ma = sum(a) / n
    mb = sum(b) / n
    da = [v - ma for v in a]
    db = [v - mb for v in b]
    num = sum(x * y for x, y in zip(da, db))
    va = sum(x * x for x in da)
    vb = sum(y * y for y in db)
    if va < 1e-30 or vb < 1e-30:
        return float("nan")
    return num / math.sqrt(va * vb)


def cosine(a, b) -> float:
    num = sum(x * y for x, y in zip(a, b))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na < 1e-30 or nb < 1e-30:
        return float("nan")
    return num / (na * nb)


def col(mat, j):
    return [row[j] for row in mat]


def restrict(mat, idxs):
    return [[row[j] for j in idxs] for row in mat]


def mean_cols(mat, idxs):
    if not idxs:
        return [0.0] * len(mat)
    k = len(idxs)
    return [[sum(row[j] for j in idxs) / k] for row in mat]


def seed_bins():
    out = []
    for x in range(READOUT_X):
        x0 = x * VOXEL_X
        x1 = (x + 1) * VOXEL_X
        if x1 <= SEED_X0 or x0 >= SEED_X1:
            continue
        for y in range(READOUT_Y):
            y0 = y * VOXEL_Y
            y1 = (y + 1) * VOXEL_Y
            if y1 <= SEED_Y0 or y0 >= SEED_Y1:
                continue
            out.append(x * READOUT_Y + y)
    return out


def scalar_multiple_of_u(series, u, tol=1e-9) -> bool:
    on = [(a, b) for a, b in zip(series, u) if abs(b) > 1e-15]
    if not on:
        return all(abs(a) <= tol for a in series)
    c = sum(a * b for a, b in on) / sum(b * b for _, b in on)
    return all(abs(a - c * b) <= tol for a, b in zip(series, u))


def matrices_identical(a, b, tol=0.0) -> bool:
    if len(a) != len(b):
        return False
    for ra, rb in zip(a, b):
        if len(ra) != len(rb):
            return False
        for xa, xb in zip(ra, rb):
            if abs(xa - xb) > tol:
                return False
    return True


def is_constant(mat, tol=1e-15) -> bool:
    vals = [v for row in mat for v in row]
    if not vals:
        return True
    lo = min(vals)
    hi = max(vals)
    return hi - lo <= tol


def max_abs(mat) -> float:
    m = 0.0
    for row in mat:
        for v in row:
            a = abs(v)
            if a > m:
                m = a
    return m


def main() -> None:
    refuse_narma()
    occ = OCC_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: FAIL**" not in car or "LANE_A_CARRIER" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js:
        raise SystemExit("audit PROTOCOL not frozen")
    if '"narma": false' not in js or '"parent_carrier_rewrite": false' not in js:
        raise SystemExit("audit must freeze narma false and not rewrite carrier")
    if '"window_mean_u": "VOID"' not in js:
        raise SystemExit("window-mean u must stay VOID")
    if '"occupied_mask_is_living_layer_gate": false' not in js:
        raise SystemExit("occupied-mask must not be a living-layer gate")

    driven_rows, _ = load_rows(RESULTS / "java_LANE_A_CARRIER_DRIVEN.csv")
    silent_rows, _ = load_rows(RESULTS / "java_LANE_A_CARRIER_SILENT.csv")
    if len(driven_rows) != 128 or len(silent_rows) != 128:
        raise SystemExit(f"expected 128 samples, got {len(driven_rows)} {len(silent_rows)}")

    y = u_of(driven_rows)
    windows = windows_of(driven_rows)
    ahl = mat_from(driven_rows, "AHL")
    rmap = mat_from(driven_rows, "R")
    lmap = mat_from(driven_rows, "L")
    silent_ahl = mat_from(silent_rows, "AHL")
    silent_r = mat_from(silent_rows, "R")
    train_idx = [i for i, w in enumerate(windows) if w in TRAIN]

    occupied = [
        j
        for j in range(N_BINS)
        if any(abs(rmap[i][j]) > OCCUPIED_EPS for i in train_idx)
    ]
    empty_train = [j for j in range(N_BINS) if j not in occupied]
    empty_frac = len(empty_train) / N_BINS
    seed = seed_bins()

    scores = {}
    scores["FIELD"] = fit_score(ahl, y, windows)
    scores["R"] = fit_score(rmap, y, windows)
    scores["L"] = fit_score(lmap, y, windows)
    scores["SILENT_R"] = fit_score(silent_r, y, windows)
    scores["DELAY_U_10"] = fit_score(delay_features(y, 10, list(range(len(y)))), y, windows)
    raw = fit_score([[v] for v in y], y, windows)
    raw["void"] = True
    scores["RAW_U"] = raw

    repro = {}
    repro_ok = True
    for name, target in STANDING_3DEC.items():
        got = round3(scores[name]["test_nrmse"])
        match = got == target
        repro[name] = {"got": scores[name]["test_nrmse"], "rounded": got, "standing": target, "match": match}
        if not match:
            repro_ok = False
    raw_ok = scores["RAW_U"]["test_nrmse"] < 1e-6
    repro["RAW_U"] = {
        "got": scores["RAW_U"]["test_nrmse"],
        "standing": "~0",
        "match": raw_ok,
    }
    if not raw_ok:
        repro_ok = False

    if occupied:
        scores["FIELD_OCC"] = fit_score(restrict(ahl, occupied), y, windows)
        scores["R_OCC"] = fit_score(restrict(rmap, occupied), y, windows)
        scores["MEAN_AHL_OCC"] = fit_score(mean_cols(ahl, occupied), y, windows)
        scores["MEAN_R_OCC"] = fit_score(mean_cols(rmap, occupied), y, windows)
        scores["MEAN_L_OCC"] = fit_score(mean_cols(lmap, occupied), y, windows)
    else:
        scores["FIELD_OCC"] = {"test_nrmse": float("nan"), "train_nrmse": float("nan"), "lambda": None}
        scores["R_OCC"] = {"test_nrmse": float("nan"), "train_nrmse": float("nan"), "lambda": None}
        scores["MEAN_AHL_OCC"] = scores["FIELD_OCC"]
        scores["MEAN_R_OCC"] = scores["R_OCC"]
        scores["MEAN_L_OCC"] = scores["R_OCC"]
    scores["MEAN_AHL_ALL"] = fit_score(mean_cols(ahl, list(range(N_BINS))), y, windows)
    scores["CENTRE_AHL"] = fit_score([[row[CENTRE_INDEX]] for row in ahl], y, windows)

    occ_idx = [i for i, w in enumerate(windows) if w in OCC_WINDOWS]
    if occupied:
        occ_mean_r = sum(sum(rmap[i][j] for j in occupied) / len(occupied) for i in occ_idx) / len(occ_idx)
    else:
        occ_mean_r = 0.0
    occupancy_alive = occ_mean_r >= 0.05

    pearson_bins = [pearson(col(ahl, j), y) for j in range(N_BINS)]
    finite_p = [p for p in pearson_bins if p == p]
    mean_ahl_all = [sum(row) / N_BINS for row in ahl]
    pearson_mean = pearson(mean_ahl_all, y)
    cosine_mean = cosine(mean_ahl_all, y)
    max_abs_pearson = max((abs(p) for p in finite_p), default=float("nan"))

    field_vs_r_identical = matrices_identical(ahl, rmap)
    field_constant = is_constant(ahl)
    every_bin_scaled_u = all(scalar_multiple_of_u(col(ahl, j), y) for j in range(N_BINS))
    field_nrmse = scores["FIELD"]["test_nrmse"]
    raw_nrmse = scores["RAW_U"]["test_nrmse"]
    field_is_u_clone = field_nrmse < 0.05
    silent_ahl_max = max_abs(silent_ahl)
    silent_r_max = max_abs(silent_r)
    silent_ok = silent_ahl_max <= SILENT_EPS and (
        silent_r_max <= SILENT_EPS or abs(scores["SILENT_R"]["test_nrmse"] - 1.0) < 1e-12
    )

    leak = {
        "field_vs_r_identical": field_vs_r_identical,
        "field_constant": field_constant,
        "every_ahl_bin_scaled_u": every_bin_scaled_u,
        "field_nrmse_near_zero": field_is_u_clone,
        "raw_u_not_void": not raw_ok,
        "silent_not_empty": not silent_ok,
        "max_abs_pearson_ahl_vs_u": max_abs_pearson,
        "pearson_spatial_mean_ahl_vs_u": pearson_mean,
        "cosine_spatial_mean_ahl_vs_u": cosine_mean,
        "silent_ahl_max_abs": silent_ahl_max,
        "silent_r_max_abs": silent_r_max,
    }
    leak_kill = (
        field_vs_r_identical
        or field_constant
        or every_bin_scaled_u
        or field_is_u_clone
        or not raw_ok
        or not silent_ok
    )

    field_full = scores["FIELD"]["test_nrmse"]
    r_full = scores["R"]["test_nrmse"]
    field_occ_n = scores["FIELD_OCC"]["test_nrmse"]
    r_occ_n = scores["R_OCC"]["test_nrmse"]
    occ_field_le_r = field_occ_n <= r_occ_n if occupied else False
    ranking_flipped = occupied and (r_occ_n < field_occ_n) and (field_full <= r_full)

    if not occupancy_alive:
        verdict = "NOT_SCORED"
    elif not repro_ok or leak_kill:
        verdict = "KILL"
    elif ranking_flipped:
        verdict = "SCOPE_NOTE"
    elif occ_field_le_r:
        verdict = "CLEAN"
    else:
        verdict = "KILL"

    summary = {
        "gate": "LaneA_CARRIER_AUDIT",
        "status_label": "LANE_A_CARRIER_AUDIT",
        "LANE_A_CARRIER_AUDIT": verdict,
        "occupancy_parent": "PASS",
        "carrier_parent": "FAIL",
        "carrier_restaged_pass": False,
        "living_layer_reopened": False,
        "occupied_mask_is_living_layer_gate": False,
        "reconstruction": True,
        "forecasting": False,
        "narma": False,
        "window_mean_u": "VOID",
        "reproduction_ok": repro_ok,
        "reproduction": repro,
        "occupancy_replay": {
            "occupied_bin_mean_R_windows_4_7": occ_mean_r,
            "alive": occupancy_alive,
            "carrier_standing_cell_mean_R": 0.661,
        },
        "n_occupied_bins": len(occupied),
        "empty_r_bin_fraction_train": empty_frac,
        "n_empty_r_bins_train": len(empty_train),
        "n_seed_rectangle_bins": len(seed),
        "occupied_indices": occupied,
        "leak": leak,
        "leak_kill": leak_kill,
        "occupied_mask": {
            "FIELD_test_nrmse": field_occ_n,
            "R_test_nrmse": r_occ_n,
            "field_le_R": occ_field_le_r,
            "ranking_flipped_vs_full": ranking_flipped,
        },
        "test_nrmse": {k: v["test_nrmse"] for k, v in scores.items()},
        "train_nrmse": {k: v.get("train_nrmse") for k, v in scores.items()},
        "lambda": {k: v.get("lambda") for k, v in scores.items()},
        "scheduler_lag_s": 0.05,
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
        "next_if_clean_or_scope_note": "LaneA_NARMA10",
    }
    (RESULTS / "lane_a_carrier_audit_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print("LANE_A_CARRIER_AUDIT maps sanity. Not living-layer reopen. Not NARMA.")
    print(
        f"reproduction={'OK' if repro_ok else 'MISMATCH'} "
        f"FIELD={scores['FIELD']['test_nrmse']:.6g} R={scores['R']['test_nrmse']:.6g} "
        f"L={scores['L']['test_nrmse']:.6g} SILENT_R={scores['SILENT_R']['test_nrmse']:.6g} "
        f"DELAY_U_10={scores['DELAY_U_10']['test_nrmse']:.6g} "
        f"RAW_U={scores['RAW_U']['test_nrmse']:.6g}"
    )
    print(
        f"occupied_bins={len(occupied)}/200 empty_train_frac={empty_frac:.3f} "
        f"seed_rectangle_bins={len(seed)}"
    )
    print(
        f"FIELD_OCC test_NRMSE={field_occ_n:.6g} R_OCC test_NRMSE={r_occ_n:.6g} "
        f"field_le_R={occ_field_le_r} ranking_flipped={ranking_flipped}"
    )
    for name in (
        "MEAN_AHL_OCC",
        "MEAN_R_OCC",
        "MEAN_L_OCC",
        "MEAN_AHL_ALL",
        "CENTRE_AHL",
    ):
        s = scores[name]
        print(
            f"  {name} test_NRMSE={s['test_nrmse']:.6g} "
            f"lambda={s.get('lambda')} LANE_A_CARRIER_AUDIT"
        )
    print(
        f"leak_kill={leak_kill} max_|pearson|_AHL_vs_u={max_abs_pearson:.6g} "
        f"pearson_mean_AHL={pearson_mean:.6g} field_identical_R={field_vs_r_identical} "
        f"field_constant={field_constant} scaled_u={every_bin_scaled_u} "
        f"silent_AHL_max={silent_ahl_max:.3g} silent_R_max={silent_r_max:.3g}"
    )
    print(
        f"occupancy_replay occupied_bin_mean_R={occ_mean_r:.6g} "
        f"{'ALIVE' if occupancy_alive else 'DEAD'}"
    )
    print(
        f"LANE_A_CARRIER_AUDIT={verdict}. Occupancy parent remains PASS. "
        f"Carrier parent remains FAIL. HybridDish Overall unchanged. "
        f"Paper 1 unchanged. Lane B not started. Not Fig. 4b. Not C1c."
    )
    if verdict in ("CLEAN", "SCOPE_NOTE"):
        print(
            "Next named extra may be LaneA_NARMA10 as a system test with a field arm. "
            "It does not reopen LANE_A_CARRIER FAIL."
        )
    elif verdict == "KILL":
        print("KILL: do not start NARMA. Do not restage carrier as PASS. Stop.")
    else:
        print("NOT_SCORED: occupancy replay DEAD. Do not raise J_max.")


if __name__ == "__main__":
    main()
