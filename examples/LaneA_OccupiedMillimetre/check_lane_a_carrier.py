#!/usr/bin/env python3
"""Lane A carrier checker. Reconstruct u. Living-layer vs field. Not NARMA."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError:
    np = None

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_CARRIER.md"
PROTOCOL_JSON = HERE / "configs" / "carrier_protocol.json"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"

N_BINS = 200
LAMBDA_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
DELTA = 0.02
TRAIN = {2, 3, 4, 5}
TEST = {6, 7}
DISCARD = {0, 1}


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("narma", "charc", "lorenz", "mackey", "fig4b")):
        raise SystemExit("Lane A carrier refuses NARMA/CHARC/Lorenz/MG/Fig4b")


def load_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if "LANE_A_CARRIER" not in names or "u" not in names:
            raise SystemExit(f"{path.name} missing LANE_A_CARRIER/u")
        return list(reader), names


def col_block(row, prefix):
    return [float(row[f"{prefix}_{i}"]) for i in range(N_BINS)]


def mat_from(rows, prefix):
    return [col_block(r, prefix) for r in rows]


def u_of(rows):
    return [float(r["u"]) for r in rows]


def windows_of(rows):
    return [int(r["window"]) for r in rows]


def pop_std(y):
    n = len(y)
    if n == 0:
        return 0.0
    mu = sum(y) / n
    var = sum((v - mu) ** 2 for v in y) / n
    return math.sqrt(var)


def nrmse(ytrue, ypred):
    n = len(ytrue)
    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(ytrue, ypred)) / n)
    denom = pop_std(ytrue)
    return rmse / denom if denom > 0 else float("inf")


def standardize_fit(X):
    n = len(X)
    p = len(X[0]) if X else 0
    mu = [0.0] * p
    sg = [1.0] * p
    if n == 0 or p == 0:
        return mu, sg
    for j in range(p):
        mu[j] = sum(row[j] for row in X) / n
        var = sum((row[j] - mu[j]) ** 2 for row in X) / n
        s = math.sqrt(var)
        sg[j] = 1.0 if s < 1e-12 else s
    return mu, sg


def apply_std(X, mu, sg):
    return [[(row[j] - mu[j]) / sg[j] for j in range(len(mu))] for row in X]


def solve_linear(A, b):
    n = len(b)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[pivot] = M[pivot], M[k]
        if abs(M[k][k]) < 1e-18:
            M[k][k] = 1e-18
        diag = M[k][k]
        for j in range(k, n + 1):
            M[k][j] /= diag
        for i in range(n):
            if i == k:
                continue
            f = M[i][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    return [M[i][n] for i in range(n)]


def ridge_fit(X, y, lam):
    if np is not None:
        Xb = np.column_stack([np.ones(len(X)), np.asarray(X, dtype=float)])
        G = Xb.T @ Xb
        if G.shape[0] > 1:
            G[1:, 1:] = G[1:, 1:] + lam * np.eye(G.shape[0] - 1)
        rhs = Xb.T @ np.asarray(y, dtype=float)
        return np.linalg.solve(G, rhs).tolist()
    n = len(X)
    p = len(X[0])
    dim = p + 1
    G = [[0.0] * dim for _ in range(dim)]
    rhs = [0.0] * dim
    for i in range(n):
        row = [1.0] + X[i]
        yi = y[i]
        for a in range(dim):
            rhs[a] += row[a] * yi
            for b in range(dim):
                G[a][b] += row[a] * row[b]
    for j in range(1, dim):
        G[j][j] += lam
    return solve_linear(G, rhs)


def ridge_predict(X, w):
    out = []
    for row in X:
        s = w[0]
        for j, v in enumerate(row):
            s += w[j + 1] * v
        out.append(s)
    return out


def delay_features(u_all, k, idxs):
    X = []
    for i in idxs:
        feat = []
        for lag in range(1, k + 1):
            j = i - lag
            feat.append(0.0 if j < 0 else u_all[j])
        X.append(feat)
    return X


def select_lambda(X, y, windows, train_idx):
    folds = sorted(TRAIN)
    best_lam = LAMBDA_GRID[-1]
    best_score = float("inf")
    val_scores = []
    for lam in LAMBDA_GRID:
        fold_scores = []
        for held in folds:
            inner = [i for i in train_idx if windows[i] != held]
            val = [i for i in train_idx if windows[i] == held]
            if not inner or not val:
                continue
            Xtr = [X[i] for i in inner]
            ytr = [y[i] for i in inner]
            Xva = [X[i] for i in val]
            yva = [y[i] for i in val]
            mu, sg = standardize_fit(Xtr)
            w = ridge_fit(apply_std(Xtr, mu, sg), ytr, lam)
            pred = ridge_predict(apply_std(Xva, mu, sg), w)
            fold_scores.append(nrmse(yva, pred))
        score = sum(fold_scores) / len(fold_scores) if fold_scores else float("inf")
        val_scores.append((lam, score))
        if score < best_score - 1e-15 or (abs(score - best_score) <= 1e-15 and lam > best_lam):
            best_score = score
            best_lam = lam
    return best_lam, val_scores


def fit_score(X, y, windows):
    train_idx = [i for i, w in enumerate(windows) if w in TRAIN]
    test_idx = [i for i, w in enumerate(windows) if w in TEST]
    lam, val = select_lambda(X, y, windows, train_idx)
    Xtr = [X[i] for i in train_idx]
    ytr = [y[i] for i in train_idx]
    Xte = [X[i] for i in test_idx]
    yte = [y[i] for i in test_idx]
    mu, sg = standardize_fit(Xtr)
    w = ridge_fit(apply_std(Xtr, mu, sg), ytr, lam)
    train_pred = ridge_predict(apply_std(Xtr, mu, sg), w)
    test_pred = ridge_predict(apply_std(Xte, mu, sg), w)
    return {
        "lambda": lam,
        "train_nrmse": nrmse(ytr, train_pred),
        "test_nrmse": nrmse(yte, test_pred),
        "n_train": len(train_idx),
        "n_test": len(test_idx),
        "n_features": len(X[0]) if X else 0,
        "val_scores": val,
    }


def main() -> None:
    refuse_narma()
    occ = OCC_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js:
        raise SystemExit("carrier PROTOCOL not frozen")
    if '"narma": false' not in js or '"forecasting": false' not in js:
        raise SystemExit("carrier must freeze narma/forecasting false")
    if '"window_mean_u": "VOID"' not in js:
        raise SystemExit("window-mean u must stay VOID")

    driven_rows, _ = load_rows(RESULTS / "java_LANE_A_CARRIER_DRIVEN.csv")
    silent_rows, _ = load_rows(RESULTS / "java_LANE_A_CARRIER_SILENT.csv")
    if len(driven_rows) != 128 or len(silent_rows) != 128:
        raise SystemExit(f"expected 128 samples, got {len(driven_rows)} {len(silent_rows)}")

    y = u_of(driven_rows)
    windows = windows_of(driven_rows)
    if set(windows) != set(range(8)):
        raise SystemExit(f"windows {sorted(set(windows))}")
    # pulse samples in a window: t_in < 75 → u=0.5; four of sixteen
    train_u = [y[i] for i, w in enumerate(windows) if w in TRAIN]
    if abs(sum(1 for v in train_u if v > 0) / len(train_u) - 4 / 16) > 1e-9:
        raise SystemExit("train u duty is not 4/16; window-mean would be VOID")

    scores = {}
    scores["FIELD"] = fit_score(mat_from(driven_rows, "AHL"), y, windows)
    scores["R"] = fit_score(mat_from(driven_rows, "R"), y, windows)
    scores["L"] = fit_score(mat_from(driven_rows, "L"), y, windows)
    scores["SILENT_R"] = fit_score(mat_from(silent_rows, "R"), y, windows)

    idxs = list(range(len(y)))
    scores["DELAY_U_1"] = fit_score(delay_features(y, 1, idxs), y, windows)
    scores["DELAY_U_10"] = fit_score(delay_features(y, 10, idxs), y, windows)

    raw = fit_score([[v] for v in y], y, windows)
    raw["void"] = True
    scores["RAW_U"] = raw

    field = scores["FIELD"]["test_nrmse"]
    r = scores["R"]["test_nrmse"]
    ell = scores["L"]["test_nrmse"]
    silent = scores["SILENT_R"]["test_nrmse"]
    best_living = min(r, ell)
    winner = "R" if r <= ell else "L"

    field_wins_both = field <= r and field <= ell
    living_beats_field = best_living <= field - DELTA
    living_beats_silent = best_living <= silent - DELTA
    if field_wins_both or not (living_beats_field and living_beats_silent):
        verdict = "FAIL"
        living = "FAIL"
    else:
        verdict = "PASS"
        living = "PASS"

    summary = {
        "gate": "LaneA_CARRIER",
        "status_label": "LANE_A_CARRIER",
        "LANE_A_CARRIER": verdict,
        "living_layer": living,
        "occupancy_parent": "PASS",
        "reconstruction": True,
        "forecasting": False,
        "narma": False,
        "window_mean_u": "VOID",
        "margin_delta": DELTA,
        "n_train": 64,
        "n_test": 32,
        "small_n": True,
        "winner_living": winner,
        "test_nrmse": {k: v["test_nrmse"] for k, v in scores.items()},
        "train_nrmse": {k: v["train_nrmse"] for k, v in scores.items()},
        "lambda": {k: v["lambda"] for k, v in scores.items()},
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
        "RAW_U_void": True,
    }
    (RESULTS / "lane_a_carrier_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print("LANE_A_CARRIER reconstruction of u (sample-level). Not forecasting. Not NARMA.")
    for name in ("FIELD", "R", "L", "SILENT_R", "DELAY_U_1", "DELAY_U_10"):
        s = scores[name]
        print(
            f"  {name} test_NRMSE={s['test_nrmse']:.6g} train_NRMSE={s['train_nrmse']:.6g} "
            f"lambda={s['lambda']:g} LANE_A_CARRIER"
        )
    print(f"  RAW_U test_NRMSE={raw['test_nrmse']:.6g} VOID answer key LANE_A_CARRIER")
    print(
        f"living_layer={living} field={field:.6g} R={r:.6g} L={ell:.6g} "
        f"silent_R={silent:.6g} delta={DELTA} LANE_A_CARRIER={verdict}"
    )
    print(
        "Occupancy parent remains PASS. HybridDish Overall unchanged. "
        "Paper 1 unchanged. Lane B not started. Not Fig. 4b. Not C1c."
    )
    if verdict != "PASS":
        print(
            "LANE_A_CARRIER FAIL living-layer vs field (allowed). "
            "Do not raise J_max. Do not rewrite occupancy. NARMA extras later "
            "are system tests and do not reopen this living-layer FAIL."
        )
        raise SystemExit(0)
    print(
        "LANE_A_CARRIER PASS living-layer: R or L beat field and silent by frozen margin. "
        "Still not NARMA. Still not paper 2 improved 0.928."
    )


if __name__ == "__main__":
    main()
