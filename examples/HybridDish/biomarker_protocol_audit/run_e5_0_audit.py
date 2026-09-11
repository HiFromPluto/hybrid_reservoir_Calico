#!/usr/bin/env python3
"""E5.0 biomarker protocol audit. No Java, no BSim.

Split is written to disk before any metric is computed.
Confirmation is scored once, last, and never used for lambda.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[3]
PKG = Path(__file__).resolve().parent
RESULTS = PKG / "results"
MAT_VENDOR = (
    REPO / "examples" / "HybridDish" / "external_analysis"
    / "synthetic_biomarker_BSim_inputs.mat"
)
MAT_SOURCE = Path(
    r"C:\Users\Ceylin\Downloads\PRC_Module8_Integrated_QS_Audit01"
    r"\PRC_Module8_Integrated_QS_Audit01\reservoir_new"
    r"\data\synthetic_biomarker_BSim_inputs.mat"
)
EXT_ROOT = Path(
    r"C:\Users\Ceylin\Downloads\PRC_Module8_Integrated_QS_Audit01"
    r"\PRC_Module8_Integrated_QS_Audit01\reservoir_new"
)
MAT_GENERATOR = EXT_ROOT / "matlab" / "MAT_biomarker_progression.m"

CHANNEL_NAMES = ["CRP", "IL6", "Ferritin", "Lactate", "Glucose"]
LAMBDA_GRID = np.array([1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6])
SPLIT_SEED = 20260817
N_DEV, N_VAL, N_CONF = 400, 200, 400
PHYS_KEYWORDS = ("physionet", "mimic", "eicu")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def index_string(idx) -> str:
    return ",".join(str(int(i)) for i in idx)


def hamilton_allocate(counts, seats):
    total = int(sum(counts))
    quotas = [c * seats / total for c in counts]
    floors = [int(q) for q in quotas]
    remainders = [q - f for q, f in zip(quotas, floors)]
    leftover = seats - sum(floors)
    order = sorted(range(len(counts)), key=lambda i: (-remainders[i], i))
    alloc = floors[:]
    for k in range(leftover):
        alloc[order[k]] += 1
    return alloc


def standardize_fit(X):
    mu = X.mean(axis=0)
    sg = X.std(axis=0, ddof=0)
    sg = np.where(sg < 1e-12, 1.0, sg)
    return (X - mu) / sg, mu, sg


def standardize_apply(X, mu, sg):
    return (X - mu) / sg


def ridge_fit(X, y, lam):
    n = X.shape[1]
    Xb = np.column_stack([np.ones(X.shape[0]), X])
    G = Xb.T @ Xb
    G[1:, 1:] = G[1:, 1:] + lam * np.eye(n)
    w = np.linalg.solve(G, Xb.T @ y)
    return w


def ridge_predict(X, w):
    Xb = np.column_stack([np.ones(X.shape[0]), X])
    return Xb @ w


def nrmse_pop(ytrue, ypred):
    ytrue = np.asarray(ytrue, dtype=float).ravel()
    ypred = np.asarray(ypred, dtype=float).ravel()
    rmse = np.sqrt(np.mean((ytrue - ypred) ** 2))
    denom = ytrue.std(ddof=0)
    if denom > 0:
        return float(rmse / denom)
    return float("inf")


def binary_auc(y_bin, scores):
    y_bin = np.asarray(y_bin).ravel().astype(int)
    scores = np.asarray(scores, dtype=float).ravel()
    n_pos = int(y_bin.sum())
    n_neg = int(len(y_bin) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    return float(roc_auc_score(y_bin, scores))


def macro_auc(y, scores):
    y = np.asarray(y).ravel().astype(int)
    scores = np.asarray(scores, dtype=float)
    if scores.ndim == 1:
        scores = np.column_stack([-scores, scores])
    aucs = []
    for c in range(scores.shape[1]):
        aucs.append(binary_auc(y == c, scores[:, c]))
    aucs = [a for a in aucs if np.isfinite(a)]
    if not aucs:
        return float("nan")
    return float(np.mean(aucs))


def accuracy(y, pred):
    y = np.asarray(y).ravel()
    pred = np.asarray(pred).ravel()
    return float(np.mean(y == pred))


def patient_slopes(Xnorm):
    n, t, b = Xnorm.shape
    tt = np.arange(t, dtype=float)
    t0 = tt - tt.mean()
    denom = np.dot(t0, t0)
    slopes = np.empty((n, b), dtype=float)
    for i in range(n):
        for c in range(b):
            y = Xnorm[i, :, c]
            slopes[i, c] = np.dot(t0, y - y.mean()) / denom
    return slopes


def fit_ridge_ovr(Xtr, ytr, Xva, yva):
    classes = np.array([0, 1, 2])
    Xtr_z, mu, sg = standardize_fit(Xtr)
    Xva_z = standardize_apply(Xva, mu, sg)
    best_auc = -np.inf
    best_lam = LAMBDA_GRID[-1]
    best_ws = None
    for lam in LAMBDA_GRID:
        ws = []
        scores = np.zeros((Xva_z.shape[0], 3))
        for c in classes:
            ybin = (ytr == c).astype(float)
            w = ridge_fit(Xtr_z, ybin, lam)
            ws.append(w)
            scores[:, c] = ridge_predict(Xva_z, w)
        auc = macro_auc(yva, scores)
        if auc > best_auc + 1e-15 or (
            abs(auc - best_auc) <= 1e-15 and lam > best_lam
        ):
            best_auc = auc
            best_lam = lam
            best_ws = ws
    return {
        "lambda": float(best_lam),
        "val_auc": float(best_auc),
        "weights": best_ws,
        "mu": mu,
        "sg": sg,
    }


def predict_ridge_ovr(model, X):
    Xz = standardize_apply(X, model["mu"], model["sg"])
    scores = np.column_stack(
        [ridge_predict(Xz, w) for w in model["weights"]]
    )
    pred = np.argmax(scores, axis=1)
    return pred, scores


def fit_multinomial(Xtr, ytr, Xva, yva):
    Xtr_z, mu, sg = standardize_fit(Xtr)
    Xva_z = standardize_apply(Xva, mu, sg)
    best_auc = -np.inf
    best_lam = LAMBDA_GRID[-1]
    best_clf = None
    for lam in LAMBDA_GRID:
        clf = LogisticRegression(
            C=1.0 / lam,
            solver="lbfgs",
            max_iter=4000,
            random_state=0,
        )
        clf.fit(Xtr_z, ytr)
        scores = clf.predict_proba(Xva_z)
        auc = macro_auc(yva, scores)
        if auc > best_auc + 1e-15 or (
            abs(auc - best_auc) <= 1e-15 and lam > best_lam
        ):
            best_auc = auc
            best_lam = lam
            best_clf = clf
    return {
        "lambda": float(best_lam),
        "val_auc": float(best_auc),
        "clf": best_clf,
        "mu": mu,
        "sg": sg,
    }


def predict_multinomial(model, X):
    Xz = standardize_apply(X, model["mu"], model["sg"])
    scores = model["clf"].predict_proba(Xz)
    pred = model["clf"].predict(Xz)
    return pred, scores


def class_metrics(y, pred, scores):
    return {
        "accuracy": accuracy(y, pred),
        "macro_auc": macro_auc(y, scores),
        "n": int(len(y)),
        "y_counts": {str(c): int(np.sum(y == c)) for c in range(3)},
    }


def delay10_rows(Xnorm, patients, t_lo, t_hi):
    """50-D delay-10 of 5 channels; target = channel 0 at t+1."""
    feats = []
    targets = []
    pids = []
    t_idx = []
    for p in patients:
        for t in range(t_lo, t_hi + 1):
            row = []
            for b in range(5):
                for k in range(10):
                    row.append(Xnorm[p, t - k, b])
            feats.append(row)
            targets.append(Xnorm[p, t + 1, 0])
            pids.append(p)
            t_idx.append(t)
    return (
        np.asarray(feats, dtype=float),
        np.asarray(targets, dtype=float),
        np.asarray(pids, dtype=int),
        np.asarray(t_idx, dtype=int),
    )


def persist_rows(Xnorm, patients, t_lo, t_hi):
    feats = []
    targets = []
    pids = []
    for p in patients:
        for t in range(t_lo, t_hi + 1):
            feats.append(Xnorm[p, t, 0])
            targets.append(Xnorm[p, t + 1, 0])
            pids.append(p)
    return (
        np.asarray(feats, dtype=float),
        np.asarray(targets, dtype=float),
        np.asarray(pids, dtype=int),
    )


def patient_pooled_nrmse(ytrue, ypred, pids):
    per = []
    skipped = 0
    for p in np.unique(pids):
        m = pids == p
        val = nrmse_pop(ytrue[m], ypred[m])
        if np.isfinite(val):
            per.append(val)
        else:
            skipped += 1
    pooled = float(np.mean(per)) if per else float("inf")
    return pooled, per, skipped


def concat_nrmse(ytrue, ypred):
    return nrmse_pop(ytrue, ypred)


def fit_ridge_forecast(Xtr, ytr, ptr, Xva, yva, pva):
    Xtr_z, mu, sg = standardize_fit(Xtr)
    Xva_z = standardize_apply(Xva, mu, sg)
    best = np.inf
    best_lam = LAMBDA_GRID[-1]
    best_w = None
    for lam in LAMBDA_GRID:
        w = ridge_fit(Xtr_z, ytr, lam)
        pred = ridge_predict(Xva_z, w)
        pooled, _, _ = patient_pooled_nrmse(yva, pred, pva)
        if pooled < best - 1e-15 or (
            abs(pooled - best) <= 1e-15 and lam > best_lam
        ):
            best = pooled
            best_lam = lam
            best_w = w
    return {
        "lambda": float(best_lam),
        "val_nrmse": float(best),
        "weights": best_w,
        "mu": mu,
        "sg": sg,
    }


def search_physionet_files(roots):
    hits = []
    for root in roots:
        if not Path(root).exists():
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [
                d
                for d in dirnames
                if d not in {".git", "__pycache__", "dist", "node_modules"}
            ]
            for name in filenames:
                low = name.lower()
                if any(k in low for k in PHYS_KEYWORDS):
                    hits.append(str(Path(dirpath) / name))
    return hits


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, (float, np.floating)):
        if not np.isfinite(x):
            return "inf"
        return f"{x:.{nd}f}"
    return str(x)


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    assert MAT_VENDOR.exists(), f"STOP: vendored .mat absent: {MAT_VENDOR}"

    mat_sha = sha256_file(MAT_VENDOR)
    src_sha = sha256_file(MAT_SOURCE) if MAT_SOURCE.exists() else "ABSENT"
    gen_sha = sha256_file(MAT_GENERATOR) if MAT_GENERATOR.exists() else "ABSENT"
    assert src_sha == mat_sha, "STOP: vendor SHA-256 does not match source"

    mat = loadmat(str(MAT_VENDOR), squeeze_me=False)
    X = np.asarray(mat["X"], dtype=float)
    Xnorm = np.asarray(mat["Xnorm"], dtype=float)
    Y = np.asarray(mat["Y"]).ravel().astype(int)
    Q = np.asarray(mat["Q"], dtype=float)
    n, t_windows, b = X.shape
    assert (n, t_windows, b) == (1000, 48, 5)
    assert Xnorm.shape == X.shape
    assert Y.shape == (1000,)
    missing_x = int(np.isnan(X).sum() + np.isnan(Xnorm).sum() + np.isnan(Q).sum())

    y_counts = {c: int(np.sum(Y == c)) for c in range(3)}
    blocked = (
        np.all(Y[:333] == 0)
        and np.all(Y[333:666] == 1)
        and np.all(Y[666:] == 2)
    )
    naive_counts = {
        "development_0_399": {c: int(np.sum(Y[:400] == c)) for c in range(3)},
        "validation_400_599": {c: int(np.sum(Y[400:600] == c)) for c in range(3)},
        "confirmation_600_999": {c: int(np.sum(Y[600:] == c)) for c in range(3)},
    }

    # ----- Module 1: split BEFORE any metric -----
    counts = [y_counts[0], y_counts[1], y_counts[2]]
    alloc_dev = hamilton_allocate(counts, N_DEV)
    alloc_val = hamilton_allocate(counts, N_VAL)
    rng = np.random.RandomState(SPLIT_SEED)
    assignment = {c: {"dev": [], "val": [], "conf": []} for c in range(3)}
    for c in range(3):
        idx = np.where(Y == c)[0].astype(int).copy()
        rng.shuffle(idx)
        n_d = alloc_dev[c]
        n_v = alloc_val[c]
        assignment[c]["dev"] = [int(i) for i in idx[:n_d]]
        assignment[c]["val"] = [int(i) for i in idx[n_d : n_d + n_v]]
        assignment[c]["conf"] = [int(i) for i in idx[n_d + n_v :]]

    development = np.sort(
        np.concatenate([assignment[c]["dev"] for c in range(3)])
    ).astype(int)
    validation = np.sort(
        np.concatenate([assignment[c]["val"] for c in range(3)])
    ).astype(int)
    confirmation = np.sort(
        np.concatenate([assignment[c]["conf"] for c in range(3)])
    ).astype(int)
    assert len(development) == N_DEV
    assert len(validation) == N_VAL
    assert len(confirmation) == N_CONF
    union = np.concatenate([development, validation, confirmation])
    assert len(np.unique(union)) == 1000
    assert np.all(np.sort(union) == np.arange(1000))

    def ycount(idx):
        return {str(c): int(np.sum(Y[idx] == c)) for c in range(3)}

    scout_dev_4_per_class = []
    for c in range(3):
        scout_dev_4_per_class.extend(assignment[c]["dev"][:4])
    scout_conf_3_3_2 = (
        assignment[0]["conf"][:3]
        + assignment[1]["conf"][:3]
        + assignment[2]["conf"][:2]
    )
    scout_dev_first12 = development[:12].tolist()
    scout_conf_first8 = confirmation[:8].tolist()

    split_doc = {
        "audit": "E5.0",
        "frozen_before_metrics": True,
        "split_unit": "patient_index",
        "index_base": 0,
        "rng_seed": SPLIT_SEED,
        "rng": "numpy.random.RandomState",
        "rule": "stratified_hamilton_40_20_40",
        "mat_class_blocked": bool(blocked),
        "naive_0_399_not_stratified": bool(blocked),
        "naive_y_counts": naive_counts,
        "hamilton_seats": {
            "development": alloc_dev,
            "validation": alloc_val,
            "confirmation": [
                counts[c] - alloc_dev[c] - alloc_val[c] for c in range(3)
            ],
        },
        "development": development.tolist(),
        "validation": validation.tolist(),
        "confirmation": confirmation.tolist(),
        "assignment_order_by_class": assignment,
        "y_counts": {
            "full": {str(k): v for k, v in y_counts.items()},
            "development": ycount(development),
            "validation": ycount(validation),
            "confirmation": ycount(confirmation),
        },
        "index_string_sha256": {
            "development": sha256_text(index_string(development)),
            "validation": sha256_text(index_string(validation)),
            "confirmation": sha256_text(index_string(confirmation)),
        },
        "scout_draws_predeclared": {
            "development_4_per_Y_class": [int(i) for i in scout_dev_4_per_class],
            "development_first_12_sorted": [int(i) for i in scout_dev_first12],
            "confirmation_3_3_2": [int(i) for i in scout_conf_3_3_2],
            "confirmation_first_8_sorted": [int(i) for i in scout_conf_first8],
            "note": (
                "If Task A survives: use 4-per-class development and "
                "3-3-2 confirmation. Else: first 12 sorted development "
                "and first 8 sorted confirmation. Not run in E5.0."
            ),
        },
        "mat_sha256": mat_sha,
    }
    split_path = RESULTS / "patient_split.json"
    split_path.write_text(
        json.dumps(split_doc, indent=2) + "\n", encoding="utf-8"
    )

    # ----- Module 0 inventory (no predictive scores) -----
    phys_hits = search_physionet_files([REPO, EXT_ROOT])
    channel_stats = []
    for c, name in enumerate(CHANNEL_NAMES):
        a = X[:, :, c]
        channel_stats.append(
            {
                "channel_index": c,
                "name": name,
                "min": float(a.min()),
                "max": float(a.max()),
                "mean": float(a.mean()),
                "std": float(a.std(ddof=0)),
                "missing": int(np.isnan(a).sum()),
            }
        )

    # ----- Module 2: features from Xnorm, fits on development only -----
    mean_x = Xnorm.mean(axis=1)
    last_x = Xnorm[:, -1, :]
    slope_x = patient_slopes(Xnorm)
    mean_slope = slope_x.mean(axis=1, keepdims=True)

    feature_sets = {
        "MEAN_XNORM": mean_x,
        "LAST_XNORM": last_x,
        "SLOPE_XNORM": slope_x,
        "MEAN_SLOPE": mean_slope,
    }

    y_dev = Y[development]
    y_val = Y[validation]
    y_conf = Y[confirmation]
    majority = int(np.bincount(y_dev, minlength=3).argmax())

    task_a_rows = []
    task_a_detail = {}

    def eval_constant(split_name, idx, y_split):
        pred = np.full(len(idx), majority, dtype=int)
        scores = np.zeros((len(idx), 3))
        scores[:, majority] = 1.0
        m = class_metrics(y_split, pred, scores)
        task_a_rows.append(
            {
                "task": "A",
                "feature": "MAJORITY_CLASS",
                "classifier": "constant",
                "split": split_name,
                "lambda": "",
                "accuracy": m["accuracy"],
                "macro_auc": m["macro_auc"],
                "nrmse_patient_pooled": "",
                "n": m["n"],
            }
        )
        return m

    maj = {
        "development": eval_constant("development", development, y_dev),
        "validation": eval_constant("validation", validation, y_val),
        "confirmation": eval_constant("confirmation", confirmation, y_conf),
    }
    task_a_detail["MAJORITY_CLASS"] = {
        "majority_class": majority,
        "metrics": maj,
        "lambda": None,
        "classifier": "constant",
    }

    for feat_name, feat in feature_sets.items():
        Xtr = feat[development]
        Xva = feat[validation]
        Xcf = feat[confirmation]
        for clf_name, fit_fn, pred_fn in (
            ("ridge_ovr", fit_ridge_ovr, predict_ridge_ovr),
            ("multinomial_logistic", fit_multinomial, predict_multinomial),
        ):
            model = fit_fn(Xtr, y_dev, Xva, y_val)
            split_metrics = {}
            for split_name, Xs, ys in (
                ("development", Xtr, y_dev),
                ("validation", Xva, y_val),
                ("confirmation", Xcf, y_conf),
            ):
                pred, scores = pred_fn(model, Xs)
                m = class_metrics(ys, pred, scores)
                split_metrics[split_name] = m
                task_a_rows.append(
                    {
                        "task": "A",
                        "feature": feat_name,
                        "classifier": clf_name,
                        "split": split_name,
                        "lambda": model["lambda"],
                        "accuracy": m["accuracy"],
                        "macro_auc": m["macro_auc"],
                        "nrmse_patient_pooled": "",
                        "n": m["n"],
                    }
                )
            task_a_detail[f"{feat_name}:{clf_name}"] = {
                "lambda": model["lambda"],
                "metrics": split_metrics,
            }

    def a_void_hit(feat_name):
        hits = []
        for clf_name in ("ridge_ovr", "multinomial_logistic"):
            m = task_a_detail[f"{feat_name}:{clf_name}"]["metrics"][
                "confirmation"
            ]
            if m["accuracy"] >= 0.80 or m["macro_auc"] >= 0.90:
                hits.append(
                    {
                        "feature": feat_name,
                        "classifier": clf_name,
                        "accuracy": m["accuracy"],
                        "macro_auc": m["macro_auc"],
                    }
                )
        return hits

    void_a_hits = a_void_hit("SLOPE_XNORM") + a_void_hit("MEAN_XNORM")
    task_a_void = bool(void_a_hits)

    # Task B
    T = t_windows
    persist_all = {}
    task_b_rows = []
    for split_name, idx in (
        ("development", development),
        ("validation", validation),
        ("confirmation", confirmation),
    ):
        xp, yp, pp = persist_rows(Xnorm, idx, 0, T - 2)
        pooled, per, skipped = patient_pooled_nrmse(yp, xp, pp)
        persist_all[split_name] = {
            "nrmse_patient_pooled": pooled,
            "nrmse_concat_diagnostic": concat_nrmse(yp, xp),
            "n_patients": int(len(idx)),
            "n_pairs": int(len(yp)),
            "skipped_zero_std": skipped,
            "mean_per_patient": float(np.mean(per)),
        }
        task_b_rows.append(
            {
                "task": "B",
                "feature": "PERSISTENCE",
                "classifier": "x_t",
                "split": split_name,
                "lambda": "",
                "accuracy": "",
                "macro_auc": "",
                "nrmse_patient_pooled": pooled,
                "n": int(len(yp)),
            }
        )

    # train-mean intercept on persistence windows
    _, ytr_p, _ = persist_rows(Xnorm, development, 0, T - 2)
    train_mean_p = float(ytr_p.mean())
    intercept_all = {}
    for split_name, idx in (
        ("development", development),
        ("validation", validation),
        ("confirmation", confirmation),
    ):
        _, yp, pp = persist_rows(Xnorm, idx, 0, T - 2)
        pred = np.full_like(yp, train_mean_p)
        pooled, _, skipped = patient_pooled_nrmse(yp, pred, pp)
        intercept_all[split_name] = {
            "nrmse_patient_pooled": pooled,
            "nrmse_concat_diagnostic": concat_nrmse(yp, pred),
            "train_mean": train_mean_p,
            "skipped_zero_std": skipped,
        }
        task_b_rows.append(
            {
                "task": "B",
                "feature": "TRAIN_MEAN_INTERCEPT",
                "classifier": "constant",
                "split": split_name,
                "lambda": "",
                "accuracy": "",
                "macro_auc": "",
                "nrmse_patient_pooled": pooled,
                "n": int(len(yp)),
            }
        )

    Xtr_d, ytr_d, ptr_d, _ = delay10_rows(Xnorm, development, 9, T - 2)
    Xva_d, yva_d, pva_d, _ = delay10_rows(Xnorm, validation, 9, T - 2)
    Xcf_d, ycf_d, pcf_d, _ = delay10_rows(Xnorm, confirmation, 9, T - 2)
    delay_model = fit_ridge_forecast(Xtr_d, ytr_d, ptr_d, Xva_d, yva_d, pva_d)
    delay_all = {}
    for split_name, Xs, ys, ps in (
        ("development", Xtr_d, ytr_d, ptr_d),
        ("validation", Xva_d, yva_d, pva_d),
        ("confirmation", Xcf_d, ycf_d, pcf_d),
    ):
        Xz = standardize_apply(Xs, delay_model["mu"], delay_model["sg"])
        pred = ridge_predict(Xz, delay_model["weights"])
        pooled, _, skipped = patient_pooled_nrmse(ys, pred, ps)
        delay_all[split_name] = {
            "nrmse_patient_pooled": pooled,
            "nrmse_concat_diagnostic": concat_nrmse(ys, pred),
            "skipped_zero_std": skipped,
            "n_pairs": int(len(ys)),
        }
        task_b_rows.append(
            {
                "task": "B",
                "feature": "DELAY10_5CH",
                "classifier": "ridge",
                "split": split_name,
                "lambda": delay_model["lambda"],
                "accuracy": "",
                "macro_auc": "",
                "nrmse_patient_pooled": pooled,
                "n": int(len(ys)),
            }
        )

    persist_aligned = {}
    for split_name, idx in (
        ("development", development),
        ("validation", validation),
        ("confirmation", confirmation),
    ):
        xp, yp, pp = persist_rows(Xnorm, idx, 9, T - 2)
        pooled, _, skipped = patient_pooled_nrmse(yp, xp, pp)
        persist_aligned[split_name] = {
            "nrmse_patient_pooled": pooled,
            "nrmse_concat_diagnostic": concat_nrmse(yp, xp),
            "skipped_zero_std": skipped,
        }
        task_b_rows.append(
            {
                "task": "B",
                "feature": "PERSISTENCE_ALIGNED_T9",
                "classifier": "x_t",
                "split": split_name,
                "lambda": "",
                "accuracy": "",
                "macro_auc": "",
                "nrmse_patient_pooled": pooled,
                "n": int(len(yp)),
            }
        )

    conf_persist = persist_all["confirmation"]["nrmse_patient_pooled"]
    conf_delay = delay_all["confirmation"]["nrmse_patient_pooled"]
    task_b_void = bool(conf_persist <= 0.35 or conf_delay <= 0.35)

    if task_a_void and task_b_void:
        overall = "TASK_VOID"
        e51 = False
        surviving = "none"
    elif task_a_void:
        overall = "PROTOCOL_FROZEN_NO_BSIM"
        e51 = True
        surviving = "B"
    elif task_b_void:
        overall = "PROTOCOL_FROZEN_NO_BSIM"
        e51 = True
        surviving = "A"
    else:
        overall = "PROTOCOL_FROZEN_NO_BSIM"
        e51 = True
        surviving = "A+B"

    if surviving in ("A", "A+B"):
        scout_dev = [int(i) for i in scout_dev_4_per_class]
        scout_conf = [int(i) for i in scout_conf_3_3_2]
        scout_rule = "4_per_Y_class_development; confirmation_3_3_2"
    else:
        scout_dev = [int(i) for i in scout_dev_first12]
        scout_conf = [int(i) for i in scout_conf_first8]
        scout_rule = "first_12_sorted_development; first_8_sorted_confirmation"

    # ----- write inventory -----
    inv = []
    inv.append("# E5.0 biomarker inventory")
    inv.append("")
    inv.append("Audit date: 2026-08-17. No BSim. Not a clinical prediction task.")
    inv.append("")
    inv.append("COHORT_STATUS: SYNTHETIC_COHORT")
    inv.append("")
    inv.append("This is a synthetic generator cohort. It is not a clinical")
    inv.append("dataset. No PhysioNet, MIMIC, or eICU records were used.")
    inv.append("")
    inv.append("## Vendored file")
    inv.append("")
    inv.append("| Item | Value |")
    inv.append("|---|---|")
    inv.append(
        f"| Source | `{MAT_SOURCE}` |"
    )
    inv.append(
        "| Vendor | `examples/HybridDish/external_analysis/"
        "synthetic_biomarker_BSim_inputs.mat` |"
    )
    inv.append(f"| Bytes | {MAT_VENDOR.stat().st_size} |")
    inv.append(f"| SHA-256 | `{mat_sha}` |")
    inv.append(f"| Source SHA-256 | `{src_sha}` |")
    inv.append("| Copy | byte-for-byte; not regenerated |")
    inv.append(
        f"| Generator SHA-256 | `{gen_sha}` "
        "(`MAT_biomarker_progression.m`) |"
    )
    inv.append("")
    inv.append("## Shape")
    inv.append("")
    inv.append("| Quantity | Value |")
    inv.append("|---|---|")
    inv.append(f"| N patients | {n} |")
    inv.append(f"| T windows | {t_windows} |")
    inv.append(f"| B channels | {b} |")
    inv.append("| Arrays | `X`, `Xnorm`, `Y`, `Q` |")
    inv.append(
        f"| Y counts | class0={y_counts[0]}, class1={y_counts[1]}, "
        f"class2={y_counts[2]} |"
    )
    inv.append(
        f"| Class blocking | {blocked} "
        "(0..332 class0, 333..665 class1, 666..999 class2) |"
    )
    inv.append(f"| Missingness | {missing_x} NaNs in X/Xnorm/Q; expect none |")
    inv.append("")
    inv.append("## Channel names and original units (before min-max)")
    inv.append("")
    inv.append("`MAT_biomarker_progression.m` names the five channels and")
    inv.append("sets generator baselines. It does **not** declare units.")
    inv.append("Observed ranges below are from vendored `X`, not from a")
    inv.append("lab instrument. Common clinical unit names are listed only")
    inv.append("as a scale resemblance; they are not measured units.")
    inv.append("")
    inv.append("| b | Name | Generator baseline | Floor | Observed min | Observed max | Units in MATLAB |")
    inv.append("|---|---|---|---|---|---|---|")
    floors = ["10", "5", "100", "0.8", "70"]
    bases = [
        "80 + 30*randn",
        "70 + 25*randn",
        "700 + 200*randn",
        "3.0 + 0.8*randn",
        "160 + 30*randn",
    ]
    resemble = [
        "CRP often mg/L",
        "IL-6 often pg/mL",
        "ferritin often ng/mL",
        "lactate often mmol/L",
        "glucose often mg/dL",
    ]
    for i, name in enumerate(CHANNEL_NAMES):
        s = channel_stats[i]
        inv.append(
            f"| {i} | {name} | `{bases[i]}` | {floors[i]} | "
            f"{s['min']:.6g} | {s['max']:.6g} | UNDECLARED "
            f"({resemble[i]}; not a measurement) |"
        )
    inv.append("")
    inv.append("`Xnorm` is a global min-max of each channel over all N×T")
    inv.append("values, including patients later placed in confirmation:")
    inv.append("")
    inv.append("```")
    inv.append("Xnorm(:,:,b) = (X(:,:,b) - minVal) ./ (maxVal - minVal + eps);")
    inv.append("```")
    inv.append("")
    inv.append("This audit uses stored `Xnorm` as specified. The global")
    inv.append("scale is a generator property, not a living-dish feature.")
    inv.append("`Q` is a leaky-filtered map of `Xnorm` (`alpha = 0.05`) and")
    inv.append("is not used in Module 2.")
    inv.append("")
    inv.append("## How Y was generated")
    inv.append("")
    inv.append("Quoted from `MAT_biomarker_progression.m`:")
    inv.append("")
    inv.append("```")
    inv.append("% Labels:")
    inv.append("% 0 = improving")
    inv.append("% 1 = stable")
    inv.append("% 2 = deteriorating")
    inv.append("")
    inv.append("    % Assign class")
    inv.append("    if i <= N/3")
    inv.append("        label = 0;          % improving")
    inv.append("    elseif i <= 2*N/3")
    inv.append("        label = 1;          % stable")
    inv.append("    else")
    inv.append("        label = 2;          % deteriorating")
    inv.append("    end")
    inv.append("")
    inv.append("    Y(i) = label;")
    inv.append("```")
    inv.append("")
    inv.append("Within-patient traces are class-conditional trends on every")
    inv.append("channel, plus independent Gaussian noise:")
    inv.append("")
    inv.append("```")
    inv.append("            case 0   % Improving")
    inv.append("                trend = startValue * (1 - 0.45*time);")
    inv.append("            case 1   % Stable")
    inv.append("                trend = startValue * (1 + 0.05*sin(2*pi*time));")
    inv.append("            case 2   % Deteriorating")
    inv.append("                trend = startValue * (1 + 0.45*time);")
    inv.append("        noise = 0.08 * startValue * randn(1,T);")
    inv.append("        X(i,:,b) = trend + noise;")
    inv.append("```")
    inv.append("")
    inv.append("`rng(1)` is set once at the top of the generator. `N=1000`,")
    inv.append("`T=48`, `B=5`, `time = linspace(0,1,T)`. `Y` is the trend")
    inv.append("class, not an external clinical endpoint.")
    inv.append("")
    inv.append("## PhysioNet / MIMIC / eICU")
    inv.append("")
    if phys_hits:
        inv.append("UNEXPECTED hits:")
        for h in phys_hits:
            inv.append(f"- `{h}`")
    else:
        inv.append("None. No PhysioNet, MIMIC, or eICU filename exists in")
        inv.append("`bsim_clean` or the `reservoir_new` external tree.")
    inv.append("")
    inv.append("## Track B relationship")
    inv.append("")
    inv.append("Different object. Track B is a 40-patient × 5-window")
    inv.append("synthetic HIGH/LOW classification on five HybridDish sites")
    inv.append("(AHL, acid, attA, attB, rep) with a field-only gate. Overall")
    inv.append("FAIL stays: the plumes classify the patient")
    inv.append("(`examples/HybridDish/README.md`;")
    inv.append("`examples/BSimReservoirPlanTrackB/PROTOCOL.md`).")
    inv.append("This E5.0 tensor is 1000 × 48 × 5 from")
    inv.append("`MAT_biomarker_progression.m`. Do not reopen Track B. Do")
    inv.append("not copy the Track B site table.")
    inv.append("")
    inv.append("## HybridDish-illegal maps")
    inv.append("")
    inv.append("| Map | Why illegal here |")
    inv.append("|---|---|")
    inv.append(
        "| glucose field | Channel 5 / `Q(:,:,5)` is labelled Nutrient/"
        "Glucose in the generator. HybridDish does not have a glucose PDE. |"
    )
    inv.append(
        "| acid-as-lactate analog | `Q(:,:,4)` maps Lactate to "
        "`k_str` stress/toxicity. Acid in this dish is a mortality "
        "channel, not an analog lactate wire. |"
    )
    inv.append(
        "| 5-AC Track B layout | Track B positions and rates are a "
        "closed FAIL object. They are not this cohort’s map. |"
    )
    inv.append("")
    inv.append("Official living readout remains 408 (`R`, `L`, tagged")
    inv.append("deaths). Any extra chemical would need its own field-only")
    inv.append("baseline. This audit does not freeze a 5-AC layout.")
    inv.append("")
    (RESULTS / "BIOMARKER_INVENTORY.md").write_text(
        "\n".join(inv) + "\n", encoding="utf-8"
    )

    # ----- CSV -----
    csv_path = RESULTS / "direct_input_baselines.csv"
    fieldnames = [
        "task",
        "feature",
        "classifier",
        "split",
        "lambda",
        "accuracy",
        "macro_auc",
        "nrmse_patient_pooled",
        "n",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in task_a_rows + task_b_rows:
            out = {}
            for k in fieldnames:
                v = row[k]
                if v == "":
                    out[k] = ""
                elif isinstance(v, float):
                    out[k] = f"{v:.10g}"
                else:
                    out[k] = v
            w.writerow(out)

    # ----- DIRECT_INPUT_BASELINES.md -----
    bmd = []
    bmd.append("# E5.0 direct-input baselines (no BSim)")
    bmd.append("")
    bmd.append("Fits on development only. Lambda on validation only.")
    bmd.append("Confirmation scored once, last. Split hashes were written")
    bmd.append("to `patient_split.json` before these numbers.")
    bmd.append("")
    bmd.append(f"Split seed `{SPLIT_SEED}`. Chance accuracy = 1/3.")
    bmd.append(
        f"Development majority class = {majority} "
        f"(counts {ycount(development)})."
    )
    bmd.append("")
    bmd.append("Ridge: intercept unregularized; population standardization")
    bmd.append("on development; grid `{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`.")
    bmd.append("Classification selects highest validation macro AUC, ties")
    bmd.append("larger lambda. Forecast selects lowest validation")
    bmd.append("patient-pooled NRMSE, ties larger lambda.")
    bmd.append("")
    bmd.append("## Task A — patient-level 3-class Y")
    bmd.append("")
    bmd.append("Void if confirmation accuracy ≥ 0.80 or macro AUC ≥ 0.90")
    bmd.append("on `SLOPE_XNORM` or `MEAN_XNORM` alone.")
    bmd.append("")
    bmd.append(
        f"**TASK_A_VOID: {'yes' if task_a_void else 'no'}.**"
    )
    if void_a_hits:
        bmd.append("")
        bmd.append("Voiding confirmation scores:")
        bmd.append("")
        bmd.append("| Feature | Classifier | Accuracy | Macro AUC |")
        bmd.append("|---|---|---|---|")
        for h in void_a_hits:
            bmd.append(
                f"| {h['feature']} | {h['classifier']} | "
                f"{h['accuracy']:.4f} | {h['macro_auc']:.4f} |"
            )
    bmd.append("")
    bmd.append("| Feature | Classifier | λ | Dev acc | Dev AUC | Val acc | Val AUC | Conf acc | Conf AUC |")
    bmd.append("|---|---|---|---|---|---|---|---|---|")
    order_a = [
        ("MAJORITY_CLASS", "constant"),
        ("MEAN_XNORM", "ridge_ovr"),
        ("MEAN_XNORM", "multinomial_logistic"),
        ("LAST_XNORM", "ridge_ovr"),
        ("LAST_XNORM", "multinomial_logistic"),
        ("SLOPE_XNORM", "ridge_ovr"),
        ("SLOPE_XNORM", "multinomial_logistic"),
        ("MEAN_SLOPE", "ridge_ovr"),
        ("MEAN_SLOPE", "multinomial_logistic"),
    ]
    lookup = {
        (r["feature"], r["classifier"], r["split"]): r for r in task_a_rows
    }
    for feat, clf in order_a:
        lam = lookup[(feat, clf, "confirmation")].get("lambda", "")
        lam_s = f"{lam:g}" if lam != "" else "—"
        def cell(split, key):
            return fmt(lookup[(feat, clf, split)][key])
        bmd.append(
            f"| {feat} | {clf} | {lam_s} | "
            f"{cell('development','accuracy')} | {cell('development','macro_auc')} | "
            f"{cell('validation','accuracy')} | {cell('validation','macro_auc')} | "
            f"{cell('confirmation','accuracy')} | {cell('confirmation','macro_auc')} |"
        )
    bmd.append("")
    bmd.append("`Y` is the generator’s trend class. A linear slope of the")
    bmd.append("five `Xnorm` channels is a subject-legal map of the raw")
    bmd.append("channels, not a living-dish computation. Multinomial")
    bmd.append("logistic AUC uses class probabilities. Ridge one-vs-rest")
    bmd.append("argmax can collapse to the majority class at large λ even")
    bmd.append("when ranking AUC is above chance; the void uses either")
    bmd.append("allowed classifier. `LAST_XNORM` is reported but is not a")
    bmd.append("void criterion.")
    bmd.append("")
    bmd.append("## Task B — horizon-1 forecast of channel 0 (CRP / Xnorm)")
    bmd.append("")
    bmd.append("Void if confirmation persistence NRMSE ≤ 0.35 or delay-10")
    bmd.append("NRMSE ≤ 0.35. Metric is patient-pooled NRMSE (mean of")
    bmd.append("per-patient NRMSE). Concat NRMSE is a diagnostic only.")
    bmd.append("")
    bmd.append(
        f"**TASK_B_VOID: {'yes' if task_b_void else 'no'}.**"
    )
    bmd.append("")
    bmd.append(
        f"Confirmation persistence NRMSE = {conf_persist:.4f}. "
        f"Confirmation delay-10 NRMSE = {conf_delay:.4f}. "
        f"Delay-10 λ = {delay_model['lambda']:g}."
    )
    bmd.append("")
    bmd.append("| Baseline | Dev pooled | Val pooled | Conf pooled | Conf concat (diagnostic) |")
    bmd.append("|---|---|---|---|---|")
    bmd.append(
        f"| persistence `t=0..T-2` | "
        f"{persist_all['development']['nrmse_patient_pooled']:.4f} | "
        f"{persist_all['validation']['nrmse_patient_pooled']:.4f} | "
        f"{persist_all['confirmation']['nrmse_patient_pooled']:.4f} | "
        f"{persist_all['confirmation']['nrmse_concat_diagnostic']:.4f} |"
    )
    bmd.append(
        f"| train-mean intercept | "
        f"{intercept_all['development']['nrmse_patient_pooled']:.4f} | "
        f"{intercept_all['validation']['nrmse_patient_pooled']:.4f} | "
        f"{intercept_all['confirmation']['nrmse_patient_pooled']:.4f} | "
        f"{intercept_all['confirmation']['nrmse_concat_diagnostic']:.4f} |"
    )
    bmd.append(
        f"| persistence aligned `t=9..T-2` | "
        f"{persist_aligned['development']['nrmse_patient_pooled']:.4f} | "
        f"{persist_aligned['validation']['nrmse_patient_pooled']:.4f} | "
        f"{persist_aligned['confirmation']['nrmse_patient_pooled']:.4f} | "
        f"{persist_aligned['confirmation']['nrmse_concat_diagnostic']:.4f} |"
    )
    bmd.append(
        f"| linear delay-10 (50-D ridge) | "
        f"{delay_all['development']['nrmse_patient_pooled']:.4f} | "
        f"{delay_all['validation']['nrmse_patient_pooled']:.4f} | "
        f"{delay_all['confirmation']['nrmse_patient_pooled']:.4f} | "
        f"{delay_all['confirmation']['nrmse_concat_diagnostic']:.4f} |"
    )
    bmd.append("")
    bmd.append("No delay crosses a patient bound. Delay-10 uses 38 pairs")
    bmd.append("per patient (`t=9..46`). Persistence void uses 47 pairs")
    bmd.append("(`t=0..46`). Concat NRMSE is smaller because between-patient")
    bmd.append("baseline spread inflates the denominator; it is **not** the")
    bmd.append("gate. `reservoir_new` forecast scores were not copied.")
    bmd.append("")
    (RESULTS / "DIRECT_INPUT_BASELINES.md").write_text(
        "\n".join(bmd) + "\n", encoding="utf-8"
    )

    # ----- RESET_ISOLATED_DESIGN.md -----
    rmd = []
    rmd.append("# E5.0 isolated vs sequential design freeze")
    rmd.append("")
    rmd.append("Paper freeze only. No living Java. No BSim in this audit.")
    rmd.append("HybridDish times, not `reservoir_new` 6 s windows.")
    rmd.append("")
    rmd.append("## Timescales")
    rmd.append("")
    rmd.append("| Quantity | Value |")
    rmd.append("|---|---|")
    rmd.append("| analysis window | 300 s |")
    rmd.append("| pulse (if a later map uses the claim dish) | 75 s |")
    rmd.append("| `tau_L` | 1500 s |")
    rmd.append("| warmup | 18000 s (same as the claim dish) |")
    rmd.append("| T analysis windows per patient | 48 |")
    rmd.append(
        "| E0 surrogate, `L` to 1% | on the order of 25–30 zero windows "
        "(`E0_PREFLIGHT_REPORT.md` §8: `L` 1% at 7395 s; residual after "
        "30 zero windows 0.003) |"
    )
    rmd.append("| 5 × 6 s reset | **forbidden** as a HybridDish protocol |")
    rmd.append("")
    rmd.append("Ideal exponential `L` carryover after 5 / 10 / 15 zero")
    rmd.append("windows is about 37% / 14% / 5% (planning estimates, not")
    rmd.append("acceptance evidence). Old `reservoir_new` reset-5 crashed")
    rmd.append("after 58/100 patients; reset-10/15 did not run; five 6 s")
    rmd.append("windows retained ~97–99% state norm; smoke NRMSE ~1517 is")
    rmd.append("unusable. Port the method, not those scores.")
    rmd.append("")
    rmd.append("## Primary design: ISOLATED_PATIENT")
    rmd.append("")
    rmd.append("One patient per BSim. Fresh init. Same warmup as the claim")
    rmd.append("dish (18000 s), then T=48 analysis windows. No sequential")
    rmd.append("batching. Patient identity is the run. Windows from patient")
    rmd.append("`p` never share a dish with patient `q`.")
    rmd.append("")
    rmd.append("## Sequential batching: NOT authorized")
    rmd.append("")
    rmd.append("A later living reset study may compare sequential vs")
    rmd.append("isolated on the **same hashed patients** from")
    rmd.append("`patient_split.json`. That study is not this audit.")
    rmd.append("Required residuals then (not now):")
    rmd.append("")
    rmd.append("- AHL, `R`, `L`, occupancy, population")
    rmd.append("- reset rows excluded from readout fitting and delayed features")
    rmd.append("- no material patient-position drift")
    rmd.append("- paired isolated reference for every sequential patient")
    rmd.append("")
    rmd.append("A 100-patient sequential dish is forbidden here.")
    rmd.append("")
    rmd.append("## Chemical map: constraints only, not a freeze")
    rmd.append("")
    rmd.append("Do **not** freeze a 5-AC HybridDish layout in this audit.")
    rmd.append("")
    rmd.append("- no glucose field")
    rmd.append("- no acid-as-lactate analog wire")
    rmd.append("- no Track B site table")
    rmd.append("- official living readout remains 408 (`R`, `L`, tagged deaths)")
    rmd.append("- any extra chemical needs its own field-only baseline")
    rmd.append("")
    rmd.append("A later E5.1 prompt may pick a legal one-AHL or two-chemical")
    rmd.append("map only if a task survives Module 2. Surviving Task B is a")
    rmd.append("horizon-1 CRP forecast. A later map may encode CRP (and at")
    rmd.append("most one other legal chemical) through the claim-dish AHL")
    rmd.append("pathway. It may not inject glucose, wire acid as lactate,")
    rmd.append("or copy the Track B site table. The five-channel delay-10")
    rmd.append("baseline remains an analysis comparison, not a dish layout.")
    rmd.append("")
    rmd.append("## Affordable isolated scout (not run)")
    rmd.append("")
    if overall == "TASK_VOID":
        rmd.append(
            "Overall `TASK_VOID`. No living scout is authorized. The"
        )
        rmd.append("predeclared draw is recorded only so it is not reinvented.")
    else:
        rmd.append(
            f"Surviving task: {surviving}. Scout authorized for a later "
            "E5.1 prompt only, seed 111, isolated patients, no 5-AC map."
        )
        rmd.append("Not run in this audit.")
    rmd.append("")
    rmd.append(f"Draw rule: `{scout_rule}`.")
    rmd.append("")
    if surviving == "B":
        rmd.append("Task A voided, so the scout is the first 12 sorted")
        rmd.append("development indices and first 8 sorted confirmation")
        rmd.append("indices, not 4-per-class. Under the blocked generator")
        rmd.append("labels those low indices are class 0. Task B does not")
        rmd.append("use `Y`.")
        rmd.append("")
    rmd.append("| Arm | n | Patient indices (0-based) |")
    rmd.append("|---|---|---|")
    rmd.append(
        f"| development scout | 12 | {','.join(str(i) for i in scout_dev)} |"
    )
    rmd.append(
        f"| confirmation scout | 8 | {','.join(str(i) for i in scout_conf)} |"
    )
    rmd.append("| bacterial seed | 111 only | — |")
    rmd.append("")
    rmd.append("SHA-256 of concatenated scout development indices:")
    rmd.append(f"`{sha256_text(index_string(scout_dev))}`")
    rmd.append("")
    rmd.append("SHA-256 of concatenated scout confirmation indices:")
    rmd.append(f"`{sha256_text(index_string(scout_conf))}`")
    rmd.append("")
    rmd.append("Do not run these patients in E5.0.")
    rmd.append("")
    (RESULTS / "RESET_ISOLATED_DESIGN.md").write_text(
        "\n".join(rmd) + "\n", encoding="utf-8"
    )

    # ----- verdict -----
    v = []
    v.append("# E5.0 verdict")
    v.append("")
    v.append("Date: 2026-08-17. Analysis only. No Java. No BSim.")
    v.append("Track B Overall FAIL is not rewritten.")
    v.append("")
    v.append("| Item | Value |")
    v.append("|---|---|")
    v.append("| COHORT_STATUS | SYNTHETIC_COHORT |")
    v.append(f"| TASK_A_VOID | {'yes' if task_a_void else 'no'} |")
    v.append(f"| TASK_B_VOID | {'yes' if task_b_void else 'no'} |")
    v.append(f"| OVERALL | {overall} |")
    v.append(
        f"| E5.1 living BSim authorized | "
        f"{'yes' if e51 else 'no'} |"
    )
    v.append(f"| Surviving task | {surviving} |")
    v.append("")
    v.append("## Why")
    v.append("")
    if task_a_void:
        v.append(
            "Task A is void because a subject-legal map of the raw"
        )
        v.append(
            "channels already classifies the synthetic trend label `Y`."
        )
        for h in void_a_hits:
            v.append(
                f"Confirmation `{h['feature']}` / `{h['classifier']}` "
                f"accuracy {h['accuracy']:.4f}, macro AUC {h['macro_auc']:.4f}."
            )
        v.append(
            "`Y` is generated as improving / stable / deteriorating trends "
            "in `MAT_biomarker_progression.m`. This is not a clinical endpoint."
        )
    else:
        v.append(
            "Task A did not meet the predeclared void on confirmation "
            "`MEAN_XNORM` or `SLOPE_XNORM`."
        )
    v.append("")
    if task_b_void:
        v.append(
            f"Task B is void: confirmation persistence NRMSE {conf_persist:.4f} "
            f"and/or delay-10 NRMSE {conf_delay:.4f} ≤ 0.35."
        )
    else:
        v.append(
            f"Task B is not void: confirmation persistence NRMSE "
            f"{conf_persist:.4f}, delay-10 NRMSE {conf_delay:.4f}, "
            "both > 0.35 (patient-pooled)."
        )
    v.append("")
    if overall == "TASK_VOID":
        v.append(
            "Both tasks void. No living HybridDish biomarker protocol is"
        )
        v.append(
            "frozen. E5.1 is not authorized. Do not invent a new clinical"
        )
        v.append("target to escape the void.")
    else:
        v.append(
            f"Keep task {surviving}. Isolated-patient design is frozen in"
        )
        v.append(
            "`RESET_ISOLATED_DESIGN.md`. Sequential batching is not"
        )
        v.append(
            "authorized. No 5-AC map is frozen. E5.1 may pick a legal"
        )
        v.append(
            "one-AHL or two-chemical map only in a later prompt, seed 111"
        )
        v.append("scout, hashed patients above. Not run here.")
    v.append("")
    v.append("Confirmation was not used to choose the task. Kinetics,")
    v.append("layout, flow, clamp, and mortality were not retuned.")
    v.append("")
    (RESULTS / "E5_0_VERDICT.md").write_text("\n".join(v) + "\n", encoding="utf-8")

    print("MAT_SHA256", mat_sha)
    print("TASK_A_VOID", task_a_void)
    print("TASK_B_VOID", task_b_void)
    print("OVERALL", overall)
    print("E5.1", e51)
    print("surviving", surviving)
    print("conf_persist", conf_persist)
    print("conf_delay", conf_delay)
    print("majority", majority)
    for h in void_a_hits:
        print("void_hit", h)
    print("scout_dev", scout_dev)
    print("scout_conf", scout_conf)
    print("y_dev", ycount(development))
    print("y_val", ycount(validation))
    print("y_conf", ycount(confirmation))
    print("split hashes", split_doc["index_string_sha256"])
    print("DONE")


if __name__ == "__main__":
    main()
