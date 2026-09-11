#!/usr/bin/env python3
"""Score Lane A seeds 101/202/303. Does not rewrite parent standings."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys_path_extra = (
    HERE,
    HERE.parent / "BSimReservoirPlanWaveform2c",
)
import sys

for p in sys_path_extra:
    sys.path.insert(0, str(p))

from check_lane_a_carrier import (  # noqa: E402
    TEST as CARRIER_TEST,
    TRAIN as CARRIER_TRAIN,
    apply_std,
    col_block,
    nrmse,
    ridge_fit,
    ridge_predict,
    standardize_fit,
    u_of,
    windows_of,
)
from check_lane_a_carrier import fit_score as carrier_fit_score  # noqa: E402
from check_lane_a_narma10 import (  # noqa: E402
    N_BINS,
    SAMPLES_PER_WINDOW,
    concat,
    fit_eval,
    load_u,
    narma10,
    sha256_u,
)
from check_lane_a_mg import load_target as load_mg_target  # noqa: E402
from check_lane_a_lorenz import ar10, persist_features  # noqa: E402
from check_lane_a_lorenz import load_target as load_lorenz_target  # noqa: E402
import check_waveform2c as W2C  # noqa: E402

RESULTS = HERE / "results"
OUT = RESULTS / "seed_replicate"
PROTOCOL = HERE / "PROTOCOL_SEED_REPLICATE.md"
PROTOCOL_JSON = HERE / "configs" / "seed_replicate_protocol.json"
STANDING = ROOT / "examples" / "PocketDish" / "LANE_A_SEED_REPLICATE_STANDING.md"

SEEDS = (101, 202, 303)
DELTA = 0.02
BRIER_DELTA = 0.02
NARMA_SHA = "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b"
MG_SHA = "e06810ea70f61bbe4baf61939f93b1e1cae9d92d9eee8d757865c4935dbac780"
LORENZ_SHA = "69b861f256fb39527c037b994885d119ccc8ab2873a3d2bca59f7a8802c44aff"
WAVE_SHA = "0979090c051c513dce14632e83ae9bf025c4b666e497dc879ae1db4bba30d465"

AXIS = {
    "carrier": "field",
    "mg": "field",
    "narma10": "living",
    "lorenz": "living",
    "waveform": "living",
}


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
        out.append([sum(block[t][j] for t in range(SAMPLES_PER_WINDOW)) / SAMPLES_PER_WINDOW for j in range(p)])
    return out


def driven_path(seed: int, name: str) -> Path:
    if seed == 101:
        return RESULTS / name
    return RESULTS / f"seed_{seed}" / name


def load_rows(path: Path, label: str):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if label not in names:
            raise SystemExit(f"{path} missing {label}")
        return list(reader)


def fit_eval_trace(X_all, y_all):
    pack = fit_eval(X_all, y_all)
    X = X_all[40:]
    y = y_all[40:]
    Xtr, ytr = X[:110], y[:110]
    Xte, yte = X[110:], y[110:]
    from check_lane_a_narma10 import select_lambda

    lam, _ = select_lambda(X, y)
    mu, sg = standardize_fit(Xtr)
    w = ridge_fit(apply_std(Xtr, mu, sg), ytr, lam)
    pred = ridge_predict(apply_std(Xte, mu, sg), w)
    pack["test_y"] = list(yte)
    pack["test_pred"] = [float(v) for v in pred]
    pack["test_n"] = list(range(150, 200))
    return pack


def carrier_trace(rows, prefix):
    y = u_of(rows)
    windows = windows_of(rows)
    X = [col_block(r, prefix) for r in rows]
    pack = carrier_fit_score(X, y, windows)
    train_idx = [i for i, w in enumerate(windows) if w in CARRIER_TRAIN]
    test_idx = [i for i, w in enumerate(windows) if w in CARRIER_TEST]
    Xtr = [X[i] for i in train_idx]
    ytr = [y[i] for i in train_idx]
    Xte = [X[i] for i in test_idx]
    yte = [y[i] for i in test_idx]
    mu, sg = standardize_fit(Xtr)
    w = ridge_fit(apply_std(Xtr, mu, sg), ytr, pack["lambda"])
    pred = ridge_predict(apply_std(Xte, mu, sg), w)
    pack["test_y"] = yte
    pack["test_pred"] = [float(v) for v in pred]
    pack["test_n"] = test_idx
    return pack


def softmax(scores):
    z = scores - scores.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def waveform_pack(X_windows, y_window):
    X_windows = np.asarray(X_windows, dtype=float)
    y_window = np.asarray(y_window, dtype=int)
    metrics = W2C.evaluate_window_readout(X_windows, y_window, list(range(W2C.NUM_WINDOWS)))
    X = X_windows[W2C.WASHOUT :]
    y = y_window[W2C.WASHOUT :]
    lam = metrics["lambda"]
    Xtr, ytr = X[: W2C.TRAIN], y[: W2C.TRAIN]
    Xte, yte = X[W2C.TRAIN :], y[W2C.TRAIN :]
    (Xtr_z, Xte_z), _, _ = W2C.standardize(Xtr, Xte)
    test_scores, _ = W2C.ovr_scores(Xtr_z, ytr, Xte_z, lam)
    pooled = []
    yb = []
    for block in range(W2C.TEST_BLOCKS):
        sl = slice(block * 8, (block + 1) * 8)
        pooled.append(test_scores[sl].mean(axis=0))
        labels = np.unique(yte[sl])
        if len(labels) != 1:
            raise SystemExit("waveform test block is not a constant label")
        yb.append(int(labels[0]))
    pooled = np.vstack(pooled)
    probs = softmax(pooled)
    onehot = np.eye(3)[yb]
    brier = float(np.mean(np.sum((probs - onehot) ** 2, axis=1)))
    pred = [int(v) for v in np.argmax(pooled, axis=1)]
    confusion = np.zeros((3, 3), dtype=int)
    for a, b in zip(yb, pred):
        confusion[a, b] += 1
    return {
        "block_macro_ovr_auc": float(metrics["block_macro_ovr_auc"]),
        "window_macro_ovr_auc": float(metrics["test_window_macro_ovr_auc"]),
        "block_accuracy": float(metrics["block_accuracy"]),
        "block_brier": brier,
        "lambda": float(lam),
        "block_true": yb,
        "block_pred": pred,
        "block_probs": probs.tolist(),
        "confusion": confusion.tolist(),
    }


def write_trace(path: Path, n, y, **preds) -> None:
    keys = ["n", "y"] + list(preds)
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=keys)
        w.writeheader()
        for i, ni in enumerate(n):
            row = {"n": ni, "y": y[i]}
            for k, series in preds.items():
                row[k] = series[i]
            w.writerow(row)


def score_nrmse_task(seed, driven_name, label, y, silent_rows, silent_label):
    driven = load_rows(driven_path(seed, driven_name), label)
    ahl = window_means(driven, "AHL", 200)
    rmap = window_means(driven, "R", 200)
    lmap = window_means(driven, "L", 200)
    s_r = window_means(silent_rows, "R", 200)
    s_l = window_means(silent_rows, "L", 200)
    rl = fit_eval_trace(concat(rmap, lmap), y)
    field = fit_eval_trace(ahl, y)
    silent = fit_eval(concat(s_r, s_l), y)
    return {
        "RL": rl["test_nrmse"],
        "FIELD": field["test_nrmse"],
        "SILENT_RL": silent["test_nrmse"],
        "trace": {
            "n": rl["test_n"],
            "y": rl["test_y"],
            "RL": rl["test_pred"],
            "FIELD": field["test_pred"],
        },
    }


def require_frozen() -> None:
    text = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in text or '"frozen_before_traces": true' not in js:
        raise SystemExit("seed-replicate PROTOCOL not frozen")
    if '"waveform_primary": "block_brier"' not in js:
        raise SystemExit("waveform primary must be block_brier")


def write_standing(summary: dict) -> None:
    lines = [
        "# Lane A seed-replicate standing",
        "",
        f"**Status: {summary['LANE_A_SEED_REPLICATE']}**",
        f"**Date:** {summary['date']}",
        "**Label:** `LANE_A_SEED_REPLICATE`",
        "**Does not rewrite:** occupancy, carrier, audit, NARMA, MG, Lorenz, Waveform parent standings.",
        "",
        "Protocol (frozen before 202/303 traces):",
        "[`../LaneA_OccupiedMillimetre/PROTOCOL_SEED_REPLICATE.md`](../LaneA_OccupiedMillimetre/PROTOCOL_SEED_REPLICATE.md).",
        "",
        "## Verdict",
        "",
        f"- Sign matches: **{summary['n_sign_match']}/5** (need ≥4).",
        f"- Mean signed NRMSE margin: **{summary['mean_nrmse_margin']:.4f}** (need > 0.02).",
        f"- Waveform axis uses **Brier**, not saturated AUC.",
        "",
        "## Per-task means (seeds 101/202/303)",
        "",
        "| Task | Living mean | Field | Expected | Match |",
        "|---|---|---|---|---|",
    ]
    for key, title in (
        ("carrier", "Carrier"),
        ("mg", "MG"),
        ("narma10", "NARMA-10"),
        ("lorenz", "Lorenz"),
        ("waveform", "Waveform"),
    ):
        row = summary["tasks"][key]
        lines.append(
            f"| {title} | {row['living_mean']:.4f} | {row['field']:.4f} | "
            f"{row['expected']} | {row['match']} |"
        )
    lines.extend(
        [
            "",
            "Field is seed-invariant. Living is the mean of three placement seeds.",
            "Waveform living/field columns are **Brier** (lower better). Others are NRMSE.",
            "Parent seed-101 numbers are not rewritten.",
            "",
        ]
    )
    STANDING.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    require_frozen()
    OUT.mkdir(parents=True, exist_ok=True)

    u_narma = load_u(HERE / "input_u_narma200.txt")
    if sha256_u(u_narma) != NARMA_SHA:
        raise SystemExit("NARMA u hash drifted")
    y_narma = narma10(u_narma)[1:]
    _u_mg, _x_mg, y_mg = load_mg_target(HERE / "mg_target.csv")
    _u_lor, _x_lor, y_lor = load_lorenz_target(HERE / "lorenz_target.csv")
    if sha256_u(load_u(HERE / "input_u_mg200.txt")) != MG_SHA:
        raise SystemExit("MG u hash drifted")
    if sha256_u(load_u(HERE / "input_u_lorenz_k10.txt")) != LORENZ_SHA:
        raise SystemExit("Lorenz u hash drifted")

    labels = list(csv.DictReader((HERE / "waveform2c_labels.csv").open(encoding="utf-8"), delimiter=";"))
    y_wave = [int(r["class"]) for r in labels]
    u_wave = load_u(HERE / "input_u_waveform2c_400.txt")
    if sha256_u(u_wave) != WAVE_SHA:
        raise SystemExit("Waveform u hash drifted")

    silent_narma = load_rows(RESULTS / "java_LANE_A_NARMA10_SILENT.csv", "LANE_A_NARMA10")
    silent_carrier = load_rows(RESULTS / "java_LANE_A_CARRIER_SILENT.csv", "LANE_A_CARRIER")
    silent_wave = load_rows(RESULTS / "java_LANE_A_WAVEFORM_SILENT.csv", "LANE_A_WAVEFORM")
    s_r_w = window_means(silent_wave, "R", 400)
    s_l_w = window_means(silent_wave, "L", 400)

    per_seed = {seed: {} for seed in SEEDS}

    for seed in SEEDS:
        carrier_rows = load_rows(driven_path(seed, "java_LANE_A_CARRIER_DRIVEN.csv"), "LANE_A_CARRIER")
        c_field = carrier_trace(carrier_rows, "AHL")
        c_r = carrier_trace(carrier_rows, "R")
        c_l = carrier_trace(carrier_rows, "L")
        per_seed[seed]["carrier"] = {
            "R": c_r["test_nrmse"],
            "L": c_l["test_nrmse"],
            "FIELD": c_field["test_nrmse"],
            "trace": {
                "n": c_r["test_n"],
                "y": c_r["test_y"],
                "R": c_r["test_pred"],
                "L": c_l["test_pred"],
                "FIELD": c_field["test_pred"],
            },
        }
        write_trace(
            OUT / f"trace_carrier_seed{seed}.csv",
            c_r["test_n"],
            c_r["test_y"],
            R=c_r["test_pred"],
            L=c_l["test_pred"],
            FIELD=c_field["test_pred"],
        )

        nar = score_nrmse_task(
            seed, "java_LANE_A_NARMA10_DRIVEN.csv", "LANE_A_NARMA10", y_narma, silent_narma, "LANE_A_NARMA10"
        )
        per_seed[seed]["narma10"] = nar
        write_trace(OUT / f"trace_narma10_seed{seed}.csv", nar["trace"]["n"], nar["trace"]["y"],
                    RL=nar["trace"]["RL"], FIELD=nar["trace"]["FIELD"])

        mg = score_nrmse_task(seed, "java_LANE_A_MG_DRIVEN.csv", "LANE_A_MG", y_mg, silent_narma, "LANE_A_NARMA10")
        per_seed[seed]["mg"] = mg
        write_trace(OUT / f"trace_mg_seed{seed}.csv", mg["trace"]["n"], mg["trace"]["y"],
                    RL=mg["trace"]["RL"], FIELD=mg["trace"]["FIELD"])

        lor = score_nrmse_task(
            seed, "java_LANE_A_LORENZ_DRIVEN.csv", "LANE_A_LORENZ", y_lor, silent_narma, "LANE_A_NARMA10"
        )
        per_seed[seed]["lorenz"] = lor
        write_trace(OUT / f"trace_lorenz_seed{seed}.csv", lor["trace"]["n"], lor["trace"]["y"],
                    RL=lor["trace"]["RL"], FIELD=lor["trace"]["FIELD"])

        driven_w = load_rows(driven_path(seed, "java_LANE_A_WAVEFORM_DRIVEN.csv"), "LANE_A_WAVEFORM")
        ahl = window_means(driven_w, "AHL", 400)
        rmap = window_means(driven_w, "R", 400)
        lmap = window_means(driven_w, "L", 400)
        w_rl = waveform_pack(concat(rmap, lmap), y_wave)
        w_field = waveform_pack(ahl, y_wave)
        w_silent = waveform_pack(concat(s_r_w, s_l_w), y_wave)
        per_seed[seed]["waveform"] = {
            "RL_brier": w_rl["block_brier"],
            "FIELD_brier": w_field["block_brier"],
            "SILENT_brier": w_silent["block_brier"],
            "RL_auc": w_rl["block_macro_ovr_auc"],
            "FIELD_auc": w_field["block_macro_ovr_auc"],
            "RL_window_auc": w_rl["window_macro_ovr_auc"],
            "FIELD_window_auc": w_field["window_macro_ovr_auc"],
            "confusion": w_rl["confusion"],
        }
        write_trace(
            OUT / f"trace_waveform_seed{seed}.csv",
            list(range(34, 50)),
            w_rl["block_true"],
            pred=w_rl["block_pred"],
            brier_contrib=[float(np.sum((np.array(p) - np.eye(3)[t]) ** 2)) for p, t in zip(w_rl["block_probs"], w_rl["block_true"])],
        )

    def mean_of(task, key):
        return float(np.mean([per_seed[s][task][key] for s in SEEDS]))

    tasks = {}
    # Carrier: living = R (closer than L); field wins if field <= R and field <= L
    r_mean = mean_of("carrier", "R")
    l_mean = mean_of("carrier", "L")
    f_car = mean_of("carrier", "FIELD")
    car_match = f_car <= r_mean and f_car <= l_mean
    tasks["carrier"] = {
        "living_mean": r_mean,
        "L_mean": l_mean,
        "field": f_car,
        "expected": "field",
        "match": bool(car_match),
        "signed_margin": r_mean - f_car,
        "per_seed": {str(s): {k: per_seed[s]["carrier"][k] for k in ("R", "L", "FIELD")} for s in SEEDS},
    }

    mg_rl = mean_of("mg", "RL")
    mg_f = mean_of("mg", "FIELD")
    tasks["mg"] = {
        "living_mean": mg_rl,
        "field": mg_f,
        "expected": "field",
        "match": bool(mg_f <= mg_rl),
        "signed_margin": mg_rl - mg_f,
        "per_seed": {str(s): {k: per_seed[s]["mg"][k] for k in ("RL", "FIELD", "SILENT_RL")} for s in SEEDS},
    }

    na_rl = mean_of("narma10", "RL")
    na_f = mean_of("narma10", "FIELD")
    tasks["narma10"] = {
        "living_mean": na_rl,
        "field": na_f,
        "expected": "living",
        "match": bool(na_rl <= na_f - DELTA),
        "signed_margin": na_f - na_rl,
        "per_seed": {str(s): {k: per_seed[s]["narma10"][k] for k in ("RL", "FIELD", "SILENT_RL")} for s in SEEDS},
    }

    lo_rl = mean_of("lorenz", "RL")
    lo_f = mean_of("lorenz", "FIELD")
    tasks["lorenz"] = {
        "living_mean": lo_rl,
        "field": lo_f,
        "expected": "living",
        "match": bool(lo_rl <= lo_f - DELTA),
        "signed_margin": lo_f - lo_rl,
        "per_seed": {str(s): {k: per_seed[s]["lorenz"][k] for k in ("RL", "FIELD", "SILENT_RL")} for s in SEEDS},
    }

    w_rl = mean_of("waveform", "RL_brier")
    w_f = mean_of("waveform", "FIELD_brier")
    tasks["waveform"] = {
        "living_mean": w_rl,
        "field": w_f,
        "expected": "living",
        "match": bool(w_rl <= w_f - BRIER_DELTA),
        "legacy_auc_rl": mean_of("waveform", "RL_auc"),
        "legacy_auc_field": mean_of("waveform", "FIELD_auc"),
        "window_auc_rl": mean_of("waveform", "RL_window_auc"),
        "window_auc_field": mean_of("waveform", "FIELD_window_auc"),
        "per_seed": {
            str(s): {k: per_seed[s]["waveform"][k] for k in (
                "RL_brier", "FIELD_brier", "RL_auc", "FIELD_auc", "RL_window_auc", "FIELD_window_auc"
            )}
            for s in SEEDS
        },
    }

    n_match = sum(1 for t in tasks.values() if t["match"])
    mean_margin = float(np.mean([
        tasks["carrier"]["signed_margin"],
        tasks["mg"]["signed_margin"],
        tasks["narma10"]["signed_margin"],
        tasks["lorenz"]["signed_margin"],
    ]))
    holds = n_match >= 4 and mean_margin > DELTA
    verdict = "AXIS_HOLDS" if holds else "AXIS_SPLITS"

    from datetime import date

    summary = {
        "gate": "LaneA_SEED_REPLICATE",
        "LANE_A_SEED_REPLICATE": verdict,
        "date": str(date.today()),
        "seeds": list(SEEDS),
        "n_sign_match": n_match,
        "mean_nrmse_margin": mean_margin,
        "parent_rewrite": False,
        "waveform_primary": "block_brier",
        "waveform_auc_not_axis": True,
        "tasks": tasks,
        "axis": AXIS,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_standing(summary)

    print("LANE_A_SEED_REPLICATE " + verdict)
    print(f"sign_match={n_match}/5 mean_nrmse_margin={mean_margin:.4f}")
    for key in ("carrier", "mg", "narma10", "lorenz", "waveform"):
        t = tasks[key]
        print(
            f"  {key} living={t['living_mean']:.4f} field={t['field']:.4f} "
            f"expected={t['expected']} match={t['match']}"
        )
        if key == "waveform":
            print(
                f"    legacy_block_AUC RL={t['legacy_auc_rl']:.4f} field={t['legacy_auc_field']:.4f} "
                f"window_AUC RL={t['window_auc_rl']:.4f} field={t['window_auc_field']:.4f}"
            )
    print("Parents unchanged. Waveform axis is Brier, not AUC 1.000.")
    print("standing " + str(STANDING))


if __name__ == "__main__":
    main()
