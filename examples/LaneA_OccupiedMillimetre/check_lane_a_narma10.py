#!/usr/bin/env python3
"""Lane A NARMA-10 system checker. Not Narma10b. Does not reopen carrier."""

from __future__ import annotations

import csv
import hashlib
import json
import random
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

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_NARMA10.md"
PROTOCOL_JSON = HERE / "configs" / "narma10_protocol.json"
U_FILE = HERE / "input_u_narma200.txt"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"
AUDIT_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_AUDIT_STANDING.md"

WASHOUT = 40
TRAIN = 110
INNER_VAL = 22
NUM_WINDOWS = 200
SAMPLES_PER_WINDOW = 16
LAMBDA_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
DELTA = 0.02
U_SHA = "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b"
NARMA10B_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
OCCUPIED = (
    63, 64, 65, 66, 73, 74, 75, 76, 83, 84, 85, 86, 93, 94, 95, 96,
    103, 104, 105, 106, 113, 114, 115, 116, 123, 124, 125, 126, 133, 134, 135, 136,
)
J_MAX = 128000000.0
PULSE_S = 75.0
N0_REL = 1e-11
LEDGER_REL = 1e-3
ALIVE_MIN = 0.05
RNG_SEED = 20260828


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("charc", "lorenz", "mackey", "fig4b")):
        raise SystemExit("Lane A NARMA-10 refuses CHARC/Lorenz/MG/Fig4b")


def sha256_u(u) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def load_u(path: Path):
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(float(line))
    return values


def narma10(u):
    y = [0.0] * (len(u) + 1)
    for t, u_t in enumerate(u):
        acc = sum(y[t - i] if t - i >= 0 else 0.0 for i in range(10))
        u_lag = u[t - 9] if t >= 9 else 0.0
        y[t + 1] = 0.3 * y[t] + 0.05 * y[t] * acc + 1.5 * u_lag * u_t + 0.1
    return y


def delay_u10(u):
    X = []
    for n in range(len(u)):
        feat = []
        for lag in range(10):
            j = n - lag
            feat.append(0.0 if j < 0 else u[j])
        X.append(feat)
    return X


def window_means(rows, prefix):
    by_w = {}
    for row in rows:
        w = int(row["window"])
        by_w.setdefault(w, []).append([float(row[f"{prefix}_{i}"]) for i in range(N_BINS)])
    out = []
    for w in range(NUM_WINDOWS):
        block = by_w.get(w)
        if not block or len(block) != SAMPLES_PER_WINDOW:
            raise SystemExit(f"window {w} has {0 if not block else len(block)} samples")
        p = len(block[0])
        mean = [sum(block[t][j] for t in range(SAMPLES_PER_WINDOW)) / SAMPLES_PER_WINDOW for j in range(p)]
        out.append(mean)
    return out


def concat(a, b):
    return [x + y for x, y in zip(a, b)]


def restrict(mat, idxs):
    return [[row[j] for j in idxs] for row in mat]


def select_lambda(X, y):
    """Post-washout rows: inner 0..87, val 88..109 only. Test unused."""
    inner_end = TRAIN - INNER_VAL
    X_inner = X[:inner_end]
    y_inner = y[:inner_end]
    X_val = X[inner_end:TRAIN]
    y_val = y[inner_end:TRAIN]
    best_lam = LAMBDA_GRID[-1]
    best_score = float("inf")
    val_scores = []
    for lam in LAMBDA_GRID:
        mu, sg = standardize_fit(X_inner)
        w = ridge_fit(apply_std(X_inner, mu, sg), y_inner, lam)
        pred = ridge_predict(apply_std(X_val, mu, sg), w)
        score = nrmse(y_val, pred)
        val_scores.append((lam, score))
        if score < best_score - 1e-15 or (abs(score - best_score) <= 1e-15 and lam > best_lam):
            best_score = score
            best_lam = lam
    return best_lam, val_scores


def fit_eval(X_all, y_all):
    X = X_all[WASHOUT:]
    y = y_all[WASHOUT:]
    if len(X) != TRAIN + 50:
        raise SystemExit(f"post-washout rows {len(X)}")
    lam, val = select_lambda(X, y)
    Xtr, ytr = X[:TRAIN], y[:TRAIN]
    Xte, yte = X[TRAIN:], y[TRAIN:]
    mu, sg = standardize_fit(Xtr)
    w = ridge_fit(apply_std(Xtr, mu, sg), ytr, lam)
    return {
        "lambda": lam,
        "train_nrmse": nrmse(ytr, ridge_predict(apply_std(Xtr, mu, sg), w)),
        "test_nrmse": nrmse(yte, ridge_predict(apply_std(Xte, mu, sg), w)),
        "n_train": TRAIN,
        "n_test": 50,
        "n_features": len(X[0]) if X else 0,
        "val_scores": val,
    }


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


def main() -> None:
    refuse_forbidden()
    occ = OCC_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    aud = AUDIT_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: FAIL**" not in car or "LANE_A_CARRIER" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    if "**Status: SCOPE_NOTE**" not in aud:
        raise SystemExit("audit standing must remain SCOPE_NOTE")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js:
        raise SystemExit("NARMA-10 PROTOCOL not frozen")
    if '"narma10b_rewrite": false' not in js or '"parent_carrier_rewrite": false' not in js:
        raise SystemExit("must freeze no Narma10b/carrier rewrite")
    if U_SHA not in js or '"occupied_mask_is_gate": false' not in js:
        raise SystemExit("PROTOCOL json missing u hash or occupied-mask freeze")

    u = load_u(U_FILE)
    if len(u) != NUM_WINDOWS:
        raise SystemExit(f"u has {len(u)}")
    if any(v < 0.0 or v > 0.5 for v in u):
        raise SystemExit("u escaped [0, 0.5]; Uniform[0,1] forbidden")
    digest = sha256_u(u)
    if digest != U_SHA:
        raise SystemExit(f"u sha256 drifted: {digest}")
    if digest == NARMA10B_SHA:
        raise SystemExit("copied Narma10b u hash")
    rng = random.Random(RNG_SEED)
    regen = [rng.uniform(0.0, 0.5) for _ in range(NUM_WINDOWS)]
    if sha256_u(regen) != digest:
        raise SystemExit("Random(20260828) does not reproduce frozen u")
    y_full = narma10(u)
    y = y_full[1:]

    driven = load_rows(RESULTS / "java_LANE_A_NARMA10_DRIVEN.csv")
    silent = load_rows(RESULTS / "java_LANE_A_NARMA10_SILENT.csv")
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
    last = driven[-1]
    expected_mass = J_MAX * sum(u) * PULSE_S
    commanded = float(last["commanded"])
    source = float(last["source_added"])
    remaining = float(last["remaining"])
    residual = float(last["residual"])
    cmd_rel = abs(commanded - source) / max(abs(commanded), abs(source), expected_mass, 1e-15)
    n0_rel = abs(residual) / max(abs(commanded), abs(source), abs(remaining), 1e-15)

    if not occupancy_alive:
        summary = {
            "gate": "LaneA_NARMA10",
            "status_label": "LANE_A_NARMA10",
            "LANE_A_NARMA10": "NOT_SCORED",
            "system": "NOT_SCORED",
            "occupancy_parent": "PASS",
            "carrier_parent": "FAIL",
            "audit_parent": "SCOPE_NOTE",
            "test_mean_R": test_r,
            "train_mean_R": train_r,
            "narma10b_rewrite": False,
        }
        (RESULTS / "lane_a_narma10_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )
        print("LANE_A_NARMA10=NOT_SCORED occupancy DEAD on this u. Do not raise J_max.")
        raise SystemExit(0)

    if cmd_rel > LEDGER_REL or n0_rel > N0_REL:
        raise SystemExit(f"ledger failed cmd_rel={cmd_rel} n0_rel={n0_rel}")
    if silent_r >= ALIVE_MIN:
        raise SystemExit(f"silent occupied like driven mean_R={silent_r}")

    ahl = window_means(driven, "AHL")
    rmap = window_means(driven, "R")
    lmap = window_means(driven, "L")
    s_r = window_means(silent, "R")
    s_l = window_means(silent, "L")

    scores = {
        "RL": fit_eval(concat(rmap, lmap), y),
        "R": fit_eval(rmap, y),
        "L": fit_eval(lmap, y),
        "FIELD": fit_eval(ahl, y),
        "SILENT_RL": fit_eval(concat(s_r, s_l), y),
        "DELAY_U_10": fit_eval(delay_u10(u), y),
        "FIELD_OCC": fit_eval(restrict(ahl, OCCUPIED), y),
        "R_OCC": fit_eval(restrict(rmap, OCCUPIED), y),
        "RL_OCC": fit_eval(concat(restrict(rmap, OCCUPIED), restrict(lmap, OCCUPIED)), y),
    }
    raw = fit_eval([[v] for v in y], y)
    raw["void"] = True
    scores["RAW_Y"] = raw

    r_n = scores["R"]["test_nrmse"]
    l_n = scores["L"]["test_nrmse"]
    rl_n = scores["RL"]["test_nrmse"]
    field_n = scores["FIELD"]["test_nrmse"]
    silent_n = scores["SILENT_RL"]["test_nrmse"]
    delay_n = scores["DELAY_U_10"]["test_nrmse"]
    best_living = min(r_n, l_n, rl_n)
    system_pass = best_living <= silent_n - DELTA
    living_note = best_living <= field_n - DELTA
    field_wins_living = field_n <= r_n and field_n <= l_n and field_n <= rl_n
    verdict = "PASS" if system_pass else "FAIL"

    summary = {
        "gate": "LaneA_NARMA10",
        "status_label": "LANE_A_NARMA10",
        "LANE_A_NARMA10": verdict,
        "system": verdict,
        "occupancy_parent": "PASS",
        "carrier_parent": "FAIL",
        "carrier_restaged_pass": False,
        "audit_parent": "SCOPE_NOTE",
        "narma10b_rewrite": False,
        "forecasting": True,
        "reconstruction": False,
        "u_sha256": digest,
        "margin_delta": DELTA,
        "test_mean_R": test_r,
        "train_mean_R": train_r,
        "silent_mean_R": silent_r,
        "occupancy": "ALIVE",
        "n0_rel": n0_rel,
        "cmd_rel": cmd_rel,
        "commanded": commanded,
        "expected_mass": expected_mass,
        "best_living": best_living,
        "living_layer_note": living_note,
        "field_wins_or_ties_living": field_wins_living,
        "delay_u_10_ceiling": delay_n,
        "occupied_mask_is_gate": False,
        "test_nrmse": {k: v["test_nrmse"] for k, v in scores.items()},
        "train_nrmse": {k: v["train_nrmse"] for k, v in scores.items()},
        "lambda": {k: v["lambda"] for k, v in scores.items()},
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
        "next": "MG_Lorenz_waveform_only_after_this_standing",
    }
    (RESULTS / "lane_a_narma10_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print("LANE_A_NARMA10 system vs silent. Field reported. Not Narma10b 0.928. Not carrier PASS.")
    print(f"u_sha256={digest} occupancy test_mean_R={test_r:.6g} train_mean_R={train_r:.6g} ALIVE")
    print(f"ledger n0_rel={n0_rel:.3e} cmd_rel={cmd_rel:.3e} silent_R={silent_r:.6g}")
    for name in ("RL", "R", "L", "FIELD", "SILENT_RL", "DELAY_U_10", "FIELD_OCC", "R_OCC", "RL_OCC"):
        s = scores[name]
        print(
            f"  {name} test_NRMSE={s['test_nrmse']:.6g} train_NRMSE={s['train_nrmse']:.6g} "
            f"lambda={s['lambda']:g} LANE_A_NARMA10"
        )
    print(f"  RAW_Y test_NRMSE={raw['test_nrmse']:.6g} VOID answer key LANE_A_NARMA10")
    print(
        f"system={verdict} best_living={best_living:.6g} silent={silent_n:.6g} "
        f"field={field_n:.6g} delay10={delay_n:.6g} living_layer_note={living_note} "
        f"field_wins_living={field_wins_living} LANE_A_NARMA10={verdict}"
    )
    print(
        "Occupancy parent remains PASS. Carrier parent remains FAIL. "
        "Audit remains SCOPE_NOTE. HybridDish Overall unchanged. "
        "Paper 1 unchanged. Lane B not started. Not Fig. 4b. Not C1c."
    )
    if verdict == "PASS":
        print(
            "LANE_A_NARMA10 PASS system: living beats silent by frozen margin. "
            "Not 0.928. Not cells-beat-the-plume unless living_layer_note."
        )
    else:
        print(
            "LANE_A_NARMA10 FAIL system: living did not beat silent by delta. "
            "Occupancy pulse-train PASS unchanged. Do not raise J_max."
        )


if __name__ == "__main__":
    main()
