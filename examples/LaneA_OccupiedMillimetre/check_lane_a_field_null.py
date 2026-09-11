#!/usr/bin/env python3
"""Lane A field-null XOR checker. u-only before Java. Full ridge after maps. No CHARC."""

from __future__ import annotations

import argparse
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
    ridge_fit,
    ridge_predict,
    standardize_fit,
)
from check_lane_a_narma10 import (  # noqa: E402
    J_MAX,
    LEDGER_REL,
    N0_REL,
    PULSE_S,
    concat,
    sha256_u,
)

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_FIELD_NULL.md"
PROTOCOL_JSON = HERE / "configs" / "field_null_protocol.json"
U_FILE = HERE / "input_u_field_null_400.txt"
LABELS = HERE / "field_null_labels.csv"
U_ONLY_JSON = RESULTS / "lane_a_field_null_u_only.json"
SUMMARY_JSON = RESULTS / "lane_a_field_null_summary.json"
DRIVEN_CSV = RESULTS / "java_LANE_A_FIELD_NULL_DRIVEN.csv"
SILENT_CSV = RESULTS / "java_LANE_A_FIELD_NULL_SILENT.csv"
NARMA_SILENT = RESULTS / "java_LANE_A_NARMA10_SILENT.csv"
WAVEFORM_SILENT = RESULTS / "java_LANE_A_WAVEFORM_SILENT.csv"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"
NARMA_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_NARMA10_STANDING.md"
WAVEFORM_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_WAVEFORM_STANDING.md"
KR_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_KR_GR_STANDING.md"
METRIC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_METRIC_LOCK_STANDING.md"

U_SHA = "3025394aa72732228f4db2a23ba33c84b4a2b2b6b3a215d098d7d3b6f8072d53"
LABEL_SHA = "4b64b64466155d6f12a0a02fd30401cef234587bdee90ad756a90c830e811495"
NUM_WINDOWS = 400
SAMPLES_PER_WINDOW = 16
TRAIN_LO, TRAIN_HI = 64, 271
VAL_LO, VAL_HI = 208, 271
TEST_LO, TEST_HI = 272, 399
OCC_LO, OCC_HI = 272, 399
LAMBDA_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
VOID_LEAK = 0.70
XOR_MIN = 0.90
DELTA = 0.05
ALIVE_MIN = 0.05
SILENT_EPS = 1e-12
HILL_K = 1.6
HILL_N = 2.0


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("charc", "fig4b", "two-way", "twoway", "0.928", "gamma")):
        raise SystemExit("Lane A field-null refuses CHARC / Fig4b / two-way / 0.928 / gamma")


def load_u(path: Path):
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "," in line:
            raise SystemExit("u line contains a comma; one u per line")
        values.append(float(line))
    return values


def load_labels(path: Path):
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    bits = [int(r["bit"]) for r in rows]
    prev = [int(r["bit_prev"]) for r in rows]
    y = [int(r["xor"]) for r in rows]
    u = [float(r["u"]) for r in rows]
    return bits, prev, y, u


def sha256_labels(y) -> str:
    payload = ",".join(str(int(v)) for v in y)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def mann_whitney_auc(y, scores) -> float:
    pairs = sorted(zip(scores, y), key=lambda t: t[0])
    n1 = sum(y)
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        avg = 0.5 * (i + 1 + j)
        for k in range(i, j):
            ranks[k] = avg
        i = j
    rank_pos = sum(r for r, (_, lab) in zip(ranks, pairs) if lab)
    return (rank_pos - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def delay_u(u, n, k):
    feat = []
    for lag in range(k):
        j = n - lag
        feat.append(0.0 if j < 0 else u[j])
    return feat


def moments_of_scalar(u_n):
    return [u_n, 0.0, u_n, u_n]


def static_hill(rows):
    out = []
    kn = HILL_K ** HILL_N
    for row in rows:
        feat = []
        for v in row:
            c = max(float(v), 0.0)
            feat.append((c ** HILL_N) / (kn + c ** HILL_N + 1e-30))
        out.append(feat)
    return out


def select_lambda(X, y, train_idx):
    Xtr = [X[i] for i in train_idx]
    ytr = [y[i] for i in train_idx]
    inner_idx = [i for i, n in enumerate(train_idx) if n < VAL_LO]
    val_loc = [i for i, n in enumerate(train_idx) if VAL_LO <= n <= VAL_HI]
    if not inner_idx or not val_loc:
        raise SystemExit("inner/val split empty")
    Xinner = [Xtr[i] for i in inner_idx]
    yinner = [ytr[i] for i in inner_idx]
    mu, sg = standardize_fit(Xinner)
    zin = apply_std(Xinner, mu, sg)
    zval = apply_std([Xtr[i] for i in val_loc], mu, sg)
    yval = [ytr[i] for i in val_loc]
    best_lam = LAMBDA_GRID[-1]
    best_auc = -1.0
    for lam in LAMBDA_GRID:
        w = ridge_fit(zin, yinner, lam)
        auc = mann_whitney_auc(yval, ridge_predict(zval, w))
        if auc > best_auc + 1e-12 or (abs(auc - best_auc) < 1e-12 and lam > best_lam):
            best_auc = auc
            best_lam = lam
    return best_lam


def fit_test_auc(X, y, train_idx, test_idx, lam):
    xtr = [X[i] for i in train_idx]
    ytr = [float(y[i]) for i in train_idx]
    xte = [X[i] for i in test_idx]
    yte = [y[i] for i in test_idx]
    mu, sg = standardize_fit(xtr)
    w = ridge_fit(apply_std(xtr, mu, sg), ytr, lam)
    scores = ridge_predict(apply_std(xte, mu, sg), w)
    return mann_whitney_auc(yte, scores), lam


def score_arm(X, y, train_idx, test_idx):
    lam = select_lambda(X, y, train_idx)
    auc, lam = fit_test_auc(X, y, train_idx, test_idx, lam)
    return {"test_auc": auc, "lambda": lam}


def u_only(u, y):
    train_idx = list(range(TRAIN_LO, TRAIN_HI + 1))
    test_idx = list(range(TEST_LO, TEST_HI + 1))
    feats = {
        "MEAN_ONLY": [[u[n]] for n in range(NUM_WINDOWS)],
        "MOMENTS": [moments_of_scalar(u[n]) for n in range(NUM_WINDOWS)],
        "DELAY_U_2": [delay_u(u, n, 2) for n in range(NUM_WINDOWS)],
        "DELAY_U_10": [delay_u(u, n, 10) for n in range(NUM_WINDOWS)],
        "XOR_OF_U": [[float(y[n])] for n in range(NUM_WINDOWS)],
    }
    return {name: score_arm(X, y, train_idx, test_idx) for name, X in feats.items()}


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
        mean = [
            sum(block[t][j] for t in range(SAMPLES_PER_WINDOW)) / SAMPLES_PER_WINDOW
            for j in range(p)
        ]
        out.append(mean)
    return out


def load_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if "LANE_A_FIELD_NULL_TASK" not in names or "u" not in names:
            raise SystemExit(f"{path.name} missing LANE_A_FIELD_NULL_TASK/u")
        return list(reader)


def band_mean(rows, key, lo, hi):
    vals = [float(r[key]) for r in rows if lo <= int(r["window"]) <= hi]
    return sum(vals) / len(vals) if vals else float("nan")


def max_abs_prefix(rows, prefix):
    m = 0.0
    for row in rows:
        for i in range(N_BINS):
            a = abs(float(row[f"{prefix}_{i}"]))
            if a > m:
                m = a
    return m


def pearson_r(a, b):
    n = len(a)
    ma = sum(a) / n
    mb = sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = sum((x - ma) ** 2 for x in a) ** 0.5
    db = sum((y - mb) ** 2 for y in b) ** 0.5
    if da < 1e-15 or db < 1e-15:
        return float("nan")
    return num / (da * db)


def require_parents() -> None:
    occ = OCC_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    nar = NARMA_STANDING.read_text(encoding="utf-8")
    wav = WAVEFORM_STANDING.read_text(encoding="utf-8")
    kr = KR_STANDING.read_text(encoding="utf-8")
    met = METRIC_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: FAIL**" not in car or "LANE_A_CARRIER" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    if "**Status: PASS**" not in nar or "LANE_A_NARMA10" not in nar:
        raise SystemExit("NARMA standing must remain system PASS")
    if "**Status: PASS**" not in wav or "LANE_A_WAVEFORM" not in wav:
        raise SystemExit("waveform standing must remain system PASS")
    if "**Status: SCORED**" not in kr or "LANE_A_KR_GR" not in kr:
        raise SystemExit("KR/GR standing must remain SCORED")
    if "NOT CHARC" not in kr:
        raise SystemExit("KR/GR must remain not CHARC")
    if "**Status: PASS**" not in met or "LANE_A_METRIC_LOCK" not in met:
        raise SystemExit("metric-lock standing must remain PASS")


def require_protocol() -> None:
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("PROTOCOL_FIELD_NULL is not frozen")
    if js.get("lane_a_charc_sweep_started") is not False:
        raise SystemExit("CHARC sweep must stay closed")
    if js.get("u_sha256") != U_SHA:
        raise SystemExit("PROTOCOL json u sha drifted")
    if js.get("extra_ACs") is not False or js.get("two_way") is not False:
        raise SystemExit("must not add ACs or two-way")
    if js.get("raise_J_max") is not False or js.get("second_Hill_K") is not False:
        raise SystemExit("must not raise J_max or add a second K")


def main() -> None:
    refuse_forbidden()
    parser = argparse.ArgumentParser()
    parser.add_argument("--u-only", action="store_true")
    args = parser.parse_args()
    require_parents()
    require_protocol()

    u = load_u(U_FILE)
    bits, prev, y, u_lab = load_labels(LABELS)
    if len(u) != NUM_WINDOWS or len(y) != NUM_WINDOWS:
        raise SystemExit("u/labels length")
    if any(abs(a - b) > 1e-12 for a, b in zip(u, u_lab)):
        raise SystemExit("labels u column drifted")
    for n in range(NUM_WINDOWS):
        expect = 0 if n == 0 else bits[n] ^ prev[n]
        if y[n] != expect:
            raise SystemExit(f"xor label mismatch at n={n}")

    u_sha = sha256_u(u)
    y_sha = sha256_labels(y)
    if u_sha != U_SHA or y_sha != LABEL_SHA:
        raise SystemExit("SHA mismatch vs frozen PROTOCOL")

    RESULTS.mkdir(parents=True, exist_ok=True)
    table_u = u_only(u, y)
    leak_names = ("MEAN_ONLY", "MOMENTS", "DELAY_U_2", "DELAY_U_10")
    void = any(table_u[n]["test_auc"] >= VOID_LEAK for n in leak_names)
    key_bad = table_u["XOR_OF_U"]["test_auc"] < XOR_MIN
    u_verdict = "VOID" if void or key_bad else "U_ONLY_PASS"
    U_ONLY_JSON.write_text(
        json.dumps(
            {
                "verdict": u_verdict,
                "u_sha256": u_sha,
                "label_sha256": y_sha,
                "test_windows": [TEST_LO, TEST_HI],
                "arms": table_u,
                "void_if_leak_ge": VOID_LEAK,
                "java": False,
                "charc": False,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print("LANE_A_FIELD_NULL_TASK u-only")
    for name in (*leak_names, "XOR_OF_U"):
        row = table_u[name]
        print(f"  {name} test_AUC={row['test_auc']:.4f} lambda={row['lambda']}")
    if void or key_bad:
        print(f"LANE_A_FIELD_NULL_TASK={u_verdict}")
        raise SystemExit("VOID / STOP: linear u leak or XOR_OF_U is not the answer key; no Java")
    print("LANE_A_FIELD_NULL_TASK=U_ONLY_PASS")
    if args.u_only:
        return

    if not DRIVEN_CSV.exists() or not SILENT_CSV.exists():
        raise SystemExit("field-null maps missing; run ant lane-a-field-null after u-only PASS")
    if SILENT_CSV.resolve() == NARMA_SILENT.resolve() or SILENT_CSV.resolve() == WAVEFORM_SILENT.resolve():
        raise SystemExit("must not reuse NARMA or Waveform silent")

    driven = load_rows(DRIVEN_CSV)
    silent = load_rows(SILENT_CSV)
    if len(driven) != NUM_WINDOWS * SAMPLES_PER_WINDOW or len(silent) != NUM_WINDOWS * SAMPLES_PER_WINDOW:
        raise SystemExit(f"sample counts driven={len(driven)} silent={len(silent)}")
    silent_max = max(max_abs_prefix(silent, "AHL"), max_abs_prefix(silent, "R"), max_abs_prefix(silent, "L"))
    if silent_max > SILENT_EPS:
        raise SystemExit(f"400-window silent maps not empty max_abs={silent_max}")

    csv_u = []
    mean_r_win = []
    for w in range(NUM_WINDOWS):
        block = [r for r in driven if int(r["window"]) == w]
        on = [float(r["u"]) for r in block if float(r["u"]) > 0]
        csv_u.append(on[0] if on else 0.0)
        mean_r_win.append(sum(float(r["mean_R"]) for r in block) / len(block))
    if any(abs(a - b) > 1e-12 for a, b in zip(csv_u, u)):
        raise SystemExit("CSV pulse amplitudes do not match frozen u")

    test_r = band_mean(driven, "mean_R", OCC_LO, OCC_HI)
    train_r = band_mean(driven, "mean_R", TRAIN_LO, TRAIN_HI)
    r_ru = pearson_r(mean_r_win, u)
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
        SUMMARY_JSON.write_text(
            json.dumps(
                {
                    "gate": "LaneA_FIELD_NULL",
                    "LANE_A_FIELD_NULL_TASK": "NOT_SCORED",
                    "test_mean_R": test_r,
                    "r_meanR_u": r_ru,
                    "charc": False,
                    "raise_J_max": False,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print("LANE_A_FIELD_NULL_TASK=NOT_SCORED occupancy DEAD on this u. Do not raise J_max.")
        raise SystemExit(0)
    if cmd_rel > LEDGER_REL or n0_rel > N0_REL:
        raise SystemExit(f"ledger failed cmd_rel={cmd_rel} n0_rel={n0_rel}")

    ahl = window_means(driven, "AHL")
    rmap = window_means(driven, "R")
    lmap = window_means(driven, "L")
    s_r = window_means(silent, "R")
    s_l = window_means(silent, "L")
    train_idx = list(range(TRAIN_LO, TRAIN_HI + 1))
    test_idx = list(range(TEST_LO, TEST_HI + 1))

    scores = {
        "RL": score_arm(concat(rmap, lmap), y, train_idx, test_idx),
        "R": score_arm(rmap, y, train_idx, test_idx),
        "L": score_arm(lmap, y, train_idx, test_idx),
        "FIELD": score_arm(ahl, y, train_idx, test_idx),
        "SILENT_RL": score_arm(concat(s_r, s_l), y, train_idx, test_idx),
        "DELAY_U_2": table_u["DELAY_U_2"],
        "DELAY_U_10": table_u["DELAY_U_10"],
        "MOMENTS": table_u["MOMENTS"],
        "STATIC_HILL": score_arm(static_hill(ahl), y, train_idx, test_idx),
        "XOR_OF_U": table_u["XOR_OF_U"],
    }

    def auc(name):
        return scores[name]["test_auc"]

    rl, field, silent_auc = auc("RL"), auc("FIELD"), auc("SILENT_RL")
    mom, d2, d10 = auc("MOMENTS"), auc("DELAY_U_2"), auc("DELAY_U_10")
    hill = auc("STATIC_HILL")
    field_null = field < 0.5 + DELTA
    living_lift = rl >= 0.5 + DELTA
    vs_silent = rl >= silent_auc + DELTA
    leak_ok = mom < VOID_LEAK and d2 < VOID_LEAK and d10 < VOID_LEAK
    if not field_null:
        verdict = "FAIL field-null"
    elif not living_lift:
        verdict = "FAIL living lift"
    elif occupancy_alive and field_null and living_lift and vs_silent and leak_ok:
        verdict = "PASS"
    else:
        verdict = "FAIL living lift"

    hill_ties_rl = abs(hill - rl) <= 0.01 or hill >= rl - 1e-12
    summary = {
        "gate": "LaneA_FIELD_NULL",
        "status_label": "LANE_A_FIELD_NULL_TASK",
        "LANE_A_FIELD_NULL_TASK": verdict,
        "occupancy_parent": "PASS",
        "carrier_parent": "FAIL",
        "narma_parent": "PASS",
        "waveform_parent": "PASS",
        "kr_gr_parent": "SCORED",
        "metric_lock_parent": "PASS",
        "charc": False,
        "raise_J_max": False,
        "u_sha256": u_sha,
        "label_sha256": y_sha,
        "margin_delta_auc": DELTA,
        "test_mean_R": test_r,
        "train_mean_R": train_r,
        "r_meanR_u": r_ru,
        "occupancy": "ALIVE",
        "n0_rel": n0_rel,
        "cmd_rel": cmd_rel,
        "commanded": commanded,
        "expected_mass": expected_mass,
        "silent_400": True,
        "silent_max_abs": silent_max,
        "static_hill_ties_rl": hill_ties_rl,
        "test_auc": {k: v["test_auc"] for k, v in scores.items()},
        "lambda": {k: v["lambda"] for k, v in scores.items()},
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("LANE_A_FIELD_NULL_TASK successive-bit XOR. Field-null identity. Not CHARC.")
    print(f"occupancy test_mean_R={test_r:.6g} train_mean_R={train_r:.6g} r(mean_R,u)={r_ru:.4g} ALIVE")
    print(f"ledger n0_rel={n0_rel:.3e} cmd_rel={cmd_rel:.3e} silent_400 max_abs={silent_max:.3g}")
    for name in ("RL", "R", "L", "FIELD", "SILENT_RL", "DELAY_U_2", "DELAY_U_10", "MOMENTS", "STATIC_HILL", "XOR_OF_U"):
        s = scores[name]
        print(f"  {name} test_AUC={s['test_auc']:.4f} lambda={s['lambda']:g} LANE_A_FIELD_NULL_TASK")
    if hill_ties_rl:
        print("STATIC_HILL ties RL: the lift is Hill of the snapshot, not L memory.")
    print(
        f"LANE_A_FIELD_NULL_TASK={verdict} RL={rl:.4f} FIELD={field:.4f} SILENT={silent_auc:.4f} "
        f"MOMENTS={mom:.4f} DELAY_U_2={d2:.4f} DELAY_U_10={d10:.4f} STATIC_HILL={hill:.4f}"
    )
    print(
        "Occupancy parent remains PASS. Carrier parent remains FAIL. "
        "NARMA 0.895 unchanged. KR/GR SCORED not CHARC. Not Fig. 4b. Not extra ACs."
    )


if __name__ == "__main__":
    main()
