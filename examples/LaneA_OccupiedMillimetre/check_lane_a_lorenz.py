#!/usr/bin/env python3
"""Lane A Lorenz L3 k=10 AUTO_X system checker. Not L3 0.826. Does not reopen carrier."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from check_lane_a_narma10 import (  # noqa: E402
    ALIVE_MIN,
    DELTA,
    J_MAX,
    LEDGER_REL,
    N0_REL,
    N_BINS,
    NUM_WINDOWS,
    OCCUPIED,
    PULSE_S,
    band_mean,
    concat,
    delay_u10,
    fit_eval,
    restrict,
    sha256_u,
    window_means,
)

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_LORENZ.md"
PROTOCOL_JSON = HERE / "configs" / "lorenz_protocol.json"
U_FILE = HERE / "input_u_lorenz_k10.txt"
TARGET = HERE / "lorenz_target.csv"
SILENT = RESULTS / "java_LANE_A_NARMA10_SILENT.csv"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"
AUDIT_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_AUDIT_STANDING.md"
NARMA_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_NARMA10_STANDING.md"
MG_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_MG_STANDING.md"

U_SHA = "69b861f256fb39527c037b994885d119ccc8ab2873a3d2bca59f7a8802c44aff"
XMIN = -16.849584157354
XMAX = 15.478808766204
SILENT_EPS = 1e-12


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("charc", "waveform", "fig4b", "cross_y", "k=1", "k=50")):
        raise SystemExit("Lane A Lorenz refuses CHARC/waveform/Fig4b/CROSS_Y/k=1/k=50")


def load_u(path: Path):
    text = path.read_text(encoding="utf-8")
    if text.count("\n") < NUM_WINDOWS:
        raise SystemExit("u file looks like a one-line comma dump; write one u per line")
    values = []
    for line in text.splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            if "," in line:
                raise SystemExit("u line contains a comma; one u per line")
            values.append(float(line))
    return values


def load_target(path: Path):
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    u = [float(r["u"]) for r in rows]
    x = [float(r["x"]) for r in rows]
    x_next = [float(r["x_next"]) for r in rows]
    return u, x, x_next


def load_rows(path: Path, label: str):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if label not in names or "u" not in names:
            raise SystemExit(f"{path.name} missing {label}/u")
        return list(reader)


def persist_features(x):
    return [[v] for v in x]


def ar10(x):
    X = []
    for n in range(len(x)):
        feat = []
        for lag in range(10):
            j = n - lag
            feat.append(0.0 if j < 0 else x[j])
        X.append(feat)
    return X


def max_abs_prefix(rows, prefix):
    m = 0.0
    for row in rows:
        for i in range(N_BINS):
            a = abs(float(row[f"{prefix}_{i}"]))
            if a > m:
                m = a
    return m


def main() -> None:
    refuse_forbidden()
    occ = OCC_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    aud = AUDIT_STANDING.read_text(encoding="utf-8")
    nar = NARMA_STANDING.read_text(encoding="utf-8")
    mg = MG_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: FAIL**" not in car or "LANE_A_CARRIER" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    if "**Status: SCOPE_NOTE**" not in aud:
        raise SystemExit("audit standing must remain SCOPE_NOTE")
    if "**Status: PASS**" not in nar or "LANE_A_NARMA10" not in nar:
        raise SystemExit("NARMA standing must remain system PASS")
    if "**Status: PASS**" not in mg or "LANE_A_MG" not in mg:
        raise SystemExit("MG standing must remain system PASS")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js:
        raise SystemExit("Lorenz PROTOCOL not frozen")
    if '"l3_nrmse_copy": false' not in js or '"parent_narma_rewrite": false' not in js:
        raise SystemExit("must freeze no L3 NRMSE copy / NARMA rewrite")
    if '"parent_mg_rewrite": false' not in js or '"cross_y": false' not in js:
        raise SystemExit("must freeze no MG rewrite / no CROSS_Y")
    if U_SHA not in js or '"occupied_mask_is_gate": false' not in js:
        raise SystemExit("PROTOCOL json missing u hash or occupied-mask freeze")
    if "reuse_LANE_A_NARMA10_SILENT" not in js:
        raise SystemExit("silent reuse of NARMA empty dish must be frozen")
    if '"lorenz_k": 10' not in js or '"charc": false' not in js:
        raise SystemExit("must freeze k=10 and not CHARC")
    if '"ar_m": 10' not in js or '"ar_m_search_on_test": false' not in js:
        raise SystemExit("AR m=10 must be frozen")

    u_file = load_u(U_FILE)
    u_tgt, x, x_next = load_target(TARGET)
    digest = sha256_u(u_file)
    if digest != U_SHA:
        raise SystemExit(f"u sha256 drifted: {digest}")
    if any(abs(a - b) > 1e-12 for a, b in zip(u_file, u_tgt)):
        raise SystemExit("input_u_lorenz_k10.txt does not match lorenz_target.csv")
    xmin, xmax = min(x), max(x)
    if abs(xmin - XMIN) > 1e-12 or abs(xmax - XMAX) > 1e-12:
        raise SystemExit(f"xmin/xmax drifted {xmin} {xmax}")
    if any(v < 0.0 or v > 0.5 for v in u_file):
        raise SystemExit("u escaped [0, 0.5]")

    driven = load_rows(RESULTS / "java_LANE_A_LORENZ_DRIVEN.csv", "LANE_A_LORENZ")
    silent = load_rows(SILENT, "LANE_A_NARMA10")
    if len(driven) != NUM_WINDOWS * 16:
        raise SystemExit(f"driven samples {len(driven)}")
    if len(silent) != NUM_WINDOWS * 16:
        raise SystemExit(f"silent samples {len(silent)}")
    silent_max = max(
        max_abs_prefix(silent, "AHL"),
        max_abs_prefix(silent, "R"),
        max_abs_prefix(silent, "L"),
    )
    if silent_max > SILENT_EPS:
        raise SystemExit(f"reused silent maps not empty max_abs={silent_max}")

    csv_u = []
    for w in range(NUM_WINDOWS):
        block = [r for r in driven if int(r["window"]) == w]
        on = [float(r["u"]) for r in block if float(r["u"]) > 0]
        csv_u.append(on[0] if on else 0.0)
    if any(abs(a - b) > 1e-12 for a, b in zip(csv_u, u_file)):
        raise SystemExit("CSV pulse amplitudes do not match frozen u")

    test_r = band_mean(driven, "mean_R", 150, 199)
    train_r = band_mean(driven, "mean_R", 40, 149)
    occupancy_alive = test_r >= ALIVE_MIN
    last = driven[-1]
    expected_mass = J_MAX * sum(u_file) * PULSE_S
    commanded = float(last["commanded"])
    source = float(last["source_added"])
    remaining = float(last["remaining"])
    residual = float(last["residual"])
    cmd_rel = abs(commanded - source) / max(abs(commanded), abs(source), expected_mass, 1e-15)
    n0_rel = abs(residual) / max(abs(commanded), abs(source), abs(remaining), 1e-15)

    if not occupancy_alive:
        summary = {
            "gate": "LaneA_LORENZ",
            "status_label": "LANE_A_LORENZ",
            "LANE_A_LORENZ": "NOT_SCORED",
            "system": "NOT_SCORED",
            "occupancy_parent": "PASS",
            "carrier_parent": "FAIL",
            "narma_parent": "PASS",
            "mg_parent": "PASS",
            "test_mean_R": test_r,
            "l3_nrmse_copy": False,
        }
        (RESULTS / "lane_a_lorenz_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print("LANE_A_LORENZ=NOT_SCORED occupancy DEAD on this u. Do not raise J_max.")
        raise SystemExit(0)

    if cmd_rel > LEDGER_REL or n0_rel > N0_REL:
        raise SystemExit(f"ledger failed cmd_rel={cmd_rel} n0_rel={n0_rel}")

    ahl = window_means(driven, "AHL")
    rmap = window_means(driven, "R")
    lmap = window_means(driven, "L")
    s_r = window_means(silent, "R")
    s_l = window_means(silent, "L")

    scores = {
        "RL": fit_eval(concat(rmap, lmap), x_next),
        "R": fit_eval(rmap, x_next),
        "L": fit_eval(lmap, x_next),
        "FIELD": fit_eval(ahl, x_next),
        "SILENT_RL": fit_eval(concat(s_r, s_l), x_next),
        "PERSIST": fit_eval(persist_features(x), x_next),
        "AR10": fit_eval(ar10(x), x_next),
        "DELAY_U_10": fit_eval(delay_u10(u_file), x_next),
        "FIELD_OCC": fit_eval(restrict(ahl, OCCUPIED), x_next),
        "RL_OCC": fit_eval(concat(restrict(rmap, OCCUPIED), restrict(lmap, OCCUPIED)), x_next),
    }
    raw = fit_eval([[v] for v in x_next], x_next)
    raw["void"] = True
    scores["RAW_X"] = raw

    r_n = scores["R"]["test_nrmse"]
    l_n = scores["L"]["test_nrmse"]
    rl_n = scores["RL"]["test_nrmse"]
    field_n = scores["FIELD"]["test_nrmse"]
    silent_n = scores["SILENT_RL"]["test_nrmse"]
    persist_n = scores["PERSIST"]["test_nrmse"]
    ar_n = scores["AR10"]["test_nrmse"]
    delay_n = scores["DELAY_U_10"]["test_nrmse"]
    best_living = min(r_n, l_n, rl_n)
    system_pass = best_living <= silent_n - DELTA
    living_note = best_living <= field_n - DELTA and best_living <= ar_n - DELTA
    field_or_ar_wins = field_n <= best_living or ar_n <= best_living
    verdict = "PASS" if system_pass else "FAIL"

    summary = {
        "gate": "LaneA_LORENZ",
        "status_label": "LANE_A_LORENZ",
        "LANE_A_LORENZ": verdict,
        "system": verdict,
        "occupancy_parent": "PASS",
        "carrier_parent": "FAIL",
        "carrier_restaged_pass": False,
        "audit_parent": "SCOPE_NOTE",
        "narma_parent": "PASS",
        "narma_note_generalized": False,
        "mg_parent": "PASS",
        "mg_rewrite": False,
        "l3_u_taken": True,
        "l3_nrmse_copy": False,
        "lorenz_k": 10,
        "lorenz_target": "AUTO_X",
        "cross_y": False,
        "charc": False,
        "forecasting": True,
        "u_sha256": digest,
        "margin_delta": DELTA,
        "test_mean_R": test_r,
        "train_mean_R": train_r,
        "occupancy": "ALIVE",
        "n0_rel": n0_rel,
        "cmd_rel": cmd_rel,
        "commanded": commanded,
        "expected_mass": expected_mass,
        "silent_reused": True,
        "silent_max_abs": silent_max,
        "best_living": best_living,
        "living_layer_note": living_note,
        "field_or_ar_wins_living": field_or_ar_wins,
        "persist_ceiling": persist_n,
        "ar10_ceiling": ar_n,
        "delay_u_10_ceiling": delay_n,
        "occupied_mask_is_gate": False,
        "test_nrmse": {k: v["test_nrmse"] for k, v in scores.items()},
        "train_nrmse": {k: v["train_nrmse"] for k, v in scores.items()},
        "lambda": {k: v["lambda"] for k, v in scores.items()},
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
        "next": "LaneA_WAVEFORM_only_after_this_standing",
    }
    (RESULTS / "lane_a_lorenz_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("LANE_A_LORENZ system vs silent. Field, persist, AR10 reported. TAKEN L3 u. Not 0.826. Not 0.928.")
    print(f"u_sha256={digest} occupancy test_mean_R={test_r:.6g} train_mean_R={train_r:.6g} ALIVE")
    print(f"ledger n0_rel={n0_rel:.3e} cmd_rel={cmd_rel:.3e} silent_reused max_abs={silent_max:.3g}")
    for name in ("RL", "R", "L", "FIELD", "SILENT_RL", "PERSIST", "AR10", "DELAY_U_10", "FIELD_OCC", "RL_OCC"):
        s = scores[name]
        print(
            f"  {name} test_NRMSE={s['test_nrmse']:.6g} train_NRMSE={s['train_nrmse']:.6g} "
            f"lambda={s['lambda']:g} LANE_A_LORENZ"
        )
    print(f"  RAW_X test_NRMSE={raw['test_nrmse']:.6g} VOID answer key LANE_A_LORENZ")
    print(
        f"system={verdict} best_living={best_living:.6g} silent={silent_n:.6g} "
        f"field={field_n:.6g} persist={persist_n:.6g} ar10={ar_n:.6g} delay10={delay_n:.6g} "
        f"living_layer_note={living_note} field_or_ar_wins={field_or_ar_wins} LANE_A_LORENZ={verdict}"
    )
    print(
        "Occupancy parent remains PASS. Carrier parent remains FAIL. "
        "Audit remains SCOPE_NOTE. NARMA remains system PASS. MG remains system PASS. "
        "HybridDish Overall unchanged. L3 scout unchanged. "
        "Paper 1 unchanged. Lane B not started. Not Fig. 4b. Not C1c. Not CHARC."
    )
    if verdict == "PASS":
        print(
            "LANE_A_LORENZ PASS system: living beats silent by frozen margin. "
            "Not L3 0.826. Not cells-beat-the-clock unless living_layer_note."
        )
    else:
        print(
            "LANE_A_LORENZ FAIL system: living did not beat silent by delta. "
            "NARMA/MG system PASS unchanged. Do not raise J_max."
        )


if __name__ == "__main__":
    main()
