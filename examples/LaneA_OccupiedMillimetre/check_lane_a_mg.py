#!/usr/bin/env python3
"""Lane A Mackey-Glass system checker. Not BenchA 0.260. Does not reopen carrier."""

from __future__ import annotations

import csv
import hashlib
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
PROTOCOL = HERE / "PROTOCOL_MG.md"
PROTOCOL_JSON = HERE / "configs" / "mg_protocol.json"
U_FILE = HERE / "input_u_mg200.txt"
TARGET = HERE / "mg_target.csv"
SILENT = RESULTS / "java_LANE_A_NARMA10_SILENT.csv"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"
AUDIT_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_AUDIT_STANDING.md"
NARMA_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_NARMA10_STANDING.md"

U_SHA = "e06810ea70f61bbe4baf61939f93b1e1cae9d92d9eee8d757865c4935dbac780"
XMIN = 0.411423008739756
XMAX = 1.2952960679588346
SILENT_EPS = 1e-12


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("charc", "lorenz", "fig4b", "waveform")):
        raise SystemExit("Lane A MG refuses CHARC/Lorenz/Fig4b/waveform")


def load_u(path: Path):
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
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


def mackey_glass_recorded():
    beta, gamma, n_pow, tau = 0.2, 0.1, 10, 17
    x_neg, transient, n_rec = 1.2, 1000, 201
    n_total = transient + n_rec
    x = [0.0] * n_total
    x[0] = x_neg
    for t in range(n_total - 1):
        xt = x[t - tau] if t - tau >= 0 else x_neg
        x[t + 1] = x[t] + beta * xt / (1.0 + xt ** n_pow) - gamma * x[t]
    recorded = x[transient : transient + n_rec]
    x_drive = recorded[:NUM_WINDOWS]
    x_next = recorded[1 : NUM_WINDOWS + 1]
    xmin, xmax = min(x_drive), max(x_drive)
    u = [0.5 * (v - xmin) / (xmax - xmin) for v in x_drive]
    return u, x_drive, x_next, xmin, xmax


def persist_features(x):
    return [[v] for v in x]


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
    if "**Status: PASS**" not in occ or "LANE_A_OCCUPIED_MILLIMETRE" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: FAIL**" not in car or "LANE_A_CARRIER" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    if "**Status: SCOPE_NOTE**" not in aud:
        raise SystemExit("audit standing must remain SCOPE_NOTE")
    if "**Status: PASS**" not in nar or "LANE_A_NARMA10" not in nar:
        raise SystemExit("NARMA standing must remain system PASS")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js:
        raise SystemExit("MG PROTOCOL not frozen")
    if '"bencha_nrmse_copy": false' not in js or '"parent_narma_rewrite": false' not in js:
        raise SystemExit("must freeze no BenchA NRMSE copy / NARMA rewrite")
    if U_SHA not in js or '"occupied_mask_is_gate": false' not in js:
        raise SystemExit("PROTOCOL json missing u hash or occupied-mask freeze")
    if "reuse_LANE_A_NARMA10_SILENT" not in js:
        raise SystemExit("silent reuse of NARMA empty dish must be frozen")

    u_file = load_u(U_FILE)
    u_tgt, x, x_next = load_target(TARGET)
    u_gen, x_gen, x_next_gen, xmin, xmax = mackey_glass_recorded()
    if abs(xmin - XMIN) > 1e-12 or abs(xmax - XMAX) > 1e-12:
        raise SystemExit(f"xmin/xmax drifted {xmin} {xmax}")
    digest = sha256_u(u_file)
    if digest != U_SHA:
        raise SystemExit(f"u sha256 drifted: {digest}")
    if any(abs(a - b) > 1e-12 for a, b in zip(u_file, u_tgt)):
        raise SystemExit("input_u_mg200.txt does not match mg_target.csv")
    if any(abs(a - b) > 1e-12 for a, b in zip(u_file, u_gen)):
        raise SystemExit("frozen u does not match regenerated MG map")
    if any(abs(a - b) > 1e-12 for a, b in zip(x_next, x_next_gen)):
        raise SystemExit("x_next does not match regenerated MG map")
    if any(v < 0.0 or v > 0.5 for v in u_file):
        raise SystemExit("u escaped [0, 0.5]")

    driven = load_rows(RESULTS / "java_LANE_A_MG_DRIVEN.csv", "LANE_A_MG")
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
            "gate": "LaneA_MG",
            "status_label": "LANE_A_MG",
            "LANE_A_MG": "NOT_SCORED",
            "system": "NOT_SCORED",
            "occupancy_parent": "PASS",
            "carrier_parent": "FAIL",
            "narma_parent": "PASS",
            "test_mean_R": test_r,
            "bencha_nrmse_copy": False,
        }
        (RESULTS / "lane_a_mg_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print("LANE_A_MG=NOT_SCORED occupancy DEAD on this u. Do not raise J_max.")
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
    delay_n = scores["DELAY_U_10"]["test_nrmse"]
    best_living = min(r_n, l_n, rl_n)
    system_pass = best_living <= silent_n - DELTA
    living_note = best_living <= field_n - DELTA
    field_wins_living = field_n <= r_n and field_n <= l_n and field_n <= rl_n
    verdict = "PASS" if system_pass else "FAIL"

    summary = {
        "gate": "LaneA_MG",
        "status_label": "LANE_A_MG",
        "LANE_A_MG": verdict,
        "system": verdict,
        "occupancy_parent": "PASS",
        "carrier_parent": "FAIL",
        "carrier_restaged_pass": False,
        "audit_parent": "SCOPE_NOTE",
        "narma_parent": "PASS",
        "narma_note_generalized": False,
        "bencha_u_taken": True,
        "bencha_nrmse_copy": False,
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
        "field_wins_or_ties_living": field_wins_living,
        "persist_ceiling": persist_n,
        "delay_u_10_ceiling": delay_n,
        "occupied_mask_is_gate": False,
        "test_nrmse": {k: v["test_nrmse"] for k, v in scores.items()},
        "train_nrmse": {k: v["train_nrmse"] for k, v in scores.items()},
        "lambda": {k: v["lambda"] for k, v in scores.items()},
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
        "next": "LaneA_Lorenz_only_after_this_standing",
    }
    (RESULTS / "lane_a_mg_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("LANE_A_MG system vs silent. Field and persist reported. TAKEN BenchA u. Not 0.260. Not 0.928.")
    print(f"u_sha256={digest} occupancy test_mean_R={test_r:.6g} train_mean_R={train_r:.6g} ALIVE")
    print(f"ledger n0_rel={n0_rel:.3e} cmd_rel={cmd_rel:.3e} silent_reused max_abs={silent_max:.3g}")
    for name in ("RL", "R", "L", "FIELD", "SILENT_RL", "PERSIST", "DELAY_U_10", "FIELD_OCC", "RL_OCC"):
        s = scores[name]
        print(
            f"  {name} test_NRMSE={s['test_nrmse']:.6g} train_NRMSE={s['train_nrmse']:.6g} "
            f"lambda={s['lambda']:g} LANE_A_MG"
        )
    print(f"  RAW_X test_NRMSE={raw['test_nrmse']:.6g} VOID answer key LANE_A_MG")
    print(
        f"system={verdict} best_living={best_living:.6g} silent={silent_n:.6g} "
        f"field={field_n:.6g} persist={persist_n:.6g} delay10={delay_n:.6g} "
        f"living_layer_note={living_note} field_wins_living={field_wins_living} LANE_A_MG={verdict}"
    )
    print(
        "Occupancy parent remains PASS. Carrier parent remains FAIL. "
        "Audit remains SCOPE_NOTE. NARMA remains system PASS. "
        "HybridDish Overall unchanged. BenchA GATE_EVIDENCE unchanged. "
        "Paper 1 unchanged. Lane B not started. Not Fig. 4b. Not C1c."
    )
    if verdict == "PASS":
        print(
            "LANE_A_MG PASS system: living beats silent by frozen margin. "
            "Not BenchA 0.260. Not cells-beat-the-plume unless living_layer_note."
        )
    else:
        print(
            "LANE_A_MG FAIL system: living did not beat silent by delta. "
            "NARMA system PASS unchanged. Do not raise J_max."
        )


if __name__ == "__main__":
    main()
