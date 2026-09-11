#!/usr/bin/env python3
"""Lane A Waveform2c system checker. AUC, not NRMSE. Not 0.835. U-only before Java."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "BSimReservoirPlanWaveform2c"))

from check_lane_a_narma10 import (  # noqa: E402
    J_MAX,
    LEDGER_REL,
    N0_REL,
    N_BINS,
    OCCUPIED,
    PULSE_S,
    concat,
    restrict,
    sha256_u,
)
import check_waveform2c as W2C  # noqa: E402

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_WAVEFORM.md"
PROTOCOL_JSON = HERE / "configs" / "waveform_protocol.json"
U_FILE = HERE / "input_u_waveform2c_400.txt"
LABELS = HERE / "waveform2c_labels.csv"
U_ONLY_JSON = RESULTS / "lane_a_waveform_u_only.json"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"
AUDIT_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_AUDIT_STANDING.md"
NARMA_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_NARMA10_STANDING.md"
MG_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_MG_STANDING.md"
LORENZ_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_LORENZ_STANDING.md"

U_SHA = "0979090c051c513dce14632e83ae9bf025c4b666e497dc879ae1db4bba30d465"
CLASS_SHA = "465b3ad7132bd1e709f346fa4be2c946d28130ce398b47f41b9c2859baf70add"
DELTA = 0.05
ALIVE_MIN = 0.05
SILENT_EPS = 1e-12
OCC_LO, OCC_HI = 272, 399
NUM_WINDOWS = 400
VOID_LEAK = 0.70
RAW_MIN = 0.90


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("charc", "fig4b", "narma10b", "sine", "square", "triangle")):
        raise SystemExit("Lane A Waveform refuses CHARC/Fig4b/Narma10b/sine-square-triangle")


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


def load_labels(path: Path):
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    y = [int(r["class"]) for r in rows]
    phase = [int(r["phase"]) for r in rows]
    u = [float(r["u"]) for r in rows]
    return y, phase, u


def load_rows(path: Path, label: str):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if label not in names or "u" not in names:
            raise SystemExit(f"{path.name} missing {label}/u")
        return list(reader)


SAMPLES_PER_WINDOW = 16


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


def require_parents() -> None:
    occ = OCC_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    aud = AUDIT_STANDING.read_text(encoding="utf-8")
    nar = NARMA_STANDING.read_text(encoding="utf-8")
    mg = MG_STANDING.read_text(encoding="utf-8")
    lor = LORENZ_STANDING.read_text(encoding="utf-8")
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
    if "**Status: PASS**" not in lor or "LANE_A_LORENZ" not in lor:
        raise SystemExit("Lorenz standing must remain system PASS")


def require_protocol() -> None:
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js:
        raise SystemExit("Waveform PROTOCOL not frozen")
    if "PREHASH_GATES_PASS" not in proto:
        raise SystemExit("pre-hash gates not recorded in PROTOCOL")
    if '"waveform2c_auc_copy": false' not in js or '"sine_square_triangle": false' not in js:
        raise SystemExit("must freeze no Waveform2c AUC copy / no sine-square-triangle")
    if '"reuse_LANE_A_NARMA10_SILENT": false' not in js:
        raise SystemExit("must not reuse 200-window NARMA silent")
    if U_SHA not in js or '"occupied_mask_is_gate": false' not in js:
        raise SystemExit("PROTOCOL json missing u hash or occupied-mask freeze")
    if '"u_only_before_java": true' not in js or '"num_windows": 400' not in js:
        raise SystemExit("must freeze u-only before Java and 400 windows")
    if '"margin_delta_auc": 0.05' not in js:
        raise SystemExit("AUC margin must be 0.05")


def cyclic_shifts(seq):
    seq = tuple(seq)
    return {tuple(seq[p:] + seq[:p]) for p in range(len(seq))}


def prehash_gates(u, y_window, phase) -> None:
    templates = W2C.TEMPLATES
    orbits = [cyclic_shifts(t) for t in templates]
    for i in range(3):
        for j in range(i + 1, 3):
            if orbits[i] & orbits[j]:
                raise SystemExit("orbit intersection nonempty")
        if round(templates[i][0], 12) != W2C.MID:
            raise SystemExit("template does not start at MID")
    if len(set(templates)) != 3:
        raise SystemExit("canonical 8-vectors not pairwise unequal")
    bu = W2C.u_blocks(u)
    ref = {
        "mean": float(bu[0].mean()),
        "var": float(bu[0].var(ddof=0)),
        "mn": float(bu[0].min()),
        "mx": float(bu[0].max()),
    }
    for b in range(W2C.NUM_BLOCKS):
        if abs(float(bu[b].mean()) - ref["mean"]) > 1e-12:
            raise SystemExit("block mean mismatch")
        if abs(float(bu[b].var(ddof=0)) - ref["var"]) > 1e-12:
            raise SystemExit("block variance mismatch")
        if abs(float(bu[b].min()) - ref["mn"]) > 1e-12 or abs(float(bu[b].max()) - ref["mx"]) > 1e-12:
            raise SystemExit("block min/max mismatch")
    if any(p != 0 for p in phase):
        raise SystemExit("phase column is not identically 0")
    _ = y_window


def moments4(block_u):
    mean = block_u.mean(axis=1)
    var = block_u.var(axis=1, ddof=0)
    mn = block_u.min(axis=1)
    mx = block_u.max(axis=1)
    return np.column_stack([mean, var, mn, mx])


def u_only_audit(u, y_window):
    y_block = W2C.block_labels(np.asarray(y_window, dtype=int))
    bu = W2C.u_blocks(u)
    feats = W2C.block_moments_from_u(bu)
    feats["MOMENTS"] = moments4(bu)
    leak_names = ("MEAN_ONLY", "POWER_ONLY", "VARIANCE_ONLY", "MOMENTS", "START_U")
    rows = []
    leak = False
    raw_ok = True
    for name in leak_names + ("RAW_U_8",):
        metrics = W2C.evaluate_block_readout(feats[name], y_block)
        auc = metrics["block_macro_ovr_auc"]
        if name == "RAW_U_8":
            passed = auc >= RAW_MIN
            raw_ok = passed
        else:
            passed = auc < VOID_LEAK
            if not passed:
                leak = True
        rows.append({
            "baseline": name,
            "n_features": metrics["n_features"],
            "lambda": metrics["lambda"],
            "block_macro_ovr_auc": auc,
            "block_accuracy": metrics["block_accuracy"],
            "pass": passed,
        })
    overall = (not leak) and raw_ok
    table = {r["baseline"]: r for r in rows}
    return {
        "overall_pass": overall,
        "void": leak or not raw_ok,
        "leak": leak,
        "raw_ok": raw_ok,
        "rows": rows,
        "MEAN_ONLY": table["MEAN_ONLY"]["block_macro_ovr_auc"],
        "POWER_ONLY": table["POWER_ONLY"]["block_macro_ovr_auc"],
        "VARIANCE_ONLY": table["VARIANCE_ONLY"]["block_macro_ovr_auc"],
        "MOMENTS": table["MOMENTS"]["block_macro_ovr_auc"],
        "START_U": table["START_U"]["block_macro_ovr_auc"],
        "RAW_U_8": table["RAW_U_8"]["block_macro_ovr_auc"],
    }


def write_u_only(audit) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    payload = {
        "gate": "LaneA_WAVEFORM",
        "status_label": "LANE_A_WAVEFORM",
        "u_only_pass": audit["overall_pass"],
        "void": audit["void"],
        "waveform2c_auc_copy": False,
        "u_sha256": U_SHA,
        **{k: audit[k] for k in (
            "MEAN_ONLY", "POWER_ONLY", "VARIANCE_ONLY", "MOMENTS", "START_U", "RAW_U_8"
        )},
        "rows": audit["rows"],
        "next": "java_only_if_u_only_pass",
    }
    U_ONLY_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def pearson_r(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.std() < 1e-15 or b.std() < 1e-15:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def arm_window_auc(X_windows, y_window):
    metrics = W2C.evaluate_window_readout(np.asarray(X_windows, dtype=float), np.asarray(y_window, dtype=int), list(range(NUM_WINDOWS)))
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--u-only", action="store_true")
    args = parser.parse_args()
    refuse_forbidden()
    require_parents()
    require_protocol()

    u = load_u(U_FILE)
    y_window, phase, u_lab = load_labels(LABELS)
    if len(u) != NUM_WINDOWS:
        raise SystemExit(f"u has {len(u)}")
    digest = sha256_u(u)
    if digest != U_SHA:
        raise SystemExit(f"u sha256 drifted: {digest}")
    class_sha = W2C.sha256_ascii_ints([y_window[i] for i in range(0, NUM_WINDOWS, 8)])
    if class_sha != CLASS_SHA:
        raise SystemExit(f"class vector sha drifted: {class_sha}")
    if any(abs(a - b) > 1e-12 for a, b in zip(u, u_lab)):
        raise SystemExit("u file does not match labels")
    if any(v < 0.0 or v > 0.5 for v in u):
        raise SystemExit("u escaped [0, 0.5]")
    prehash_gates(u, y_window, phase)

    audit = u_only_audit(u, y_window)
    write_u_only(audit)
    print("LANE_A_WAVEFORM u-only audit. TAKEN Waveform2c u. Not 0.835. Not sine/square/triangle.")
    print(f"u_sha256={digest}")
    for name in ("MEAN_ONLY", "POWER_ONLY", "VARIANCE_ONLY", "MOMENTS", "START_U", "RAW_U_8"):
        row = next(r for r in audit["rows"] if r["baseline"] == name)
        print(
            f"  {name} block_macro_OVR_AUC={row['block_macro_ovr_auc']:.4f} "
            f"acc={row['block_accuracy']:.4f} lambda={row['lambda']:g} pass={row['pass']}"
        )
    if audit["void"]:
        print("LANE_A_WAVEFORM=VOID u-only leak or RAW_U_8 < 0.90. STOP. No Java.")
        raise SystemExit(0)
    print("u-only PASS. MOMENTS/MEAN/POWER/VARIANCE/START_U < 0.70. RAW_U_8 >= 0.90.")
    if args.u_only:
        print("LANE_A_WAVEFORM u-only done. Java may start.")
        return

    if not U_ONLY_JSON.exists():
        raise SystemExit("u-only json missing; run --u-only first")
    prior = json.loads(U_ONLY_JSON.read_text(encoding="utf-8"))
    if not prior.get("u_only_pass") or prior.get("void"):
        raise SystemExit("u-only did not PASS; Java maps must not be scored")

    driven_path = RESULTS / "java_LANE_A_WAVEFORM_DRIVEN.csv"
    silent_path = RESULTS / "java_LANE_A_WAVEFORM_SILENT.csv"
    if not driven_path.exists() or not silent_path.exists():
        raise SystemExit("waveform maps missing; run ant lane-a-waveform after u-only PASS")
    if silent_path.resolve() == (RESULTS / "java_LANE_A_NARMA10_SILENT.csv").resolve():
        raise SystemExit("must not reuse NARMA 200-window silent")

    driven = load_rows(driven_path, "LANE_A_WAVEFORM")
    silent = load_rows(silent_path, "LANE_A_WAVEFORM")
    if len(driven) != NUM_WINDOWS * 16 or len(silent) != NUM_WINDOWS * 16:
        raise SystemExit(f"sample counts driven={len(driven)} silent={len(silent)}")
    silent_max = max(
        max_abs_prefix(silent, "AHL"),
        max_abs_prefix(silent, "R"),
        max_abs_prefix(silent, "L"),
    )
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
    train_r = band_mean(driven, "mean_R", 64, 271)
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
        summary = {
            "gate": "LaneA_WAVEFORM",
            "LANE_A_WAVEFORM": "NOT_SCORED",
            "system": "NOT_SCORED",
            "test_mean_R": test_r,
            "r_meanR_u": r_ru,
            "waveform2c_auc_copy": False,
        }
        (RESULTS / "lane_a_waveform_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print("LANE_A_WAVEFORM=NOT_SCORED occupancy DEAD on this u. Do not raise J_max.")
        raise SystemExit(0)
    if cmd_rel > LEDGER_REL or n0_rel > N0_REL:
        raise SystemExit(f"ledger failed cmd_rel={cmd_rel} n0_rel={n0_rel}")

    ahl = window_means(driven, "AHL")
    rmap = window_means(driven, "R")
    lmap = window_means(driven, "L")
    s_r = window_means(silent, "R")
    s_l = window_means(silent, "L")
    y = np.asarray(y_window, dtype=int)
    bu = W2C.u_blocks(u)
    y_block = W2C.block_labels(y)

    scores = {
        "RL": arm_window_auc(concat(rmap, lmap), y),
        "R": arm_window_auc(rmap, y),
        "L": arm_window_auc(lmap, y),
        "FIELD": arm_window_auc(ahl, y),
        "SILENT_RL": arm_window_auc(concat(s_r, s_l), y),
        "MOMENTS": W2C.evaluate_block_readout(moments4(bu), y_block),
        "RAW_U_8": W2C.evaluate_block_readout(bu, y_block),
        "FIELD_OCC": arm_window_auc(restrict(ahl, OCCUPIED), y),
        "RL_OCC": arm_window_auc(concat(restrict(rmap, OCCUPIED), restrict(lmap, OCCUPIED)), y),
    }
    scores["RAW_U_8"]["void"] = True

    def auc(name):
        return scores[name]["block_macro_ovr_auc"]

    rl, silent_auc, mom, field = auc("RL"), auc("SILENT_RL"), auc("MOMENTS"), auc("FIELD")
    system_pass = rl >= silent_auc + DELTA and rl >= mom + DELTA
    living_note = rl >= field + DELTA
    field_ge_driven = field >= rl
    verdict = "PASS" if system_pass else "FAIL"

    summary = {
        "gate": "LaneA_WAVEFORM",
        "status_label": "LANE_A_WAVEFORM",
        "LANE_A_WAVEFORM": verdict,
        "system": verdict,
        "occupancy_parent": "PASS",
        "carrier_parent": "FAIL",
        "narma_parent": "PASS",
        "mg_parent": "PASS",
        "lorenz_parent": "PASS",
        "narma_note_generalized": False,
        "lorenz_note_generalized": False,
        "waveform2c_u_taken": True,
        "waveform2c_auc_copy": False,
        "waveform1_rewrite": False,
        "charc": False,
        "u_sha256": digest,
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
        "living_layer_note": living_note,
        "field_ge_driven": field_ge_driven,
        "block_macro_ovr_auc": {k: v["block_macro_ovr_auc"] for k, v in scores.items()},
        "block_accuracy": {k: v.get("block_accuracy", v.get("block_accuracy")) for k, v in scores.items()},
        "lambda": {k: v["lambda"] for k, v in scores.items()},
        "u_only": {k: audit[k] for k in ("MEAN_ONLY", "POWER_ONLY", "VARIANCE_ONLY", "MOMENTS", "START_U", "RAW_U_8")},
        "hybriddish_overall_rewrite": False,
        "paper1_rewrite": False,
        "lane_b_started": False,
        "path_step_2_closed": True,
    }
    (RESULTS / "lane_a_waveform_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("LANE_A_WAVEFORM system vs silent and MOMENTS. Field reported. TAKEN Waveform2c u. Not 0.835.")
    print(f"occupancy test_mean_R={test_r:.6g} train_mean_R={train_r:.6g} r(mean_R,u)={r_ru:.4g} ALIVE")
    print(f"ledger n0_rel={n0_rel:.3e} cmd_rel={cmd_rel:.3e} silent_400 max_abs={silent_max:.3g}")
    for name in ("RL", "R", "L", "FIELD", "SILENT_RL", "MOMENTS", "RAW_U_8", "FIELD_OCC", "RL_OCC"):
        s = scores[name]
        acc = s.get("block_accuracy", float("nan"))
        print(
            f"  {name} block_macro_OVR_AUC={s['block_macro_ovr_auc']:.4f} "
            f"acc={acc:.4f} lambda={s['lambda']:g} LANE_A_WAVEFORM"
        )
    print(
        f"system={verdict} RL={rl:.4f} silent={silent_auc:.4f} moments={mom:.4f} field={field:.4f} "
        f"living_layer_note={living_note} field_ge_driven={field_ge_driven} LANE_A_WAVEFORM={verdict}"
    )
    print(
        "Occupancy parent remains PASS. Carrier parent remains FAIL. "
        "NARMA/MG/Lorenz unchanged. Waveform2c GATE_EVIDENCE unchanged. "
        "Not Fig. 4b. Not C1c. Not CHARC. Path step 2 may close."
    )
    if verdict == "PASS":
        print("LANE_A_WAVEFORM PASS system: RL beats silent and MOMENTS by frozen AUC margin. Not 0.835.")
    else:
        print("LANE_A_WAVEFORM FAIL system: RL did not beat silent and MOMENTS by delta. Do not raise J_max.")


if __name__ == "__main__":
    main()
