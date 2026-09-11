#!/usr/bin/env python3
"""Evaluate frozen Monod M3 identity gates from the on/off production CSVs."""

import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
M3_JAVA = ROOT / "BSimMonodM3.java"
M2_JAVA = REPO / "examples" / "BSimMonodM2" / "BSimMonodM2.java"
M1_EVIDENCE = REPO / "examples" / "BSimMonodM1" / "results" / "GATE_EVIDENCE.md"
M2_EVIDENCE = REPO / "examples" / "BSimMonodM2" / "results" / "GATE_EVIDENCE.md"
STAGE6_JAVA = REPO / "examples" / "BSimReservoirPlanStage6" / "BSimReservoirPlanStage6.java"

K_S_UM = 0.18
MU_MAX = math.log(2.0) / 1800.0
K_UPTAKE = 1.2e3
K_SOURCE_ON = 1.0e5
T_COMPLETE = 7200.0
LATE0, LATE1 = 6000.0, 7200.0
REQUIRED_COLUMNS = [
    "t_s",
    "G_uM_mean",
    "G_uM_near",
    "G_uM_far",
    "N",
    "births_near_cum",
    "births_far_cum",
    "arm",
]
SHARED_CONSTS = ["K_S_UM", "MU_MAX", "GROWTH_RATE_SA", "K_UPTAKE", "MOL_PER_UM3_PER_UM"]

RUNS = {
    "on": {
        "path": ROOT / "results" / "m3on_seed101" / "m3_timeseries.csv",
        "arm": "on",
    },
    "off": {
        "path": ROOT / "results" / "m3off_seed101" / "m3_timeseries.csv",
        "arm": "off",
    },
}


def load_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        header = reader.fieldnames or []
        rows = []
        for raw in reader:
            if not raw or raw.get("t_s") in (None, ""):
                continue
            rows.append({
                "t": float(raw["t_s"]),
                "g_mean": float(raw["G_uM_mean"]),
                "g_near": float(raw["G_uM_near"]),
                "g_far": float(raw["G_uM_far"]),
                "n": int(float(raw["N"])),
                "births_near": int(float(raw["births_near_cum"])),
                "births_far": int(float(raw["births_far_cum"])),
                "arm": raw.get("arm", ""),
            })
    return header, rows


def window(rows, t0, t1):
    return [row for row in rows if t0 <= row["t"] <= t1]


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def strip_java_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//.*", "", text)
    text = re.sub(r'"[^"]*"', "", text)
    return text


def parse_constants(path, names):
    text = path.read_text(encoding="utf-8")
    found = {}
    for name in names:
        match = re.search(rf"static final double {name}\s*=\s*([^;]+);", text)
        if match:
            found[name] = re.sub(r"\s+", "", match.group(1))
    return found


def static_gates():
    m2 = parse_constants(M2_JAVA, SHARED_CONSTS)
    m3 = parse_constants(M3_JAVA, SHARED_CONSTS + ["K_SOURCE_ON", "G0_UM", "GLU_DECAY_RATE"])
    m3_code = strip_java_comments(M3_JAVA.read_text(encoding="utf-8"))
    stage6_code = strip_java_comments(STAGE6_JAVA.read_text(encoding="utf-8"))
    m1_ev = M1_EVIDENCE.read_text(encoding="utf-8") if M1_EVIDENCE.is_file() else ""
    m2_ev = M2_EVIDENCE.read_text(encoding="utf-8") if M2_EVIDENCE.is_file() else ""

    mismatches = [
        name for name in SHARED_CONSTS
        if name not in m2 or name not in m3 or m2[name] != m3[name]
    ]
    k_source_ok = m3.get("K_SOURCE_ON") == "1.0e5"
    g0_ks = m3.get("G0_UM") == "0.18"
    decay_zero = m3.get("GLU_DECAY_RATE") == "0.0"
    has_replenisher = "GlucoseReplenisher" in m3_code or "k_supply" in m3_code
    has_ac = "ArtificialCell" in m3_code
    has_goal = "setGoal" in m3_code
    has_source = "addQuantity(source" in m3_code
    stage6_glucose = bool(re.search(r"glucose|Monod|K_G|K_S", stage6_code, re.I))
    m1_pass = "Overall: PASS" in m1_ev
    m2_pass = "Overall: PASS" in m2_ev

    return {
        "mismatches": mismatches,
        "k_source_ok": k_source_ok,
        "g0_ks": g0_ks,
        "decay_zero": decay_zero,
        "has_replenisher": has_replenisher,
        "has_ac": has_ac,
        "has_goal": has_goal,
        "has_source": has_source,
        "stage6_glucose": stage6_glucose,
        "m1_pass": m1_pass,
        "m2_pass": m2_pass,
        "gate5": (
            not mismatches and k_source_ok and g0_ks and decay_zero
            and has_source and not has_replenisher and not has_ac and not has_goal
        ),
        "gate6": m1_pass and m2_pass and not stage6_glucose,
    }


def summarize(path, expected_arm):
    header, rows = load_rows(path)
    missing_cols = [name for name in REQUIRED_COLUMNS if name not in header]
    last_t = rows[-1]["t"] if rows else float("nan")
    late = window(rows, LATE0, LATE1)
    return {
        "path": str(path),
        "n_rows": len(rows),
        "last_t": last_t,
        "complete": last_t == T_COMPLETE,
        "missing_columns": missing_cols,
        "labels_ok": all(row["arm"] == expected_arm for row in rows) if rows else False,
        "n_initial": rows[0]["n"] if rows else 0,
        "n_final": rows[-1]["n"] if rows else 0,
        "g_near_late": mean([row["g_near"] for row in late]),
        "g_far_late": mean([row["g_far"] for row in late]),
        "g_mean_late": mean([row["g_mean"] for row in late]),
        "births_near": rows[-1]["births_near"] if rows else 0,
        "births_far": rows[-1]["births_far"] if rows else 0,
        "g0_mean": rows[0]["g_mean"] if rows else float("nan"),
    }


def evaluate():
    static = static_gates()
    missing = [spec["path"] for spec in RUNS.values() if not spec["path"].is_file()]
    per = {}
    if not missing:
        per["on"] = summarize(RUNS["on"]["path"], "on")
        per["off"] = summarize(RUNS["off"]["path"], "off")

    have = "on" in per and "off" in per
    on_ratio_ok = False
    off_ratio = float("nan")
    if have:
        off_den = per["off"]["g_far_late"]
        off_ratio = (
            per["off"]["g_near_late"] / off_den if off_den not in (0.0, 0) else float("nan")
        )
        # If both late means are exactly 0, the silent field is uniform empty: ratio ≡ 1.
        if per["off"]["g_near_late"] == 0.0 and per["off"]["g_far_late"] == 0.0:
            off_ratio = 1.0

    g1 = have and per["on"]["g_near_late"] > per["on"]["g_far_late"]
    g2 = have and math.isfinite(off_ratio) and 0.5 <= off_ratio <= 2.0
    g3 = have and per["on"]["g_near_late"] > per["off"]["g_near_late"]
    g4 = have and per["on"]["n_final"] > per["on"]["n_initial"] and (
        per["off"]["n_final"] > per["off"]["n_initial"]
    )
    complete = have and per["on"]["complete"] and per["off"]["complete"]
    missing_cols = have and (per["on"]["missing_columns"] or per["off"]["missing_columns"])
    all_pass = (
        complete and not missing_cols and not missing
        and g1 and g2 and g3 and g4 and static["gate5"] and static["gate6"]
    )
    return {
        "missing_files": [str(p) for p in missing],
        "per": per,
        "off_ratio": off_ratio,
        "gate1": g1,
        "gate2": g2,
        "gate3": g3,
        "gate4": g4,
        "gate5": static["gate5"],
        "gate6": static["gate6"],
        "complete": complete,
        "static": static,
        "all_pass": all_pass,
    }


def main():
    r = evaluate()
    print(f"K_s = {K_S_UM} uM")
    print(f"mu_max = {MU_MAX:.12e} /s")
    print(f"K_UPTAKE = {K_UPTAKE:.6g} molecules/s")
    print(f"K_SOURCE_ON = {K_SOURCE_ON:.6g} molecules/s (engineering identity rate)")
    if r["missing_files"]:
        print("Missing CSVs:")
        for path in r["missing_files"]:
            print(f"  {path}")
    for arm in ("on", "off"):
        run = r["per"].get(arm)
        if not run:
            print(f"M3{arm}: NO CSV")
            continue
        print(f"M3{arm}  CSV {run['path']}")
        print(f"  rows={run['n_rows']} last_t={run['last_t']} complete={run['complete']}")
        if run["missing_columns"]:
            print(f"  missing columns: {run['missing_columns']}")
        print(f"  G(0)_mean={run['g0_mean']:.6e}")
        print(f"  late G_near={run['g_near_late']:.6e}  G_far={run['g_far_late']:.6e}"
              f"  G_mean={run['g_mean_late']:.6e}")
        print(f"  N_initial={run['n_initial']} N_final={run['n_final']}")
        print(f"  births_near={run['births_near']} births_far={run['births_far']} (diagnostic)")

    print(f"Gate 1 M3on late G_near > G_far: {'PASS' if r['gate1'] else 'FAIL'}")
    if "on" in r["per"]:
        print(f"  G_near={r['per']['on']['g_near_late']:.6e}  G_far={r['per']['on']['g_far_late']:.6e}")
    if not r["gate1"]:
        print("  Point source does not beat uptake/diffusion on this box. FAIL, stop.")
        print("  Do not raise K_SOURCE. Do not add chemotaxis or a replenisher.")
    print(f"Gate 2 M3off late G_near/G_far in [0.5, 2.0]: {'PASS' if r['gate2'] else 'FAIL'}")
    print(f"  ratio={r['off_ratio']:.6f}")
    print(f"Gate 3 M3on late G_near > M3off late G_near: {'PASS' if r['gate3'] else 'FAIL'}")
    if "on" in r["per"] and "off" in r["per"]:
        print(f"  on={r['per']['on']['g_near_late']:.6e}  off={r['per']['off']['g_near_late']:.6e}")
    print(f"Gate 4 N_final > N_initial both arms: {'PASS' if r['gate4'] else 'FAIL'}")
    s = r["static"]
    print(f"Gate 5 K_UPTAKE/K_s/mu_max from M2; no replenisher/AC/setGoal:"
          f" {'PASS' if r['gate5'] else 'FAIL'}")
    print(f"  mismatches={s['mismatches']} K_SOURCE={s['k_source_ok']} G0=Ks={s['g0_ks']}"
          f" decay0={s['decay_zero']} source={s['has_source']}")
    print(f"  replenisher={s['has_replenisher']} AC={s['has_ac']} setGoal={s['has_goal']}")
    print(f"Gate 6 M1 PASS, M2 PASS, Plan Stage 6 has no glucose:"
          f" {'PASS' if r['gate6'] else 'FAIL'}")
    print(f"  m1_pass={s['m1_pass']} m2_pass={s['m2_pass']} stage6_glucose={s['stage6_glucose']}")
    print(f"Overall: {'PASS' if r['all_pass'] else 'FAIL'}")
    return 0 if r["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
