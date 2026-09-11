#!/usr/bin/env python3
"""Evaluate frozen Danino D1 identity gates from the production CSV."""

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
CSV_PATH = ROOT / "results" / "d1_seed101" / "d1_timeseries.csv"
STAGE10 = REPO / "examples" / "BSimReservoirStage10" / "BSimReservoirStage10.java"
D1_JAVA = ROOT / "BSimDaninoD1.java"
STAGE3_ARCHIVE = REPO / "examples" / "BSimReservoirPlanStage3" / "ARCHIVED_NEGATIVE_RESULT.md"

BATH_UM = 0.05
BATH_TOL = 1e-6
T_PULSE_ON = 600.0
T_PULSE_OFF = 2400.0
AHL_IN_LO = 0.001
AHL_IN_HI = 10.0
REQUIRED_COLUMNS = [
    "t_s",
    "AHL_ext_uM_mean",
    "LuxI_mean",
    "AHL_in_mean",
    "AiiA_mean",
    "LA_mean",
    "N",
]
QS_NAMES = [
    "TIME_ADJ",
    "QS_DELTA1",
    "QS_DELTA2",
    "QS_G",
    "QS_KP2",
    "QS_KR1OFF",
    "QS_KR1ON",
    "QS_KCAT_AIIA",
    "QS_T_A",
    "QS_T_LA",
    "QS_A0LI",
    "QS_A0AA",
    "QS_KPLI",
    "QS_KPAA",
    "QS_KMLA",
    "QS_KMAA",
    "QS_LTOT",
    "QS_N",
    "CELL_WALL_DIFF",
]


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
                "ahl_ext": float(raw["AHL_ext_uM_mean"]),
                "luxi": float(raw["LuxI_mean"]),
                "ahl_in": float(raw["AHL_in_mean"]),
                "aiia": float(raw["AiiA_mean"]),
                "la": float(raw["LA_mean"]),
                "n": int(float(raw["N"])),
            })
    return header, rows


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def window(rows, t0, t1):
    return [row for row in rows if t0 <= row["t"] <= t1]


def parse_qs_constants(path):
    text = path.read_text(encoding="utf-8")
    found = {}
    for name in QS_NAMES:
        match = re.search(
            rf"static final double {name}\s*=\s*([^;]+);",
            text,
        )
        if match:
            found[name] = re.sub(r"\s+", "", match.group(1))
    return found


def gate6():
    archive = STAGE3_ARCHIVE.read_text(encoding="utf-8")
    stage3_fail = (
        "ARCHIVED FAIL" in archive
        or re.search(r"\bFAIL\b", archive) is not None
    )
    stage10 = parse_qs_constants(STAGE10)
    d1 = parse_qs_constants(D1_JAVA)
    missing = [name for name in QS_NAMES if name not in stage10 or name not in d1]
    mismatches = [
        name for name in QS_NAMES
        if name in stage10 and name in d1 and stage10[name] != d1[name]
    ]
    d1_text = D1_JAVA.read_text(encoding="utf-8")
    uses_602 = "MOL_PER_UM3_PER_UM = 602.0" in d1_text
    code = re.sub(r"/\*.*?\*/", "", d1_text, flags=re.S)
    code = re.sub(r"//.*", "", code)
    code = re.sub(r'"[^"]*"', "", code)
    uses_bad = bool(re.search(r"1e-?15", code))
    return {
        "stage3_fail": stage3_fail,
        "qs_match": not missing and not mismatches,
        "missing": missing,
        "mismatches": mismatches,
        "uses_602": uses_602,
        "uses_1e15": uses_bad,
        "pass": stage3_fail and not missing and not mismatches and uses_602 and not uses_bad,
    }


def evaluate(path=CSV_PATH):
    header, rows = load_rows(path)
    missing_cols = [name for name in REQUIRED_COLUMNS if name not in header]
    last_t = rows[-1]["t"] if rows else float("nan")
    complete = last_t == 5400.0
    pulse = [row for row in rows if T_PULSE_ON <= row["t"] < T_PULSE_OFF]
    late_pulse = window(rows, 1800.0, 2400.0)
    pre = window(rows, 0.0, 600.0)
    wash = window(rows, 4800.0, 5400.0)

    ext_err = max((abs(row["ahl_ext"] - BATH_UM) for row in pulse), default=float("inf"))
    ahl_in_late = mean([row["ahl_in"] for row in late_pulse])
    la_late = mean([row["la"] for row in late_pulse])
    la_pre = mean([row["la"] for row in pre])
    ahl_in_wash = mean([row["ahl_in"] for row in wash])
    n_ok = all(row["n"] == 50 for row in rows) and len(rows) > 0

    g1 = bool(pulse) and ext_err <= BATH_TOL
    g2 = AHL_IN_LO <= ahl_in_late <= AHL_IN_HI
    g3 = la_late > la_pre
    g4 = ahl_in_wash < ahl_in_late
    g5 = n_ok
    g6 = gate6()

    results = {
        "path": str(path),
        "n_rows": len(rows),
        "last_t": last_t,
        "complete": complete,
        "missing_columns": missing_cols,
        "gate1_pulse_ext_err": ext_err,
        "gate1": g1,
        "gate2_ahl_in_late": ahl_in_late,
        "gate2": g2,
        "gate3_la_late": la_late,
        "gate3_la_pre": la_pre,
        "gate3": g3,
        "gate4_ahl_in_wash": ahl_in_wash,
        "gate4_ahl_in_late": ahl_in_late,
        "gate4": g4,
        "gate5": g5,
        "gate6": g6["pass"],
        "gate6_detail": g6,
        "all_pass": complete and not missing_cols and g1 and g2 and g3 and g4 and g5 and g6["pass"],
    }
    return results


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else CSV_PATH
    r = evaluate(path)
    print(f"CSV {r['path']}")
    print(f"  rows={r['n_rows']} last_t={r['last_t']} complete={r['complete']}")
    if r["missing_columns"]:
        print(f"  missing columns: {r['missing_columns']}")
    print(f"Gate 1 bath ON AHL_ext=0.05 +/- 1e-6: {'PASS' if r['gate1'] else 'FAIL'}"
          f"  max_|err|={r['gate1_pulse_ext_err']:.3e}")
    print(f"Gate 2 mean AHL_in [1800,2400] in [0.001,10]: {'PASS' if r['gate2'] else 'FAIL'}"
          f"  AHL_in={r['gate2_ahl_in_late']:.6e} uM")
    print(f"Gate 3 mean LA [1800,2400] > mean LA [0,600]: {'PASS' if r['gate3'] else 'FAIL'}")
    print(f"  LA_pulse={r['gate3_la_late']:.6e}  LA_pre={r['gate3_la_pre']:.6e}")
    print(f"Gate 4 mean AHL_in [4800,5400] < [1800,2400]: {'PASS' if r['gate4'] else 'FAIL'}")
    print(f"  AHL_in_wash={r['gate4_ahl_in_wash']:.6e}  AHL_in_pulse={r['gate4_ahl_in_late']:.6e}")
    print(f"Gate 5 N=50 every row: {'PASS' if r['gate5'] else 'FAIL'}")
    g6 = r["gate6_detail"]
    print(f"Gate 6 Stage3 FAIL + QS_* unchanged: {'PASS' if r['gate6'] else 'FAIL'}")
    print(f"  stage3_archive_fail={g6['stage3_fail']} qs_match={g6['qs_match']}"
          f" uses_602={g6['uses_602']} uses_1e15={g6['uses_1e15']}")
    if g6["missing"] or g6["mismatches"]:
        print(f"  missing={g6['missing']} mismatches={g6['mismatches']}")
    if r["gate2"] is False:
        print("STOP: unit coupling still wrong. Do not retune QS_KPLI, QS_KP2, or bath.")
    elif r["gate2"] and (not r["gate3"]) and (not r["gate4"]):
        print("ODE does not respond to a 0.05 µM bath on this timescale")
        print("STOP: do not start D2.")
    elif r["gate2"] and (not r["gate3"] or not r["gate4"]):
        print("STOP: identity gate 3/4 incomplete. Do not retune QS_* or start D2.")
    print(f"Overall: {'PASS' if r['all_pass'] else 'FAIL'}")
    return 0 if r["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
