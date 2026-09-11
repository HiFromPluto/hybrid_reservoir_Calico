#!/usr/bin/env python3
"""Evaluate frozen Danino D1g identity gates from the production CSV."""

import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
CSV_PATH = ROOT / "results" / "d1g_seed101" / "d1g_timeseries.csv"
PARAMS_PATH = ROOT / "results" / "d1g_seed101" / "d1g_params.csv"
D1G_JAVA = ROOT / "BSimDaninoD1g.java"
D1_JAVA = REPO / "examples" / "BSimDaninoD1" / "BSimDaninoD1.java"
D1_EVIDENCE = REPO / "examples" / "BSimDaninoD1" / "results" / "GATE_EVIDENCE.md"
STAGE10 = REPO / "examples" / "BSimReservoirStage10" / "BSimReservoirStage10.java"
STAGE3_ARCHIVE = REPO / "examples" / "BSimReservoirPlanStage3" / "ARCHIVED_NEGATIVE_RESULT.md"

BATH_UM = 0.05
BATH_TOL = 1e-6
T_PULSE_ON = 600.0
T_PULSE_OFF = 2400.0
AHL_IN_LO = 0.001
AHL_IN_HI = 10.0
MU_EXPECTED = math.log(2.0) / 1800.0
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


def load_params(path):
    values = {}
    if not path.exists():
        return values
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def window(rows, t0, t1):
    return [row for row in rows if t0 <= row["t"] <= t1]


def strip_code(text):
    code = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    code = re.sub(r"//.*", "", code)
    return re.sub(r'"[^"]*"', "", code)


def parse_qs_constants(path):
    text = path.read_text(encoding="utf-8")
    found = {}
    for name in QS_NAMES:
        match = re.search(rf"static final double {name}\s*=\s*([^;]+);", text)
        if match:
            found[name] = re.sub(r"\s+", "", match.group(1))
    return found


def ics_are_zero(text):
    if re.search(r"new double\[\]\s*\{\s*0\.05", text):
        return False
    return bool(re.search(
        r"(ICS\s*=\s*\{0(?:\.0)?,\s*0(?:\.0)?,\s*0(?:\.0)?,\s*0(?:\.0)?\}"
        r"|new double\[\]\s*\{\s*0(?:\.0)?,\s*0(?:\.0)?,\s*0(?:\.0)?,\s*0(?:\.0)?\s*\})",
        text,
    ))


def gate6():
    archive = STAGE3_ARCHIVE.read_text(encoding="utf-8")
    d1_evidence = D1_EVIDENCE.read_text(encoding="utf-8") if D1_EVIDENCE.exists() else ""
    d1_java = D1_JAVA.read_text(encoding="utf-8") if D1_JAVA.exists() else ""
    d1g_text = D1G_JAVA.read_text(encoding="utf-8")
    d1g_code = strip_code(d1g_text)
    params = load_params(PARAMS_PATH)

    stage3_fail = "ARCHIVED FAIL" in archive
    d1_fail = "Overall: FAIL" in d1_evidence
    d1_ics_untouched = "new double[]{0.05, 0.05, 0.05, 0.05}" in d1_java

    stage10 = parse_qs_constants(STAGE10)
    d1 = parse_qs_constants(D1_JAVA)
    d1g = parse_qs_constants(D1G_JAVA)
    missing = [name for name in QS_NAMES if name not in d1g or name not in d1 or name not in stage10]
    mismatches = [
        name for name in QS_NAMES
        if name in stage10 and name in d1g and stage10[name] != d1g[name]
    ]
    mismatches += [
        name for name in QS_NAMES
        if name in d1 and name in d1g and d1[name] != d1g[name] and name not in mismatches
    ]

    uses_602 = "MOL_PER_UM3_PER_UM = 602.0" in d1g_text
    uses_1e15 = bool(re.search(r"1e-?15", d1g_code))
    mu_present = "Math.log(2.0) / 1800.0" in d1g_text.replace(" ", "") or \
        "Math.log(2.0)/1800.0" in d1g_text.replace(" ", "")
    mu_on_proteins = all(f"- MU * y[{i}]" in d1g_code or f"-MU*y[{i}]" in d1g_code.replace(" ", "")
                         for i in (0, 2, 3))
    mu_on_ahl = bool(re.search(r"- *MU *\* *y\[1\]", d1g_code))
    zero_ics = ics_are_zero(d1g_text)
    sidecar_mu = abs(float(params.get("mu_per_s", "nan")) - MU_EXPECTED) < 1e-12 if "mu_per_s" in params else False
    sidecar_ics = params.get("ics", "") in ("0,0,0,0", "0.0,0.0,0.0,0.0")

    return {
        "stage3_fail": stage3_fail,
        "d1_fail": d1_fail,
        "d1_ics_untouched": d1_ics_untouched,
        "qs_match": not missing and not mismatches,
        "missing": missing,
        "mismatches": mismatches,
        "uses_602": uses_602,
        "uses_1e15": uses_1e15,
        "mu_present": mu_present,
        "mu_on_proteins": mu_on_proteins,
        "mu_on_ahl": mu_on_ahl,
        "zero_ics": zero_ics,
        "sidecar_mu": sidecar_mu,
        "sidecar_ics": sidecar_ics,
        "pass": (
            stage3_fail and d1_fail and d1_ics_untouched and not missing
            and not mismatches and uses_602 and not uses_1e15 and mu_present
            and mu_on_proteins and not mu_on_ahl and zero_ics
            and sidecar_mu and sidecar_ics
        ),
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

    return {
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
    print(f"Gate 6 D1 FAIL + Stage3 FAIL + zero ICs + mu + 602: {'PASS' if r['gate6'] else 'FAIL'}")
    print(f"  d1_fail={g6['d1_fail']} stage3_fail={g6['stage3_fail']}"
          f" qs_match={g6['qs_match']} uses_602={g6['uses_602']}"
          f" uses_1e15={g6['uses_1e15']}")
    print(f"  zero_ics={g6['zero_ics']} mu_present={g6['mu_present']}"
          f" mu_on_proteins={g6['mu_on_proteins']} mu_on_ahl={g6['mu_on_ahl']}")
    print(f"  sidecar_mu={g6['sidecar_mu']} sidecar_ics={g6['sidecar_ics']}"
          f" d1_ics_untouched={g6['d1_ics_untouched']}")
    if g6["missing"] or g6["mismatches"]:
        print(f"  missing={g6['missing']} mismatches={g6['mismatches']}")
    if not r["gate3"]:
        print("D1g FAIL on gate 3. Do not retune QS_KPLI, QS_KP2, QS_KR1ON, mu, or bath.")
        print("STOP: do not start D2.")
    print(f"Overall: {'PASS' if r['all_pass'] else 'FAIL'}")
    if r["all_pass"]:
        print("PASS. Stop. Do not merge into BSimReservoirPlanStage6.")
    return 0 if r["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
