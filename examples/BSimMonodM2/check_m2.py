#!/usr/bin/env python3
"""Evaluate frozen Monod M2 identity gates from the production CSV."""

import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
M2_JAVA = ROOT / "BSimMonodM2.java"
M1_JAVA = REPO / "examples" / "BSimMonodM1" / "BSimMonodM1.java"
M1_EVIDENCE = REPO / "examples" / "BSimMonodM1" / "results" / "GATE_EVIDENCE.md"
STAGE6_JAVA = REPO / "examples" / "BSimReservoirPlanStage6" / "BSimReservoirPlanStage6.java"
CSV_PATH = ROOT / "results" / "m2a_seed101" / "m2_timeseries.csv"

K_S_UM = 0.18
MU_MAX = math.log(2.0) / 1800.0
K_UPTAKE = 1.2e3
G0_UM = 1.80
G0_TOL = 0.01
T_COMPLETE = 7200.0
T_MU_START = 600.0
T_MU_END = 7200.0
M1C_MU_FRAC = 0.880
REQUIRED_COLUMNS = ["t_s", "G_uM_mean", "G_uM_min", "N", "run_label"]
M1_CONST_NAMES = ["K_S_UM", "MU_MAX", "GROWTH_RATE_SA", "MOL_PER_UM3_PER_UM"]


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
                "g_min": float(raw["G_uM_min"]),
                "n": int(float(raw["N"])),
                "label": raw.get("run_label", ""),
            })
    return header, rows


def nearest_row(rows, t):
    return min(rows, key=lambda row: abs(row["t"] - t))


def window(rows, t0, t1):
    return [row for row in rows if t0 <= row["t"] <= t1]


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def mu_hat(rows):
    start = nearest_row(rows, T_MU_START)
    end = nearest_row(rows, T_MU_END)
    dt = end["t"] - start["t"]
    if start["n"] <= 0 or end["n"] <= 0 or dt <= 0:
        return float("nan"), start, end
    return (math.log(end["n"]) - math.log(start["n"])) / dt, start, end


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


def gate5_and_6_static():
    m1 = parse_constants(M1_JAVA, M1_CONST_NAMES)
    m2 = parse_constants(M2_JAVA, M1_CONST_NAMES + ["K_UPTAKE", "GLU_DECAY_RATE", "GLU_DIFFUSIVITY"])
    m2_text = M2_JAVA.read_text(encoding="utf-8")
    m2_code = strip_java_comments(m2_text)
    stage6_code = strip_java_comments(STAGE6_JAVA.read_text(encoding="utf-8"))
    m1_evidence = M1_EVIDENCE.read_text(encoding="utf-8") if M1_EVIDENCE.is_file() else ""

    mismatches = [
        name for name in M1_CONST_NAMES
        if name not in m1 or name not in m2 or m1[name] != m2[name]
    ]
    k_uptake_ok = m2.get("K_UPTAKE") == "1.2e3"
    decay_zero = m2.get("GLU_DECAY_RATE") == "0.0"
    diff_100 = m2.get("GLU_DIFFUSIVITY") == "100.0"
    has_uptake = "addQuantity" in m2_code
    bath_reset = bool(re.search(r"setConc\(commandedUm|setConc\(G0_UM \* MOL", m2_code))
    # Initialization setConc(G0) is required; a per-tick bath reset is not.
    ticker_reset = "setConc(G0_UM * MOL_PER_UM3_PER_UM)" in m2_code and m2_code.count(
        "setConc(G0_UM * MOL_PER_UM3_PER_UM)"
    ) > 1
    has_replenisher = "GlucoseReplenisher" in m2_code or "k_supply" in m2_code
    has_ac = "ArtificialCell" in m2_code
    stage6_glucose = bool(re.search(r"glucose|Monod|K_G|K_S", stage6_code, re.I))
    m1_pass = "Overall: PASS" in m1_evidence

    return {
        "mismatches": mismatches,
        "k_uptake_ok": k_uptake_ok,
        "decay_zero": decay_zero,
        "diff_100": diff_100,
        "has_uptake": has_uptake,
        "ticker_reset": ticker_reset,
        "has_replenisher": has_replenisher,
        "has_ac": has_ac,
        "stage6_glucose": stage6_glucose,
        "m1_pass": m1_pass,
        "gate5": not mismatches and k_uptake_ok and decay_zero and diff_100
                 and has_uptake and not ticker_reset
                 and not has_replenisher and not has_ac,
        "gate6": m1_pass and not stage6_glucose,
    }


def evaluate(path=CSV_PATH):
    static = gate5_and_6_static()
    if not path.is_file():
        return {
            "path": str(path),
            "missing": True,
            "static": static,
            "all_pass": False,
        }
    header, rows = load_rows(path)
    missing_cols = [name for name in REQUIRED_COLUMNS if name not in header]
    last_t = rows[-1]["t"] if rows else float("nan")
    complete = last_t == T_COMPLETE
    row0 = nearest_row(rows, 0.0) if rows else None
    row_end = nearest_row(rows, T_COMPLETE) if rows else None
    g0 = row0["g_mean"] if row0 else float("nan")
    g_end = row_end["g_mean"] if row_end else float("nan")
    early = window(rows, 0.0, 600.0)
    late = window(rows, 6000.0, 7200.0)
    g_early = mean([row["g_mean"] for row in early])
    g_late = mean([row["g_mean"] for row in late])
    n0 = rows[0]["n"] if rows else 0
    n1 = rows[-1]["n"] if rows else 0
    mu, start, end = mu_hat(rows) if rows else (float("nan"), None, None)
    frac = mu / MU_MAX if math.isfinite(mu) else float("nan")
    labels_ok = all(row["label"] == "m2a" for row in rows) if rows else False

    g1 = row0 is not None and abs(g0 - G0_UM) <= G0_TOL
    g2 = math.isfinite(g0) and math.isfinite(g_end) and g_end < 0.90 * g0
    g3 = math.isfinite(g_early) and math.isfinite(g_late) and g_late < g_early
    g4 = n1 > n0
    all_pass = (
        complete and not missing_cols and labels_ok
        and g1 and g2 and g3 and g4 and static["gate5"] and static["gate6"]
    )
    return {
        "path": str(path),
        "missing": False,
        "n_rows": len(rows),
        "last_t": last_t,
        "complete": complete,
        "missing_columns": missing_cols,
        "labels_ok": labels_ok,
        "g0": g0,
        "g_end": g_end,
        "g_early": g_early,
        "g_late": g_late,
        "n_initial": n0,
        "n_final": n1,
        "n_start_mu": start["n"] if start else 0,
        "n_end_mu": end["n"] if end else 0,
        "mu_hat": mu,
        "mu_frac": frac,
        "gate1": g1,
        "gate2": g2,
        "gate3": g3,
        "gate4": g4,
        "gate5": static["gate5"],
        "gate6": static["gate6"],
        "static": static,
        "all_pass": all_pass,
        "mu_max": MU_MAX,
        "k_s": K_S_UM,
        "k_uptake": K_UPTAKE,
    }


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else CSV_PATH
    r = evaluate(path)
    print(f"mu_max = {r.get('mu_max', MU_MAX):.12e} /s")
    print(f"K_s = {r.get('k_s', K_S_UM)} uM")
    print(f"K_UPTAKE = {r.get('k_uptake', K_UPTAKE):.6g} molecules/s "
          "(engineering identity rate, not Senn 1994)")
    if r.get("missing"):
        print(f"Missing CSV: {r['path']}")
        print("Overall: FAIL")
        return 1
    print(f"CSV {r['path']}")
    print(f"  rows={r['n_rows']} last_t={r['last_t']} complete={r['complete']}")
    if r["missing_columns"]:
        print(f"  missing columns: {r['missing_columns']}")
    print(f"  G(0)={r['g0']:.12f}  G(7200)={r['g_end']:.12f}  ratio={r['g_end']/r['g0'] if r['g0'] else float('nan'):.6f}")
    print(f"  N_initial={r['n_initial']} N_final={r['n_final']}")
    print(f"  N[600]={r['n_start_mu']} N[7200]={r['n_end_mu']}")
    print(f"  mu_hat={r['mu_hat']:.12e} /s  mu_hat/mu_max={r['mu_frac']:.6f}"
          f"  (M1c reference {M1C_MU_FRAC:.3f}, not a gate)")
    print(f"Gate 1 G_uM_mean(0)=1.80 +/- 0.01: {'PASS' if r['gate1'] else 'FAIL'}"
          f"  G0={r['g0']:.6f}")
    print(f"Gate 2 G(7200) < 0.90*G(0): {'PASS' if r['gate2'] else 'FAIL'}")
    if not r["gate2"]:
        print("  G is flat: uptake is not hitting the field. FAIL, stop.")
        print("  Do not raise K_UPTAKE. Do not add a replenisher or AC. Do not start M3.")
    print(f"Gate 3 mean G [6000,7200] < mean G [0,600]: {'PASS' if r['gate3'] else 'FAIL'}")
    print(f"  G_late={r['g_late']:.6e}  G_early={r['g_early']:.6e}")
    print(f"Gate 4 N_final > N_initial: {'PASS' if r['gate4'] else 'FAIL'}")
    s = r["static"]
    print(f"Gate 5 K_s/mu_max/SA from M1; no replenisher/AC: {'PASS' if r['gate5'] else 'FAIL'}")
    print(f"  mismatches={s['mismatches']} K_UPTAKE={s['k_uptake_ok']} decay0={s['decay_zero']}"
          f" D=100={s['diff_100']} uptake={s['has_uptake']} ticker_reset={s['ticker_reset']}")
    print(f"  replenisher={s['has_replenisher']} AC={s['has_ac']}")
    print(f"Gate 6 M1 still PASS; Plan Stage 6 has no glucose: {'PASS' if r['gate6'] else 'FAIL'}")
    print(f"  m1_pass={s['m1_pass']} stage6_glucose={s['stage6_glucose']}")
    print(f"Overall: {'PASS' if r['all_pass'] else 'FAIL'}")
    return 0 if r["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
