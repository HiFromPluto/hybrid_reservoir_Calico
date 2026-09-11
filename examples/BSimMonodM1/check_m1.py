#!/usr/bin/env python3
"""Evaluate frozen Monod M1 identity gates from the three production CSVs."""

import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
M1_JAVA = ROOT / "BSimMonodM1.java"
STAGE6_JAVA = REPO / "examples" / "BSimReservoirPlanStage6" / "BSimReservoirPlanStage6.java"
STAGE11_JAVA = REPO / "examples" / "BSimReservoirStage11" / "BSimReservoirStage11.java"

K_S_UM = 0.18
MU_MAX = math.log(2.0) / 1800.0
BATH_TOL = 1e-6
T_MU_START = 600.0
T_MU_END = 7200.0
T_COMPLETE = 7200.0
REQUIRED_COLUMNS = ["t_s", "G_uM", "N", "run_label"]

RUNS = {
    "m1a": {
        "g_um": 0.02,
        "path": ROOT / "results" / "m1a_seed101" / "m1_timeseries.csv",
        "expected_frac": 0.02 / 0.20,
    },
    "m1b": {
        "g_um": 0.18,
        "path": ROOT / "results" / "m1b_seed101" / "m1_timeseries.csv",
        "expected_frac": 0.50,
    },
    "m1c": {
        "g_um": 1.80,
        "path": ROOT / "results" / "m1c_seed101" / "m1_timeseries.csv",
        "expected_frac": 1.80 / 1.98,
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
                "g": float(raw["G_uM"]),
                "n": int(float(raw["N"])),
                "label": raw.get("run_label", ""),
            })
    return header, rows


def nearest_row(rows, t):
    return min(rows, key=lambda row: abs(row["t"] - t))


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


def gate6_static():
    m1 = M1_JAVA.read_text(encoding="utf-8")
    m1_code = strip_java_comments(m1)
    stage6 = STAGE6_JAVA.read_text(encoding="utf-8")
    stage6_code = strip_java_comments(stage6)

    has_ks = "K_S_UM = 0.18" in m1
    has_mu = "MU_MAX = Math.log(2.0) / 1800.0" in m1
    has_sa = "GROWTH_RATE_SA = 4.0 * Math.PI / 1800.0" in m1
    has_602 = "MOL_PER_UM3_PER_UM = 602.0" in m1
    uses_unitless = bool(re.search(r"K_G\s*=\s*5", m1_code)) or bool(
        re.search(r"mu_max\s*=\s*1\b", m1_code, re.I)
    )
    uses_senn_umax = "0.92" in m1_code
    has_uptake = "addQuantity" in m1_code
    has_replenisher = "GlucoseReplenisher" in m1_code or "k_supply" in m1_code
    has_ac = "ArtificialCell" in m1_code
    stage6_glucose = bool(re.search(r"glucose|Monod|K_G|K_S", stage6_code, re.I))
    stage11_exists = STAGE11_JAVA.is_file()

    return {
        "has_ks": has_ks,
        "has_mu": has_mu,
        "has_sa": has_sa,
        "has_602": has_602,
        "uses_unitless": uses_unitless,
        "uses_senn_umax": uses_senn_umax,
        "has_uptake": has_uptake,
        "has_replenisher": has_replenisher,
        "has_ac": has_ac,
        "stage6_glucose": stage6_glucose,
        "stage11_exists": stage11_exists,
        "pass_static": (
            has_ks and has_mu and has_sa and has_602
            and not uses_unitless and not uses_senn_umax
            and not has_uptake and not has_replenisher and not has_ac
            and not stage6_glucose and stage11_exists
        ),
    }


def evaluate():
    static = gate6_static()
    per_run = {}
    missing_files = []
    for label, spec in RUNS.items():
        path = spec["path"]
        if not path.is_file():
            missing_files.append(str(path))
            continue
        header, rows = load_rows(path)
        missing_cols = [name for name in REQUIRED_COLUMNS if name not in header]
        last_t = rows[-1]["t"] if rows else float("nan")
        complete = last_t == T_COMPLETE
        g_err = max((abs(row["g"] - spec["g_um"]) for row in rows), default=float("inf"))
        g_min = min((row["g"] for row in rows), default=float("nan"))
        g_max = max((row["g"] for row in rows), default=float("nan"))
        n0 = rows[0]["n"] if rows else 0
        n1 = rows[-1]["n"] if rows else 0
        labels_ok = all(row["label"] == label for row in rows) if rows else False
        mu, start, end = mu_hat(rows) if rows else (float("nan"), None, None)
        frac = mu / MU_MAX if math.isfinite(mu) else float("nan")
        per_run[label] = {
            "path": str(path),
            "n_rows": len(rows),
            "last_t": last_t,
            "complete": complete,
            "missing_columns": missing_cols,
            "g_err": g_err,
            "g_min": g_min,
            "g_max": g_max,
            "n_initial": n0,
            "n_final": n1,
            "n_start_mu": start["n"] if start else 0,
            "n_end_mu": end["n"] if end else 0,
            "t_start_mu": start["t"] if start else float("nan"),
            "t_end_mu": end["t"] if end else float("nan"),
            "labels_ok": labels_ok,
            "mu_hat": mu,
            "mu_frac": frac,
            "gate1": bool(rows) and g_err <= BATH_TOL,
            "gate2": n1 > n0,
            "g_fell": bool(rows) and g_min < spec["g_um"] - BATH_TOL,
        }

    have_all = all(label in per_run for label in RUNS)
    g1 = have_all and all(per_run[label]["gate1"] for label in RUNS)
    g2 = have_all and all(per_run[label]["gate2"] for label in RUNS)
    mus = [per_run[label]["mu_hat"] for label in ("m1a", "m1b", "m1c")] if have_all else []
    same_n = have_all and (
        per_run["m1a"]["n_final"] == per_run["m1b"]["n_final"]
        == per_run["m1c"]["n_final"]
    )
    g3 = (
        have_all
        and all(math.isfinite(mu) for mu in mus)
        and mus[0] < mus[1] < mus[2]
        and not same_n
    )
    g4 = have_all and 0.85 <= per_run["m1c"]["mu_frac"] <= 1.00
    g5 = have_all and 0.40 <= per_run["m1b"]["mu_frac"] <= 0.60
    g_no_fall = have_all and not any(per_run[label]["g_fell"] for label in RUNS)
    complete = have_all and all(per_run[label]["complete"] for label in RUNS)
    missing_cols = have_all and any(per_run[label]["missing_columns"] for label in RUNS)
    g6 = g_no_fall and static["pass_static"]
    all_pass = (
        complete and not missing_cols and not missing_files
        and g1 and g2 and g3 and g4 and g5 and g6
    )
    return {
        "missing_files": missing_files,
        "per_run": per_run,
        "gate1": g1,
        "gate2": g2,
        "gate3": g3,
        "gate4": g4,
        "gate5": g5,
        "gate6": g6,
        "same_n": same_n,
        "complete": complete,
        "static": static,
        "all_pass": all_pass,
        "mu_max": MU_MAX,
        "k_s": K_S_UM,
    }


def main():
    r = evaluate()
    print(f"mu_max = {r['mu_max']:.12e} /s")
    print(f"K_s = {r['k_s']} uM")
    if r["missing_files"]:
        print("Missing CSVs:")
        for path in r["missing_files"]:
            print(f"  {path}")
    for label, spec in RUNS.items():
        run = r["per_run"].get(label)
        if not run:
            print(f"{label}: NO CSV")
            continue
        print(f"{label} G={spec['g_um']} uM  expected mu/mu_max={spec['expected_frac']:.6f}")
        print(f"  CSV {run['path']}")
        print(f"  rows={run['n_rows']} last_t={run['last_t']} complete={run['complete']}")
        if run["missing_columns"]:
            print(f"  missing columns: {run['missing_columns']}")
        print(f"  G min/max={run['g_min']:.12f}/{run['g_max']:.12f}  max_|err|={run['g_err']:.3e}")
        print(f"  N_initial={run['n_initial']} N_final={run['n_final']}"
              f"  N[{run['t_start_mu']:.0f}]={run['n_start_mu']}"
              f"  N[{run['t_end_mu']:.0f}]={run['n_end_mu']}")
        print(f"  mu_hat={run['mu_hat']:.12e} /s  mu_hat/mu_max={run['mu_frac']:.6f}")

    print(f"Gate 1 G_uM commanded +/- 1e-6 every row: {'PASS' if r['gate1'] else 'FAIL'}")
    print(f"Gate 2 N_final > N_initial every run: {'PASS' if r['gate2'] else 'FAIL'}")
    print(f"Gate 3 mu_hat(0.02) < mu_hat(0.18) < mu_hat(1.80): {'PASS' if r['gate3'] else 'FAIL'}")
    if r["same_n"]:
        print("  N climbs the same at all three G. Monod is not wired: FAIL, stop.")
        print("  Do not add uptake or an AC. Do not retune K_s or mu_max.")
    print(f"Gate 4 mu_hat(1.80)/mu_max in [0.85, 1.00]: {'PASS' if r['gate4'] else 'FAIL'}")
    print(f"Gate 5 mu_hat(0.18)/mu_max in [0.40, 0.60]: {'PASS' if r['gate5'] else 'FAIL'}")
    s = r["static"]
    print(f"Gate 6 no uptake; Stage 11 / Plan Stage 6 untouched; no glucose in dish:"
          f" {'PASS' if r['gate6'] else 'FAIL'}")
    print(f"  uptake={s['has_uptake']} replenisher={s['has_replenisher']} AC={s['has_ac']}"
          f" stage6_glucose={s['stage6_glucose']}")
    print(f"  K_s={s['has_ks']} mu_max=ln2/1800={s['has_mu']} SA={s['has_sa']} 602={s['has_602']}"
          f" unitless={s['uses_unitless']} senn_0.92={s['uses_senn_umax']}")
    print(f"Overall: {'PASS' if r['all_pass'] else 'FAIL'}")
    return 0 if r["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
