#!/usr/bin/env python3
"""Evaluate job-2 isolated-cell gates from the production CSV."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "results" / "job2_seed101" / "job2_timeseries.csv"

LAMBDA_S_PER_HOUR = 1.0
NU_PER_SECOND = LAMBDA_S_PER_HOUR / 3600.0
KS_MM = 0.02
NUTRIENT_MM = 0.5
L0_UM = 1.0
LDIV_UM = 3.0
RADIUS_UM = 0.5
SIM_TIME_S = 9000.0
LOG_DT_S = 10.0
REL_LEN_TOL = 1e-6
LEN_AFTER_DIV_TOL_UM = 0.02
RADIUS_TOL_UM = 1e-9

REQUIRED_COLUMNS = [
    "t_s",
    "N",
    "founder_L_um",
    "founder_L_analytic_um",
    "L_min_um",
    "L_max_um",
    "radius_um",
    "nutrient_mM",
    "n_negative",
]


def monod(n: float) -> float:
    return n / (n + KS_MM)


def analytic_tdiv() -> float:
    return math.log(2.0) / (NU_PER_SECOND * monod(NUTRIENT_MM))


def load_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        header = reader.fieldnames or []
        rows = []
        for raw in reader:
            if not raw or raw.get("t_s") in (None, ""):
                continue
            rows.append({
                "t": float(raw["t_s"]),
                "n": int(float(raw["N"])),
                "founder_L": float(raw["founder_L_um"]),
                "analytic": float(raw["founder_L_analytic_um"]),
                "lmin": float(raw["L_min_um"]),
                "lmax": float(raw["L_max_um"]),
                "radius": float(raw["radius_um"]),
                "nutrient": float(raw["nutrient_mM"]),
                "n_negative": int(float(raw["n_negative"])),
            })
    return header, rows


def evaluate(path: Path = CSV_PATH):
    header, rows = load_rows(path)
    missing = [name for name in REQUIRED_COLUMNS if name not in header]
    last_t = rows[-1]["t"] if rows else float("nan")
    complete = abs(last_t - SIM_TIME_S) < 1e-6
    tdiv = analytic_tdiv()

    g1 = bool(rows) and all(row["n_negative"] == 0 for row in rows)
    g1 = g1 and all(
        row["founder_L"] >= 0 and row["lmin"] >= 0 and row["nutrient"] >= 0
        for row in rows
    )

    pre = [row for row in rows if row["n"] == 1]
    rel_err = 0.0
    for row in pre:
        rel_err = max(rel_err, abs(row["founder_L"] - row["analytic"]) / L0_UM)
    g2 = bool(pre) and rel_err <= REL_LEN_TOL

    first_n2 = next((row for row in rows if row["n"] == 2), None)
    t_n2 = first_n2["t"] if first_n2 else float("inf")
    g3 = first_n2 is not None and abs(t_n2 - tdiv) <= LOG_DT_S

    g4 = False
    if first_n2 is not None:
        g4 = (
            abs(first_n2["lmin"] - L0_UM) <= LEN_AFTER_DIV_TOL_UM
            and abs(first_n2["lmax"] - L0_UM) <= LEN_AFTER_DIV_TOL_UM
        )

    g5 = bool(rows) and all(abs(row["radius"] - RADIUS_UM) <= RADIUS_TOL_UM for row in rows)

    allowed_n = {1, 2, 4, 8}
    g6 = complete and not missing and all(row["n"] in allowed_n for row in rows)

    results = {
        "path": str(path),
        "n_rows": len(rows),
        "last_t": last_t,
        "complete": complete,
        "missing_columns": missing,
        "tdiv_s": tdiv,
        "tdiv_min": tdiv / 60.0,
        "t_n2": t_n2,
        "rel_err": rel_err,
        "gate1": g1,
        "gate2": g2,
        "gate3": g3,
        "gate4": g4,
        "gate5": g5,
        "gate6": g6,
        "lmin_at_div": first_n2["lmin"] if first_n2 else float("nan"),
        "lmax_at_div": first_n2["lmax"] if first_n2 else float("nan"),
        "n_final": rows[-1]["n"] if rows else 0,
        "all_pass": g1 and g2 and g3 and g4 and g5 and g6,
    }
    return results


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else CSV_PATH
    r = evaluate(path)
    print(f"CSV {r['path']}")
    print(f"  rows={r['n_rows']} last_t={r['last_t']} complete={r['complete']}")
    if r["missing_columns"]:
        print(f"  missing columns: {r['missing_columns']}")
    print(f"  analytic T_div={r['tdiv_s']:.3f} s ({r['tdiv_min']:.2f} min)")
    print(f"Gate 1 non-negative L, radius, nutrient: {'PASS' if r['gate1'] else 'FAIL'}")
    print(f"Gate 2 founder vs Valdez analytic before division: {'PASS' if r['gate2'] else 'FAIL'}"
          f"  max_|L-Lhat|/L0={r['rel_err']:.3e}")
    print(f"Gate 3 first N=2 within {LOG_DT_S:.0f} s of T_div: {'PASS' if r['gate3'] else 'FAIL'}"
          f"  t_N2={r['t_n2']}")
    print(f"Gate 4 daughters at L0 after first division: {'PASS' if r['gate4'] else 'FAIL'}"
          f"  Lmin={r['lmin_at_div']} Lmax={r['lmax_at_div']}")
    print(f"Gate 5 radius=0.5 um: {'PASS' if r['gate5'] else 'FAIL'}")
    print(f"Gate 6 complete binary fission N in {{1,2,4,8}}: {'PASS' if r['gate6'] else 'FAIL'}"
          f"  N_final={r['n_final']}")
    print(f"Overall: {'PASS' if r['all_pass'] else 'FAIL'}")
    if not r["all_pass"]:
        print("STOP: do not retune lambda_S or K_S.")
    return 0 if r["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
