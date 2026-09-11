#!/usr/bin/env python3
"""Evaluate frozen Job 6 standalone starvation-viability gates."""

from __future__ import annotations

import csv
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "job6_starvation"
VIABILITY = RESULTS / "viability_timeseries.csv"
SUMMARY = RESULTS / "job6_summary.csv"


def load(path: Path):
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=";"))


def gate1_units(rows):
    errors = []
    unit_rows = [r for r in rows if r["section"] == "units"]
    if len(unit_rows) != 4:
        return [f"expected 4 EV1 rows, got {len(unit_rows)}"]
    for row in unit_rows:
        mu = float(row["arg1"])
        sd = float(row["arg2"])
        predicted = float(row["value"])
        observed = float(row["reference"])
        if abs(predicted - observed) > sd + 1e-12:
            errors.append(
                f"mu={mu:.2f}: {predicted:.4f} outside {observed:.2f}+/-{sd:.2f}"
            )
        wrong_hourly = 24.0 * predicted
        if wrong_hourly / observed <= 10.0:
            errors.append(f"unit typo guard did not exceed 10x at mu={mu:.2f}")
    return errors


def gate2_viability(rows):
    errors = []
    if not rows:
        return ["empty viability CSV"]
    previous = float(rows[0]["N_rk4"])
    max_error = 0.0
    for row in rows:
        value = float(row["N_rk4"])
        exact = float(row["N_exact"])
        if not math.isfinite(value) or value < 0:
            errors.append(f"bad viability {value} at t={row['t_day']}")
            break
        if value > previous + 1e-14:
            errors.append(f"viability increased at t={row['t_day']}")
            break
        previous = value
        max_error = max(max_error, abs(value - exact))
    if max_error >= 1e-9:
        errors.append(f"max RK4 error {max_error:.3e} >= 1e-9")
    if abs(float(rows[-1]["t_day"]) - 10.0) > 1e-12:
        errors.append(f"last t={rows[-1]['t_day']} instead of 10 day")
    return errors, max_error


def one(summary, section, case):
    matches = [r for r in summary if r["section"] == section and r["case"] == case]
    if len(matches) != 1:
        raise ValueError(f"expected one {section}/{case} row, got {len(matches)}")
    return matches[0]


def gate3_uv_lag(summary):
    errors = []
    cases = [("UV_50_50", 2.3, 0.10), ("UV_30_70", 1.2, 0.25)]
    for case, observed, tol in cases:
        predicted = float(one(summary, "lag", case)["value"])
        if abs(predicted - observed) > tol:
            errors.append(f"{case}: {predicted:.3f} not within {tol:.2f} d of {observed}")
    return errors


def gate4_maintenance(summary):
    predicted = float(one(summary, "lag", "maintenance")["value"])
    rel = abs(predicted - 0.825) / 0.825
    return ([] if rel <= 0.02 else
            [f"E0/beta={predicted:.6f}, relative error {rel:.2%} > 2%"]), rel


def gate5_optimum(summary):
    errors = []
    for row in [r for r in summary if r["section"] == "optimum"]:
        numeric = float(row["value"])
        closed = float(row["reference"])
        if abs(numeric - closed) >= 2e-4:
            errors.append(
                f"cycle {row['arg1']}h/{row['arg2']}d: "
                f"|{numeric:.6f}-{closed:.6f}| >= 2e-4"
            )
    corrected = one(summary, "erratum", "corrected_3h_6d")
    old = one(summary, "erratum", "old_3h_3d")
    cmu, cg = float(corrected["value"]), float(corrected["reference"])
    omu, og = float(old["value"]), float(old["reference"])
    if not (0.84 <= cmu <= 0.89 and 0.48 <= cg <= 0.52):
        errors.append(f"corrected worked point mu={cmu:.4f}, gamma={cg:.4f}")
    if 0.84 <= omu <= 0.89 or 0.48 <= og <= 0.52:
        errors.append(f"old 3-day point unexpectedly passes: mu={omu:.4f}, gamma={og:.4f}")
    return errors, (cmu, cg, omu, og)


def gate6_mechanism(summary):
    errors = []
    alpha = float(one(summary, "mechanism", "alpha_beta_over_gamma")["value"])
    wrong = float(one(summary, "mechanism", "wrong_alphaG_substitution")["value"])
    if not (1.10 <= alpha <= 1.18):
        errors.append(f"alpha=beta/gamma={alpha:.4f} fmol/CFU")
    if wrong / 0.43 <= 5.0:
        errors.append(f"dimensionless alpha_G guard too weak: wrong gamma={wrong:.3f}")
    return errors, (alpha, wrong)


def gate7_chassis():
    errors = []
    for script in ("check_job2.py", "check_job3.py", "check_job3b.py", "check_job3c.py"):
        result = subprocess.run(
            [sys.executable, str(ROOT / script)],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode != 0:
            errors.append(f"{script} exit {result.returncode}")
    return errors


def main() -> int:
    if not VIABILITY.is_file() or not SUMMARY.is_file():
        print("FAIL: missing Job 6 output; run Job6Sims first")
        return 1
    viability = load(VIABILITY)
    summary = load(SUMMARY)

    e1 = gate1_units(summary)
    e2, max_err = gate2_viability(viability)
    e3 = gate3_uv_lag(summary)
    e4, lag_rel = gate4_maintenance(summary)
    e5, worked = gate5_optimum(summary)
    e6, mechanism = gate6_mechanism(summary)
    e7 = gate7_chassis()

    failures = []
    for label, errors in (
        ("1 units / EV1", e1),
        ("2 viability integration", e2),
        ("3 UV-killed lag", e3),
        ("4 maintenance lag", e4),
        ("5 feast-famine optimum", e5),
        ("6 mechanism dimensions", e6),
        ("7 chassis regression", e7),
    ):
        print(f"Gate {label}: {'PASS' if not errors else 'FAIL'}")
        for error in errors:
            print("   ", error)
        failures.extend(errors)

    cmu, cg, omu, og = worked
    alpha, wrong = mechanism
    print()
    print(f"  viability RK4 max error = {max_err:.3e}")
    print(f"  maintenance lag relative error = {lag_rel:.2%}")
    print(f"  corrected 3h/6d: mu*={cmu:.4f} h^-1 gamma={cg:.4f} day^-1")
    print(f"  old 3h/3d:       mu*={omu:.4f} h^-1 gamma={og:.4f} day^-1")
    print(f"  alpha=beta/gamma={alpha:.4f} fmol/CFU; alpha_G misuse -> {wrong:.3f} day^-1")
    print()
    print(f"job6 OVERALL: {'FAIL' if failures else 'PASS'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
