#!/usr/bin/env python3
"""Evaluate frozen Job 5b mm-scale spatial-AHL transport gates."""

from __future__ import annotations

import csv
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DIR = ROOT / "results" / "job5b_spatial_ahl"


def load(name):
    with (DIR / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=";"))


def gate1_units(scales):
    row = scales[0]
    d = float(row["D_um2_s"])
    errors = []
    if abs(d - 550.0) > 1e-9:
        errors.append(f"1.98 mm2/h converted to {d}, not 550 um2/s")
    if abs(d - 490.0) / 490.0 > 0.15:
        errors.append(f"D={d:.1f} differs from Leaman 490 by >15%")
    if not 379.0 <= d <= 991.0:
        errors.append(f"D={d:.1f} outside Grant 379-991 bracket")
    return errors


def by_dx(rows, dx):
    matches = [r for r in rows if abs(float(r["dx_mm"]) - dx) < 1e-9]
    if len(matches) != 1:
        raise ValueError(f"expected one dx={dx} row, got {len(matches)}")
    return matches[0]


def gate2_conservation(summary):
    errors = []
    for row in summary:
        dx = float(row["dx_mm"])
        mass = float(row["mass_rel_err"])
        minimum = float(row["min_nM"])
        mean = float(row["final_mean_nM"])
        if mass >= 1e-11:
            errors.append(f"dx={dx}: mass error {mass:.3e} >= 1e-11")
        if not math.isfinite(minimum) or minimum < -1e-12:
            errors.append(f"dx={dx}: minimum {minimum}")
        if abs(mean - 4.0) >= 1e-10:
            errors.append(f"dx={dx}: mean {mean:.12f} != 4 nM")
    return errors


def gate3_grid(summary):
    a10 = float(by_dx(summary, 0.10)["arrival_x10_h"])
    a05 = float(by_dx(summary, 0.05)["arrival_x10_h"])
    movement = abs(a10 - a05)
    errors = [] if movement < 0.05 else [
        f"arrival moved {movement:.4f} h (0.10 -> 0.05 mm), limit 0.05"
    ]
    return errors, (a10, a05, movement)


def gate4_profile():
    coarse = [float(r["C_nM"]) for r in load("profile_dx100_t10h.csv")]
    fine = [float(r["C_nM"]) for r in load("profile_dx050_t10h.csv")]
    if len(fine) != 2 * len(coarse):
        return [f"profile sizes coarse={len(coarse)} fine={len(fine)}"], math.nan
    restricted = [0.5 * (fine[2 * i] + fine[2 * i + 1])
                  for i in range(len(coarse))]
    mse = sum((a - b) ** 2 for a, b in zip(coarse, restricted)) / len(coarse)
    rmse = math.sqrt(mse)
    scale = max(restricted) - min(restricted)
    nrmse = rmse / scale
    return ([] if nrmse < 0.01 else
            [f"profile NRMSE {nrmse:.4e} >= 0.01"]), nrmse


def gate5_geometry(scales):
    row = scales[0]
    t10 = float(row["t_diff_10mm_h"])
    t60 = float(row["t_diff_60um_s"])
    errors = []
    if not 10.0 <= t10 <= 15.0:
        errors.append(f"10 mm scale {t10:.3f} h outside [10,15]")
    if not t60 < 2.0:
        errors.append(f"60 um scale {t60:.3f} s is not <2 s")
    return errors, (t10, t60)


def gate6_regression():
    errors = []
    for script in (
        "check_job2.py", "check_job3.py", "check_job3b.py",
        "check_job3c.py", "check_job6.py",
    ):
        result = subprocess.run(
            [sys.executable, str(ROOT / script)],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode:
            errors.append(f"{script} exit {result.returncode}")
    if not (ROOT / "JOB5_STANDING.md").is_file():
        errors.append("missing Job 5 partial standing")
    return errors


def main():
    required = [
        "job5b_summary.csv", "geometry_scales.csv",
        "profile_dx100_t10h.csv", "profile_dx050_t10h.csv",
    ]
    missing = [name for name in required if not (DIR / name).is_file()]
    if missing:
        print("FAIL: missing", ", ".join(missing))
        return 1

    summary = load("job5b_summary.csv")
    scales = load("geometry_scales.csv")
    e1 = gate1_units(scales)
    e2 = gate2_conservation(summary)
    e3, arrivals = gate3_grid(summary)
    e4, nrmse = gate4_profile()
    e5, times = gate5_geometry(scales)
    e6 = gate6_regression()

    failures = []
    for label, errors in (
        ("1 units / literature range", e1),
        ("2 mass / positivity", e2),
        ("3 grid convergence", e3),
        ("4 profile convergence", e4),
        ("5 geometry separation", e5),
        ("6 regression", e6),
    ):
        print(f"Gate {label}: {'PASS' if not errors else 'FAIL'}")
        for error in errors:
            print("   ", error)
        failures.extend(errors)

    a10, a05, movement = arrivals
    t10, t60 = times
    print()
    print(f"  C(10 mm)>=1.5 nM: dx=.10 {a10:.4f} h; dx=.05 {a05:.4f} h"
          f"; movement {movement:.4f} h")
    print(f"  profile NRMSE (0.10 vs restricted 0.05) = {nrmse:.3e}")
    print(f"  diffusion scales: 10 mm={t10:.3f} h; 60 um={t60:.3f} s")
    print("  scope: transport PASS does not imply spatial variation in Job 3c")
    print()
    print(f"job5b OVERALL: {'FAIL' if failures else 'PASS'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
