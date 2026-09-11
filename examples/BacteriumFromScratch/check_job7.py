#!/usr/bin/env python3
"""Evaluate frozen Job 7 unbiased run-and-tumble gates."""

from __future__ import annotations

import csv
import hashlib
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DIR = ROOT / "results" / "job7_motility"


def load(name):
    with (DIR / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=";"))


def metrics(rows):
    return {r["metric"]: float(r["value"]) for r in rows}


def gate1_distributions(values):
    errors = []
    expected = {
        "mean_run": 2.39,
        "mean_tumble": 0.38,
        "mean_speed": 16.0,
        "speed_sd": 5.78,
    }
    for name, target in expected.items():
        rel = abs(values[name] - target) / target
        if rel > 0.02:
            errors.append(f"{name}={values[name]:.5f}, error {rel:.2%} > 2%")
    return errors


def gate2_run(values):
    errors = []
    if abs(values["run_fraction"] - 0.86) > 0.01:
        errors.append(f"run fraction {values['run_fraction']:.4f}")
    target = 16.0 * 2.39
    rel = abs(values["mean_run_length"] - target) / target
    if rel > 0.02:
        errors.append(
            f"run length {values['mean_run_length']:.3f}, error {rel:.2%}"
        )
    return errors


def gate3_diffusion(values):
    theory = values["D_theory"]
    simulated = values["D_simulated"]
    errors = []
    rel = abs(simulated - theory) / theory
    if rel > 0.07:
        errors.append(f"D_sim={simulated:.2f}, D_theory={theory:.2f}, error {rel:.2%}")
    if not 160.0 <= simulated <= 220.0:
        errors.append(f"D_sim={simulated:.2f} outside [160,220] um2/s")
    return errors, rel


def gate4_crossover(rows):
    by_t = {float(r["t_s"]): float(r["msd_um2"]) for r in rows}
    short = math.log(by_t[0.5] / by_t[0.1]) / math.log(0.5 / 0.1)
    long = math.log(by_t[1000.0] / by_t[100.0]) / math.log(10.0)
    errors = []
    if short <= 1.7:
        errors.append(f"short-time MSD slope {short:.3f} <= 1.7")
    if not 0.90 <= long <= 1.10:
        errors.append(f"long-time MSD slope {long:.3f} outside [0.90,1.10]")
    return errors, (short, long)


def gate5_isotropy_and_repeat(rows):
    errors = []
    final = rows[-1]
    rms = math.sqrt(float(final["msd_um2"]))
    means = [abs(float(final[k])) for k in ("mean_x_um", "mean_y_um", "mean_z_um")]
    for axis, value in zip("xyz", means):
        if value >= 0.02 * rms:
            errors.append(f"mean {axis}={value:.3f} >= 2% RMS ({0.02*rms:.3f})")

    paths = [DIR / "job7_summary.csv", DIR / "msd.csv"]
    before = [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]
    rerun = subprocess.run(
        ["cmd", "/c", str(ROOT / "compile_and_run.cmd"), "job7"],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if rerun.returncode:
        errors.append(f"deterministic rerun exit {rerun.returncode}")
    else:
        after = [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]
        if before != after:
            errors.append("seed-101 rerun changed output bytes")
    return errors, (rms, means)


def gate6_regression():
    errors = []
    for script in (
        "check_job2.py", "check_job3.py", "check_job3b.py",
        "check_job3c.py", "check_job5b.py", "check_job6.py",
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
    return errors


def main():
    required = [DIR / "job7_summary.csv", DIR / "msd.csv"]
    if any(not p.is_file() for p in required):
        print("FAIL: missing Job 7 output; run Job7Sims first")
        return 1
    summary = load("job7_summary.csv")
    msd = load("msd.csv")
    values = metrics(summary)

    e1 = gate1_distributions(values)
    e2 = gate2_run(values)
    e3, drel = gate3_diffusion(values)
    e4, slopes = gate4_crossover(msd)
    e5, isotropy = gate5_isotropy_and_repeat(msd)
    e6 = gate6_regression()

    failures = []
    for label, errors in (
        ("1 distributions", e1),
        ("2 run fraction / length", e2),
        ("3 effective diffusion", e3),
        ("4 crossover", e4),
        ("5 isotropy / reproducibility", e5),
        ("6 regression", e6),
    ):
        print(f"Gate {label}: {'PASS' if not errors else 'FAIL'}")
        for error in errors:
            print("   ", error)
        failures.extend(errors)

    short, long = slopes
    rms, means = isotropy
    print()
    print(f"  D_sim={values['D_simulated']:.2f}, D_theory={values['D_theory']:.2f}"
          f" um2/s ({drel:.2%})")
    print(f"  MSD slopes: short={short:.3f}, long={long:.3f}")
    print(f"  final RMS={rms:.2f} um; mean xyz={means}")
    print("  scope: basal motility only; no chemotaxis/chassis coupling")
    print()
    print(f"job7 OVERALL: {'FAIL' if failures else 'PASS'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
