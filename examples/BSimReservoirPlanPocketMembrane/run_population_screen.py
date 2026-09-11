#!/usr/bin/env python3
"""Run the predeclared PocketMembrane population arm.

Not a task scout. --smoke integrates 50 s and must not be scored.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARM_CONFIGS = (
    ("LIVE_W20_MEM", "config/live_w20_mem.properties"),
)
FORBIDDEN = ("narma", "mackey-glass", "mackey_glass", "waveform_auc", "ridge_lambda")


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    for token in FORBIDDEN:
        if token in blob:
            raise SystemExit("PocketMembrane runner cannot see a NARMA target.")


def compile_job() -> None:
    cmd = ["cmd", "/c", "compile_and_run.cmd"]
    print("Compiling PocketMembrane...")
    subprocess.run(cmd, cwd=HERE, check=True)


def run_arm(arm: str, config: str, smoke: bool) -> None:
    cmd = ["cmd", "/c", "compile_and_run.cmd", config.replace("/", "\\")]
    if smoke:
        cmd.append("smoke")
    print(f"Running {arm} ({'smoke' if smoke else 'population'})...")
    subprocess.run(cmd, cwd=HERE, check=True)


def occupancy_of(arm: str) -> str:
    path = HERE / "results" / "runs" / arm / "summary.csv"
    if not path.exists():
        return "NA"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return "NA"
    return rows[0].get("Occupancy", "NA")


def spillover_of(arm: str) -> int:
    path = HERE / "results" / "runs" / arm / "summary.csv"
    if not path.exists():
        return -1
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return -1
    try:
        return int(float(rows[0].get("Spillover", "-1")))
    except (TypeError, ValueError):
        return -1


def main() -> int:
    refuse_narma()
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--skip-compile", action="store_true")
    args = parser.parse_args()
    if not args.skip_compile:
        compile_job()
    for arm, config in ARM_CONFIGS:
        run_arm(arm, config, args.smoke)
        if not args.smoke:
            occ = occupancy_of(arm)
            if occ == "DEAD":
                print("LIVE_W20_MEM occupancy DEAD. Stop. Do not raise J_max.")
                check = subprocess.run(
                    [sys.executable, str(HERE / "check_pocketmembrane.py")],
                    cwd=HERE,
                )
                return check.returncode if check.returncode != 0 else 4
            spill = spillover_of(arm)
            if spill > 0:
                print("LIVE_W20_MEM spillover > 0. Membrane leaked. Fix bounce.")
                check = subprocess.run(
                    [sys.executable, str(HERE / "check_pocketmembrane.py")],
                    cwd=HERE,
                )
                return check.returncode if check.returncode != 0 else 6
    if args.smoke:
        print("smoke complete; not population evidence")
        return 0
    check = subprocess.run(
        [sys.executable, str(HERE / "check_pocketmembrane.py")],
        cwd=HERE,
    )
    return check.returncode


if __name__ == "__main__":
    raise SystemExit(main())
