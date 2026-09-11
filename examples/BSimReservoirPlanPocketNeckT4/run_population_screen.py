#!/usr/bin/env python3
"""Run the two predeclared PocketNeck-T4 population arms.

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
    ("LIVE_W20_GROWTH", "config/live_w20_growth.properties"),
    ("LIVE_W100_GROWTH", "config/live_w100_growth.properties"),
)
FORBIDDEN = ("narma", "mackey-glass", "mackey_glass", "waveform_auc", "ridge_lambda")


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    for token in FORBIDDEN:
        if token in blob:
            raise SystemExit("PocketNeck-T4 runner cannot see a NARMA target.")


def compile_job() -> None:
    cmd = ["cmd", "/c", "compile_and_run.cmd"]
    print("Compiling PocketNeckT4...")
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


def main() -> int:
    refuse_narma()
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--skip-compile", action="store_true")
    parser.add_argument("--from-arm", default="")
    parser.add_argument("--primary-only", action="store_true")
    args = parser.parse_args()
    if not args.skip_compile:
        compile_job()
    started = not args.from_arm
    for arm, config in ARM_CONFIGS:
        if args.primary_only and arm != "LIVE_W20_GROWTH":
            continue
        if not started:
            if arm == args.from_arm:
                started = True
            else:
                continue
        run_arm(arm, config, args.smoke)
        if arm == "LIVE_W20_GROWTH" and not args.smoke:
            occ = occupancy_of(arm)
            if occ == "DEAD":
                print("LIVE_W20_GROWTH occupancy DEAD. Stop. Do not raise J_max.")
                check = subprocess.run(
                    [sys.executable, str(HERE / "check_pocketneck_t4.py")],
                    cwd=HERE,
                )
                return check.returncode if check.returncode != 0 else 4
    if args.smoke:
        print("smoke complete; not population evidence")
        return 0
    check = subprocess.run(
        [sys.executable, str(HERE / "check_pocketneck_t4.py")],
        cwd=HERE,
    )
    return check.returncode


if __name__ == "__main__":
    raise SystemExit(main())
