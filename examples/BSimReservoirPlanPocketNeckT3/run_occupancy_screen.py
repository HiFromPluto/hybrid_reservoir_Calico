#!/usr/bin/env python3
"""Run the four predeclared PocketNeck-T3 occupancy arms.

Not a task scout. --smoke integrates 50 s and must not be scored.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARM_CONFIGS = (
    ("FIELD_W100_FLUSH", "config/field_w100.properties"),
    ("FIELD_W20_L20", "config/field_w20.properties"),
    ("LIVE_W100_FLUSH", "config/live_w100.properties"),
    ("LIVE_W20_L20", "config/live_w20.properties"),
)
FORBIDDEN = ("narma", "mackey-glass", "mackey_glass", "waveform_auc", "ridge_lambda")


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    for token in FORBIDDEN:
        if token in blob:
            raise SystemExit("PocketNeck-T3 runner cannot see a NARMA target.")


def compile_job() -> None:
    cmd = ["cmd", "/c", "compile_and_run.cmd"]
    print("Compiling PocketNeckT3...")
    subprocess.run(cmd, cwd=HERE, check=True)


def run_arm(arm: str, config: str, smoke: bool) -> None:
    cmd = ["cmd", "/c", "compile_and_run.cmd", config.replace("/", "\\")]
    if smoke:
        cmd.append("smoke")
    print(f"Running {arm} ({'smoke' if smoke else 'occupancy'})...")
    subprocess.run(cmd, cwd=HERE, check=True)


def main() -> int:
    refuse_narma()
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--skip-compile", action="store_true")
    parser.add_argument("--from-arm", default="")
    args = parser.parse_args()
    if not args.skip_compile:
        compile_job()
    started = not args.from_arm
    for arm, config in ARM_CONFIGS:
        if not started:
            if arm == args.from_arm:
                started = True
            else:
                continue
        run_arm(arm, config, args.smoke)
        if arm == "FIELD_W100_FLUSH" and not args.smoke:
            summary = HERE / "results" / "runs" / arm / "summary.csv"
            if not summary.exists():
                raise SystemExit("FIELD_W100_FLUSH summary missing")
    if args.smoke:
        print("smoke complete; not occupancy evidence")
        return 0
    check = subprocess.run(
        [sys.executable, str(HERE / "check_pocketneck_t3.py")],
        cwd=HERE,
    )
    return check.returncode


if __name__ == "__main__":
    raise SystemExit(main())
