#!/usr/bin/env python3
"""Run E4.2 living BSim sequentially. Seeds 222/333 are not started."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CP = str(HERE.parent) + ";" + str(ROOT / "dist" / "build") + ";" + ";".join(
    str(ROOT / "lib" / name) for name in ("core.jar", "vecmath.jar", "objimport.jar")
)


def run_job(config: str) -> int:
    log_dir = HERE / "logs"
    log_dir.mkdir(exist_ok=True)
    stdout = log_dir / f"{Path(config).stem}.stdout.log"
    stderr = log_dir / f"{Path(config).stem}.log.err"
    print(f"start {config}", flush=True)
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        proc = subprocess.run(
            [
                "java",
                "-cp",
                CP,
                "BSimReservoirPlanE42TW.BSimReservoirPlanE42TW",
                config,
            ],
            cwd=HERE,
            stdout=out,
            stderr=err,
        )
    print(f"done {config} exit={proc.returncode}", flush=True)
    return proc.returncode


def compile_java() -> None:
    cmd = [
        "javac",
        "-cp",
        CP,
        "-d",
        str(HERE.parent),
        str(HERE / "BSimReservoirPlanE42TW.java"),
        str(HERE / "VoxelAnalyzer.java"),
    ]
    print("compile", flush=True)
    proc = subprocess.run(cmd, cwd=ROOT)
    if proc.returncode != 0:
        raise SystemExit("STOP javac failed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--compile", action="store_true")
    parser.add_argument("configs", nargs="*")
    args = parser.parse_args()
    if args.compile or args.smoke or not args.configs:
        compile_java()
    jobs = args.configs or (
        ["sim_config_e42_smoke.properties"] if args.smoke
        else ["sim_config_e42_closed_seed111.properties"]
    )
    failed = []
    for job in jobs:
        code = run_job(job)
        if code != 0:
            failed.append((job, code))
            break
    if failed:
        for job, code in failed:
            print(f"FAIL {job} {code}", file=sys.stderr)
        raise SystemExit(1)
    print(f"completed {len(jobs)} E42 job(s)")


if __name__ == "__main__":
    main()
