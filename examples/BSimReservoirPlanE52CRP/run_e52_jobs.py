#!/usr/bin/env python3
"""Run E5.2 living BSim. Confirmation starts only if living n_ALIVE >= 10."""

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
DEV = [1, 4, 6, 7, 9, 13, 16, 18, 20, 22, 24, 25]
CONF = [8, 12, 14, 15, 17, 19, 21, 23]


def compile_java() -> None:
    cmd = [
        "javac",
        "-encoding",
        "UTF-8",
        "-cp",
        CP,
        "-d",
        str(HERE.parent),
        str(HERE / "BSimReservoirPlanE52CRP.java"),
        str(HERE / "VoxelAnalyzer.java"),
    ]
    print("compile", flush=True)
    proc = subprocess.run(cmd, cwd=ROOT)
    if proc.returncode != 0:
        raise SystemExit("STOP javac failed")


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
                "BSimReservoirPlanE52CRP.BSimReservoirPlanE52CRP",
                config,
            ],
            cwd=HERE,
            stdout=out,
            stderr=err,
        )
    print(f"done {config} exit={proc.returncode}", flush=True)
    return proc.returncode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--compile", action="store_true")
    parser.add_argument("--development", action="store_true")
    parser.add_argument("--confirmation", action="store_true")
    parser.add_argument("configs", nargs="*")
    args = parser.parse_args()
    freeze = HERE / "results" / "e52_encoding_freeze.json"
    if freeze.exists():
        import json
        rec = json.loads(freeze.read_text(encoding="utf-8"))
        if rec.get("status") != "ENCODING_FROZEN":
            print("ENCODING_FAIL: living BSim is not started.", flush=True)
            raise SystemExit(0)
        print("frozen map", rec.get("winner"), flush=True)
    compile_java()
    if args.smoke:
        jobs = ["sim_config_e52_driven_p0001_seed111.properties"]
    elif args.development:
        jobs = [f"sim_config_e52_driven_p{p:04d}_seed111.properties" for p in DEV]
    elif args.confirmation:
        jobs = [f"sim_config_e52_driven_p{p:04d}_seed111.properties" for p in CONF]
    else:
        jobs = args.configs
    if args.compile and not jobs and not args.smoke and not args.development and not args.confirmation:
        print("compile only")
        return
    if not jobs:
        print("no jobs; pass --smoke, --development, --confirmation, or config names")
        raise SystemExit(2)
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
    print(f"completed {len(jobs)} E52 job(s)")


if __name__ == "__main__":
    main()
