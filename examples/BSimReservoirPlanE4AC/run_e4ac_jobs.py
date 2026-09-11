#!/usr/bin/env python3
"""Run E4.1 A1 living BSim. Seeds 222/333 are not started from this list."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CP = str(HERE.parent) + ";" + str(ROOT / "dist" / "build") + ";" + ";".join(
    str(ROOT / "lib" / name) for name in ("core.jar", "vecmath.jar", "objimport.jar")
)

SMOKE = ["sim_config_e4ac_smoke.properties"]
SCOUT = [
    "sim_config_e4ac_driven_seed111.properties",
    "sim_config_e4ac_brownian_seed111.properties",
]


def run_job(config: str) -> tuple[str, int]:
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
                "BSimReservoirPlanE4AC.BSimReservoirPlanE4AC",
                config,
            ],
            cwd=HERE,
            stdout=out,
            stderr=err,
        )
    print(f"done {config} exit={proc.returncode}", flush=True)
    return config, proc.returncode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("configs", nargs="*")
    args = parser.parse_args()
    jobs = args.configs or (SMOKE if args.smoke else SCOUT)
    failed = []
    workers = min(2, len(jobs))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_job, job) for job in jobs]
        for future in as_completed(futures):
            config, code = future.result()
            if code != 0:
                failed.append((config, code))
    if failed:
        for job, code in failed:
            print(f"FAIL {job} {code}", file=sys.stderr)
        raise SystemExit(1)
    print(f"completed {len(jobs)} E4AC job(s)")


if __name__ == "__main__":
    main()
