#!/usr/bin/env python3
"""Run the three frozen WashoutReset W=5 production BSim jobs in parallel."""

from __future__ import annotations

import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CP = str(HERE.parent) + ";" + str(ROOT / "dist" / "build") + ";" + ";".join(
    str(ROOT / "lib" / name) for name in ("core.jar", "vecmath.jar", "objimport.jar")
)
JOBS = [
    "sim_config_wash_decay_w5_seed101.properties",
    "sim_config_wash_trickle_w5_seed101.properties",
    "sim_config_wash_flush_w5_seed101.properties",
]
FOLLOW_ON = [
    "sim_config_wash_decay_w15_seed101.properties",
    "sim_config_wash_flush_w15_seed101.properties",
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
                "BSimReservoirPlanWashoutReset.BSimReservoirPlanWashoutReset",
                config,
            ],
            cwd=HERE,
            stdout=out,
            stderr=err,
        )
    print(f"done {config} exit={proc.returncode}", flush=True)
    return config, proc.returncode


def main() -> None:
    jobs = FOLLOW_ON if "--w15" in sys.argv else JOBS
    failed = []
    with ThreadPoolExecutor(max_workers=min(3, len(jobs))) as pool:
        futures = [pool.submit(run_job, job) for job in jobs]
        for future in as_completed(futures):
            config, code = future.result()
            if code != 0:
                failed.append((config, code))
    if failed:
        for config, code in failed:
            print(f"FAIL {config} {code}", file=sys.stderr)
        raise SystemExit(1)
    print(f"all {len(jobs)} jobs completed")


if __name__ == "__main__":
    main()
