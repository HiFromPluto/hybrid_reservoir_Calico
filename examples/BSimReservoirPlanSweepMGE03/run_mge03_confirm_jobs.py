#!/usr/bin/env python3
"""Seed 202/303 confirmation BSim for SweepMGE03 story-move conditions only."""

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
    "sim_config_mge03_center_f0p25_driven_seed202.properties",
    "sim_config_mge03_center_f0p25_silent_seed202.properties",
    "sim_config_mge03_center_f0p25_brownian_seed202.properties",
    "sim_config_mge03_upstream_f0p0_driven_seed202.properties",
    "sim_config_mge03_center_f0p25_driven_seed303.properties",
    "sim_config_mge03_center_f0p25_silent_seed303.properties",
    "sim_config_mge03_center_f0p25_brownian_seed303.properties",
    "sim_config_mge03_upstream_f0p0_driven_seed303.properties",
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
                "BSimReservoirPlanSweepMGE03.BSimReservoirPlanSweepMGE03",
                config,
            ],
            cwd=HERE,
            stdout=out,
            stderr=err,
        )
    print(f"done {config} exit={proc.returncode}", flush=True)
    return config, proc.returncode


def main() -> None:
    failed = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(run_job, job) for job in JOBS]
        for future in as_completed(futures):
            config, code = future.result()
            if code != 0:
                failed.append((config, code))
    if failed:
        for config, code in failed:
            print(f"FAIL {config} {code}", file=sys.stderr)
        raise SystemExit(1)
    print("all 8 confirmation jobs completed")


if __name__ == "__main__":
    main()
