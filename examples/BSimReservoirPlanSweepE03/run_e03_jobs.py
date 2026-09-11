#!/usr/bin/env python3
"""Run the 16 frozen E0.3 production BSim jobs, at most four at a time."""

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
    "sim_config_e03_center_f0p25_driven_seed111.properties",
    "sim_config_e03_center_f0p5_driven_seed111.properties",
    "sim_config_e03_center_f0p666667_driven_seed111.properties",
    "sim_config_e03_center_f1p0_driven_seed111.properties",
    "sim_config_e03_upstream_f0p0_driven_seed111.properties",
    "sim_config_e03_upstream_f0p666667_driven_seed111.properties",
    "sim_config_e03_downstream_f0p0_driven_seed111.properties",
    "sim_config_e03_downstream_f0p666667_driven_seed111.properties",
    "sim_config_e03_f0p25_silent_seed111.properties",
    "sim_config_e03_f0p5_silent_seed111.properties",
    "sim_config_e03_f0p666667_silent_seed111.properties",
    "sim_config_e03_f1p0_silent_seed111.properties",
    "sim_config_e03_f0p25_brownian_seed111.properties",
    "sim_config_e03_f0p5_brownian_seed111.properties",
    "sim_config_e03_f0p666667_brownian_seed111.properties",
    "sim_config_e03_f1p0_brownian_seed111.properties",
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
                "BSimReservoirPlanSweepE03.BSimReservoirPlanSweepE03",
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
    print("all 16 production jobs completed")


if __name__ == "__main__":
    main()
