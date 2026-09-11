#!/usr/bin/env python3
"""One-shot Lane A seed-replicate runner. Resume-safe. Do not retune."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
LOG_DIR = RESULTS / "seed_replicate"
PROTOCOL = HERE / "PROTOCOL_SEED_REPLICATE.md"
PROTOCOL_JSON = HERE / "configs" / "seed_replicate_protocol.json"

SEEDS = (202, 303)
TASKS = (
    ("lane-a-carrier", "lanea.carrier.args", "java_LANE_A_CARRIER_DRIVEN.csv", 128),
    ("lane-a-narma10", "lanea.narma10.args", "java_LANE_A_NARMA10_DRIVEN.csv", 3200),
    ("lane-a-mg", "lanea.mg.args", "java_LANE_A_MG_DRIVEN.csv", 3200),
    ("lane-a-lorenz", "lanea.lorenz.args", "java_LANE_A_LORENZ_DRIVEN.csv", 3200),
    ("lane-a-waveform", "lanea.waveform.args", "java_LANE_A_WAVEFORM_DRIVEN.csv", 6400),
)

SEED101 = {
    "java_LANE_A_CARRIER_DRIVEN.csv": 128,
    "java_LANE_A_NARMA10_DRIVEN.csv": 3200,
    "java_LANE_A_MG_DRIVEN.csv": 3200,
    "java_LANE_A_LORENZ_DRIVEN.csv": 3200,
    "java_LANE_A_WAVEFORM_DRIVEN.csv": 6400,
    "java_LANE_A_CARRIER_SILENT.csv": 128,
    "java_LANE_A_NARMA10_SILENT.csv": 3200,
    "java_LANE_A_WAVEFORM_SILENT.csv": 6400,
}


def log(handle, message: str) -> None:
    line = f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}  {message}"
    print(line, flush=True)
    handle.write(line + "\n")
    handle.flush()


def csv_complete(path: Path, n_rows: int) -> bool:
    if not path.exists():
        return False
    with path.open(encoding="utf-8") as fh:
        n = sum(1 for _ in fh) - 1
    return n >= n_rows


def find_ant() -> list[str]:
    which = shutil.which("ant")
    if which:
        return [which]
    candidates = []
    temp = os.environ.get("TEMP") or os.environ.get("TMP") or ""
    if temp:
        candidates.append(
            Path(temp) / "cursor-apache-ant-1.10.17" / "apache-ant-1.10.17" / "bin" / "ant.bat"
        )
    local_app = Path.home() / "AppData" / "Local" / "Temp"
    candidates.append(
        local_app / "cursor-apache-ant-1.10.17" / "apache-ant-1.10.17" / "bin" / "ant.bat"
    )
    home = Path.home()
    candidates.extend(home.glob("**/apache-ant-1.10.17/bin/ant.bat"))
    for path in candidates:
        if path.exists():
            return [str(path)]
    raise SystemExit("ant not found. Put Apache Ant on PATH, then rerun this script.")


def require_frozen() -> None:
    text = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in text or '"frozen_before_traces": true' not in js:
        raise SystemExit("PROTOCOL_SEED_REPLICATE is not frozen_before_traces")
    if '"parent_rewrite": false' not in js or '"seeds": [101, 202, 303]' not in js:
        raise SystemExit("seed list / parent_rewrite freeze drifted")
    if '"waveform_primary": "block_brier"' not in js:
        raise SystemExit("waveform primary must stay block_brier")


def require_seed101() -> None:
    missing = [name for name, n in SEED101.items() if not csv_complete(RESULTS / name, n)]
    if missing:
        raise SystemExit(
            "seed 101 maps missing or short: " + ", ".join(missing)
            + ". Do not rerun 101 from this extra unless the parent standing said to."
        )


def run_job(ant: list[str], target: str, prop: str, seed: int, log_fh) -> int:
    args = f"--seed {seed} --driven-only"
    cmd = ant + ["-D" + prop + "=" + args, target]
    log(log_fh, "RUN " + " ".join(cmd))
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(ROOT))
    log(log_fh, f"EXIT {target} seed={seed} code={proc.returncode} ({time.time() - t0:.0f}s)")
    return proc.returncode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--java-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    require_frozen()
    require_seed101()
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / "run.log"
    status = {"jobs": [], "started": datetime.now(timezone.utc).isoformat()}

    with log_path.open("a", encoding="utf-8") as log_fh:
        log(log_fh, "=== LaneA_SEED_REPLICATE start ===")
        if not args.check_only:
            ant = find_ant()
            log(log_fh, "ant=" + " ".join(ant))
            compile_cmd = ant + ["compile"]
            log(log_fh, "RUN " + " ".join(compile_cmd))
            if subprocess.run(compile_cmd, cwd=str(ROOT)).returncode != 0:
                raise SystemExit("ant compile failed")
            for seed in SEEDS:
                out = RESULTS / f"seed_{seed}"
                out.mkdir(parents=True, exist_ok=True)
                for target, prop, name, n_rows in TASKS:
                    dest = out / name
                    if csv_complete(dest, n_rows) and not args.force:
                        log(log_fh, f"SKIP {target} seed={seed} ({dest.name} complete)")
                        status["jobs"].append({"target": target, "seed": seed, "status": "skip"})
                        continue
                    code = run_job(ant, target, prop, seed, log_fh)
                    status["jobs"].append({"target": target, "seed": seed, "status": "ok" if code == 0 else "fail"})
                    if code != 0:
                        (LOG_DIR / "STATUS.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
                        raise SystemExit(
                            f"{target} seed={seed} failed. Fix the error, rerun this same script."
                        )
                    if not csv_complete(dest, n_rows):
                        raise SystemExit(f"{dest} is short after a 0 exit")

        if args.java_only:
            log(log_fh, "java-only; checker not run")
            (LOG_DIR / "STATUS.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
            return

        log(log_fh, "RUN checker")
        check = subprocess.run(
            [sys.executable, str(HERE / "check_lane_a_seed_replicate.py")],
            cwd=str(ROOT),
        )
        log(log_fh, f"EXIT checker code={check.returncode}")
        status["checker"] = check.returncode
        (LOG_DIR / "STATUS.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
        if check.returncode != 0:
            raise SystemExit("checker failed")
        log(log_fh, "=== LaneA_SEED_REPLICATE done ===")


if __name__ == "__main__":
    main()
