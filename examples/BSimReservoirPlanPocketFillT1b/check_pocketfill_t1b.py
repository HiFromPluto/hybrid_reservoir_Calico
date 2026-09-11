#!/usr/bin/env python3
"""NARMA-blind occupancy checker for PocketFill-T1b."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("PocketFill-T1b checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "occupancy_screen.csv"
    if not path.exists():
        raise SystemExit("run occupancy screen first")
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    by = {row["arm"]: row for row in rows}
    zero = by["VOLUMETRIC_ZERO"]["occupancy"]
    pulse = by["SEED_PULSE"]["occupancy"]
    bath = by["SEED_BATH"]["occupancy"]
    status = "NA"
    flags_path = RESULTS / "occupancy_flags.json"
    if flags_path.exists():
        import json

        status = json.loads(flags_path.read_text(encoding="utf-8")).get(
            "TRANSPORT_MODEL_STATUS", "NA"
        )
    if zero in ("ALIVE", "SATURATED"):
        raise SystemExit("VOLUMETRIC_ZERO ALIVE: T1 port broke. Do not celebrate.")
    print(f"OCCUPANCY_VOLUMETRIC_ZERO={zero}")
    print(f"OCCUPANCY_SEED_PULSE={pulse}")
    print(f"OCCUPANCY_SEED_BATH={bath}")
    print(f"TRANSPORT_MODEL_STATUS={status}")
    print("LIVING_DECISION=STOP_AFTER_T1B_SEED")
    if pulse == "DEAD":
        print("STOP. Do not retune QS_KMLA or the seed. Do not start NARMA.")
        if bath in ("ALIVE", "SATURATED"):
            print("SEED_BATH occupied under a maintained bath; leak still wins on the pulse.")
        sys.exit(0)
    print("SEED_PULSE occupancy not DEAD. NARMA is still not authorized by this job.")


if __name__ == "__main__":
    main()
