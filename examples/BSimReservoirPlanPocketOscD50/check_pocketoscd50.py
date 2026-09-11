#!/usr/bin/env python3
"""NARMA-blind occupancy checker for PocketOsc-D50."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("PocketOsc-D50 checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "occupancy_screen.csv"
    if not path.exists():
        raise SystemExit("run occupancy screen first")
    rows = {row["arm"]: row for row in csv.DictReader(path.open(encoding="utf-8"))}
    flags = json.loads((RESULTS / "occupancy_flags.json").read_text(encoding="utf-8"))
    claim = rows["D05_H165_CLOSED"]["occupancy"]
    replay = rows["D005_H10_REPLAY"]["occupancy"]
    print(f"OCCUPANCY_D05_H165_CLOSED={claim}")
    print(f"OCCUPANCY_D005_H10_REPLAY={replay}")
    print(f"TRANSPORT_MODEL_STATUS={flags.get('TRANSPORT_MODEL_STATUS')}")
    print(f"LIVING_DECISION={flags.get('LIVING_DECISION')}")
    if replay in ("ALIVE", "SATURATED"):
        raise SystemExit("REPLAY ALIVE: T1c-class port broke.")
    if claim == "DEAD":
        print("STOP. Do not ease QS_KMLA. Do not raise d. No NARMA.")
        sys.exit(0)
    print("CLAIM not DEAD. Period-vs-flow and NARMA are still not authorized by this job.")


if __name__ == "__main__":
    main()
