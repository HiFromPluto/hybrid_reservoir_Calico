#!/usr/bin/env python3
"""NARMA-blind occupancy checker for PocketFill-T1."""

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
        raise SystemExit("PocketFill-T1 checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "occupancy_screen.csv"
    if not path.exists():
        raise SystemExit("run occupancy screen first")
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    by = {row["arm"]: row for row in rows}
    empty = by["EMPTY_POINT"]["occupancy"]
    vol = by["VOLUMETRIC_LUXI"]["occupancy"]
    if empty == "ALIVE":
        raise SystemExit("EMPTY_POINT ALIVE: T0 replay is wrong. Do not celebrate.")
    print(f"EMPTY_POINT={empty} PACKED_PLUG={by['PACKED_PLUG']['occupancy']} VOLUMETRIC_LUXI={vol}")
    if vol == "DEAD":
        print("STOP. Do not retune QS_KMLA or N_pack. Do not start NARMA.")
        sys.exit(0)
    print("VOLUMETRIC occupancy not DEAD. Period-vs-flow is still not authorized by this job.")


if __name__ == "__main__":
    main()
