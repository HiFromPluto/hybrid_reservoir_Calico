#!/usr/bin/env python3
"""NARMA-blind occupancy checker for PocketFill-T1c."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def refuse_narma() -> None:
    if "narma" in " ".join(sys.argv).lower() or "mackey" in " ".join(sys.argv).lower():
        raise SystemExit("PocketFill-T1c checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "occupancy_screen.csv"
    if not path.exists():
        raise SystemExit("run occupancy screen first")
    rows = {row["arm"]: row for row in csv.DictReader(path.open(encoding="utf-8"))}
    flags = json.loads((RESULTS / "occupancy_flags.json").read_text(encoding="utf-8"))
    flush = rows["VOLUMETRIC_FLUSH"]["occupancy"]
    closed = rows["VOLUMETRIC_CLOSED"]["occupancy"]
    w20 = rows["VOLUMETRIC_W20"]["occupancy"]
    print(f"OCCUPANCY_VOLUMETRIC_FLUSH={flush}")
    print(f"OCCUPANCY_VOLUMETRIC_CLOSED={closed}")
    print(f"OCCUPANCY_VOLUMETRIC_W20={w20}")
    print(f"TRANSPORT_MODEL_STATUS={flags.get('TRANSPORT_MODEL_STATUS')}")
    print(f"LIVING_DECISION={flags.get('LIVING_DECISION')}")
    if flush in ("ALIVE", "SATURATED"):
        raise SystemExit("FLUSH ALIVE: T1 port broke.")
    if closed == "DEAD":
        print("STOP Track C. Do not ease QS_KMLA.")
        sys.exit(0)
    print("CLOSED not DEAD. NARMA is still not authorized by this job.")


if __name__ == "__main__":
    main()
