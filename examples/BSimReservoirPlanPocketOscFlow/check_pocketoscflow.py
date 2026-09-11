#!/usr/bin/env python3
"""NARMA-blind occupancy/period checker for PocketOsc-Flow."""

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
        raise SystemExit("PocketOsc-Flow checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "occupancy_screen.csv"
    if not path.exists():
        raise SystemExit("run occupancy screen first")
    rows = {row["arm"]: row for row in csv.DictReader(path.open(encoding="utf-8"))}
    flags = json.loads((RESULTS / "occupancy_flags.json").read_text(encoding="utf-8"))
    closed = rows["CLOSED_D50"]["occupancy"]
    closed_w = rows["CLOSED_D50"]["occupancy_d50_window"]
    lo = rows["FLOW_180"]["occupancy"]
    hi = rows["FLOW_296"]["occupancy"]
    print(f"OCCUPANCY_CLOSED_D50={closed}")
    print(f"OCCUPANCY_CLOSED_D50_WINDOW={closed_w}")
    print(f"OCCUPANCY_FLOW_180={lo}")
    print(f"OCCUPANCY_FLOW_296={hi}")
    print(f"TRANSPORT_MODEL_STATUS={flags.get('TRANSPORT_MODEL_STATUS')}")
    print(f"LIVING_DECISION={flags.get('LIVING_DECISION')}")
    if closed_w == "DEAD":
        raise SystemExit("CLOSED D50 window DEAD: port broke.")
    if lo == "DEAD" or hi == "DEAD":
        print("STOP. Open flow DEAD. Do not ease QS_KMLA. Do not raise d. No NARMA.")
        sys.exit(0)
    print("Open arms not DEAD. GFP maps and NARMA are still not authorized by this job.")
    print(f"LIVING_DECISION={flags.get('LIVING_DECISION')}")


if __name__ == "__main__":
    main()
