#!/usr/bin/env python3
"""NARMA-blind occupancy/period checker for PocketOsc-SI."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
COVER = ("MU_0", "MU_025", "MU_05", "MU_1", "MU_2")


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("PocketOsc-SI checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "occupancy_screen.csv"
    if not path.exists():
        raise SystemExit("run occupancy screen first")
    rows = {row["arm"]: row for row in csv.DictReader(path.open(encoding="utf-8"))}
    flags = json.loads((RESULTS / "occupancy_flags.json").read_text(encoding="utf-8"))
    for arm in COVER + ("MU_CSTR_MIN",):
        rec = rows[arm]
        print(f"{arm} mu={rec['mu']} flag={rec['flag']} T={rec['period']} mean_I={rec['mean_I']}")
    print(f"OSC_COVER={flags.get('OSC_COVER')}")
    print(f"LIVING_DECISION={flags.get('LIVING_DECISION')}")
    decision = flags.get("LIVING_DECISION")
    if decision == "STOP_AFTER_SI_BULK_DEAD":
        print("STOP. No COVER oscillation. Do not retune. No NARMA.")
        sys.exit(0)
    if decision == "STOP_AFTER_SI_NO_IDENTITY":
        print("STOP. Oscillation without T(mu) identity. Do not retune. No NARMA.")
        sys.exit(0)
    print("SI bulk identity scored. GFP maps, Fig. 2c, spatial N=200, and NARMA stay unauthorized.")


if __name__ == "__main__":
    main()
