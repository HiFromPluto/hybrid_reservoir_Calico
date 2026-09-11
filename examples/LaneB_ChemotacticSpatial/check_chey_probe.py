#!/usr/bin/env python3
"""B4_CHEY_PROBE checker. Not colony kappa. Parents B0–B3 FAIL kept."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_CHEY_PROBE.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_chey.json"

HOLD_LO, HOLD_HI = 0.85, 1.15
STEP_UP_MAX = 0.88
STEP_DOWN_MIN = 1.12
WIN_LO, WIN_HI = 5.4, 7.0


def main() -> int:
    if any(k in " ".join(sys.argv).lower() for k in ("narma", "charc", "ipc")):
        raise SystemExit("B4 checker refuses NARMA/CHARC/IPC")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("B4 PROTOCOL is not frozen")
    if js.get("status_label") != "B4_CHEY_PROBE" or js.get("narma") is not False:
        raise SystemExit("B4 must freeze probe and narma=false")

    hold = load(RESULTS / "b4_hold.csv")
    up = load(RESULTS / "b4_step_up.csv")
    down = load(RESULTS / "b4_step_down.csv")
    hold_mean = mean_y(hold, 3.0, 15.0)
    up_min = min_y(up, WIN_LO, WIN_HI)
    down_max = max_y(down, WIN_LO, WIN_HI)
    hold_ok = HOLD_LO <= hold_mean <= HOLD_HI
    up_ok = up_min <= STEP_UP_MAX
    down_ok = down_max >= STEP_DOWN_MIN
    print(f"B4.1 HOLD mean_y={hold_mean:.4f} {'PASS' if hold_ok else 'FAIL'}")
    print(f"B4.2 STEP_UP min_y={up_min:.4f} {'PASS' if up_ok else 'FAIL'}")
    print(f"B4.3 STEP_DOWN max_y={down_max:.4f} {'PASS' if down_ok else 'FAIL'}")
    if hold_ok and up_ok and down_ok:
        print("B4_CHEY_PROBE=PASS")
        return 0
    killer = "hold" if not hold_ok else ("step_up" if not up_ok else "step_down")
    print(f"B4_CHEY_PROBE=FAIL {killer}")
    return 1


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def ys(rows: list[dict[str, str]], t0: float, t1: float) -> list[float]:
    out = []
    for row in rows:
        t = float(row["t"])
        if t0 <= t <= t1:
            out.append(float(row["cheyP"]))
    return out


def mean_y(rows: list[dict[str, str]], t0: float, t1: float) -> float:
    vals = ys(rows, t0, t1)
    return sum(vals) / len(vals)


def min_y(rows: list[dict[str, str]], t0: float, t1: float) -> float:
    return min(ys(rows, t0, t1))


def max_y(rows: list[dict[str, str]], t0: float, t1: float) -> float:
    return max(ys(rows, t0, t1))


if __name__ == "__main__":
    raise SystemExit(main())
