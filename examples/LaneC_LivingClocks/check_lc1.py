#!/usr/bin/env python3
"""LC1_DEATH_VISIBLE checker. Death only. Does not rewrite LC0."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_DEATH.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_lc1.json"

N0 = 64
SEED = 101
GAMMA = 0.43
T_HORIZON_D = 4.0
DELTA = 0.20
VIS_MIN = 5
N_FULL_WINDOWS = 4


def main() -> int:
    if any(k in " ".join(sys.argv).lower() for k in ("narma", "charc", "ipc")):
        raise SystemExit("LC1 checker refuses NARMA/CHARC/IPC")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("LC1 PROTOCOL is not frozen")
    if js.get("status_label") != "LC1_DEATH_VISIBLE" or js.get("narma") is not False:
        raise SystemExit("LC1 must freeze extra and narma=false")
    if js.get("n0") != N0 or js.get("seed") != SEED:
        raise SystemExit("LC1 must freeze N0=64 seed=101")
    if float(js.get("gamma_per_d", 0.0)) != GAMMA:
        raise SystemExit("LC1 must freeze Schink gamma=0.43/d")
    if js.get("growth") is not False or js.get("paper1_toxic") is not False:
        raise SystemExit("LC1 must freeze growth=false and paper1_toxic=false")
    if js.get("n_full_windows") != N_FULL_WINDOWS or float(js.get("window_d", 0)) != 1.0:
        raise SystemExit("LC1 must freeze four 1 d windows")

    summary = json.loads((RESULTS / "lc1_summary.json").read_text(encoding="utf-8"))
    die = load(RESULTS / "lc1_die.csv")
    off = load(RESULTS / "lc1_off.csv")
    windows = load(RESULTS / "lc1_windows.csv")

    honesty = (
        summary.get("seed") == SEED
        and summary.get("narma") is False
        and summary.get("death") is True
        and summary.get("growth") is False
        and summary.get("motility") is False
        and summary.get("chemistry") is False
        and summary.get("gate") == "LC1_DEATH_VISIBLE"
        and abs(float(summary.get("gamma_per_d", 0.0)) - GAMMA) < 1e-12
    )
    immortal = (
        int(summary["off_deaths"]) == 0
        and int(summary["off_N_end"]) == N0
        and int(summary["off_N_ever"]) == N0
        and all(int(row["N"]) == N0 and int(row["N_ever"]) == N0 for row in off)
    )
    full = [
        int(row["deaths"])
        for row in windows
        if row["arm"] == "LC1_DIE" and row["full"] == "true"
    ]
    visibility = any(d >= VIS_MIN for d in full)
    clock = clock_ok(die)
    no_silent = (
        int(summary["die_N_ever"]) == N0
        and int(summary["off_N_ever"]) == N0
        and int(summary["die_deaths"]) == N0 - int(summary["die_N_end"])
        and int(summary["off_deaths"]) == 0
        and all(int(row["N_ever"]) == N0 for row in die)
    )

    print(f"LC1.1 honesty {'PASS' if honesty else 'FAIL'}")
    print(
        f"LC1.2 immortal deaths={summary['off_deaths']} "
        f"N_end={summary['off_N_end']} {'PASS' if immortal else 'FAIL'}"
    )
    max_full = max(full) if full else 0
    vis_word = "PASS" if visibility else "BOUND"
    print(f"LC1.3 visibility max_full_window_deaths={max_full} {vis_word}")
    print(f"LC1.4 clock delta={DELTA:.2f} {'PASS' if clock else 'FAIL'}")
    print(f"LC1.5 no_silent_kill {'PASS' if no_silent else 'FAIL'}")

    if honesty and immortal and visibility and clock and no_silent:
        print("LC1_DEATH_VISIBLE=PASS")
        return 0
    killer = (
        "honesty"
        if not honesty
        else (
            "immortal"
            if not immortal
            else (
                "visibility_bound"
                if not visibility
                else ("clock" if not clock else "silent_kill")
            )
        )
    )
    print(f"LC1_DEATH_VISIBLE=FAIL {killer}")
    return 1


def clock_ok(rows: list[dict[str, str]]) -> bool:
    if len(rows) != N_FULL_WINDOWS + 1:
        return False
    last_t = float(rows[-1]["t_d"])
    if abs(last_t - T_HORIZON_D) > 1e-9:
        return False
    for row in rows:
        t = float(row["t_d"])
        pred = N0 * math.exp(-GAMMA * t)
        n = int(row["N"])
        if pred <= 0:
            return False
        if abs(n / pred - 1.0) > DELTA + 1e-12:
            return False
    return True


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
