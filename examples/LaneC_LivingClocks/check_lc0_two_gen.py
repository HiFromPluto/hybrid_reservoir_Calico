#!/usr/bin/env python3
"""LC0_TWO_GENERATION checker. Two generations. Does not rewrite LC0."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_TWO_GEN.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_lc0_two_gen.json"

N0 = 64
SEED = 101
T_DIV = math.log(2.0) / ((1.0 / 3600.0) * (0.5 / 0.52))
T_HORIZON = 2.0 * T_DIV
DELTA = 0.20
VIS_MIN = 5
N_FULL_WINDOWS = 17


def main() -> int:
    if any(k in " ".join(sys.argv).lower() for k in ("narma", "charc", "ipc")):
        raise SystemExit("LC0 two-gen checker refuses NARMA/CHARC/IPC")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("LC0 two-gen PROTOCOL is not frozen")
    if js.get("status_label") != "LC0_TWO_GENERATION" or js.get("narma") is not False:
        raise SystemExit("LC0 two-gen must freeze extra and narma=false")
    if js.get("n0") != N0 or js.get("seed") != SEED:
        raise SystemExit("LC0 two-gen must freeze N0=64 seed=101")
    if js.get("birth_rate_law") != "ln2/T_div" or js.get("death") is not False:
        raise SystemExit("LC0 two-gen must freeze ln2/T and death=false")
    if js.get("n_full_windows") != N_FULL_WINDOWS:
        raise SystemExit("LC0 two-gen must freeze 17 full windows")
    if abs(float(js.get("t_horizon_s", 0.0)) - T_HORIZON) > 1e-6:
        raise SystemExit("LC0 two-gen must freeze T=2 T_div")

    summary = json.loads((RESULTS / "lc0_two_gen_summary.json").read_text(encoding="utf-8"))
    grow = load(RESULTS / "lc0_two_gen_grow.csv")
    off = load(RESULTS / "lc0_two_gen_off.csv")
    windows = load(RESULTS / "lc0_two_gen_windows.csv")

    honesty = (
        summary.get("seed") == SEED
        and summary.get("narma") is False
        and summary.get("death") is False
        and summary.get("motility") is False
        and summary.get("chemistry") is False
        and summary.get("gate") == "LC0_TWO_GENERATION"
    )
    process_off = (
        int(summary["off_births"]) == 0
        and int(summary["off_N_end"]) == N0
        and int(summary["off_N_ever"]) == N0
        and all(int(row["N"]) == N0 and int(row["N_ever"]) == N0 for row in off)
    )
    full = [
        int(row["births"])
        for row in windows
        if row["arm"] == "LC0_TWO_GROW" and row["full"] == "true"
    ]
    visibility = any(b >= VIS_MIN for b in full)
    clock = clock_ok(grow)
    no_kill = (
        int(summary["grow_N_end"]) == int(summary["grow_N_ever"])
        and int(summary["off_N_end"]) == int(summary["off_N_ever"])
    )
    cap = bool(summary.get("cap_fired"))

    print(f"LC0_TWO.1 honesty {'PASS' if honesty else 'FAIL'}")
    print(
        f"LC0_TWO.2 process-off births={summary['off_births']} "
        f"N_end={summary['off_N_end']} {'PASS' if process_off else 'FAIL'}"
    )
    max_full = max(full) if full else 0
    vis_word = "PASS" if visibility else ("SCOPE_NOTE" if cap else "BOUND")
    print(f"LC0_TWO.3 visibility max_full_window_births={max_full} {vis_word}")
    print(f"LC0_TWO.4 clock delta={DELTA:.2f} {'PASS' if clock else 'FAIL'}")
    print(
        f"LC0_TWO.5 no_silent_kill "
        f"{'PASS' if no_kill and not cap else ('SCOPE_NOTE' if cap else 'FAIL')}"
    )

    if cap:
        print("LC0_TWO_GENERATION=SCOPE_NOTE compute_cap")
        return 2
    if honesty and process_off and visibility and clock and no_kill:
        print("LC0_TWO_GENERATION=PASS")
        return 0
    killer = (
        "honesty"
        if not honesty
        else (
            "process_off"
            if not process_off
            else (
                "visibility_bound"
                if not visibility
                else ("clock" if not clock else "silent_kill")
            )
        )
    )
    print(f"LC0_TWO_GENERATION=FAIL {killer}")
    return 1


def clock_ok(rows: list[dict[str, str]]) -> bool:
    if not rows:
        return False
    last_t = float(rows[-1]["t"])
    if abs(last_t - T_HORIZON) > 1e-6:
        return False
    for row in rows:
        t = float(row["t"])
        pred = N0 * (2.0 ** (t / T_DIV))
        n = int(row["N"])
        ever = int(row["N_ever"])
        if pred <= 0:
            return False
        if abs(n / pred - 1.0) > DELTA + 1e-12:
            return False
        if abs(ever / pred - 1.0) > DELTA + 1e-12:
            return False
    return True


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
