#!/usr/bin/env python3
"""LC1B_DEATH_BINOMIAL checker. Does not rewrite LC1. Does not replay seed 101."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_DEATH_BINOMIAL.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_lc1b.json"

N0 = 64
SEED = 303
GAMMA = 0.43
T_HORIZON_D = 4.0
K = 2.0
VIS_MIN = 5
N_FULL_WINDOWS = 4


def main() -> int:
    if any(tok in " ".join(sys.argv).lower() for tok in ("narma", "charc", "ipc")):
        raise SystemExit("LC1b checker refuses NARMA/CHARC/IPC")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("LC1b PROTOCOL is not frozen")
    if js.get("status_label") != "LC1B_DEATH_BINOMIAL" or js.get("narma") is not False:
        raise SystemExit("LC1b must freeze extra and narma=false")
    if js.get("n0") != N0 or js.get("seed") != SEED:
        raise SystemExit("LC1b must freeze N0=64 seed=303")
    if js.get("seed") == 101:
        raise SystemExit("LC1b must not replay seed 101")
    if float(js.get("gamma_per_d", 0.0)) != GAMMA:
        raise SystemExit("LC1b must freeze Schink gamma=0.43/d")
    if js.get("clock_law") != "binomial_2sigma" or float(js.get("clock_k", 0.0)) != K:
        raise SystemExit("LC1b must freeze binomial_2sigma k=2")
    if js.get("growth") is not False or js.get("paper1_toxic") is not False:
        raise SystemExit("LC1b must freeze growth=false and paper1_toxic=false")

    summary = json.loads((RESULTS / "lc1b_summary.json").read_text(encoding="utf-8"))
    die = load(RESULTS / "lc1b_die.csv")
    off = load(RESULTS / "lc1b_off.csv")
    windows = load(RESULTS / "lc1b_windows.csv")

    honesty = (
        summary.get("seed") == SEED
        and summary.get("narma") is False
        and summary.get("death") is True
        and summary.get("growth") is False
        and summary.get("motility") is False
        and summary.get("chemistry") is False
        and summary.get("gate") == "LC1B_DEATH_BINOMIAL"
        and summary.get("clock_law") == "binomial_2sigma"
        and abs(float(summary.get("clock_k", 0.0)) - K) < 1e-12
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
        if row["arm"] == "LC1B_DIE" and row["full"] == "true"
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

    print(f"LC1B.1 honesty {'PASS' if honesty else 'FAIL'}")
    print(
        f"LC1B.2 immortal deaths={summary['off_deaths']} "
        f"N_end={summary['off_N_end']} {'PASS' if immortal else 'FAIL'}"
    )
    max_full = max(full) if full else 0
    vis_word = "PASS" if visibility else "BOUND"
    print(f"LC1B.3 visibility max_full_window_deaths={max_full} {vis_word}")
    print(f"LC1B.4 clock binomial k={K:.0f} {'PASS' if clock else 'FAIL'}")
    print(f"LC1B.5 no_silent_kill {'PASS' if no_silent else 'FAIL'}")

    if honesty and immortal and visibility and clock and no_silent:
        print("LC1B_DEATH_BINOMIAL=PASS")
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
    print(f"LC1B_DEATH_BINOMIAL=FAIL {killer}")
    return 1


def mu_sigma(t: float) -> tuple[float, float]:
    mu = N0 * math.exp(-GAMMA * t)
    if mu <= 0.0 or mu >= N0:
        return mu, 0.0
    return mu, math.sqrt(mu * (1.0 - mu / N0))


def clock_ok(rows: list[dict[str, str]]) -> bool:
    if len(rows) != N_FULL_WINDOWS + 1:
        return False
    last_t = float(rows[-1]["t_d"])
    if abs(last_t - T_HORIZON_D) > 1e-9:
        return False
    for row in rows:
        t = float(row["t_d"])
        n = int(row["N"])
        mu, sig = mu_sigma(t)
        if sig <= 0.0:
            if n != N0:
                return False
            continue
        if abs(n - mu) > K * sig + 1e-12:
            return False
    return True


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
