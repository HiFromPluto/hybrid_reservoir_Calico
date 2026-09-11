#!/usr/bin/env python3
"""NARMA-blind D0b gate checker. Reads results/d0b_summary.json."""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("D0b checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "d0b_summary.json"
    if not path.exists():
        raise SystemExit("run_d0b.py first")
    summary = json.loads(path.read_text(encoding="utf-8"))
    print(f"primary={summary.get('primary')}")
    print(f"ics_label={summary.get('ics_label')}")
    print(f"D0b={summary.get('D0b')}")
    print(f"reason={summary.get('reason')}")
    print(f"d1_may_start={summary.get('d1_may_start')}")
    gates = summary.get("gates", {})
    print(
        f"identity_osc={gates.get('n_identity_osc')}/{gates.get('n_identity')} "
        f"occupancy_ok={gates.get('occupancy_ok')} "
        f"consecutive_increase={gates.get('consecutive_increase')} "
        f"spearman={gates.get('spearman')} "
        f"period_class_ok={gates.get('period_class_ok')}"
    )
    print(
        f"wrong_kick_ok={gates.get('wrong_kick_ok')} "
        f"wrong_kick_flag={gates.get('wrong_kick_flag')} "
        f"extras_off_ok={gates.get('extras_mu_ge_1p2_off_ok')} "
        f"tolerance_ok={gates.get('tolerance_ok')}"
    )
    print(f"osc_identity_mu={gates.get('osc_identity_mu')}")
    print(f"osc_identity_period={gates.get('osc_identity_period')}")
    extras = summary.get("finite_interval", {}).get("extra_flags", {})
    print(f"finite_interval_extras={extras}")
    status = summary.get("D0b")
    if status == "FAIL":
        print("D0b FAIL. Do not hunt ICs. Do not retune. Do not start D1. No NARMA.")
        sys.exit(0)
    if status == "FAIL_NO_IDENTITY":
        print(
            "D0b FAIL_NO_IDENTITY. Occupied but period does not increase with mu. "
            "Do not retune. Do not start D1. No NARMA."
        )
        sys.exit(0)
    if status == "SMOKE":
        print("Smoke only. Full grid not scored.")
        sys.exit(0)
    print("D0b PASS on predeclared primary AHL kick. D1 Java circuit parity may start.")


if __name__ == "__main__":
    main()
