#!/usr/bin/env python3
"""NARMA-blind D0 gate checker. Reads results/d0_summary.json."""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("D0 checker cannot see a NARMA target.")


def main() -> None:
    refuse_narma()
    path = RESULTS / "d0_summary.json"
    if not path.exists():
        raise SystemExit("run_d0.py first")
    summary = json.loads(path.read_text(encoding="utf-8"))
    print(f"primary={summary.get('primary')}")
    print(f"ics_label={summary.get('ics_label')}")
    print(f"D0={summary.get('D0')}")
    trend = summary.get("trend", {})
    print(
        f"identity_osc={trend.get('n_identity_osc')}/{trend.get('n_identity')} "
        f"consecutive_increase={trend.get('consecutive_increase')} "
        f"spearman={trend.get('spearman')} "
        f"T(1)={trend.get('T_mu_1')} T(2)={trend.get('T_mu_2')}"
    )
    tol = summary.get("tolerance", {})
    print(f"tolerance_pass={tol.get('pass')} rel={tol.get('rel_period_change')}")
    extras = summary.get("finite_interval", {}).get("extra_flags", {})
    print(f"finite_interval_extras={extras}")
    if summary.get("D0") == "FAIL":
        print("D0 FAIL. Do not retune. Do not start D1. No NARMA.")
        sys.exit(0)
    if summary.get("D0") == "SMOKE":
        print("Smoke only. Full grid not scored.")
        sys.exit(0)
    print("D0 PASS on predeclared primary arm. Java circuit parity is a later gate.")


if __name__ == "__main__":
    main()
