#!/usr/bin/env python3
"""NARMA-blind Lane A occupancy checker. LANE_A_OCCUPIED_MILLIMETRE."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL.md"
PROTOCOL_JSON = HERE / "configs" / "protocol.json"

ALIVE_MIN = 0.05
LEDGER_REL = 1e-3
N0_REL = 1e-11
J_MAX = 128000000.0
DRIVEN_MASS = 38400000000.0
CLOSED_MASS = 128000000.0


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("narma", "charc", "fig4b", "time_adj")):
        raise SystemExit("Lane A checker refuses NARMA/CHARC/Fig4b/TIME_ADJ")


def load_csv(path: Path):
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if "LANE_A_OCCUPIED_MILLIMETRE" not in names:
            raise SystemExit(f"{path.name} missing LANE_A_OCCUPIED_MILLIMETRE")
        return list(reader)


def main() -> None:
    refuse_narma()
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js:
        raise SystemExit("Lane A PROTOCOL is not frozen_before_traces")
    if "LANE_A_OCCUPIED_MILLIMETRE" not in proto:
        raise SystemExit("PROTOCOL must name LANE_A_OCCUPIED_MILLIMETRE")
    if '"narma": false' not in js or '"J_max": 128000000.0' not in js:
        raise SystemExit("PROTOCOL must freeze narma=false and J_max")
    summary = json.loads((RESULTS / "lane_a_summary.json").read_text(encoding="utf-8"))
    if summary.get("narma") is not False:
        raise SystemExit("summary must keep narma false")
    if summary.get("hybriddish_overall_rewrite") is not False:
        raise SystemExit("summary must not rewrite HybridDish Overall")
    if summary.get("paper1_rewrite") is not False:
        raise SystemExit("summary must not rewrite paper 1")
    if summary.get("lane_b_started") is not False:
        raise SystemExit("Lane B must not be started from this standing")

    by_name = {a["name"]: a for a in summary["arms"]}
    required = [
        "LANE_A_CLOSED_CONSERVATIVE",
        "LANE_A_DRIVEN",
        "LANE_A_SILENT",
        "LANE_A_BROWNIAN",
    ]
    missing = [n for n in required if n not in by_name]
    if missing:
        raise SystemExit(f"missing arms: {missing}")

    closed = by_name["LANE_A_CLOSED_CONSERVATIVE"]
    driven = by_name["LANE_A_DRIVEN"]
    silent = by_name["LANE_A_SILENT"]
    brownian = by_name["LANE_A_BROWNIAN"]

    ok = True
    if abs(closed["cmd_mass"] - CLOSED_MASS) > 1e-3 or closed["n0_rel"] > N0_REL:
        print("FAIL closed ledger", closed)
        ok = False
    else:
        print(
            "CLOSED cmd", closed["cmd_mass"], "n0_rel", closed["n0_rel"],
            "LANE_A_OCCUPIED_MILLIMETRE",
        )

    if driven["occupancy"] != "ALIVE" or driven["mean_R"] < ALIVE_MIN:
        print("FAIL driven occupancy", driven["occupancy"], driven["mean_R"])
        ok = False
    elif driven["cmd_rel"] > LEDGER_REL or driven["n0_rel"] > N0_REL:
        print("FAIL driven ledger", driven)
        ok = False
    elif abs(driven["cmd_mass"] - DRIVEN_MASS) > 1.0:
        print("FAIL driven commanded mass", driven["cmd_mass"])
        ok = False
    else:
        print(
            "DRIVEN occupancy ALIVE mean_R", driven["mean_R"],
            "cmd", driven["cmd_mass"], "remain", driven["remaining"],
            "decay", driven["decay"], "LANE_A_OCCUPIED_MILLIMETRE",
        )

    if silent["mean_R"] >= ALIVE_MIN or not (driven["mean_R"] > silent["mean_R"]):
        print("FAIL silent mean_R", silent["mean_R"])
        ok = False
    else:
        print("SILENT mean_R", silent["mean_R"], "LANE_A_OCCUPIED_MILLIMETRE")

    if brownian["mean_R"] != 0.0 or brownian["occupancy"] != "NO_RECEIVER":
        print("FAIL brownian", brownian)
        ok = False
    else:
        print("BROWNIAN mean_R 0 (no Hill) LANE_A_OCCUPIED_MILLIMETRE")

    driven_csv = load_csv(RESULTS / "java_LANE_A_DRIVEN.csv")
    occ = [r for r in driven_csv if 4 <= int(r["window"]) <= 7]
    if len(occ) != 64:
        print("FAIL occupancy sample count", len(occ), "expected 64")
        ok = False

    if summary.get("LANE_A_OCCUPIED_MILLIMETRE") != ("PASS" if ok else "FAIL"):
        print("FAIL summary flag mismatch", summary.get("LANE_A_OCCUPIED_MILLIMETRE"), "ok", ok)
        ok = False

    if not ok:
        print(
            "LANE_A_OCCUPIED_MILLIMETRE FAIL. Do not raise J_max. "
            "Do not edit HybridDish Overall. Not Fig. 4b. Not C1c. Not paper 1."
        )
        raise SystemExit(1)
    print(
        "LANE_A_OCCUPIED_MILLIMETRE PASS occupancy-first N0 millimetre dish. "
        "HybridDish Overall unchanged. Paper 1 unchanged. Lane B not started. "
        "Not Fig. 4b. Not C1c. Not NARMA."
    )


if __name__ == "__main__":
    main()
