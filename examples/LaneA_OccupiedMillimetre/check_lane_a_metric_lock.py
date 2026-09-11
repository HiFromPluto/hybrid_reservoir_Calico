#!/usr/bin/env python3
"""Lane A metric lock. Jaeger already locked. CHARC/IPC out of scope.

Reads existing NARMA CSVs and the Jaeger summary only. No Java. No new u.
A homemade CHARC-like number is a kill.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from check_lane_a_narma10 import load_u, sha256_u  # noqa: E402

RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_METRIC.md"
PROTOCOL_JSON = HERE / "configs" / "metric_lock_protocol.json"
U_FILE = HERE / "input_u_narma200.txt"
LANEA_JAVA = HERE.parents[1] / "src" / "bsim" / "lanea"
MC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_MEMORY_CAPACITY_STANDING.md"
NARMA_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_NARMA10_STANDING.md"
OCC_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_OCCUPIED_MILLIMETRE_STANDING.md"
CARRIER_STANDING = HERE.parents[1] / "examples" / "PocketDish" / "LANE_A_CARRIER_STANDING.md"

U_SHA = "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b"
NARMA10B_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
EXPECTED_MC = {"RL": 2.069, "FIELD": 2.865, "SILENT_RL": 0.0, "DELAY_U_10": 9.0}
MC_TOL = 0.001
FORBIDDEN_SUMMARY_KEYS = (
    "C_total",
    "C_2",
    "C_3",
    "kernel_rank",
    "generalisation_rank",
    "generalization_rank",
    "CHARC",
    "charc_quality",
    "charc_radar",
)


def refuse_forbidden() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("--charc", "--ipc", "erickson", "lane-d", "laned")):
        raise SystemExit(
            "LANE_A_METRIC_LOCK refuses CHARC/IPC flags and Lane D/Erickson"
        )


def require_protocol_frozen() -> dict:
    proto = PROTOCOL.read_text(encoding="utf-8")
    js_text = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js_text:
        raise SystemExit("METRIC_LOCK PROTOCOL not frozen")
    js = json.loads(js_text)
    if js.get("object") != "LANE_A_METRIC_LOCK":
        raise SystemExit("must freeze object LANE_A_METRIC_LOCK")
    if js.get("charc") is not False or js.get("dambre_ipc") is not False:
        raise SystemExit("must freeze charc/dambre_ipc false — CHARC/IPC are out of scope")
    if js.get("kernel_rank") is not False or js.get("generalisation_rank") is not False:
        raise SystemExit("must freeze KR/GR false")
    if js.get("homemade_charc_like_is_kill") is not True:
        raise SystemExit("homemade CHARC-like must be a kill")
    if js.get("jaeger_already_locked") is not True:
        raise SystemExit("must cite Jaeger as already locked")
    if js.get("new_u") is not False or js.get("rerun_java") is not False:
        raise SystemExit("must freeze new_u=false rerun_java=false")
    if js.get("raise_J_max") is not False or js.get("J_max") != 128000000.0:
        raise SystemExit("J_max must stay 1.28e8 and unraised")
    if js.get("ieee_draft_reopened") is not False:
        raise SystemExit("IEEE draft must stay closed")
    if js.get("hybriddish_0928_paste") is not False:
        raise SystemExit("must not paste HybridDish 0.928")
    if js.get("lane_d_gamma_paste") is not False:
        raise SystemExit("must not paste Lane D gamma")
    if js.get("u_sha256") != U_SHA:
        raise SystemExit("PROTOCOL json u hash mismatch")
    if js.get("narma_0895_is_not_mc") is not True or js.get("narma_0895_is_not_charc") is not True:
        raise SystemExit("must freeze NARMA 0.895 is not MC/CHARC")
    return js


def require_parents() -> None:
    occ = OCC_STANDING.read_text(encoding="utf-8")
    nar = NARMA_STANDING.read_text(encoding="utf-8")
    car = CARRIER_STANDING.read_text(encoding="utf-8")
    mc = MC_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in occ:
        raise SystemExit("occupancy standing must remain PASS")
    if "**Status: PASS**" not in nar or "0.895" not in nar:
        raise SystemExit("NARMA standing must remain system PASS with 0.895")
    if "**Status: FAIL**" not in car:
        raise SystemExit("carrier standing must remain FAIL vs field")
    if "**Status: PASS**" not in mc or "Not** CHARC" not in mc and "Not CHARC" not in mc:
        if "Not CHARC" not in mc and "**Not** CHARC" not in mc:
            raise SystemExit("MC standing must remain PASS and not CHARC")
    if "2.069" not in mc:
        raise SystemExit("MC standing must still carry Jaeger 2.069")


def refuse_math_random() -> None:
    hits = [p.name for p in sorted(LANEA_JAVA.glob("*.java")) if "Math.random()" in p.read_text(encoding="utf-8")]
    if hits:
        raise SystemExit(f"Math.random() forbidden in src/bsim/lanea: {hits}")


def write_summary(payload: dict) -> None:
    forbidden = [k for k in FORBIDDEN_SUMMARY_KEYS if payload.get(k) not in (None, False)]
    if forbidden:
        raise SystemExit(f"summary would emit forbidden CHARC/IPC keys: {forbidden}")
    if payload.get("charc") is not False or payload.get("dambre_ipc") is not False:
        raise SystemExit("summary must keep charc=false and dambre_ipc=false")
    (RESULTS / "lane_a_metric_lock_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def main() -> int:
    refuse_forbidden()
    js = require_protocol_frozen()
    require_parents()
    refuse_math_random()

    u = load_u(U_FILE)
    digest = sha256_u(u)
    if digest != U_SHA:
        raise SystemExit(f"u sha256 drifted: {digest}")
    if digest == NARMA10B_SHA:
        raise SystemExit("copied Narma10b u hash")

    driven = RESULTS / "java_LANE_A_NARMA10_DRIVEN.csv"
    silent = RESULTS / "java_LANE_A_NARMA10_SILENT.csv"
    jaeger_path = RESULTS / "lane_a_memory_capacity_summary.json"
    missing = [p.name for p in (driven, silent, jaeger_path) if not p.is_file()]
    if missing:
        summary = {
            "gate": "LaneA_METRIC_LOCK",
            "status_label": "LANE_A_METRIC_LOCK",
            "LANE_A_METRIC_LOCK": "SCOPE_NOTE",
            "object": "LANE_A_METRIC_LOCK",
            "charc": False,
            "dambre_ipc": False,
            "jaeger_2001": True,
            "missing": missing,
            "note": "required file missing; do not rerun Java",
        }
        write_summary(summary)
        print("LANE_A_METRIC_LOCK=SCOPE_NOTE missing " + ", ".join(missing) + ". Do not rerun Java.")
        return 0

    jaeger = json.loads(jaeger_path.read_text(encoding="utf-8"))
    if jaeger.get("charc") is not False or jaeger.get("dambre_ipc") is not False:
        raise SystemExit("Jaeger summary must keep charc/dambre_ipc false")
    if jaeger.get("LANE_A_MEMORY_CAPACITY") != "PASS":
        raise SystemExit("Jaeger parent must remain PASS")
    if jaeger.get("u_sha256") != U_SHA:
        raise SystemExit("Jaeger summary SHA drifted")
    if jaeger.get("narma_nrmse_is_not_mc") is not True:
        raise SystemExit("Jaeger summary must keep NARMA NRMSE is not MC")

    recorded = jaeger.get("MC_fading") or {}
    mismatches = {}
    for arm, expected in EXPECTED_MC.items():
        got = float(recorded[arm])
        if abs(got - expected) > MC_TOL:
            mismatches[arm] = {"expected": expected, "got": got}

    field_wins = float(recorded["FIELD"]) > float(recorded["RL"])
    if js.get("field_winning_delayed_u_is_fail") is True:
        raise SystemExit("protocol must allow field winning delayed-u")

    verdict = "FAIL" if mismatches else "PASS"
    summary = {
        "gate": "LaneA_METRIC_LOCK",
        "status_label": "LANE_A_METRIC_LOCK",
        "LANE_A_METRIC_LOCK": verdict,
        "object": "LANE_A_METRIC_LOCK",
        "identity": "ruler_lock",
        "not_a_paper": True,
        "not_a_living_computer": True,
        "frozen_before_traces": True,
        "u_sha256": digest,
        "occupancy_parent": "PASS",
        "narma_parent": "PASS",
        "mc_parent": "PASS",
        "carrier_parent": "FAIL",
        "jaeger_2001": True,
        "jaeger_already_locked": True,
        "dambre_ipc": False,
        "dambre_2012": False,
        "dambre_out_reason": js["dambre_out_reason"],
        "charc": False,
        "charc_citation": js["charc_citation"],
        "charc_out_reason": js["charc_out_reason"],
        "kernel_rank": False,
        "generalisation_rank": False,
        "paper1_ipc_kr_rewrite": False,
        "narma_0895_is_not_mc": True,
        "narma_0895_is_not_charc": True,
        "hybriddish_0928_paste": False,
        "lane_d_gamma_paste": False,
        "ieee_draft_reopened": False,
        "rerun_java": False,
        "new_u": False,
        "J_max": 128000000.0,
        "motility": "OFF",
        "growth": "OFF",
        "death": "OFF",
        "MC_fading": {arm: float(recorded[arm]) for arm in EXPECTED_MC},
        "field_wins_delayed_u": field_wins,
        "field_winning_delayed_u_is_fail": False,
        "delay_u_10_is_kill": False,
        "mc_mismatches": mismatches,
        "protocol_grade": "jaeger_delay_capacity_windowed_maps",
        "forbidden_to_report": [
            "CHARC",
            "KR",
            "GR",
            "Dambre_C_total",
            "Dambre_C2",
            "Dambre_C3",
            "NARMA_0.895_as_MC",
            "NARMA_0.895_as_CHARC",
            "Paper1_IPC_KR",
            "HybridDish_0.928",
            "LaneD_gamma",
        ],
    }
    write_summary(summary)

    print("LANE_A_METRIC_LOCK ruler. Not a paper. Not a living computer.")
    print("Jaeger 2001 delay capacity: already locked. CHARC: out of scope. Dambre IPC: out of scope.")
    print(f"u_sha256={digest} frozen_before_traces=true J_max=1.28e8 motility=OFF")
    print(
        f"  RL MC_fading={recorded['RL']:.6g} FIELD={recorded['FIELD']:.6g} "
        f"SILENT_RL={recorded['SILENT_RL']:.6g} DELAY_U_10={recorded['DELAY_U_10']:.6g}"
    )
    print(f"field_wins_delayed_u={field_wins} (allowed, not a fail)")
    print("CHARC (Dale 2019 KR/GR/MC sweep) cannot be matched on these window-mean maps.")
    print("Dambre 2012 full IPC cannot be matched (160 rows, 400-D, one u).")
    print(f"LANE_A_METRIC_LOCK={verdict}")
    if verdict == "PASS":
        print(
            "Protocol-grade: Jaeger delay capacity on frozen NARMA maps. "
            "Forbidden: CHARC, KR, GR, Dambre C_total, NARMA 0.895 as MC/CHARC."
        )
    else:
        print(f"Jaeger reprint drifted: {mismatches}. Do not raise J_max.")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
