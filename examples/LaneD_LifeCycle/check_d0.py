#!/usr/bin/env python3
"""D0_CLOSED_BATH checker. Scalar closed bath. Not NARMA. Not Erickson."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_D0.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_d0.json"

N0 = 8
SEED = 101
C_S = 0.5
K_S = 0.02
LAMBDA_S = 1.0
C_CUT = 1.0e-3
DELTA_MASS = 0.01
CAP = 256
ELL0 = 1.0

ISOLATION_PATHS = [
    "examples/BacteriumFromScratch/NutrientField.java",
    "examples/BacteriumFromScratch/EcoliRodCell.java",
    "examples/BacteriumFromScratch/ChassisParameters.java",
    "examples/BacteriumFromScratch/JOB2_STANDING.md",
    "examples/BacteriumFromScratch/JOB3_STANDING.md",
    "examples/BacteriumFromScratch/JOB3B_STANDING.md",
    "examples/BacteriumFromScratch/JOB3C_STANDING.md",
    "examples/BacteriumFromScratch/check_job2.py",
    "examples/BacteriumFromScratch/check_job3.py",
    "examples/BacteriumFromScratch/check_job3b.py",
    "examples/BacteriumFromScratch/check_job3c.py",
    "examples/BacteriumFromScratch/Job3Sims.java",
    "examples/BacteriumFromScratch/Job3bSims.java",
    "examples/BacteriumFromScratch/Job3cSims.java",
    "examples/LaneC_LivingClocks",
    "src/bsim/lanec",
    "src/bsim/lanea",
    "src/bsim/laneb",
    "examples/PocketDish/manuscript_ieee",
    "examples/PocketDish/LANE_A_OCCUPIED_MILLIMETRE_STANDING.md",
    "examples/PocketDish/LANE_B_PAINT_HOLD_STANDING.md",
    "examples/PocketDish/LANE_C_LC0_REPLICATION_STANDING.md",
    "examples/PocketDish/LANE_C_LC1_DEATH_STANDING.md",
    "examples/PocketDish/LANE_C_LC1B_DEATH_BINOMIAL_STANDING.md",
]


def main() -> int:
    if any(k in " ".join(sys.argv).lower() for k in ("narma", "charc", "ipc")):
        raise SystemExit("D0 checker refuses NARMA/CHARC/IPC")
    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("D0 PROTOCOL is not frozen")
    if js.get("status_label") != "D0_CLOSED_BATH" or js.get("device") != "LANE_D_CLOSED_BATH":
        raise SystemExit("D0 must freeze D0_CLOSED_BATH / LANE_D_CLOSED_BATH")
    if js.get("not_device") != "LC0_OPEN_BATH":
        raise SystemExit("D0 must not adopt LC0_OPEN_BATH")
    if js.get("n0") != N0 or js.get("seed") != SEED:
        raise SystemExit("D0 must freeze N0=8 seed=101")
    if js.get("death") is not False or js.get("erickson") is not False:
        raise SystemExit("D0 must freeze death=false erickson=false")
    if js.get("sink_density") != "rho_cell" or js.get("elongation") != "incremental":
        raise SystemExit("D0 must freeze rho_cell incremental")
    if js.get("c_cut_mM") != C_CUT or js.get("k_s_mM") != K_S:
        raise SystemExit("D0 must freeze C_cut=0.001 Warren K_S=0.02")

    summary = json.loads((RESULTS / "d0_summary.json").read_text(encoding="utf-8"))
    grow = load(RESULTS / "d0_grow.csv")
    off = load(RESULTS / "d0_off.csv")

    honesty = (
        summary.get("seed") == SEED
        and summary.get("device") == "LANE_D_CLOSED_BATH"
        and summary.get("narma") is False
        and summary.get("death") is False
        and summary.get("motility") is False
        and summary.get("chemistry") is False
        and summary.get("pde") is False
        and summary.get("erickson") is False
        and summary.get("sink_density") == "rho_cell"
        and summary.get("elongation") == "incremental"
    )
    process_off = (
        int(summary["off_births"]) == 0
        and int(summary["off_N_end"]) == N0
        and abs(float(summary["off_C_end"]) - C_S) <= 1e-15
        and all(int(row["N"]) == N0 and abs(float(row["C_mM"]) - C_S) <= 1e-15 for row in off)
    )
    no_shrink = (
        int(summary["grow_shrink_events"]) == 0
        and all(float(row["min_ell"]) + 1e-12 >= ELL0 for row in grow)
        and not any_length_drop(grow)
    )
    residual = float(summary["grow_mass_residual"])
    ledger = residual <= DELTA_MASS + 1e-15
    monod = monod_ok(grow)
    c_end = float(summary["grow_C_end"])
    exhausted = c_end <= C_CUT + 1e-15
    cap = bool(summary.get("cap_fired"))
    t_cap = bool(summary.get("t_cap_scope"))
    no_kill = (
        int(summary["grow_N_end"]) == int(summary["grow_N_ever"])
        and int(summary["off_N_end"]) == int(summary["off_N_ever"])
    )
    isolation = isolation_ok()

    print(f"D0.1 honesty {'PASS' if honesty else 'FAIL'}")
    print(
        f"D0.2 process-off births={summary['off_births']} "
        f"N_end={summary['off_N_end']} C_end={summary['off_C_end']} "
        f"{'PASS' if process_off else 'FAIL'}"
    )
    print(f"D0.3 no_shrink events={summary['grow_shrink_events']} {'PASS' if no_shrink else 'FAIL'}")
    print(
        f"D0.4 bath_ledger residual={residual:.6e} delta={DELTA_MASS:.2f} "
        f"{'PASS' if ledger else 'FAIL'}"
    )
    print(f"D0.5 monod_identity Warren_K_S={K_S} {'PASS' if monod else 'FAIL'}")
    exh_word = "PASS" if exhausted else ("SCOPE_NOTE" if t_cap else "FAIL")
    print(f"D0.6 exhaustion C_end={c_end:.8f} stop={summary['grow_stop']} {exh_word}")
    print(
        f"D0.7 no_silent_kill "
        f"{'PASS' if no_kill and not cap else ('SCOPE_NOTE' if cap else 'FAIL')}"
    )
    print(f"D0.8 isolation {'PASS' if isolation else 'FAIL'}")
    print(
        f"F3_DIVISION_VOLUME_DROP fissions={summary['grow_births']} "
        f"drop_cum={summary['grow_f3_drop_cum']} (reported, not gated)"
    )
    print(
        f"REPORT sum_rho_cell_V={summary['grow_rho_cell_sum_V']} gCDW "
        "(not the mass gate)"
    )

    if cap or t_cap:
        print(f"D0_CLOSED_BATH=SCOPE_NOTE {'compute_cap' if cap else 't_cap'}")
        return 2
    if honesty and process_off and no_shrink and ledger and monod and exhausted and no_kill and isolation:
        print("D0_CLOSED_BATH=PASS")
        return 0
    killer = (
        "honesty"
        if not honesty
        else (
            "process_off"
            if not process_off
            else (
                "no_shrink_F4"
                if not no_shrink
                else (
                    "bath_ledger"
                    if not ledger
                    else (
                        "monod_identity"
                        if not monod
                        else (
                            "exhaustion"
                            if not exhausted
                            else ("silent_kill" if not no_kill else "isolation")
                        )
                    )
                )
            )
        )
    )
    print(f"D0_CLOSED_BATH=FAIL {killer}")
    if killer == "bath_ledger":
        print("STOP. Do not retune Y, rho_cell, C_s, or N0.")
    return 1


def monod_ok(rows: list[dict[str, str]]) -> bool:
    if not rows:
        return False
    for row in rows:
        c = float(row["C_mM"])
        printed = float(row["lambda_per_h"])
        expect = LAMBDA_S * (max(c, 0.0) / (max(c, 0.0) + K_S))
        if abs(printed - expect) > 1e-12:
            return False
    return True


def any_length_drop(rows: list[dict[str, str]]) -> bool:
    prev_min = None
    for row in rows:
        min_ell = float(row["min_ell"])
        if prev_min is not None and min_ell + 1e-12 < ELL0:
            return True
        prev_min = min_ell
    return False


def isolation_ok() -> bool:
    existing = [p for p in ISOLATION_PATHS if (ROOT / p).exists()]
    if not existing:
        return False
    try:
        proc = subprocess.run(
            ["git", "diff", "HEAD", "--exit-code", "--", *existing],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return False
    if proc.returncode == 0:
        return True
    sys.stderr.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    return False


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
