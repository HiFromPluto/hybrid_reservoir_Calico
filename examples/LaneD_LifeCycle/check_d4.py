#!/usr/bin/env python3
"""D4_THREE_PHASE checker. Sequential feast / handoff / famine. Not NARMA."""

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
PROTOCOL = HERE / "PROTOCOL_D4.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_d4.json"
D33 = ROOT / "examples" / "PocketDish" / "LANE_D_D33_PROTEOME_TABLE_STANDING.md"

N0 = 8
SEED = 101
C_S = 0.5
C_CUT = 1.0e-3
D0_T = 10455.0
D0_C = 9.977e-4
D0_N = 64
C_TOL = 5.0e-7
PHI_RB0 = 0.049
GAMMA_TR = 11.02
GAMMA0 = 0.21
DELTA = 1.0e-6

ISOLATION_PATHS = [
    "examples/LaneD_LifeCycle/PROTOCOL_D0.md",
    "examples/LaneD_LifeCycle/configs/protocol_d0.json",
    "examples/LaneD_LifeCycle/results/d0_grow.csv",
    "examples/LaneD_LifeCycle/results/d0_summary.json",
    "examples/PocketDish/LANE_D_D0_CLOSED_BATH_STANDING.md",
    "src/bsim/laned/LaneDD0Job.java",
    "src/bsim/laned/LaneDD0Identity.java",
    "src/bsim/laned/LaneDClosedBath.java",
    "src/bsim/laned/LaneDClosedBathCell.java",
    "examples/BacteriumFromScratch/NutrientField.java",
    "examples/BacteriumFromScratch/EcoliRodCell.java",
    "examples/BacteriumFromScratch/ChassisParameters.java",
    "examples/BacteriumFromScratch/StarvationViability.java",
    "examples/BacteriumFromScratch/JOB6_STANDING.md",
    "examples/BacteriumFromScratch/JOB2_STANDING.md",
    "examples/BacteriumFromScratch/JOB3_STANDING.md",
    "examples/BacteriumFromScratch/JOB3B_STANDING.md",
    "examples/BacteriumFromScratch/JOB3C_STANDING.md",
    "examples/LaneC_LivingClocks",
    "src/bsim/lanec",
    "src/bsim/lanea",
    "src/bsim/laneb",
    "examples/PocketDish/manuscript_ieee",
    "examples/PocketDish/LANE_C_LC1_DEATH_STANDING.md",
]


def main() -> int:
    if any(k in " ".join(sys.argv).lower() for k in ("narma", "charc", "ipc")):
        raise SystemExit("D4 checker refuses NARMA/CHARC/IPC")
    d33 = D33.read_text(encoding="utf-8")
    if "Status: PORT_HOLDS" not in d33 or "lifts D3.1 `T_OVERLAY`" not in d33:
        print("D4_THREE_PHASE=BLOCKED D3.3 does not lift T_OVERLAY")
        return 3

    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        return fail("honesty", "D4 PROTOCOL is not frozen")
    if js.get("status_label") != "D4_THREE_PHASE" or js.get("phi_x") != "phi_Rb":
        return fail("honesty", "must freeze D4_THREE_PHASE and phi_Rb")
    if js.get("erickson_k_m_in_monod") is not False:
        return fail("honesty", "K_M must stay out of Monod")
    if js.get("lc1b_coin_flip") is not False:
        return fail("honesty", "famine is not LC1b")

    summary = json.loads((RESULTS / "d4_summary.json").read_text(encoding="utf-8"))
    feast = load(RESULTS / "d4_seq_feast.csv")
    famine = load(RESULTS / "d4_seq_famine.csv")
    off = load(RESULTS / "d4_off.csv")

    honesty = (
        summary.get("device") == "LANE_D_CLOSED_BATH"
        and summary.get("phi_x") == "phi_Rb"
        and summary.get("narma") is False
        and summary.get("seed") == SEED
    )
    process_off = (
        int(summary["off_births"]) == 0
        and int(summary["off_N_end"]) == N0
        and abs(float(summary["off_C_end"]) - C_S) <= 1e-15
        and all(int(r["N"]) == N0 and float(r["deaths"]) == 0.0 for r in off)
    )
    no_conc = all(r["death_on"] == "false" for r in feast)
    fingerprint = (
        abs(float(summary["feast_t_end"]) - D0_T) <= 1e-6
        and abs(float(summary["feast_C_end"]) - D0_C) <= C_TOL
        and int(summary["feast_N_end"]) == D0_N
    )
    phi_rb = float(summary["phi_Rb_freeze"])
    mu_eff = float(summary["mu_eff"])
    mu_id = abs(mu_eff - GAMMA_TR * (phi_rb - PHI_RB0)) <= 1e-10
    gamma = float(summary["gamma_per_d"])
    gamma_id = abs(gamma - GAMMA0 * math.exp(mu_eff)) <= 1e-12
    n_h = int(summary["n_handoff"])
    famine_ok = True
    for row in famine:
        t = float(row["t_d"])
        n = float(row["N"])
        pred = n_h * math.exp(-gamma * t)
        if n_h <= 0 or abs(n - pred) / n_h > DELTA + 1e-15:
            famine_ok = False
    no_kill = int(summary["feast_N_end"]) == int(summary["feast_N_ever"])
    cap = bool(summary.get("cap_fired"))
    t_cap = bool(summary.get("t_cap_scope"))
    isolation = isolation_ok()

    print(
        "D4_THREE_PHASE death=OFF_while_C>C_cut device=LANE_D_CLOSED_BATH "
        f"phi_X=phi_Rb seed={SEED}"
    )
    print(f"D4.1 honesty {'PASS' if honesty else 'FAIL'}")
    print(f"D4.2 process-off {'PASS' if process_off else 'FAIL'}")
    print(f"D4.3 no_concurrency {'PASS' if no_conc else 'FAIL'}")
    print(
        f"D4.4 feast_fingerprint t={summary['feast_t_end']} "
        f"C={summary['feast_C_end']} N={summary['feast_N_end']} "
        f"{'PASS' if fingerprint else 'FAIL'}"
    )
    print(
        f"D4.5 phi_X=phi_Rb mu_eff={mu_eff:.6f} gamma={gamma:.6f} "
        f"{'PASS' if mu_id and gamma_id else 'FAIL'}"
    )
    print(f"D4.6 famine_clock N_h={n_h} {'PASS' if famine_ok else 'FAIL'}")
    print(
        f"D4.7 no_silent_kill "
        f"{'PASS' if no_kill and not cap else ('SCOPE_NOTE' if cap else 'FAIL')}"
    )
    print(f"D4.8 isolation {'PASS' if isolation else 'FAIL'}")

    if cap or t_cap:
        print(f"D4_THREE_PHASE=SCOPE_NOTE {'compute_cap' if cap else 't_cap'}")
        return 2
    if honesty and process_off and no_conc and fingerprint and mu_id and gamma_id and famine_ok and no_kill and isolation:
        print("D4_THREE_PHASE=PASS")
        return 0
    killer = (
        "honesty"
        if not honesty
        else (
            "process_off"
            if not process_off
            else (
                "concurrency"
                if not no_conc
                else (
                    "feast_fingerprint"
                    if not fingerprint
                    else (
                        "phi_x_identity"
                        if not (mu_id and gamma_id)
                        else (
                            "famine_clock"
                            if not famine_ok
                            else ("silent_kill" if not no_kill else "isolation")
                        )
                    )
                )
            )
        )
    )
    return fail(killer, "Do not raise C_s or put K_M in Monod.")


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


def fail(killer: str, detail: str) -> int:
    print(f"D4_THREE_PHASE=FAIL {killer}")
    print(detail)
    return 1


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
