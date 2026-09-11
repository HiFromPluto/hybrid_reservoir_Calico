#!/usr/bin/env python3
"""D5_GATES checker. Regression / memory / optimum / refuse. Not NARMA."""

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
PROTOCOL = HERE / "PROTOCOL_D5.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_d5.json"
D4_STANDING = ROOT / "examples" / "PocketDish" / "LANE_D_D4_THREE_PHASE_STANDING.md"

ISOLATION_PATHS = [
    "examples/BacteriumFromScratch/JOB2_STANDING.md",
    "examples/BacteriumFromScratch/JOB3_STANDING.md",
    "examples/BacteriumFromScratch/JOB3B_STANDING.md",
    "examples/BacteriumFromScratch/JOB3C_STANDING.md",
    "examples/BacteriumFromScratch/JOB6_STANDING.md",
    "examples/BacteriumFromScratch/StarvationViability.java",
    "examples/BacteriumFromScratch/Job6Sims.java",
    "examples/BacteriumFromScratch/ChassisParameters.java",
    "examples/BacteriumFromScratch/EcoliRodCell.java",
    "examples/BacteriumFromScratch/NutrientField.java",
    "examples/LaneD_LifeCycle/PROTOCOL_D0.md",
    "examples/LaneD_LifeCycle/results/d0_grow.csv",
    "src/bsim/laned/LaneDD0Job.java",
    "src/bsim/laned/LaneDClosedBath.java",
    "src/bsim/laned/LaneDD4Job.java",
    "examples/LaneC_LivingClocks",
    "src/bsim/lanec",
    "examples/PocketDish/LANE_C_LC0_REPLICATION_STANDING.md",
    "examples/PocketDish/LANE_C_LC0_TWO_GEN_STANDING.md",
    "examples/PocketDish/LANE_C_LC1_DEATH_STANDING.md",
    "examples/PocketDish/LANE_C_LC1B_DEATH_BINOMIAL_STANDING.md",
    "src/bsim/lanea",
    "src/bsim/laneb",
    "examples/PocketDish/manuscript_ieee",
]

STANDINGS_PASS = [
    ROOT / "examples/BacteriumFromScratch/JOB2_STANDING.md",
    ROOT / "examples/BacteriumFromScratch/JOB3_STANDING.md",
    ROOT / "examples/BacteriumFromScratch/JOB3B_STANDING.md",
    ROOT / "examples/BacteriumFromScratch/JOB3C_STANDING.md",
    ROOT / "examples/PocketDish/LANE_C_LC0_REPLICATION_STANDING.md",
    ROOT / "examples/PocketDish/LANE_C_LC0_TWO_GEN_STANDING.md",
    ROOT / "examples/PocketDish/LANE_C_LC1B_DEATH_BINOMIAL_STANDING.md",
]

CHECKERS = [
    (ROOT / "examples/BacteriumFromScratch/check_job2.py",
     ROOT / "examples/BacteriumFromScratch/results/job2_seed101/job2_timeseries.csv",
     ROOT / "examples/BacteriumFromScratch"),
    (ROOT / "examples/BacteriumFromScratch/check_job3.py",
     ROOT / "examples/BacteriumFromScratch/results/job3_seed101/job3_timeseries.csv",
     ROOT / "examples/BacteriumFromScratch"),
    (ROOT / "examples/LaneC_LivingClocks/check_lc0.py",
     ROOT / "examples/LaneC_LivingClocks/results/lc0_summary.json",
     ROOT / "examples/LaneC_LivingClocks"),
    (ROOT / "examples/LaneC_LivingClocks/check_lc0_two_gen.py",
     ROOT / "examples/LaneC_LivingClocks/results/lc0_two_gen_summary.json",
     ROOT / "examples/LaneC_LivingClocks"),
    (ROOT / "examples/LaneC_LivingClocks/check_lc1b.py",
     ROOT / "examples/LaneC_LivingClocks/results/lc1b_summary.json",
     ROOT / "examples/LaneC_LivingClocks"),
]


def isolation_ok() -> bool:
    existing = [p for p in ISOLATION_PATHS if (ROOT / p).exists()]
    proc = subprocess.run(
        ["git", "diff", "HEAD", "--exit-code", "--", *existing],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        return False
    return True


def standings_ok() -> bool:
    for path in STANDINGS_PASS:
        text = path.read_text(encoding="utf-8")
        if "PASS" not in text:
            return False
    lc1 = (ROOT / "examples/PocketDish/LANE_C_LC1_DEATH_STANDING.md").read_text(
        encoding="utf-8"
    )
    return "**Status: FAIL**" in lc1


def run_optional_checkers() -> tuple[bool, bool]:
    """Returns (all_present_passed, any_missing)."""
    missing = False
    for script, csv, cwd in CHECKERS:
        if not script.is_file():
            continue
        if not csv.is_file():
            missing = True
            continue
        proc = subprocess.run(
            [sys.executable, str(script)],
            cwd=cwd,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout)
            sys.stderr.write(proc.stderr)
            return False, missing
    return True, missing


def main() -> int:
    if any(k in " ".join(sys.argv).lower() for k in ("narma", "charc", "ipc")):
        raise SystemExit("D5 checker refuses NARMA/CHARC/IPC")

    proto = PROTOCOL.read_text(encoding="utf-8")
    spec = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    d4 = D4_STANDING.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or spec.get("frozen_before_traces") is not True:
        return fail("honesty", "D5 PROTOCOL is not frozen")
    if spec.get("status_label") != "D5_GATES" or spec.get("phi_x") != "phi_Rb":
        return fail("honesty", "must freeze D5_GATES and phi_Rb")
    if spec.get("invent_ev1_rows") is not False or spec.get("lc1_restage") is not False:
        return fail("honesty", "must not invent EV1 or restage LC1")
    if "**Status: PASS**" not in d4 or "D4_THREE_PHASE" not in d4:
        return fail("honesty", "D5 requires D4 PASS")

    summary = json.loads((RESULTS / "d5_summary.json").read_text(encoding="utf-8"))
    memory = load(RESULTS / "d5_memory.csv")
    feast = load(RESULTS / "d4_seq_feast.csv")

    honesty = (
        summary.get("device") == "LANE_D_CLOSED_BATH"
        and summary.get("phi_x") == "phi_Rb"
        and summary.get("narma") is False
        and summary.get("invent_ev1_rows") is False
    )

    iso = isolation_ok()
    stands = standings_ok()
    checks_ok, checks_missing = run_optional_checkers()
    regress = iso and stands and checks_ok

    delta = float(spec["delta_gamma_rel"])
    mem_ok = True
    for row in memory:
        mu = float(row["mu_per_h"])
        gamma = float(row["gamma_law"])
        meas = float(row["gamma_meas"])
        phi = float(row["phi_Rb"])
        phi_id = spec["phi_rb0"] + mu / spec["gamma_tr_per_h"]
        if abs(phi - phi_id) > 1e-12:
            mem_ok = False
        law = spec["gamma0_per_d"] * math.exp(spec["death_slope_h"] * mu)
        if abs(gamma - law) > 1e-12:
            mem_ok = False
        if abs(meas - gamma) / max(gamma, 1e-15) > delta:
            mem_ok = False

    ev1_ok = True
    typo_ok = True
    for mu, obs, sd in zip(
        spec["ev1_mu_per_h"], spec["ev1_gamma_per_d"], spec["ev1_sd_per_d"]
    ):
        pred = spec["gamma0_per_d"] * math.exp(spec["death_slope_h"] * mu)
        if abs(pred - obs) > sd + 1e-12:
            ev1_ok = False
        if (24.0 * pred) / obs <= 10.0:
            typo_ok = False

    mu_star = float(summary["mu_star"])
    mu_num = float(summary["mu_star_numeric"])
    g_star = float(summary["gamma_star"])
    mu_old = float(summary["mu_old"])
    g_old = float(summary["gamma_old"])
    opt_ok = (
        abs(mu_star - spec["mu_star_locked"]) <= spec["delta_mu_star"]
        and abs(mu_num - mu_star) <= spec["delta_mu_star"]
        and 0.48 <= g_star <= 0.52
        and not (0.84 <= mu_old <= 0.89 or 0.48 <= g_old <= 0.52)
    )

    refuse = (
        bool(summary.get("refuse_replete_threw"))
        and bool(summary.get("refuse_cut_allowed"))
        and all(r["death_on"] == "false" for r in feast)
    )

    print(
        "D5_GATES death=OFF_while_C>C_cut device=LANE_D_CLOSED_BATH "
        "phi_X=phi_Rb invent_ev1=false"
    )
    print(
        f"D5.1 regression git={iso} standings={stands} "
        f"checkers={'PASS' if checks_ok else 'FAIL'}"
        f"{' SCOPE_NOTE_csv_missing' if checks_missing else ''} "
        f"{'PASS' if regress else 'FAIL'}"
    )
    print(f"D5.2 memory EV1={ev1_ok} typo_guard={typo_ok} sweep={mem_ok} "
          f"{'PASS' if mem_ok and ev1_ok and typo_ok else 'FAIL'}")
    print(
        f"D5.3 optimum mu*={mu_star:.6f} gamma={g_star:.6f} "
        f"old={mu_old:.4f}/{g_old:.4f} {'PASS' if opt_ok else 'FAIL'}"
    )
    print(f"D5.4 refuse {'PASS' if refuse else 'FAIL'}")

    if honesty and regress and mem_ok and ev1_ok and typo_ok and opt_ok and refuse:
        print("D5_GATES=PASS")
        return 0
    killer = (
        "honesty" if not honesty
        else "regression" if not regress
        else "memory" if not (mem_ok and ev1_ok and typo_ok)
        else "optimum" if not opt_ok
        else "refuse"
    )
    return fail(killer, "Do not retune gamma0, EV1 SDs, or C_cut.")


def fail(killer: str, detail: str) -> int:
    print(f"D5_GATES=FAIL {killer}")
    print(detail)
    return 1


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
