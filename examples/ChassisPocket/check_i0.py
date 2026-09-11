#!/usr/bin/env python3
"""Evaluate ChassisPocket I0 packed-monolayer fill gates.

Thresholds are frozen in PROTOCOL.md before the run existed. Do not retune a
threshold after seeing output. If a gate fails, the failure is the finding.

AHL / Hill / NARMA are OFF. There is no mean_R gate.
"""

from __future__ import annotations

import csv
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "results" / "i0_seed101" / "colony_timeseries.csv"
CHASSIS = ROOT.parent / "BacteriumFromScratch"

SEED = 101
T_END_S = 23400.0
LOG_DT_S = 10.0
N_TARGET = 256
OVERLAP_CAP_UM = 0.25
MIN_CENTRE_UM = 0.75
OFFPLANE_MAX = 0.05
R_RDISC_MIN = 0.7
R_RDISC_MAX = 1.3
D_CENTRES_MAX_UM = 1.2

REQUIRED = [
    "t_s", "seed", "N", "R_um", "R_disc_um", "R_Rdisc", "offplane",
    "delta_cc_max_um", "d_centers_min_um", "spill",
]

# Frozen 2026-08-24 before I0 code. PROTOCOL.md I0.5.
FROZEN_SHA256 = {
    "results/job2_seed101/job2_timeseries.csv":
        "e5ad8c7809d6b1712b47d4f231de59d71698aab888980e644c7604392df746ae",
    "results/job3_seed101/twobody_timeseries.csv":
        "b3a41af567834cc3076056027f9608f46bd541f8d191d88c8119c690b90eb5c9",
    "results/job3_seed101/cluster_timeseries.csv":
        "bceef8bc7a417a069552a82ebeb54efbd3f2ec8df66e9e67760292e512941774",
    "results/job3_seed101/isolated_packing_timeseries.csv":
        "cb21d5ecce63b59813f33ef54af8388ce70fdaf26608356afb37d72a06c6cf9a",
    "results/job3b_seed101/depletion_timeseries.csv":
        "351f2f2ffa5d8dc23b71dd0697eab7dc3056bdd3fcde017925537c6ff311d829",
    "results/job3b_seed101/packing_cluster_timeseries.csv":
        "e6fe563fa6c1b02db5c72c0a97d960c825e2a6658f96348050d48b2ece577a1c",
    "results/job3b_seed101/uniform_timeseries.csv":
        "ebf19ffa47e980f8b04b827ce9bcfd02bb190c500958c2ea717ec55369e74bbc",
    "results/job3c_seed101/colony_timeseries.csv":
        "f1d83f2e87ee36d3217568a78eb3ec9438c94886b712f0fa76e8c66640dc20cc",
}

JOB3B = (
    "results/job3b_seed101/depletion_timeseries.csv",
    "results/job3b_seed101/packing_cluster_timeseries.csv",
    "results/job3b_seed101/uniform_timeseries.csv",
)


def load():
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=";")
        if reader.fieldnames is None:
            raise ValueError(f"missing header: {CSV_PATH}")
        missing = [c for c in REQUIRED if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"missing columns: {missing}")
        return list(reader)


def gate1_growth(rows):
    errs = []
    n_final = max(int(r["N"]) for r in rows)
    if n_final < N_TARGET:
        errs.append(f"N never reached {N_TARGET} (max {n_final})")
    prev = None
    for r in rows:
        n = int(r["N"])
        if prev is not None:
            if n < prev:
                errs.append(f"N decreased {prev}->{n} at t={r['t_s']}")
                break
            if n - prev > prev:
                errs.append(f"N grew {prev}->{n} at t={r['t_s']}: not binary")
                break
        prev = n
    last_t = float(rows[-1]["t_s"])
    if abs(last_t - T_END_S) > LOG_DT_S:
        errs.append(f"last t={last_t}, expected {T_END_S}")
    return errs, n_final


def gate2_packing(rows):
    errs = []
    for r in rows:
        if int(r["N"]) < 8:
            continue
        if float(r["delta_cc_max_um"]) > OVERLAP_CAP_UM:
            errs.append(
                f"overlap {r['delta_cc_max_um']} > {OVERLAP_CAP_UM} at t={r['t_s']}"
            )
            break
    for r in rows:
        if int(r["N"]) < 8:
            continue
        if float(r["d_centers_min_um"]) < MIN_CENTRE_UM:
            errs.append(
                f"centres {r['d_centers_min_um']} < {MIN_CENTRE_UM} at t={r['t_s']}"
            )
            break
    return errs


def first_at_target(rows):
    for r in rows:
        if int(r["N"]) >= N_TARGET:
            return r
    return None


def gate3_morphology(target):
    if target is None:
        return ["never reached N=256, gate not evaluable"]
    errs = []
    off = float(target["offplane"])
    if off > OFFPLANE_MAX:
        errs.append(f"offplane={off} > {OFFPLANE_MAX}")
    rr = float(target["R_Rdisc"])
    if rr < R_RDISC_MIN or rr > R_RDISC_MAX:
        errs.append(f"R/R_disc={rr:.3f} not in [{R_RDISC_MIN}, {R_RDISC_MAX}]")
        if rr > 2.0:
            errs.append("filament growth — check DivisionMode.SYMMETRY_BROKEN")
    dmin = float(target["d_centers_min_um"])
    if dmin > D_CENTRES_MAX_UM:
        errs.append(
            f"d_centres_min={dmin:.3f} > {D_CENTRES_MAX_UM} um — "
            "lateral Hertzian never engaged"
        )
    return errs


def gate4_fill(rows, n_final):
    errs = []
    n_end = int(rows[-1]["N"])
    if n_end != n_final:
        errs.append(f"N_end={n_end} != N_final={n_final} (emptied)")
    if n_end == 0:
        errs.append("N_end=0 EMPTY — T4/T5 class failure")
    for r in rows:
        if int(r["spill"]) != 0:
            errs.append(f"spill={r['spill']} at t={r['t_s']} — wall Hertzian leaked")
            break
    return errs


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gate5_chassis():
    errs = []
    for rel, expected in FROZEN_SHA256.items():
        path = CHASSIS / rel
        if not path.is_file():
            errs.append(f"missing chassis CSV {rel}")
            continue
        got = sha256(path)
        if got != expected:
            kind = "Job 3b byte-identity" if rel in JOB3B else "chassis CSV"
            errs.append(f"{kind} moved: {rel}\n    expected {expected}\n    got      {got}")
    for script in ("check_job2.py", "check_job3.py", "check_job3b.py", "check_job3c.py"):
        proc = subprocess.run(
            [sys.executable, str(CHASSIS / script)],
            cwd=str(CHASSIS),
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            errs.append(f"{script} exited {proc.returncode}")
            tail = (proc.stdout or proc.stderr or "").strip().splitlines()
            if tail:
                errs.append("    " + tail[-1])
    return errs


def main() -> int:
    if not CSV_PATH.is_file():
        print("FAIL: missing", CSV_PATH)
        return 1
    try:
        rows = load()
    except ValueError as exc:
        print("FAIL:", exc)
        return 1
    if not rows:
        print("FAIL: empty CSV")
        return 1
    if {int(r["seed"]) for r in rows} != {SEED}:
        print(f"FAIL: seed not uniformly {SEED}")
        return 1

    target = first_at_target(rows)
    e1, n_final = gate1_growth(rows)
    e2 = gate2_packing(rows)
    e3 = gate3_morphology(target)
    e4 = gate4_fill(rows, n_final)
    e5 = gate5_chassis()

    failures = []
    for label, errs in (
        ("I0.1 growth", e1),
        ("I0.2 packing", e2),
        ("I0.3 morphology", e3),
        ("I0.4 fill", e4),
        ("I0.5 chassis regression", e5),
    ):
        status = "PASS" if not errs else "FAIL"
        print(f"Gate {label}: {status}")
        for e in errs:
            print(f"    {e}")
        failures.extend(errs)

    print("Gate I0.6 honesty: scored by standing memo (AHL/Hill/NARMA OFF)")
    print()
    print(f"  N_final = {n_final}   N_end = {int(rows[-1]['N'])}   "
          f"t_end = {rows[-1]['t_s']} s")
    print(f"  spill_max = {max(int(r['spill']) for r in rows)}")
    if target is not None:
        print(f"  at N>=256: t = {target['t_s']} s, R = {float(target['R_um']):.2f} um")
        print(f"    R/R_disc = {float(target['R_Rdisc']):.3f}"
              f"   offplane = {float(target['offplane']):.3f}"
              f"   d_centres_min = {float(target['d_centers_min_um']):.3f} um"
              f"   delta_cc_max = {float(target['delta_cc_max_um']):.4f} um")
    print()
    print("  AHL=OFF  Hill=OFF  NARMA=OFF  motility=OFF  death=OFF  nutrient PDE=OFF")
    print()
    if failures:
        print("I0 OVERALL: FAIL")
        return 1
    print("I0 OVERALL: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
