#!/usr/bin/env python3
"""Evaluate ChassisPocket I0b DoorWeir gates (W=20 µm).

Thresholds are frozen in PROTOCOL.md before the run existed. Do not retune a
threshold, W, or founder after seeing output. If a gate fails, the failure
is the finding.

AHL / Hill / NARMA are OFF. There is no mean_R gate.
spill_cum > 0 is allowed. T5 failed because N_end = 0.
"""

from __future__ import annotations

import csv
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "results" / "i0b_seed101" / "colony_timeseries.csv"
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
W_UM = 20.0
FOUNDER_Y = 80.0

REQUIRED = [
    "t_s", "seed", "N", "R_um", "R_disc_um", "R_Rdisc", "offplane",
    "delta_cc_max_um", "d_centers_min_um", "spill_tick", "spill_cum",
    "N_ever", "door_contact", "y_max_um", "wall_leak",
]

# Same freeze as I0.5 / PROTOCOL.md. Recorded 2026-08-24 before I0 code.
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
    n_ever_final = max(int(r["N_ever"]) for r in rows)
    if n_ever_final < N_TARGET:
        errs.append(f"N_ever never reached {N_TARGET} (max {n_ever_final})")
    prev = None
    for r in rows:
        n_ever = int(r["N_ever"])
        if prev is not None:
            if n_ever < prev:
                errs.append(f"N_ever decreased {prev}->{n_ever} at t={r['t_s']}")
                break
            if n_ever - prev > prev and prev > 0:
                errs.append(f"N_ever grew {prev}->{n_ever} at t={r['t_s']}: not binary")
                break
        prev = n_ever
    last_t = float(rows[-1]["t_s"])
    if abs(last_t - T_END_S) > LOG_DT_S:
        errs.append(f"last t={last_t}, expected {T_END_S}")
    return errs, n_ever_final


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
        return ["never reached N=256 in garage, gate not evaluable"]
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


def gate4_fill(rows):
    errs = []
    n_end = int(rows[-1]["N"])
    n_max = max(int(r["N"]) for r in rows)
    if n_end == 0:
        errs.append("N_end=0 EMPTY — T4/T5 class failure")
    if n_end < N_TARGET:
        errs.append(
            f"N_end={n_end} < {N_TARGET} (drained through W={W_UM:.0f}; N_max={n_max})"
        )
    if any(int(r["wall_leak"]) != 0 for r in rows):
        errs.append("wall_leak > 0 — closed-face Hertzian leaked")
    if not any(int(r["door_contact"]) == 1 for r in rows):
        errs.append(
            "door_contact never fired: I0b did not ask the door question "
            f"(founder_y={FOUNDER_Y})"
        )
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
    e1, n_ever_final = gate1_growth(rows)
    e2 = gate2_packing(rows)
    e3 = gate3_morphology(target)
    e4 = gate4_fill(rows)
    e5 = gate5_chassis()

    failures = []
    for label, errs in (
        ("I0b.1 growth", e1),
        ("I0b.2 packing", e2),
        ("I0b.3 morphology", e3),
        ("I0b.4 fill", e4),
        ("I0b.5 chassis regression", e5),
    ):
        status = "PASS" if not errs else "FAIL"
        print(f"Gate {label}: {status}")
        for e in errs:
            print(f"    {e}")
        failures.extend(errs)

    print("Gate I0b.6 honesty: scored by standing memo (AHL/Hill/NARMA OFF, W=20 frozen)")
    print()
    last = rows[-1]
    print(f"  N_end = {int(last['N'])}   N_ever = {n_ever_final}   "
          f"spill_cum = {int(last['spill_cum'])}   t_end = {last['t_s']} s")
    print(f"  door_contact_any = {max(int(r['door_contact']) for r in rows)}   "
          f"y_max_end = {float(last['y_max_um']):.2f} um   "
          f"wall_leak_max = {max(int(r['wall_leak']) for r in rows)}")
    print(f"  W = {W_UM:.0f} um   founder_y = {FOUNDER_Y:.0f} um")
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
        print("I0b OVERALL: FAIL")
        return 1
    print("I0b OVERALL: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
