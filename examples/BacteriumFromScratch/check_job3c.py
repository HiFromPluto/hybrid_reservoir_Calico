#!/usr/bin/env python3
"""Evaluate Job 3c colony-dish gates.

Thresholds are frozen in PROTOCOL.md before the run existed. Do not retune a
threshold after seeing output. If a gate fails, the failure is the finding.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "results" / "job3c_seed101" / "colony_timeseries.csv"

SEED = 101
CS_MM = 0.5
KS_MM = 0.02
T_END_S = 23400.0
LOG_DT_S = 10.0

# --- frozen gate thresholds (PROTOCOL.md Job 3c) ---
N_TARGET = 256
OVERLAP_CAP_UM = 0.25
MIN_CENTRE_UM = 0.75            # w0 - overlap cap
DEPLETION_MIN_MM = 1.55e-2      # Gate 3
FLUX_TOL_REL = 1.0e-3           # Gate 5, derived in the acceptance study

# --- Gate 4b morphology band (PROTOCOL.md, frozen) ---
OFFPLANE_MAX = 0.0              # b_z = 2r enforces this; regression guard
R_RDISC_MAX = 2.0               # >2 means filament growth
D_CENTRES_MAX_UM = 1.60         # must approach w0=1.0, not 2*L0=2.0

REQUIRED = [
    "t_s", "seed", "N", "R_um", "R_disc_um", "R_Rdisc", "offplane",
    "N_centre_mM", "N_edge_mM", "field_min_mM",
    "field_mean_mM", "delta_cc_max_um", "d_centers_min_um", "flux_rel_err",
    "gs_iter", "gs_converged", "n_clamp", "marked_um3",
    "inner_n", "inner_rate", "inner_sem", "outer_n", "outer_rate", "outer_sem",
]


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
            # binary fission: each cell yields at most one daughter per tick
            if n - prev > prev:
                errs.append(f"N grew {prev}->{n} at t={r['t_s']}: not binary")
                break
        prev = n
    # NOTE: no power-of-two requirement. Desynchronised division is Gate 4's
    # signal, not a defect -- see PROTOCOL.md.
    return errs, n_final


def gate2_packing(rows):
    errs = []
    for r in rows:
        if float(r["delta_cc_max_um"]) > OVERLAP_CAP_UM:
            errs.append(f"overlap {r['delta_cc_max_um']} > {OVERLAP_CAP_UM} at t={r['t_s']}")
            break
    for r in rows:
        if int(r["N"]) < 2:
            continue
        if float(r["d_centers_min_um"]) < MIN_CENTRE_UM:
            errs.append(f"centres {r['d_centers_min_um']} < {MIN_CENTRE_UM} at t={r['t_s']}")
            break
    return errs


def first_at_target(rows):
    for r in rows:
        if int(r["N"]) >= N_TARGET:
            return r
    return None


def gate3_depletion(row):
    errs = []
    if row is None:
        return ["never reached N=256, gate not evaluable"]
    nc = float(row["N_centre_mM"])
    ne = float(row["N_edge_mM"])
    dn = CS_MM - nc
    if not math.isfinite(nc) or nc < 0.0:
        errs.append(f"bad N_centre={nc}")
    if dn < DEPLETION_MIN_MM:
        errs.append(f"centre depletion {dn:.4e} < {DEPLETION_MIN_MM:.4e} mM")
    if not (nc < ne):
        errs.append(f"centre {nc:.6f} not below edge {ne:.6f}")
    if float(row["field_min_mM"]) < 0.0:
        errs.append(f"negative field_min at t={row['t_s']}")
    if int(row["n_clamp"]) != 0:
        errs.append(f"n_clamp={row['n_clamp']} at t={row['t_s']}")
    return errs


def gate4_gradient(row):
    """Interior elongation rate below edge, separation beyond combined SEM."""
    if row is None:
        return ["never reached N=256, gate not evaluable"], None
    errs = []
    ni, no = int(row["inner_n"]), int(row["outer_n"])
    if ni < 2 or no < 2:
        return [f"insufficient cells to bin (inner={ni}, outer={no})"], None
    ri, si = float(row["inner_rate"]), float(row["inner_sem"])
    ro, so = float(row["outer_rate"]), float(row["outer_sem"])
    sep = ro - ri
    combined = si + so
    if sep <= 0.0:
        errs.append(f"interior rate {ri:.6e} not below edge {ro:.6e}")
    elif sep <= combined:
        errs.append(f"separation {sep:.3e} within combined SEM {combined:.3e}")
    return errs, (ri, si, ro, so, sep, combined)


def gate4b_morphology(rows, target):
    """Colony must be a monolayer disc: not a filament, not a stack.

    The filament run failed this badly (R/R_disc = 4.94, d_centres = 1.98 um)
    while still reporting PASS on Gates 3 and 4 -- which is exactly why this
    gate exists.
    """
    errs = []
    if target is None:
        return ["never reached N=256, gate not evaluable"]
    if float(target["offplane"]) > OFFPLANE_MAX:
        errs.append(f"offplane={target['offplane']} > {OFFPLANE_MAX} "
                    f"(chamber height not enforcing the monolayer)")
    rr = float(target["R_Rdisc"])
    if rr > R_RDISC_MAX:
        errs.append(f"R/R_disc={rr:.2f} > {R_RDISC_MAX} -- filament growth, "
                    f"check DivisionMode")
    dmin = float(target["d_centers_min_um"])
    if dmin > D_CENTRES_MAX_UM:
        errs.append(f"d_centres_min={dmin:.3f} > {D_CENTRES_MAX_UM} um -- "
                    f"lateral Hertzian never engaged")
    return errs


def gate5_solver(rows):
    errs = []
    for r in rows:
        if float(r["t_s"]) <= 0.0:
            continue
        if int(r["gs_converged"]) != 1:
            errs.append(f"gs_converged=0 at t={r['t_s']}")
            break
    for r in rows:
        if float(r["t_s"]) <= 0.0:
            continue
        fx = float(r["flux_rel_err"])
        if not math.isfinite(fx) or fx > FLUX_TOL_REL:
            errs.append(f"flux_rel_err={fx:.3e} > {FLUX_TOL_REL:.0e} at t={r['t_s']}")
            break
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
    last_t = float(rows[-1]["t_s"])
    if abs(last_t - T_END_S) > LOG_DT_S:
        print(f"FAIL: last t={last_t}, expected {T_END_S}")
        return 1

    target = first_at_target(rows)
    failures = []

    e1, n_final = gate1_growth(rows)
    e2 = gate2_packing(rows)
    e3 = gate3_depletion(target)
    e4, g4 = gate4_gradient(target)
    e4b = gate4b_morphology(rows, target)
    e5 = gate5_solver(rows)

    for label, errs in (("1 growth", e1), ("2 packing", e2), ("3 depletion", e3),
                        ("4 radial gradient", e4), ("4b morphology", e4b),
                        ("5 solver", e5)):
        status = "PASS" if not errs else "FAIL"
        print(f"Gate {label}: {status}")
        for e in errs:
            print(f"    {e}")
        failures.extend(errs)

    print()
    print(f"  N_final = {n_final}")
    if target is not None:
        nc = float(target["N_centre_mM"])
        print(f"  at N>=256: t = {target['t_s']} s, R = {float(target['R_um']):.2f} um")
        print(f"    morphology: R/R_disc = {float(target['R_Rdisc']):.2f}"
              f"   offplane = {float(target['offplane']):.2f}"
              f"   d_centres_min = {float(target['d_centers_min_um']):.3f} um")
        print(f"    N_centre = {nc:.6f} mM   depletion = {CS_MM - nc:.4e} mM"
              f"   ({(CS_MM - nc) / DEPLETION_MIN_MM:.2f}x threshold)")
        print(f"    N_edge   = {float(target['N_edge_mM']):.6f} mM")
    if g4 is not None:
        ri, si, ro, so, sep, comb = g4
        print(f"    inner rate {ri:.6e} +/- {si:.1e}   outer {ro:.6e} +/- {so:.1e}")
        print(f"    separation {sep:.3e}  vs combined SEM {comb:.3e}"
              f"   ({sep / comb if comb > 0 else float('nan'):.2f}x)")
        print(f"    RELATIVE separation {sep / ro * 100:.3f}% -- Gate 4 is a SIGN test:")
        print(f"      SEM collapses in a deterministic model, so any nonzero")
        print(f"      gradient clears it. PASS != biologically large.")

    print()
    if failures:
        print("job3c OVERALL: FAIL")
        return 1
    print("job3c OVERALL: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
