#!/usr/bin/env python3
"""Evaluate job-3b nutrient-field gates. Do not retune lambda_S or K_S on FAIL."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
JOB3B_DIR = ROOT / "results" / "job3b_seed101"
UNIFORM_CSV = JOB3B_DIR / "uniform_timeseries.csv"
PACKING_CSV = JOB3B_DIR / "packing_cluster_timeseries.csv"
DEPLETION_CSV = JOB3B_DIR / "depletion_timeseries.csv"

SEED = 101
LAMBDA_S_PER_HOUR = 1.0
NU_PER_SECOND = LAMBDA_S_PER_HOUR / 3600.0
KS_MM = 0.02
CS_MM = 0.5
L0_UM = 1.0
LDIV_UM = 3.0
RADIUS_UM = 0.5
LOG_DT_S = 10.0
REL_LEN_TOL = 1e-6
LEN_AFTER_DIV_TOL_UM = 0.02
RADIUS_TOL_UM = 1e-9

UNIFORM_T_END = 9000.0
PACKING_T_END = 6000.0
DEPLETION_T_END = 6000.0
CLUSTER_OVERLAP_CAP = 0.25
CLUSTER_MIN_CENTRE = 0.75  # w0 - overlap cap
CLUSTER_ALLOWED_N = {2, 4, 8}
FLUX_TOL_REL = 0.005
FLUX_WARMUP_S = 0.0  # 3b.3: conv_tol 1e-9 resolves the start transient; no rows exempt


def monod(n: float) -> float:
    return n / (n + KS_MM)


def analytic_tdiv(n: float = CS_MM) -> float:
    return math.log(2.0) / (NU_PER_SECOND * monod(n))


def load_rows(path: Path, required: list[str]):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        if reader.fieldnames is None:
            raise ValueError(f"missing header: {path}")
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} missing columns: {missing}")
        return list(reader)


def check_uniform(rows: list[dict]) -> list[str]:
    errs: list[str] = []
    if not rows:
        return ["uniform: empty CSV"]
    last = rows[-1]
    t_last = float(last["t_s"])
    if abs(t_last - UNIFORM_T_END) > LOG_DT_S:
        errs.append(f"uniform: last t={t_last}, expected {UNIFORM_T_END}")

    pre_div = [r for r in rows if int(r["N"]) == 1]
    for r in pre_div:
        l = float(r["founder_L_um"])
        la = float(r["founder_L_analytic_um"])
        if abs(l - la) / L0_UM > REL_LEN_TOL:
            errs.append(f"uniform: |L-analytic|/L0={abs(l-la)/L0_UM} at t={r['t_s']}")
            break

    div_rows = [r for r in rows if int(r["N"]) == 2]
    if not div_rows:
        errs.append("uniform: never divided")
    else:
        t_div = float(div_rows[0]["t_s"])
        if abs(t_div - analytic_tdiv()) > LOG_DT_S:
            errs.append(f"uniform: N=2 at t={t_div}, analytic T_div={analytic_tdiv():.1f}")

    for r in rows:
        if int(r["n_negative"]) != 0:
            errs.append(f"uniform: n_negative={r['n_negative']} at t={r['t_s']}")
            break
        nloc = float(r["N_local_mM"])
        if abs(nloc - CS_MM) > 1e-6:
            errs.append(f"uniform: N_local={nloc} != C_s at t={r['t_s']}")
            break
        if abs(float(r["radius_um"]) - RADIUS_UM) > RADIUS_TOL_UM:
            errs.append(f"uniform: bad radius at t={r['t_s']}")
            break

    n_set = {int(r["N"]) for r in rows}
    if not n_set.issubset({1, 2, 4, 8}):
        errs.append(f"uniform: N not binary fission: {sorted(n_set)}")
    return errs


def check_packing(rows: list[dict]) -> list[str]:
    errs: list[str] = []
    if not rows:
        return ["packing: empty CSV"]
    last = rows[-1]
    if abs(float(last["t_s"]) - PACKING_T_END) > LOG_DT_S:
        errs.append(f"packing: last t={last['t_s']}, expected {PACKING_T_END}")

    n_seen = {int(r["N"]) for r in rows}
    for target in (2, 4, 8):
        if target not in n_seen:
            errs.append(f"packing: never saw N={target}")

    for r in rows:
        if int(r["n_negative"]) != 0:
            errs.append(f"packing: n_negative at t={r['t_s']}")
            break
        if float(r["delta_cc_max_um"]) > CLUSTER_OVERLAP_CAP:
            errs.append(f"packing: overlap {r['delta_cc_max_um']} at t={r['t_s']}")
            break
        if float(r["d_centers_min_um"]) < CLUSTER_MIN_CENTRE:
            errs.append(f"packing: centres too close at t={r['t_s']}")
            break
        nloc = float(r["N_local_mM"])
        if abs(nloc - CS_MM) > 1e-4:
            errs.append(f"packing: N_local={nloc} != C_s at t={r['t_s']}")
            break

    if int(last["N"]) not in CLUSTER_ALLOWED_N:
        errs.append(f"packing: final N={last['N']}")
    return errs


def check_depletion(rows: list[dict]) -> list[str]:
    errs: list[str] = []
    if not rows:
        return ["depletion: empty CSV"]
    last = rows[-1]
    if abs(float(last["t_s"]) - DEPLETION_T_END) > LOG_DT_S:
        errs.append(f"depletion: last t={last['t_s']}, expected {DEPLETION_T_END}")

    # Must show depletion at some point after t=0
    depleted = False
    slowed = False
    monod_bath = monod(CS_MM)
    for r in rows:
        if float(r["t_s"]) <= 0.0:
            continue
        if int(r.get("gs_converged", "0")) != 1:
            errs.append(f"depletion: GS did not converge at t={r['t_s']}")
            return errs
        nloc = float(r["N_local_mM"])
        if not math.isfinite(nloc) or nloc < 0.0:
            errs.append(f"depletion: bad N_local={nloc} at t={r['t_s']}")
            return errs
        fmin = float(r["field_min_mM"])
        fmean = float(r["field_mean_mM"])
        if not math.isfinite(fmin) or not math.isfinite(fmean):
            errs.append(f"depletion: non-finite field stats at t={r['t_s']}")
            return errs
        if float(r["t_s"]) > 0.0:
            mloc = monod(nloc)
            if (nloc < CS_MM - 1e-8 or fmin < CS_MM - 1e-8) and mloc < monod_bath - 1e-9:
                depleted = True
            l = float(r["founder_L_um"])
            la = float(r["founder_L_analytic_uniform_mM"])
            if l + 1e-9 < la:
                slowed = True

    if not depleted:
        errs.append("depletion: never saw N_local < C_s with Monod < bath")
    if not slowed:
        errs.append("depletion: founder length never below uniform-bath analytic")

    flux_vals = []
    for r in rows:
        t = float(r["t_s"])
        if t <= 0.0 or t < FLUX_WARMUP_S:
            continue
        flux = float(r["flux_rel_err"])
        if not math.isfinite(flux):
            errs.append(f"depletion: non-finite flux_rel_err at t={t}")
            return errs
        flux_vals.append(flux)
        if flux > FLUX_TOL_REL:
            errs.append(f"depletion: flux_rel_err={flux} at t={t} (after warmup)")

    if not flux_vals:
        errs.append("depletion: no flux rows after warmup")

    # Coupled growth: at least one division
    max_n = max(int(r["N"]) for r in rows)
    if max_n < 4:
        errs.append(f"depletion: no division under spatial field (max N={max_n})")

    div_rows = [r for r in rows if int(r["N"]) >= 4]
    for r in div_rows:
        if float(r["delta_cc_max_um"]) > CLUSTER_OVERLAP_CAP:
            errs.append(f"depletion: overlap {r['delta_cc_max_um']} at t={r['t_s']}")
            break
        if float(r["d_centers_min_um"]) < CLUSTER_MIN_CENTRE:
            errs.append(f"depletion: centres too close at t={r['t_s']}")
            break

    return errs


def main() -> int:
    failures: list[str] = []
    for path in (UNIFORM_CSV, PACKING_CSV, DEPLETION_CSV):
        if not path.is_file():
            failures.append(f"missing {path}")
            continue

    if failures:
        for f in failures:
            print("FAIL:", f)
        return 1

    uniform_req = [
        "t_s", "seed", "N", "founder_L_um", "founder_L_analytic_um",
        "N_local_mM", "field_min_mM", "field_mean_mM", "monod_local",
        "radius_um", "n_negative",
    ]
    packing_req = [
        "t_s", "seed", "N", "delta_cc_max_um", "d_centers_min_um",
        "L_min_um", "L_max_um", "radius_um", "N_local_mM",
        "field_min_mM", "field_mean_mM", "monod_local", "n_negative", "F_pack_max",
    ]
    depletion_req = [
        "t_s", "seed", "N", "founder_L_um", "founder_L_analytic_uniform_mM",
        "N_local_mM", "field_min_mM", "field_mean_mM", "monod_local",
        "delta_cc_max_um", "d_centers_min_um", "flux_rel_err", "gs_iter",
        "gs_converged", "n_clamp", "F_pack_max",
    ]

    try:
        u_rows = load_rows(UNIFORM_CSV, uniform_req)
        p_rows = load_rows(PACKING_CSV, packing_req)
        d_rows = load_rows(DEPLETION_CSV, depletion_req)
    except ValueError as exc:
        print("FAIL:", exc)
        return 1

    for label, errs in (
        ("uniform", check_uniform(u_rows)),
        ("packing", check_packing(p_rows)),
        ("depletion", check_depletion(d_rows)),
    ):
        for e in errs:
            print(f"FAIL [{label}]:", e)
            failures.append(e)

    if failures:
        return 1

    print("PASS job3b")
    print(f"  seed={SEED}")
    print(f"  uniform t_end={UNIFORM_T_END} s")
    print(f"  packing t_end={PACKING_T_END} s")
    print(f"  depletion t_end={DEPLETION_T_END} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
