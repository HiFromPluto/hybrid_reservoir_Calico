#!/usr/bin/env python3
"""Evaluate job-3 Hertzian packing gates. Do not retune k_cc on FAIL."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
JOB3_DIR = ROOT / "results" / "job3_seed101"
TWO_BODY_CSV = JOB3_DIR / "twobody_timeseries.csv"
CLUSTER_CSV = JOB3_DIR / "cluster_timeseries.csv"
ISOLATED_CSV = JOB3_DIR / "isolated_packing_timeseries.csv"

SEED = 101
W0_UM = 1.0
L0_UM = 1.0
RADIUS_UM = 0.5
NUTRIENT_MM = 0.5
LOG_DT_S = 10.0

TWO_BODY_T_END = 7200.0
TWO_BODY_DELTA0 = 0.20
TWO_BODY_GAP_RESIDUAL = 0.02
TWO_BODY_MIN_DISP = 0.05
TWO_BODY_MAX_DISP = 0.50
TWO_BODY_MAX_GAP = 0.10

CLUSTER_T_END = 6000.0
CLUSTER_OVERLAP_CAP = 0.25
CLUSTER_ALLOWED_N = {2, 4, 8}

ISOLATED_T_END = 100.0
ISOLATED_F_PACK_MAX = 1.0e-8
REL_LEN_TOL = 1e-6
RADIUS_TOL_UM = 1e-9

TWO_BODY_COLUMNS = [
    "t_s", "seed", "delta_cc_um", "gap_um", "disp_a_um", "disp_b_um",
    "F_cc_mag", "y_a", "y_b",
]
CLUSTER_COLUMNS = [
    "t_s", "seed", "N", "delta_cc_max_um", "d_centers_min_um",
    "L_min_um", "L_max_um", "radius_um", "nutrient_mM", "n_negative",
    "F_pack_max",
]
ISOLATED_COLUMNS = [
    "t_s", "seed", "N", "founder_L_um", "founder_L_analytic_um",
    "F_pack_max", "radius_um", "nutrient_mM", "n_negative",
]


def load_rows(path: Path, columns: list[str]):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        header = reader.fieldnames or []
        rows = []
        for raw in reader:
            if not raw or raw.get("t_s") in (None, ""):
                continue
            rows.append(raw)
    missing = [name for name in columns if name not in header]
    return header, rows, missing


def complete(rows, t_end: float) -> bool:
    if not rows:
        return False
    return abs(float(rows[-1]["t_s"]) - t_end) < 1e-6


def seed_ok(rows) -> bool:
    return bool(rows) and all(int(float(row["seed"])) == SEED for row in rows)


def eval_two_body():
    header, rows, missing = load_rows(TWO_BODY_CSV, TWO_BODY_COLUMNS)
    ok_complete = complete(rows, TWO_BODY_T_END) and not missing and seed_ok(rows)
    first = rows[0] if rows else None
    last = rows[-1] if rows else None
    delta0 = float(first["delta_cc_um"]) if first else float("nan")
    gap0 = float(first["gap_um"]) if first else float("nan")
    gap_end = float(last["gap_um"]) if last else float("nan")
    disp_a = float(last["disp_a_um"]) if last else 0.0
    disp_b = float(last["disp_b_um"]) if last else 0.0
    seeded_overlap = first is not None and delta0 > 0.05
    separated = (
        last is not None
        and gap_end >= -TWO_BODY_GAP_RESIDUAL
        and gap_end <= TWO_BODY_MAX_GAP
    )
    moved = (
        TWO_BODY_MIN_DISP <= disp_a <= TWO_BODY_MAX_DISP
        and TWO_BODY_MIN_DISP <= disp_b <= TWO_BODY_MAX_DISP
    )
    passed = ok_complete and seeded_overlap and separated and moved
    return {
        "missing": missing,
        "complete": ok_complete,
        "delta0": delta0,
        "gap0": gap0,
        "gap_end": gap_end,
        "disp_a": disp_a,
        "disp_b": disp_b,
        "seeded_overlap": seeded_overlap,
        "separated": separated,
        "moved": moved,
        "pass": passed,
        "last_t": float(last["t_s"]) if last else float("nan"),
    }


def eval_cluster():
    header, rows, missing = load_rows(CLUSTER_CSV, CLUSTER_COLUMNS)
    ok_complete = complete(rows, CLUSTER_T_END) and not missing and seed_ok(rows)
    ns = [int(float(row["N"])) for row in rows]
    n_ok = bool(ns) and all(n in CLUSTER_ALLOWED_N for n in ns)
    n_final = ns[-1] if ns else 0
    two_divs = n_final >= 8 and 2 in ns and 4 in ns and 8 in ns
    max_delta = max((float(row["delta_cc_max_um"]) for row in rows), default=float("inf"))
    min_centres = min((float(row["d_centers_min_um"]) for row in rows), default=0.0)
    overlap_ok = max_delta <= CLUSTER_OVERLAP_CAP + 1e-12
    centres_ok = min_centres >= (W0_UM - CLUSTER_OVERLAP_CAP) - 1e-12
    n_neg = all(int(float(row["n_negative"])) == 0 for row in rows) if rows else False
    radius_ok = all(abs(float(row["radius_um"]) - RADIUS_UM) <= RADIUS_TOL_UM for row in rows) if rows else False
    nutrient_ok = all(abs(float(row["nutrient_mM"]) - NUTRIENT_MM) <= 1e-12 for row in rows) if rows else False
    def finite_row(row) -> bool:
        for key in ("delta_cc_max_um", "d_centers_min_um", "F_pack_max", "L_min_um"):
            try:
                if not math.isfinite(float(row[key])):
                    return False
            except (TypeError, ValueError):
                return False
        return True

    finite = bool(rows) and all(finite_row(row) for row in rows)
    last_centres = float(rows[-1]["d_centers_min_um"]) if rows else float("inf")
    in_box = finite and last_centres < 8.0
    passed = (
        ok_complete and n_ok and two_divs and overlap_ok and centres_ok
        and n_neg and radius_ok and nutrient_ok and finite and in_box
    )
    return {
        "missing": missing,
        "complete": ok_complete,
        "n_final": n_final,
        "n_ok": n_ok,
        "two_divs": two_divs,
        "max_delta": max_delta,
        "min_centres": min_centres,
        "overlap_ok": overlap_ok,
        "centres_ok": centres_ok,
        "pass": passed,
        "last_t": float(rows[-1]["t_s"]) if rows else float("nan"),
    }


def eval_isolated():
    header, rows, missing = load_rows(ISOLATED_CSV, ISOLATED_COLUMNS)
    ok_complete = complete(rows, ISOLATED_T_END) and not missing and seed_ok(rows)
    n1 = bool(rows) and all(int(float(row["N"])) == 1 for row in rows)
    fmax = max((float(row["F_pack_max"]) for row in rows), default=float("inf"))
    f_ok = fmax <= ISOLATED_F_PACK_MAX
    rel_err = 0.0
    for row in rows:
        rel_err = max(
            rel_err,
            abs(float(row["founder_L_um"]) - float(row["founder_L_analytic_um"])) / L0_UM,
        )
    len_ok = bool(rows) and rel_err <= REL_LEN_TOL
    n_neg = all(int(float(row["n_negative"])) == 0 for row in rows) if rows else False
    passed = ok_complete and n1 and f_ok and len_ok and n_neg
    return {
        "missing": missing,
        "complete": ok_complete,
        "fmax": fmax,
        "f_ok": f_ok,
        "rel_err": rel_err,
        "n1": n1,
        "pass": passed,
        "last_t": float(rows[-1]["t_s"]) if rows else float("nan"),
    }


def main() -> int:
    two = eval_two_body()
    cl = eval_cluster()
    iso = eval_isolated()

    print(f"Two-body CSV {TWO_BODY_CSV}")
    print(f"  last_t={two['last_t']} complete={two['complete']} missing={two['missing']}")
    print(f"  seed overlap delta0={two['delta0']:.4f} um  gap0={two['gap0']:.4f}")
    print(f"Gate 1a seeded delta_cc>0: {'PASS' if two['seeded_overlap'] else 'FAIL'}")
    print(f"Gate 1b gap in [-{TWO_BODY_GAP_RESIDUAL}, {TWO_BODY_MAX_GAP}] um: {'PASS' if two['separated'] else 'FAIL'}"
          f"  gap_end={two['gap_end']:.4f}")
    print(f"Gate 1c displacement in [{TWO_BODY_MIN_DISP}, {TWO_BODY_MAX_DISP}] um: {'PASS' if two['moved'] else 'FAIL'}"
          f"  disp_a={two['disp_a']:.4f} disp_b={two['disp_b']:.4f}")
    print(f"Two-body overall: {'PASS' if two['pass'] else 'FAIL'}")

    print(f"Cluster CSV {CLUSTER_CSV}")
    print(f"  last_t={cl['last_t']} complete={cl['complete']} missing={cl['missing']}")
    print(f"  N_final={cl['n_final']} max_delta={cl['max_delta']:.4f} min_centres={cl['min_centres']:.4f}")
    print(f"Gate 2a binary fission N in {{2,4,8}} and N>=8: {'PASS' if cl['n_ok'] and cl['two_divs'] else 'FAIL'}")
    print(f"Gate 2b max delta_cc <= {CLUSTER_OVERLAP_CAP} um: {'PASS' if cl['overlap_ok'] else 'FAIL'}")
    print(f"Gate 2c centres >= w0 - cap: {'PASS' if cl['centres_ok'] else 'FAIL'}")
    print(f"Cluster overall: {'PASS' if cl['pass'] else 'FAIL'}")

    print(f"Isolated packing CSV {ISOLATED_CSV}")
    print(f"  last_t={iso['last_t']} complete={iso['complete']} F_pack_max={iso['fmax']:.3e}")
    print(f"Gate 3a N=1 packing force ~0: {'PASS' if iso['n1'] and iso['f_ok'] else 'FAIL'}")
    print(f"Gate 3b Job 2 length analytic: {'PASS' if iso['rel_err'] <= REL_LEN_TOL else 'FAIL'}"
          f"  rel_err={iso['rel_err']:.3e}")
    print(f"Isolated overall: {'PASS' if iso['pass'] else 'FAIL'}")

    all_pass = two["pass"] and cl["pass"] and iso["pass"]
    print(f"Overall: {'PASS' if all_pass else 'FAIL'}")
    if not two["pass"] or not cl["pass"]:
        print("STOP: do not raise k_cc, do not enlarge radius, do not start QS.")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
