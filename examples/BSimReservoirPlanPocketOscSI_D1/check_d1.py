#!/usr/bin/env python3
"""NARMA-blind D1 checker. Object B NOT_FIG4B. Reads Java CSVs vs D0b npz."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D0B = ROOT / "examples" / "BSimReservoirPlanPocketOscSI_D0b"
RESULTS = HERE / "results"
FIXTURES = HERE / "fixtures"
CLAIM = ROOT / "examples" / "PocketDish" / "CLAIM_FREEZE_DANINO_SI.md"

if str(D0B) not in sys.path:
    sys.path.insert(0, str(D0B))

from period_check import extract_period, refuse_narma  # noqa: E402


PROTOCOL = {
    "t_discard": 180.0,
    "sample_dt": 0.5,
    "peak_min_spacing": 25.0,
    "peak_rel_prominence": 0.20,
    "min_peaks": 4,
    "persist_last_fraction": 0.30,
    "amp_last_fraction": 0.40,
    "amp_persist_ratio": 0.40,
    "rel_amplitude_min": 0.25,
}

TRAJ_RTOL = 0.02
TRAJ_ATOL = 1e-4
PERIOD_REL = 0.02

ARMS = [
    ("primary_mu_0p40", "OSC", 63.583333333333336),
    ("primary_osc", "OSC", 56.30769230769231),
    ("si_basal_perturb_mu_0p40", "NO_PERIOD", None),
    ("failed_prior_mu_1p5", "NO_PERIOD", None),
]


def load_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    t = []
    y = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            t.append(float(row["t"]))
            y.append([float(row["A"]), float(row["I"]), float(row["Hi"]), float(row["He"])])
    return np.array(t), np.array(y)


def main() -> None:
    refuse_narma()
    if not CLAIM.exists():
        raise SystemExit("CLAIM_FREEZE_DANINO_SI.md missing")
    text = CLAIM.read_text(encoding="utf-8")
    if "DaninoSI_OccupiedDF" not in text or "NOT_FIG4B" not in text:
        raise SystemExit("claim freeze must name Object B NOT_FIG4B")
    summary_path = RESULTS / "d1_summary.json"
    if not summary_path.exists():
        raise SystemExit("run ant d1-parity first")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    print(f"object={summary.get('object')} status={summary.get('object_status')}")
    print(f"D0={summary.get('D0')} D0b={summary.get('D0b')} D1={summary.get('D1')} NOT_FIG4B")

    ok = True
    for name, required, py_t in ARMS:
        csv_path = RESULTS / f"java_{name}.csv"
        npz_path = FIXTURES / f"{name}.npz"
        t_j, y_j = load_csv(csv_path)
        py = np.load(npz_path)
        t_p = py["t"]
        y_p = py["Y"]
        if t_j.shape != t_p.shape:
            print(f"{name} grid mismatch")
            ok = False
            continue
        allowed = TRAJ_ATOL + TRAJ_RTOL * np.abs(y_p)
        err = np.abs(y_j - y_p)
        viol = int(np.sum(err > allowed + 1e-15))
        max_norm = float(np.max(err / np.maximum(allowed, 1e-30)))
        peaks = extract_period(t_j, y_j[:, 1], PROTOCOL)
        flag_ok = peaks["flag"] == required
        period_ok = True
        rel = float("nan")
        if required == "OSC":
            tj = float(peaks["period"])
            rel = abs(tj - py_t) / py_t
            period_ok = rel < PERIOD_REL
        arm_ok = viol == 0 and flag_ok and period_ok
        ok = ok and arm_ok
        print(
            f"  {name} NOT_FIG4B py_flag={required} java={peaks['flag']} "
            f"viol={viol} max_norm={max_norm:.4g} relT={rel} gate={'PASS' if arm_ok else 'FAIL'}"
        )

    csv15 = RESULTS / "java_primary_mu_1p5.csv"
    t15, y15 = load_csv(csv15)
    p15 = extract_period(t15, y15[:, 1], PROTOCOL)
    off_ok = p15["flag"] == "NO_PERIOD"
    ok = ok and off_ok
    print(f"  primary_mu_1p5 NOT_FIG4B java={p15['flag']} gate={'PASS' if off_ok else 'FAIL'}")

    if summary.get("D1") != ("PASS" if ok else "FAIL") and summary.get("D1") == "PASS" and not ok:
        print("Java summary PASS disagrees with Python re-check.")
        ok = False
    if not ok:
        print("D1 FAIL NOT_FIG4B. Do not retune. Do not claim Fig. 4b. Do not start C0/C1.")
        sys.exit(1)
    if summary.get("D1") == "FAIL":
        print("D1 FAIL NOT_FIG4B (Java summary). Do not retune. Do not start C0/C1.")
        sys.exit(1)
    print("D1 PASS on Object B fixtures. NOT_FIG4B. D0/D0b identity standings unchanged.")


if __name__ == "__main__":
    main()
