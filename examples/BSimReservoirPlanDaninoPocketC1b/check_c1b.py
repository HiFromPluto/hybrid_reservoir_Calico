#!/usr/bin/env python3
"""NARMA-blind C1b checker. Island C0-equivalent Object B. NOT_FIG4B."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D0B = ROOT / "examples" / "BSimReservoirPlanPocketOscSI_D0b"
D1 = ROOT / "examples" / "BSimReservoirPlanPocketOscSI_D1"
RESULTS = HERE / "results"
CLAIM = ROOT / "examples" / "PocketDish" / "CLAIM_FREEZE_DANINO_SI.md"
C1_STANDING = ROOT / "examples" / "PocketDish" / "C1_PACKED_SPATIAL_STANDING.md"
PROTOCOL = HERE / "PROTOCOL.md"

if str(D0B) not in sys.path:
    sys.path.insert(0, str(D0B))

from period_check import extract_period, refuse_narma  # noqa: E402

PROTOCOL_PEAKS = {
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

PERIOD_REF = 63.583333333333336
PERIOD_REL = 0.02
TRAJ_RTOL = 0.02
TRAJ_ATOL = 1e-4
LEDGER_REL = 1e-3


def load_csv(path: Path):
    t = []
    y = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "NOT_FIG4B" not in (reader.fieldnames or []):
            raise SystemExit(f"{path.name} missing NOT_FIG4B")
        for row in reader:
            t.append(float(row["t"]))
            y.append([float(row["A"]), float(row["I"]), float(row["Hi"]), float(row["He"])])
    return np.array(t), np.array(y)


def main() -> None:
    refuse_narma()
    if "DaninoSI_OccupiedDF" not in CLAIM.read_text(encoding="utf-8"):
        raise SystemExit("claim freeze missing Object B")
    standing = C1_STANDING.read_text(encoding="utf-8")
    if "**Status: FAIL**" not in standing:
        raise SystemExit("C1 standing must remain FAIL")
    proto = PROTOCOL.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or "C1b_ISLAND" not in proto:
        raise SystemExit("C1b PROTOCOL not frozen")
    summary = json.loads((RESULTS / "c1b_summary.json").read_text(encoding="utf-8"))
    print(f"C1={summary.get('C1')} C1b={summary.get('C1b')} NOT_FIG4B")
    print(f"open_chip_still_off={summary.get('open_chip_still_off')} a0_automatic={summary.get('a0_automatic')}")

    ok = True
    t, y = load_csv(RESULTS / "java_C1b_ISLAND_D05_MU040.csv")
    py = np.load(D1 / "fixtures" / "primary_mu_0p40.npz")
    allowed = TRAJ_ATOL + TRAJ_RTOL * np.abs(py["Y"])
    viol = int(np.sum(np.abs(y - py["Y"]) > allowed + 1e-15))
    peaks = extract_period(t, y[:, 1], PROTOCOL_PEAKS)
    rel = abs(float(peaks["period"]) - PERIOD_REF) / PERIOD_REF if peaks["flag"] == "OSC" else float("nan")
    ident_ok = viol == 0 and peaks["flag"] == "OSC" and rel < PERIOD_REL
    ok = ok and ident_ok
    print(
        f"  C1b_ISLAND_D05_MU040 NOT_FIG4B java={peaks['flag']} T={peaks.get('period')} "
        f"relT={rel} viol={viol} gate={'PASS' if ident_ok else 'FAIL'}"
    )

    for name, required in (
        ("C1b_ISLAND_WRONG_KICK", "NO_PERIOD"),
        ("C1b_ISLAND_MU150", "NO_PERIOD"),
    ):
        tt, yy = load_csv(RESULTS / f"java_{name}.csv")
        p = extract_period(tt, yy[:, 1], PROTOCOL_PEAKS)
        arm_ok = p["flag"] == required
        ok = ok and arm_ok
        print(f"  {name} NOT_FIG4B java={p['flag']} gate={'PASS' if arm_ok else 'FAIL'}")

    arms = {a["name"]: a for a in summary.get("arms", [])}
    rel_l = float(arms["C1b_ISLAND_MU0_LEDGER"]["ledger_rel"])
    ledger_ok = rel_l <= LEDGER_REL and bool(arms["C1b_ISLAND_MU0_LEDGER"]["pass"])
    ok = ok and ledger_ok
    print(f"  C1b_ISLAND_MU0_LEDGER NOT_FIG4B ledger_rel={rel_l:.4g} gate={'PASS' if ledger_ok else 'FAIL'}")

    open_ok = arms["C1_OPEN_DILUTE"]["flag"] == "NO_PERIOD"
    ok = ok and open_ok
    print(f"  C1_OPEN_DILUTE NOT_FIG4B report-only {arms['C1_OPEN_DILUTE']['flag']} gate={'PASS' if open_ok else 'FAIL'}")

    extra = arms.get("C1b_ISLAND_D1_800", {})
    print(
        f"  C1b_ISLAND_D1_800 NOT_FIG4B not identity flag={extra.get('flag')} "
        f"T={extra.get('period')} (report-only)"
    )

    if summary.get("C1") != "FAIL":
        print("C1 standing/summary must remain FAIL")
        ok = False
    if not ok or summary.get("C1b") != "PASS":
        print("C1b FAIL NOT_FIG4B. C1 FAIL unchanged. Do not start A0/W0.")
        sys.exit(1)
    print(
        "C1b PASS on island C0-equivalent occupancy. NOT_FIG4B. "
        "C1 open chip remains FAIL. A0 is not automatic. W0 later."
    )


if __name__ == "__main__":
    main()
