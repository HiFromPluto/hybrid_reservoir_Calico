#!/usr/bin/env python3
"""NARMA-blind C0 checker. Object B well-mixed coupling. NOT_FIG4B."""

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
D1 = ROOT / "examples" / "BSimReservoirPlanPocketOscSI_D1"
RESULTS = HERE / "results"
CLAIM = ROOT / "examples" / "PocketDish" / "CLAIM_FREEZE_DANINO_SI.md"
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

TRAJ_RTOL = 0.02
TRAJ_ATOL = 1e-4
PERIOD_REL = 0.02
LEDGER_REL = 1e-3


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
    proto = PROTOCOL.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or "NOT_FIG4B" not in proto:
        raise SystemExit("C0 PROTOCOL must be frozen and NOT_FIG4B")
    summary_path = RESULTS / "c0_summary.json"
    if not summary_path.exists():
        raise SystemExit("run ant c0-wellmixed first")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    print(f"object={summary.get('object')} status={summary.get('object_status')}")
    print(
        f"D0={summary.get('D0')} D0b={summary.get('D0b')} D1={summary.get('D1')} "
        f"P0={summary.get('P0')} C0={summary.get('C0')} NOT_FIG4B"
    )
    print(
        f"N={summary.get('n_cells')} v_cell={summary.get('v_cell_um3')} "
        f"V_e={summary.get('V_e_um3')} d={summary.get('d')} NOT_FIG4B"
    )

    ok = True
    t_j, y_j = load_csv(RESULTS / "java_C0_IDENTICAL_D05_MU040.csv")
    py = np.load(D1 / "fixtures" / "primary_mu_0p40.npz")
    allowed = TRAJ_ATOL + TRAJ_RTOL * np.abs(py["Y"])
    err = np.abs(y_j - py["Y"])
    viol = int(np.sum(err > allowed + 1e-15))
    max_norm = float(np.max(err / np.maximum(allowed, 1e-30)))
    peaks = extract_period(t_j, y_j[:, 1], PROTOCOL_PEAKS)
    py_t = 63.583333333333336
    rel = abs(float(peaks["period"]) - py_t) / py_t if peaks["flag"] == "OSC" else float("nan")
    arm_ok = viol == 0 and peaks["flag"] == "OSC" and rel < PERIOD_REL
    ok = ok and arm_ok
    print(
        f"  C0_IDENTICAL_D05_MU040 NOT_FIG4B java={peaks['flag']} "
        f"viol={viol} max_norm={max_norm:.4g} relT={rel} gate={'PASS' if arm_ok else 'FAIL'}"
    )

    t32, y32 = load_csv(RESULTS / "java_C0_IDENTICAL_D05_MU032.csv")
    p32 = extract_period(t32, y32[:, 1], PROTOCOL_PEAKS)
    py32 = 56.30769230769231
    rel32 = abs(float(p32["period"]) - py32) / py32 if p32["flag"] == "OSC" else float("nan")
    ok32 = p32["flag"] == "OSC" and rel32 < PERIOD_REL
    ok = ok and ok32
    print(f"  C0_IDENTICAL_D05_MU032 NOT_FIG4B java={p32['flag']} relT={rel32} gate={'PASS' if ok32 else 'FAIL'}")

    t15, y15 = load_csv(RESULTS / "java_C0_IDENTICAL_D05_MU150.csv")
    p15 = extract_period(t15, y15[:, 1], PROTOCOL_PEAKS)
    ok15 = p15["flag"] == "NO_PERIOD"
    ok = ok and ok15
    print(f"  C0_IDENTICAL_D05_MU150 NOT_FIG4B java={p15['flag']} gate={'PASS' if ok15 else 'FAIL'}")

    tb, yb = load_csv(RESULTS / "java_C0_WRONG_KICK_BASAL.csv")
    pb = extract_period(tb, yb[:, 1], PROTOCOL_PEAKS)
    okb = pb["flag"] == "NO_PERIOD"
    ok = ok and okb
    print(f"  C0_WRONG_KICK_BASAL NOT_FIG4B java={pb['flag']} gate={'PASS' if okb else 'FAIL'}")

    arms = {a["name"]: a for a in summary.get("arms", [])}
    for name in ("C0_MU0_LEDGER", "C0_MU0_CLOSED_MEMBRANE"):
        rel_l = float(arms[name]["ledger_rel"])
        ledger_ok = rel_l <= LEDGER_REL and bool(arms[name]["pass"])
        ok = ok and ledger_ok
        print(f"  {name} NOT_FIG4B ledger_rel={rel_l:.4g} gate={'PASS' if ledger_ok else 'FAIL'}")

    hetero = arms["C0_HETERO_ONE_KICK"]
    he_move = float(hetero["he_move"])
    hetero_ok = he_move > 1e-6 and bool(hetero["pass"])
    ok = ok and hetero_ok
    print(
        f"  C0_HETERO_ONE_KICK NOT_FIG4B not Fig4e he_move={he_move:.4g} "
        f"gate={'PASS' if hetero_ok else 'FAIL'}"
    )

    if summary.get("C0") == "PASS" and not ok:
        print("Java summary PASS disagrees with Python re-check.")
        ok = False
    if not ok or summary.get("C0") != "PASS":
        print("C0 FAIL NOT_FIG4B. Do not retune Object B. Do not start C1. No NARMA.")
        sys.exit(1)
    print("C0 PASS on Object B well-mixed coupling. NOT_FIG4B. C1 may start. D0/D0b unchanged.")


if __name__ == "__main__":
    main()
