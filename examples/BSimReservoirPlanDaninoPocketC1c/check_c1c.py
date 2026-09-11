#!/usr/bin/env python3
"""NARMA-blind C1c checker. Filled-pocket Object B. NOT_FIG4B."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D0B = ROOT / "examples" / "BSimReservoirPlanPocketOscSI_D0b"
RESULTS = HERE / "results"
CLAIM = ROOT / "examples" / "PocketDish" / "CLAIM_FREEZE_DANINO_SI.md"
C1_STANDING = ROOT / "examples" / "PocketDish" / "C1_PACKED_SPATIAL_STANDING.md"
C1B_STANDING = ROOT / "examples" / "PocketDish" / "C1B_ISLAND_STANDING.md"
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
PERIOD_REL_BAND = 0.20
LEDGER_REL = 1e-3
HE_RANGE_MIN = 1e-6


def load_mean_csv(path: Path):
    t = []
    i_mean = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "NOT_FIG4B" not in (reader.fieldnames or []):
            raise SystemExit(f"{path.name} missing NOT_FIG4B")
        for row in reader:
            t.append(float(row["t"]))
            i_mean.append(float(row["I_mean"]))
    import numpy as np

    return np.array(t), np.array(i_mean)


def main() -> None:
    refuse_narma()
    if "DaninoSI_OccupiedDF" not in CLAIM.read_text(encoding="utf-8"):
        raise SystemExit("claim freeze missing Object B")
    standing = C1_STANDING.read_text(encoding="utf-8")
    if "**Status: FAIL**" not in standing:
        raise SystemExit("C1 standing must remain FAIL")
    c1b = C1B_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in c1b:
        raise SystemExit("C1b island standing must remain PASS")
    proto = PROTOCOL.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or "C1c_FILLED_POCKET" not in proto:
        raise SystemExit("C1c PROTOCOL not frozen")
    if "TIME_ADJ" not in proto:
        raise SystemExit("C1c PROTOCOL must freeze TIME_ADJ unused")
    summary = json.loads((RESULTS / "c1c_summary.json").read_text(encoding="utf-8"))
    print(f"C1={summary.get('C1')} C1b={summary.get('C1b')} C1c={summary.get('C1c')} NOT_FIG4B")
    print(
        f"open_chip_still_off={summary.get('open_chip_still_off')} "
        f"a0_on_open_dilute={summary.get('a0_on_open_dilute')} "
        f"a0_may_start_on_c1c_filled_pocket={summary.get('a0_may_start_on_c1c_filled_pocket')}"
    )

    ok = True
    t, i_mean = load_mean_csv(RESULTS / "java_C1c_FILLED_POCKET_D05_MU040.csv")
    peaks = extract_period(t, i_mean, PROTOCOL_PEAKS)
    rel = (
        abs(float(peaks["period"]) - PERIOD_REF) / PERIOD_REF
        if peaks["flag"] == "OSC"
        else float("nan")
    )
    ident_ok = peaks["flag"] == "OSC" and rel <= PERIOD_REL_BAND
    ok = ok and ident_ok
    print(
        f"  C1c_FILLED_POCKET_D05_MU040 NOT_FIG4B java={peaks['flag']} "
        f"T={peaks.get('period')} relT={rel} gate={'PASS' if ident_ok else 'FAIL'}"
    )

    arms = {a["name"]: a for a in summary.get("arms", [])}
    he = float(arms["C1c_FILLED_POCKET_D05_MU040"]["he_pocket_range"])
    he_ok = he >= HE_RANGE_MIN
    ok = ok and he_ok
    print(f"  He_space_used NOT_FIG4B pocket_range={he:.4g} gate={'PASS' if he_ok else 'FAIL'}")

    for name, required in (
        ("C1c_WRONG_KICK", "NO_PERIOD"),
        ("C1c_MU150", "NO_PERIOD"),
    ):
        tt, yy = load_mean_csv(RESULTS / f"java_{name}.csv")
        p = extract_period(tt, yy, PROTOCOL_PEAKS)
        arm_ok = p["flag"] == required
        ok = ok and arm_ok
        print(f"  {name} NOT_FIG4B java={p['flag']} gate={'PASS' if arm_ok else 'FAIL'}")

    rel_l = float(arms["C1c_MU0_LEDGER"]["ledger_rel"])
    ledger_ok = rel_l <= LEDGER_REL and bool(arms["C1c_MU0_LEDGER"]["pass"])
    ok = ok and ledger_ok
    print(f"  C1c_MU0_LEDGER NOT_FIG4B ledger_rel={rel_l:.4g} gate={'PASS' if ledger_ok else 'FAIL'}")

    open_ok = arms["C1_OPEN_DILUTE"]["flag"] == "NO_PERIOD"
    ok = ok and open_ok
    print(
        f"  C1_OPEN_DILUTE NOT_FIG4B report-only {arms['C1_OPEN_DILUTE']['flag']} "
        f"gate={'PASS' if open_ok else 'FAIL'}"
    )

    extra = arms.get("C1c_BUS_CONNECTED", {})
    print(
        f"  C1c_BUS_CONNECTED NOT_FIG4B not identity flag={extra.get('flag')} "
        f"(report-only; predeclared NO_PERIOD)"
    )

    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(t, i_mean, color="black", lw=1.2)
        ax.set_xlabel("SI-scaled time (NOT_FIG4B)")
        ax.set_ylabel("colony-mean I (LuxI)")
        ax.set_title("C1c_FILLED_POCKET_D05_MU040 colony-mean I  NOT_FIG4B  not Fig. 4b")
        fig.tight_layout()
        out = RESULTS / "c1c_filled_pocket_d05_mu040_I_mean_NOT_FIG4B.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        print(f"  plot {out.name} NOT_FIG4B")
    except Exception as exc:
        print(f"  plot skipped ({exc}) NOT_FIG4B")

    if summary.get("C1") != "FAIL":
        print("C1 standing/summary must remain FAIL")
        ok = False
    if not ok or summary.get("C1c") != "PASS":
        print(
            "C1c FAIL NOT_FIG4B. Do not retune Object B / D1_spatial / mu / d. "
            "Do not start A0/W0. Do not occupy C1_OPEN_DILUTE."
        )
        sys.exit(1)
    print(
        "C1c PASS on filled-pocket occupancy. NOT_FIG4B. "
        "C1 open chip remains FAIL. A0 may start only as an AC on C1c_FILLED_POCKET, "
        "still NOT_FIG4B, not on C1_OPEN_DILUTE. W0 later."
    )


if __name__ == "__main__":
    main()
