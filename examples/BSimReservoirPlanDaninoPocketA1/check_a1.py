#!/usr/bin/env python3
"""NARMA-blind A1 checker. A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE. NOT_FIG4B."""

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
A0_STANDING = ROOT / "examples" / "PocketDish" / "A0_IDEAL_SOURCE_STANDING.md"
C1C_STANDING = ROOT / "examples" / "PocketDish" / "C1C_FILLED_POCKET_STANDING.md"
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
A0_CONTRAST = 0.10
J_LEAK = 0.1
M0 = 5000.0


def load_csv(path: Path):
    rows = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if "NOT_FIG4B" not in names or "A1_BOUNDED_TRANSDUCER" not in names:
            raise SystemExit(f"{path.name} missing NOT_FIG4B / A1_BOUNDED_TRANSDUCER")
        if "HYPOTHETICAL_DESIGN_ENVELOPE" not in names:
            raise SystemExit(f"{path.name} missing HYPOTHETICAL_DESIGN_ENVELOPE")
        for row in reader:
            rows.append(row)
    return rows


def col(rows, key):
    import numpy as np

    return np.array([float(r[key]) for r in rows])


def main() -> None:
    refuse_narma()
    if "DaninoSI_OccupiedDF" not in CLAIM.read_text(encoding="utf-8"):
        raise SystemExit("claim freeze missing Object B")
    if "**Status: FAIL**" not in C1_STANDING.read_text(encoding="utf-8"):
        raise SystemExit("C1 standing must remain FAIL")
    if "**Status: PASS**" not in A0_STANDING.read_text(encoding="utf-8"):
        raise SystemExit("A0 standing must remain PASS")
    if "**Status: PASS**" not in C1C_STANDING.read_text(encoding="utf-8"):
        raise SystemExit("C1c standing must remain PASS")
    proto = PROTOCOL.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or "HYPOTHETICAL_DESIGN_ENVELOPE" not in proto:
        raise SystemExit("A1 PROTOCOL not frozen")
    summary = json.loads((RESULTS / "a1_summary.json").read_text(encoding="utf-8"))
    print(
        f"A1={summary.get('A1')} device={summary.get('device')} "
        f"A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B"
    )
    print(
        f"calibrated={summary.get('calibrated')} a1c_may_start={summary.get('a1c_may_start')} "
        f"a2_may_start={summary.get('a2_may_start')} a1_on_open_dilute={summary.get('a1_on_open_dilute')}"
    )
    if summary.get("device") != "C1c_FILLED_POCKET":
        raise SystemExit("A1 device must be C1c_FILLED_POCKET")
    if summary.get("calibrated") is not False or summary.get("a2_may_start") is not False:
        raise SystemExit("A1 must not be calibrated; A2 must stay closed")

    ok = True
    arms = {a["name"]: a for a in summary.get("arms", [])}

    for name in ("A1_FIELD_STEP_U", "A1_FIELD_PULSE_U"):
        rows = load_csv(RESULTS / f"java_{name}.csv")
        j_s = col(rows, "J_S")
        j_a0 = col(rows, "J_A0")
        contrast = float((abs(j_s - j_a0) / 825.0).max())
        mass_rel = float(arms[name]["mass_rel"])
        delayed = float(j_s[1]) < float(j_a0[1])
        arm_ok = contrast >= A0_CONTRAST and delayed and mass_rel <= LEDGER_REL and bool(arms[name]["pass"])
        ok = ok and arm_ok
        print(
            f"  {name} A1_BOUNDED_TRANSDUCER NOT_FIG4B A0_contrast={contrast:.4g} "
            f"delayed={delayed} mass_rel={mass_rel:.4g} gate={'PASS' if arm_ok else 'FAIL'}"
        )

    leak = load_csv(RESULTS / "java_A1_U0_LEAK.csv")
    mean_j = float(col(leak, "J_S")[1:].mean())
    leak_ok = abs(mean_j - J_LEAK) <= 0.01 and bool(arms["A1_U0_LEAK"]["pass"])
    ok = ok and leak_ok
    print(
        f"  A1_U0_LEAK A1_BOUNDED_TRANSDUCER NOT_FIG4B meanJ={mean_j:.6g} "
        f"J_leak={J_LEAK} gate={'PASS' if leak_ok else 'FAIL'}"
    )

    ex = load_csv(RESULTS / "java_A1_PAYLOAD_EXHAUST.csv")
    m_end = float(col(ex, "M")[-1])
    j_end = float(col(ex, "J_S")[-1])
    j_max = float(col(ex, "J_S").max())
    ex_ok = (m_end / M0) <= 0.10 and j_end <= 0.15 * j_max and bool(arms["A1_PAYLOAD_EXHAUST"]["pass"])
    ok = ok and ex_ok
    print(
        f"  A1_PAYLOAD_EXHAUST A1_BOUNDED_TRANSDUCER NOT_FIG4B M_frac={m_end / M0:.4g} "
        f"J_end={j_end:.6g} gate={'PASS' if ex_ok else 'FAIL'}"
    )

    rows_j0 = load_csv(RESULTS / "java_A1_J0_BACTERIA.csv")
    t = col(rows_j0, "t")
    i_mean = col(rows_j0, "I_mean")
    peaks = extract_period(t, i_mean, PROTOCOL_PEAKS)
    rel = (
        abs(float(peaks["period"]) - PERIOD_REF) / PERIOD_REF
        if peaks["flag"] == "OSC"
        else float("nan")
    )
    j0_ok = peaks["flag"] == "OSC" and rel <= PERIOD_REL_BAND and bool(arms["A1_J0_BACTERIA"]["pass"])
    ok = ok and j0_ok
    print(
        f"  A1_J0_BACTERIA A1_BOUNDED_TRANSDUCER NOT_FIG4B java={peaks['flag']} "
        f"T={peaks.get('period')} relT={rel} vs C1c T=61.25 gate={'PASS' if j0_ok else 'FAIL'}"
    )

    rows_p = load_csv(RESULTS / "java_A1_PULSE_BACTERIA.csv")
    import numpy as np

    t_p = col(rows_p, "t")
    idx = int(np.argmax(t_p >= 1.0 - 1e-12))
    he_rise = bool(col(rows_p, "He_ac")[idx] > col(rows_j0, "He_ac")[idx])
    pulse_ok = he_rise and bool(arms["A1_PULSE_BACTERIA"]["pass"])
    ok = ok and pulse_ok
    print(
        f"  A1_PULSE_BACTERIA A1_BOUNDED_TRANSDUCER NOT_FIG4B He_ac_rise={he_rise} "
        f"I_mean={float(col(rows_p, 'I_mean').mean()):.6g} "
        f"released={arms['A1_PULSE_BACTERIA']['released']} gate={'PASS' if pulse_ok else 'FAIL'}"
    )

    open_ok = arms["C1_OPEN_DILUTE"]["flag"] == "FORBIDDEN"
    ok = ok and open_ok
    print(
        f"  C1_OPEN_DILUTE A1_BOUNDED_TRANSDUCER NOT_FIG4B forbidden "
        f"gate={'PASS' if open_ok else 'FAIL'}"
    )

    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))
        pulse = load_csv(RESULTS / "java_A1_FIELD_PULSE_U.csv")
        ax.plot(col(pulse, "t"), col(pulse, "J_A0"), color="0.5", lw=1.0, label="A0 J")
        ax.plot(col(pulse, "t"), col(pulse, "J_S"), color="tab:blue", lw=1.2, label="A1 J_S")
        ax.set_xlabel("SI-scaled time (A1_BOUNDED_TRANSDUCER NOT_FIG4B)")
        ax.set_ylabel("J (quantity / SI-time)")
        ax.set_title("A1 vs A0 pulse  HYPOTHETICAL_DESIGN_ENVELOPE  NOT_FIG4B")
        ax.legend()
        fig.tight_layout()
        out = RESULTS / "a1_pulse_J_vs_A0_NOT_FIG4B_A1_BOUNDED_TRANSDUCER.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        print(f"  plot {out.name} A1_BOUNDED_TRANSDUCER NOT_FIG4B")
    except Exception as exc:
        print(f"  plot skipped ({exc}) A1_BOUNDED_TRANSDUCER NOT_FIG4B")

    if not ok or summary.get("A1") != "PASS":
        print(
            "A1 FAIL A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B. "
            "Do not retune. Do not start A1C/A2/W0."
        )
        sys.exit(1)
    print(
        "A1 PASS A1_BOUNDED_TRANSDUCER HYPOTHETICAL_DESIGN_ENVELOPE NOT_FIG4B on C1c_FILLED_POCKET. "
        "Not calibrated. A1C may not start. A2 may not start. W0 later."
    )


if __name__ == "__main__":
    main()
