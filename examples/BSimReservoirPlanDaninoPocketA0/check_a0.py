#!/usr/bin/env python3
"""NARMA-blind A0 checker. A0_IDEAL_SOURCE on C1c_FILLED_POCKET. NOT_FIG4B."""

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
FIELD_HE_REL = 1e-9
DELTA_M_PULSE = 825.0


def load_csv(path: Path):
    rows = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        if "NOT_FIG4B" not in names or "A0_IDEAL_SOURCE" not in names:
            raise SystemExit(f"{path.name} missing NOT_FIG4B / A0_IDEAL_SOURCE")
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
    if "**Status: PASS**" not in C1B_STANDING.read_text(encoding="utf-8"):
        raise SystemExit("C1b island standing must remain PASS")
    c1c = C1C_STANDING.read_text(encoding="utf-8")
    if "**Status: PASS**" not in c1c or "C1c_FILLED_POCKET" not in c1c:
        raise SystemExit("C1c standing must remain filled-pocket PASS")
    proto = PROTOCOL.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or "A0_IDEAL_SOURCE" not in proto:
        raise SystemExit("A0 PROTOCOL not frozen")
    if "C1c_FILLED_POCKET" not in proto or "NOT_FIG4B" not in proto:
        raise SystemExit("A0 PROTOCOL must freeze device and NOT_FIG4B")
    summary = json.loads((RESULTS / "a0_summary.json").read_text(encoding="utf-8"))
    print(
        f"A0={summary.get('A0')} device={summary.get('device')} "
        f"A0_IDEAL_SOURCE NOT_FIG4B"
    )
    print(
        f"C1={summary.get('C1')} C1b={summary.get('C1b')} C1c={summary.get('C1c')} "
        f"a0_on_open_dilute={summary.get('a0_on_open_dilute')} "
        f"a1_may_start={summary.get('a1_may_start')}"
    )
    if summary.get("device") != "C1c_FILLED_POCKET":
        raise SystemExit("A0 device must be C1c_FILLED_POCKET")
    if summary.get("a0_on_open_dilute") is not False:
        raise SystemExit("A0 must not run on C1_OPEN_DILUTE")

    ok = True
    arms = {a["name"]: a for a in summary.get("arms", [])}

    for name in ("A0_FIELD_IMPULSE", "A0_FIELD_STEP"):
        a0 = load_csv(RESULTS / f"java_{name}.csv")
        n0 = load_csv(RESULTS / f"java_{name}_N0_ORACLE.csv")
        he_ac = col(a0, "He_ac")
        he_far = col(a0, "He_far")
        o_ac = col(n0, "He_ac")
        o_far = col(n0, "He_far")
        rel = float(
            max(
                (abs(he_ac - o_ac) / (abs(o_ac) + 1e-15)).max(),
                (abs(he_far - o_far) / (abs(o_far) + 1e-15)).max(),
            )
        )
        mass_rel = float(arms[name]["mass_rel"])
        arm_ok = rel <= FIELD_HE_REL and mass_rel <= LEDGER_REL and bool(arms[name]["pass"])
        ok = ok and arm_ok
        print(
            f"  {name} A0_IDEAL_SOURCE NOT_FIG4B He_vs_N0_rel={rel:.4g} "
            f"mass_rel={mass_rel:.4g} gate={'PASS' if arm_ok else 'FAIL'}"
        )

    t, i_mean = None, None
    rows_j0 = load_csv(RESULTS / "java_A0_J0_BACTERIA.csv")
    import numpy as np

    t = col(rows_j0, "t")
    i_mean = col(rows_j0, "I_mean")
    peaks = extract_period(t, i_mean, PROTOCOL_PEAKS)
    rel = (
        abs(float(peaks["period"]) - PERIOD_REF) / PERIOD_REF
        if peaks["flag"] == "OSC"
        else float("nan")
    )
    j0_ok = peaks["flag"] == "OSC" and rel <= PERIOD_REL_BAND and bool(arms["A0_J0_BACTERIA"]["pass"])
    ok = ok and j0_ok
    print(
        f"  A0_J0_BACTERIA A0_IDEAL_SOURCE NOT_FIG4B java={peaks['flag']} "
        f"T={peaks.get('period')} relT={rel} vs C1c T=61.25 / Object B T=63.58 "
        f"gate={'PASS' if j0_ok else 'FAIL'}"
    )

    rows_p = load_csv(RESULTS / "java_A0_J_PULSE_BACTERIA.csv")
    he_ac_p = col(rows_p, "He_ac")
    he_ac_0 = col(rows_j0, "He_ac")
    t_p = col(rows_p, "t")
    idx = int(np.argmax(t_p >= 1.0 - 1e-12))
    he_rise = bool(he_ac_p[idx] > he_ac_0[idx])
    cmd = float(arms["A0_J_PULSE_BACTERIA"]["cmd_mass"])
    mass_ok = abs(cmd - DELTA_M_PULSE) <= 1e-6 and float(arms["A0_J_PULSE_BACTERIA"]["mass_rel"]) <= LEDGER_REL
    pulse_ok = he_rise and mass_ok and bool(arms["A0_J_PULSE_BACTERIA"]["pass"])
    ok = ok and pulse_ok
    i_mean_p = col(rows_p, "I_mean")
    i_ac_p = col(rows_p, "I_ac")
    print(
        f"  A0_J_PULSE_BACTERIA A0_IDEAL_SOURCE NOT_FIG4B He_ac_rise={he_rise} "
        f"I_mean={float(i_mean_p.mean()):.6g} I_ac={float(i_ac_p.mean()):.6g} "
        f"cmd={cmd} gate={'PASS' if pulse_ok else 'FAIL'}"
    )

    open_ok = arms["C1_OPEN_DILUTE"]["flag"] == "FORBIDDEN"
    ok = ok and open_ok
    print(
        f"  C1_OPEN_DILUTE A0_IDEAL_SOURCE NOT_FIG4B forbidden "
        f"flag={arms['C1_OPEN_DILUTE']['flag']} gate={'PASS' if open_ok else 'FAIL'}"
    )

    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(t, i_mean, color="black", lw=1.2, label="A0_J0_BACTERIA I_mean")
        ax.plot(t_p, i_mean_p, color="tab:red", lw=1.0, alpha=0.8, label="A0_J_PULSE_BACTERIA I_mean")
        ax.set_xlabel("SI-scaled time (A0_IDEAL_SOURCE NOT_FIG4B)")
        ax.set_ylabel("colony-mean I (LuxI)")
        ax.set_title("A0 on C1c_FILLED_POCKET  A0_IDEAL_SOURCE  NOT_FIG4B  not Fig. 4b")
        ax.legend()
        fig.tight_layout()
        out = RESULTS / "a0_j0_vs_pulse_I_mean_NOT_FIG4B_A0_IDEAL_SOURCE.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        print(f"  plot {out.name} A0_IDEAL_SOURCE NOT_FIG4B")

        imp = load_csv(RESULTS / "java_A0_FIELD_IMPULSE.csv")
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(col(imp, "t"), col(imp, "He_ac"), color="tab:blue", label="He AC voxel")
        ax.plot(col(imp, "t"), col(imp, "He_far"), color="tab:orange", label="He far voxel")
        ax.set_xlabel("SI-scaled time (A0_IDEAL_SOURCE NOT_FIG4B)")
        ax.set_ylabel("He")
        ax.set_title("A0_FIELD_IMPULSE  A0_IDEAL_SOURCE  NOT_FIG4B  not Fig. 4b")
        ax.legend()
        fig.tight_layout()
        out = RESULTS / "a0_field_impulse_He_NOT_FIG4B_A0_IDEAL_SOURCE.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        print(f"  plot {out.name} A0_IDEAL_SOURCE NOT_FIG4B")
    except Exception as exc:
        print(f"  plot skipped ({exc}) A0_IDEAL_SOURCE NOT_FIG4B")

    if summary.get("C1") != "FAIL":
        print("C1 standing/summary must remain FAIL")
        ok = False
    if not ok or summary.get("A0") != "PASS":
        print(
            "A0 FAIL A0_IDEAL_SOURCE NOT_FIG4B. Do not retune Object B / J_max. "
            "Do not start A1/A2/W0. Do not place the AC on C1_OPEN_DILUTE."
        )
        sys.exit(1)
    print(
        "A0 PASS A0_IDEAL_SOURCE NOT_FIG4B on C1c_FILLED_POCKET. "
        "C1 open chip remains FAIL. A1 may start only on C1c_FILLED_POCKET, "
        "still NOT_FIG4B, still not calibrated. Not Fig. 4b. W0 later."
    )


if __name__ == "__main__":
    main()
