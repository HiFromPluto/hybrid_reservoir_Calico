#!/usr/bin/env python3
"""NARMA-blind C1 checker. Packed spatial Object B. NOT_FIG4B."""

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
SYNC_FRAC_OSC = 0.90
SYNC_WINDOW = 5.0
SYNC_FRAC_WINDOW = 0.90
HE_RANGE_MIN = 1e-6


def load_mean_csv(path: Path):
    t = []
    i_mean = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "NOT_FIG4B" not in fieldnames:
            raise SystemExit(f"{path.name} missing NOT_FIG4B column")
        for row in reader:
            t.append(float(row["t"]))
            i_mean.append(float(row["I_mean"]))
    import numpy as np

    return np.array(t), np.array(i_mean)


def load_sync_csv(path: Path):
    flags = []
    last = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            flags.append(row["flag"])
            raw = row["last_peak_t"]
            last.append(float("nan") if raw == "nan" else float(raw))
    return flags, last


def main() -> None:
    refuse_narma()
    if not CLAIM.exists():
        raise SystemExit("CLAIM_FREEZE_DANINO_SI.md missing")
    text = CLAIM.read_text(encoding="utf-8")
    if "DaninoSI_OccupiedDF" not in text or "NOT_FIG4B" not in text:
        raise SystemExit("claim freeze must name Object B NOT_FIG4B")
    proto = PROTOCOL.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or "NOT_FIG4B" not in proto:
        raise SystemExit("C1 PROTOCOL must be frozen and NOT_FIG4B")
    if "TIME_ADJ" not in proto or "SYNC_PEAK_WINDOW" not in proto:
        raise SystemExit("C1 PROTOCOL must freeze TIME_ADJ and SYNC_PEAK_WINDOW")
    summary_path = RESULTS / "c1_summary.json"
    if not summary_path.exists():
        raise SystemExit("run ant c1-packed first")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    print(f"object={summary.get('object')} status={summary.get('object_status')}")
    print(
        f"D0={summary.get('D0')} D0b={summary.get('D0b')} D1={summary.get('D1')} "
        f"P0={summary.get('P0')} C0={summary.get('C0')} C1={summary.get('C1')} NOT_FIG4B"
    )
    print(
        f"N={summary.get('n_cells')} v_cell={summary.get('v_cell_um3')} "
        f"Ve_local={summary.get('Ve_local_um3')} d={summary.get('d')} "
        f"D1_spatial={summary.get('D1_spatial')} {summary.get('D1_spatial_label')} NOT_FIG4B"
    )

    ok = True
    t, i_mean = load_mean_csv(RESULTS / "java_C1_HELD_PACK_D05_MU040.csv")
    peaks = extract_period(t, i_mean, PROTOCOL_PEAKS)
    rel = (
        abs(float(peaks["period"]) - PERIOD_REF) / PERIOD_REF
        if peaks["flag"] == "OSC"
        else float("nan")
    )
    ident_ok = peaks["flag"] == "OSC" and rel <= PERIOD_REL_BAND
    ok = ok and ident_ok
    print(
        f"  C1_HELD_PACK_D05_MU040 NOT_FIG4B java={peaks['flag']} "
        f"T={peaks.get('period')} relT={rel} gate={'PASS' if ident_ok else 'FAIL'}"
    )

    flags, last = load_sync_csv(RESULTS / "sync_C1_HELD_PACK_D05_MU040.csv")
    n = len(flags)
    osc_last = [x for f, x in zip(flags, last) if f == "OSC" and x == x]
    frac_osc = sum(1 for f in flags if f == "OSC") / n
    osc_sorted = sorted(osc_last)
    median = osc_sorted[len(osc_sorted) // 2] if osc_sorted else float("nan")
    frac_win = (
        sum(1 for x in osc_last if abs(x - median) <= SYNC_WINDOW) / len(osc_last)
        if osc_last
        else 0.0
    )
    sync_ok = frac_osc >= SYNC_FRAC_OSC and frac_win >= SYNC_FRAC_WINDOW
    ok = ok and sync_ok
    print(
        f"  SYNC_PEAK_WINDOW NOT_FIG4B frac_osc={frac_osc:.3f} "
        f"frac_win={frac_win:.3f} gate={'PASS' if sync_ok else 'FAIL'}"
    )

    arms = {a["name"]: a for a in summary.get("arms", [])}
    he = float(arms["C1_HELD_PACK_D05_MU040"]["he_fluid_range"])
    he_ok = he >= HE_RANGE_MIN
    ok = ok and he_ok
    print(
        f"  He_nonuniform NOT_FIG4B fluid_range={he:.4g} "
        f"patch_range={arms['C1_HELD_PACK_D05_MU040']['he_patch_range']} "
        f"gate={'PASS' if he_ok else 'FAIL'}"
    )

    t15, i15 = load_mean_csv(RESULTS / "java_C1_HELD_PACK_MU150.csv")
    p15 = extract_period(t15, i15, PROTOCOL_PEAKS)
    ok15 = p15["flag"] == "NO_PERIOD"
    ok = ok and ok15
    print(f"  C1_HELD_PACK_MU150 NOT_FIG4B java={p15['flag']} gate={'PASS' if ok15 else 'FAIL'}")

    tb, ib = load_mean_csv(RESULTS / "java_C1_HELD_PACK_WRONG_KICK.csv")
    pb = extract_period(tb, ib, PROTOCOL_PEAKS)
    okb = pb["flag"] == "NO_PERIOD"
    ok = ok and okb
    print(f"  C1_HELD_PACK_WRONG_KICK NOT_FIG4B java={pb['flag']} gate={'PASS' if okb else 'FAIL'}")

    rel_l = float(arms["C1_HELD_PACK_MU0_LEDGER"]["ledger_rel"])
    ledger_ok = rel_l <= LEDGER_REL and bool(arms["C1_HELD_PACK_MU0_LEDGER"]["pass"])
    ok = ok and ledger_ok
    print(
        f"  C1_HELD_PACK_MU0_LEDGER NOT_FIG4B ledger_rel={rel_l:.4g} "
        f"gate={'PASS' if ledger_ok else 'FAIL'}"
    )

    tlo, ilo = load_mean_csv(RESULTS / "java_C1_HELD_PACK_MU040_LORES.csv")
    plo = extract_period(tlo, ilo, PROTOCOL_PEAKS)
    lores_ok = plo["flag"] == "OSC" and peaks["flag"] == "OSC"
    ok = ok and lores_ok
    print(
        f"  C1_HELD_PACK_MU040_LORES NOT_FIG4B java={plo['flag']} "
        f"identity={peaks['flag']} OSC_call_unchanged={lores_ok} "
        f"gate={'PASS' if lores_ok else 'FAIL'}"
    )

    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(t, i_mean, color="black", lw=1.2)
        ax.set_xlabel("SI-scaled time (NOT_FIG4B)")
        ax.set_ylabel("colony-mean I (LuxI)")
        ax.set_title("C1_HELD_PACK_D05_MU040 colony-mean I  NOT_FIG4B  not Fig. 4b")
        fig.tight_layout()
        out = RESULTS / "c1_held_pack_d05_mu040_I_mean_NOT_FIG4B.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        print(f"  plot {out.name} NOT_FIG4B")
    except Exception as exc:
        print(f"  plot skipped ({exc}) NOT_FIG4B")

    if summary.get("C1") == "PASS" and not ok:
        print("Java summary PASS disagrees with Python re-check.")
        ok = False
    if not ok or summary.get("C1") != "PASS":
        print("C1 FAIL NOT_FIG4B. Do not retune Object B / D1_spatial / d. Do not start A0/W0. No NARMA.")
        sys.exit(1)
    print(
        "C1 PASS on Object B packed spatial dynamics. NOT_FIG4B. "
        "A0 may start. W0 later. D0/D0b unchanged."
    )


if __name__ == "__main__":
    main()
