#!/usr/bin/env python3
"""B1_PAINT_HOLD checker. LANE_B_CHEMOTACTIC_SPATIAL. Parent B0 FAIL kept.

Prints B1_PAINT_HOLD=PASS or FAIL with the gate that killed it.
Does not raise J_max. Does not start NARMA. Does not rewrite B0.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_PAINT_HOLD.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_b1.json"
B0_PROTOCOL = HERE / "PROTOCOL.md"

INV_E = 1.0 / math.e
T_OFF = 60.0
T_MIX = 180.0
T_HOLD = 275.3
KAPPA_MOTILE_MIN = 0.12
KAPPA_MINUS = 0.08
KAPPA_SILENT_MAX = 0.08
KAPPA_FIELD_HOLD_MAX = 0.15
LEDGER_REL = 1e-3


def refuse_capacity() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("narma", "charc", "ipc", "fig4b")):
        raise SystemExit("B1 checker refuses NARMA/CHARC/IPC/Fig4b")


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def col(rows: list[dict[str, str]], name: str) -> list[float]:
    out = []
    for row in rows:
        raw = row[name]
        out.append(float("nan") if raw in ("NaN", "nan", "") else float(raw))
    return out


def interpolate(t: list[float], y: list[float], tq: float) -> float:
    if not t:
        return float("nan")
    if tq <= t[0]:
        return y[0]
    for i in range(1, len(t)):
        if t[i] < t[i - 1]:
            break
        if tq <= t[i]:
            den = t[i] - t[i - 1] + 1e-18
            f = (tq - t[i - 1]) / den
            return y[i - 1] + f * (y[i] - y[i - 1])
    return y[-1]


def tau_mix(t: list[float], corr: list[float]) -> float:
    target = INV_E
    for i in range(1, len(t)):
        if i > 1 and t[i] <= t[i - 1]:
            break
        c0, c1 = corr[i - 1], corr[i]
        if not math.isfinite(c0) or not math.isfinite(c1):
            continue
        if c0 >= target and c1 <= target:
            f = (c0 - target) / (c0 - c1 + 1e-18)
            return t[i - 1] + f * (t[i] - t[i - 1])
    return float("nan")


def main() -> int:
    refuse_capacity()
    proto = PROTOCOL.read_text(encoding="utf-8")
    js_text = PROTOCOL_JSON.read_text(encoding="utf-8")
    js = json.loads(js_text)
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("B1 PROTOCOL is not frozen_before_traces")
    if js.get("status_label") != "B1_PAINT_HOLD" or js.get("narma") is not False:
        raise SystemExit("B1 protocol must freeze B1_PAINT_HOLD and narma=false")
    if js.get("t_off_s") != T_OFF or js.get("J_max") != 20000.0:
        raise SystemExit("B1 protocol must freeze t_off=60 and J_max=20000")
    if js.get("T_hold_s") != T_HOLD:
        raise SystemExit(f"B1 protocol must freeze scout T_hold={T_HOLD}")
    if js.get("parent_fail") != "B0_MOTILITY_MEMORY":
        raise SystemExit("B1 must point at parent B0 FAIL")
    b0 = B0_PROTOCOL.read_text(encoding="utf-8")
    if "B0_MOTILITY_MEMORY" not in b0 or "**10 s**" not in b0:
        raise SystemExit("B0 PROTOCOL.md must remain the 10 s parent")

    mix_path = RESULTS / "b1_mix.csv"
    if not mix_path.exists():
        print("B1_PAINT_HOLD=FAIL mix (missing b1_mix.csv)")
        return 1
    mix = load_rows(mix_path)
    t_mix = col(mix, "t")
    corr = col(mix, "auto_corr")
    tau = tau_mix(t_mix, corr)
    mix_ok = math.isfinite(tau) and 0.0 < tau < T_MIX
    print(f"B1.1 stripe mix tau_mix={tau:.6g} s {('PASS' if mix_ok else 'FAIL')}")
    if not mix_ok:
        print("B1_PAINT_HOLD=FAIL mix")
        return 1

    needed = [
        "b1_paint_motile.csv",
        "b1_paint_off.csv",
        "b1_paint_field.csv",
        "b1_paint_silent.csv",
    ]
    missing = [n for n in needed if not (RESULTS / n).exists()]
    if missing:
        print(f"B1_PAINT_HOLD=FAIL occupation (missing {missing})")
        return 1

    motile = load_rows(RESULTS / "b1_paint_motile.csv")
    off = load_rows(RESULTS / "b1_paint_off.csv")
    field = load_rows(RESULTS / "b1_paint_field.csv")
    silent = load_rows(RESULTS / "b1_paint_silent.csv")

    k_mot = interpolate(col(motile, "t"), col(motile, "kappa_cell"), T_OFF)
    k_off = interpolate(col(off, "t"), col(off, "kappa_cell"), T_OFF)
    k_sil = interpolate(col(silent, "t"), col(silent, "kappa_cell"), T_OFF)
    k_fld = interpolate(col(field, "t"), col(field, "kappa_field"), T_OFF)
    occ_ok = (
        math.isfinite(k_mot)
        and math.isfinite(k_off)
        and math.isfinite(k_sil)
        and abs(k_mot) >= KAPPA_MOTILE_MIN
        and abs(k_mot) - abs(k_off) >= KAPPA_MINUS
        and abs(k_sil) < KAPPA_SILENT_MAX
    )
    print(
        f"B1.2 occupation t={T_OFF:g} kappa_motile={k_mot:.6g} kappa_off={k_off:.6g} "
        f"kappa_silent={k_sil:.6g} kappa_field={k_fld:.6g} {('PASS' if occ_ok else 'FAIL')}"
    )
    if not occ_ok:
        print("B1.3 hold NOT_SCORED (occupation FAIL)")
        print("B1_PAINT_HOLD=FAIL occupation")
        return 1

    k_mot_h = interpolate(col(motile, "t"), col(motile, "kappa_cell"), T_HOLD)
    k_off_h = interpolate(col(off, "t"), col(off, "kappa_cell"), T_HOLD)
    k_fld_h = interpolate(col(field, "t"), col(field, "kappa_field"), T_HOLD)
    hold_ok = (
        math.isfinite(k_mot_h)
        and math.isfinite(k_off_h)
        and math.isfinite(k_fld_h)
        and abs(k_mot_h) - abs(k_fld_h) >= KAPPA_MINUS
        and abs(k_mot_h) - abs(k_off_h) >= KAPPA_MINUS
    )
    fade_ok = math.isfinite(k_fld_h) and abs(k_fld_h) < KAPPA_FIELD_HOLD_MAX
    print(
        f"B1.3 hold T_hold={T_HOLD:g} kappa_motile={k_mot_h:.6g} kappa_off={k_off_h:.6g} "
        f"kappa_field={k_fld_h:.6g} {('PASS' if hold_ok else 'FAIL')}"
    )
    if not fade_ok:
        print(
            f"SCOPE_NOTE scout vs Java: |kappa_field|(T_hold)={k_fld_h:.6g} not < 0.15. "
            "Do not retune D."
        )
    if not hold_ok:
        print("B1_PAINT_HOLD=FAIL hold")
        return 1

    summary = json.loads((RESULTS / "b1_summary.json").read_text(encoding="utf-8"))
    ledger_ok = summary.get("B1.4_ledger") is True
    if summary.get("narma") is not False or summary.get("b0_overall_rewrite") is not False:
        print("B1_PAINT_HOLD=FAIL ledger (summary rewrites NARMA/B0)")
        return 1
    print(f"B1.4 N0 ledger {'PASS' if ledger_ok else 'FAIL'}")
    if not ledger_ok:
        print("B1_PAINT_HOLD=FAIL ledger")
        return 1

    pulse_l = summary.get("mean_abs_L_pulse", float("nan"))
    print(f"REPORT mean_|L|_pulse={pulse_l} (not a gate)")
    print("B1_PAINT_HOLD=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
