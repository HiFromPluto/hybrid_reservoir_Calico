#!/usr/bin/env python3
"""B2_STRIPE_OCCUPY checker. LANE_B_CHEMOTACTIC_SPATIAL. Parent B0/B1 FAILs kept.

Prints B2_STRIPE_OCCUPY=PASS or FAIL with the gate that killed it.
Does not raise L_slab. Does not drop D. Does not start NARMA. Does not rewrite B0/B1.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "PROTOCOL_STRIPE_OCCUPY.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_b2.json"
B0_PROTOCOL = HERE / "PROTOCOL.md"
B1_PROTOCOL = HERE / "PROTOCOL_PAINT_HOLD.md"

T_OCC = 1.1
T_HOLD = 3.1
KAPPA_MOTILE_MIN = 0.12
KAPPA_MINUS = 0.08
KAPPA_SILENT_MAX = 0.08
KAPPA_FIELD_HOLD_MAX = 0.15


def refuse_capacity() -> None:
    blob = " ".join(sys.argv).lower()
    if any(k in blob for k in ("narma", "charc", "ipc", "fig4b")):
        raise SystemExit("B2 checker refuses NARMA/CHARC/IPC/Fig4b")


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


def main() -> int:
    refuse_capacity()
    proto = PROTOCOL.read_text(encoding="utf-8")
    js_text = PROTOCOL_JSON.read_text(encoding="utf-8")
    js = json.loads(js_text)
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        raise SystemExit("B2 PROTOCOL is not frozen_before_traces")
    if js.get("status_label") != "B2_STRIPE_OCCUPY" or js.get("narma") is not False:
        raise SystemExit("B2 protocol must freeze B2_STRIPE_OCCUPY and narma=false")
    if js.get("t_occ_s") != T_OCC or js.get("T_hold_s") != T_HOLD:
        raise SystemExit(f"B2 protocol must freeze scout t_occ={T_OCC} T_hold={T_HOLD}")
    if js.get("D_att_um2_s") != 800.0:
        raise SystemExit("B2 protocol must freeze D=800")
    if js.get("parent_fail") != "B0_MOTILITY_MEMORY" or js.get("parent_fail_b1") != "B1_PAINT_HOLD":
        raise SystemExit("B2 must point at parent B0 and B1 FAILs")
    b0 = B0_PROTOCOL.read_text(encoding="utf-8")
    if "B0_MOTILITY_MEMORY" not in b0 or "**10 s**" not in b0:
        raise SystemExit("B0 PROTOCOL.md must remain the 10 s parent")
    b1 = B1_PROTOCOL.read_text(encoding="utf-8")
    if "B1_PAINT_HOLD" not in b1 or "**275.3 s**" not in b1:
        raise SystemExit("B1 PROTOCOL_PAINT_HOLD.md must remain the 275.3 s parent")

    needed = [
        "b2_stripe_motile.csv",
        "b2_stripe_off.csv",
        "b2_stripe_field.csv",
        "b2_stripe_silent.csv",
    ]
    missing = [n for n in needed if not (RESULTS / n).exists()]
    if missing:
        print(f"B2_STRIPE_OCCUPY=FAIL occupation (missing {missing})")
        return 1

    motile = load_rows(RESULTS / "b2_stripe_motile.csv")
    off = load_rows(RESULTS / "b2_stripe_off.csv")
    field = load_rows(RESULTS / "b2_stripe_field.csv")
    silent = load_rows(RESULTS / "b2_stripe_silent.csv")

    k_mot = interpolate(col(motile, "t"), col(motile, "kappa_cell"), T_OCC)
    k_off = interpolate(col(off, "t"), col(off, "kappa_cell"), T_OCC)
    k_sil = interpolate(col(silent, "t"), col(silent, "kappa_cell"), T_OCC)
    k_fld = interpolate(col(field, "t"), col(field, "kappa_field"), T_OCC)
    occ_ok = (
        math.isfinite(k_mot)
        and math.isfinite(k_off)
        and math.isfinite(k_sil)
        and abs(k_mot) >= KAPPA_MOTILE_MIN
        and abs(k_mot) - abs(k_off) >= KAPPA_MINUS
        and abs(k_sil) < KAPPA_SILENT_MAX
    )
    print(
        f"B2.1 occupation t={T_OCC:g} kappa_motile={k_mot:.6g} kappa_off={k_off:.6g} "
        f"kappa_silent={k_sil:.6g} kappa_field={k_fld:.6g} {('PASS' if occ_ok else 'FAIL')}"
    )
    if not occ_ok:
        print("B2.2 hold NOT_SCORED (occupation FAIL)")
        print("B2_STRIPE_OCCUPY=FAIL occupation")
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
        f"B2.2 hold T_hold={T_HOLD:g} kappa_motile={k_mot_h:.6g} kappa_off={k_off_h:.6g} "
        f"kappa_field={k_fld_h:.6g} {('PASS' if hold_ok else 'FAIL')}"
    )
    if not fade_ok:
        print(
            f"SCOPE_NOTE scout vs Java: |kappa_lambda_field|(T_hold)={k_fld_h:.6g} not < 0.15. "
            "Do not retune D."
        )
    if not hold_ok:
        print("B2_STRIPE_OCCUPY=FAIL hold")
        return 1

    summary = json.loads((RESULTS / "b2_summary.json").read_text(encoding="utf-8"))
    ledger_ok = summary.get("B2.3_ledger") is True
    if summary.get("narma") is not False or summary.get("b0_overall_rewrite") is not False:
        print("B2_STRIPE_OCCUPY=FAIL ledger (summary rewrites NARMA/B0)")
        return 1
    print(f"B2.3 N0 ledger {'PASS' if ledger_ok else 'FAIL'}")
    if not ledger_ok:
        print("B2_STRIPE_OCCUPY=FAIL ledger")
        return 1

    dL = summary.get("delta_L_ell_t_occ", float("nan"))
    print(f"REPORT DeltaL(ell) at t_occ={dL} (not a gate)")
    print("B2_STRIPE_OCCUPY=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
