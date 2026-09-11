#!/usr/bin/env python3
"""D1_MEASURED_DECLINE. Metrics on frozen d0_grow.csv. Not a gate. Not Erickson."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
GROW = RESULTS / "d0_grow.csv"
SUMMARY = RESULTS / "d0_summary.json"
PROTOCOL = HERE / "PROTOCOL_D1.md"
PROTOCOL_JSON = HERE / "configs" / "protocol_d1.json"

K_S = 0.02
C_CUT = 1.0e-3
D0_T_END = 10455.0
D0_C_END = 9.977e-4
D0_N_END = 64
C_END_TOL = 5.0e-7

WINDOWS = (
    ("W_19K", "c_lt", 19.0 * K_S, 0.816, "lam<0.95"),
    ("W_9K", "c_lt", 9.0 * K_S, 0.444, "lam<0.90"),
    ("W_3K", "c_lt", 3.0 * K_S, 0.163, "lam<0.75"),
    ("W_K", "c_lt", 1.0 * K_S, 0.057, "lam<0.50"),
    ("W_POST3", "n_eq", 64.0, None, "after last fission; D0 structure, not a threshold"),
)

ISOLATION_PATHS = [
    "examples/LaneD_LifeCycle/PROTOCOL_D0.md",
    "examples/LaneD_LifeCycle/configs/protocol_d0.json",
    "examples/LaneD_LifeCycle/results/d0_grow.csv",
    "examples/LaneD_LifeCycle/results/d0_summary.json",
    "examples/LaneD_LifeCycle/results/d0_off.csv",
    "examples/PocketDish/LANE_D_D0_CLOSED_BATH_STANDING.md",
    "src/bsim/laned",
    "examples/BacteriumFromScratch/NutrientField.java",
    "examples/BacteriumFromScratch/EcoliRodCell.java",
    "examples/BacteriumFromScratch/ChassisParameters.java",
    "examples/BacteriumFromScratch/JOB2_STANDING.md",
    "examples/BacteriumFromScratch/JOB3_STANDING.md",
    "examples/BacteriumFromScratch/JOB3B_STANDING.md",
    "examples/BacteriumFromScratch/JOB3C_STANDING.md",
    "examples/BacteriumFromScratch/JOB6_STANDING.md",
    "examples/LaneC_LivingClocks",
    "src/bsim/lanec",
    "src/bsim/lanea",
    "src/bsim/laneb",
    "examples/PocketDish/manuscript_ieee",
]


def main() -> int:
    if any(k in " ".join(sys.argv).lower() for k in ("narma", "charc", "ipc")):
        raise SystemExit("D1 checker refuses NARMA/CHARC/IPC")
    if any(k in " ".join(sys.argv).lower() for k in ("erickson", "fcr", "d2", "famine")):
        raise SystemExit("D1 checker refuses Erickson/FCR/D2/famine")

    proto = PROTOCOL.read_text(encoding="utf-8")
    js = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    if "frozen_before_traces" not in proto or js.get("frozen_before_traces") is not True:
        return incomplete("honesty", "D1 PROTOCOL is not frozen")
    if js.get("status_label") != "D1_MEASURED_DECLINE":
        return incomplete("honesty", "status_label must be D1_MEASURED_DECLINE")
    if js.get("device") != "LANE_D_CLOSED_BATH" or js.get("not_device") != "LC0_OPEN_BATH":
        return incomplete("honesty", "device must be LANE_D_CLOSED_BATH")
    if js.get("death") is not False or js.get("erickson") is not False:
        return incomplete("honesty", "death and erickson must be false")
    if js.get("k_s_mM") != K_S or js.get("resample_1s") is not False:
        return incomplete("honesty", "Warren K_S=0.02; no 1 s resample")
    if js.get("new_java") is not False or js.get("new_trace") is not False:
        return incomplete("honesty", "no new Java or new C(t)")
    if js.get("report_not_gate") is not True:
        return incomplete("honesty", "D1 is report-not-gate")

    print(
        "D1_MEASURED_DECLINE death=OFF device=LANE_D_CLOSED_BATH "
        f"K_S={K_S} mM new_java=false new_trace=false "
        "frozen_before_traces=true"
    )

    if not GROW.is_file():
        return incomplete(
            "reuse",
            "d0_grow.csv missing. Rerun D0 unchanged. Do not raise C_s or N0.",
        )
    if not SUMMARY.is_file():
        return incomplete("reuse", "d0_summary.json missing. Do not fix D0.")

    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    rows = load(GROW)
    if not rows:
        return incomplete("reuse", "d0_grow.csv is empty")

    last = rows[-1]
    t_end = float(summary["grow_t_end"])
    c_end = float(summary["grow_C_end"])
    n_end = int(summary["grow_N_end"])
    csv_t = float(last["t"])
    csv_c = float(last["C_mM"])
    csv_n = int(last["N"])

    fingerprint = (
        abs(t_end - D0_T_END) <= 1.0e-6
        and abs(csv_t - D0_T_END) <= 1.0e-6
        and n_end == D0_N_END
        and csv_n == D0_N_END
        and abs(c_end - D0_C_END) <= C_END_TOL
        and abs(csv_c - D0_C_END) <= C_END_TOL
    )
    if not fingerprint:
        print(
            f"D1.2 reuse FAIL t_end={t_end} C_end={c_end} N_end={n_end} "
            f"csv=({csv_t},{csv_c},{csv_n}) "
            f"expected=({D0_T_END},{D0_C_END},{D0_N_END})"
        )
        return incomplete("reuse", "end-row fingerprint drifted. Stop. Do not fix D0.")

    end_i = first_index(rows, lambda r: float(r["C_mM"]) <= C_CUT)
    if end_i is None:
        return incomplete("completeness", "no row with C <= C_cut")

    reports = []
    for wid, kind, thresh, budget, label in WINDOWS:
        if kind == "c_lt":
            start_i = first_index(rows, lambda r, th=thresh: float(r["C_mM"]) < th)
        else:
            start_i = first_index(rows, lambda r: int(r["N"]) == int(thresh))
        if start_i is None or start_i > end_i:
            return incomplete("completeness", f"{wid} start row missing")
        reports.append(measure(wid, rows, start_i, end_i, budget, label))

    isolation = isolation_ok()

    print(f"D1.1 honesty {'PASS' if True else 'FAIL'}")
    print(
        f"D1.2 reuse t_end={t_end:.0f} C_end={c_end:.6e} N_end={n_end} PASS"
    )
    print(f"D1.3 completeness windows={len(reports)} PASS")
    print(f"D1.4 isolation {'PASS' if isolation else 'FAIL'}")

    print(
        "window  t0_s      C(t0)_mM     lam(t0)   t1_s      "
        "D_lambda  count_dbl  vol_ratio  budget_est"
    )
    for rep in reports:
        budget = "na" if rep["budget_estimate"] is None else f"{rep['budget_estimate']:.3f}"
        print(
            f"{rep['id']:<7} {rep['t0_s']:8.1f}  {rep['C_t0_mM']:.6e}  "
            f"{rep['lambda_t0_per_h']:.6f}  {rep['t1_s']:8.1f}  "
            f"{rep['D_lambda']:.6f}  {rep['count_doublings']:.6f}  "
            f"{rep['volume_ratio_note']:.6f}  {budget}"
        )
        print(
            f"  {rep['id']} note: F3 volume-ratio is not proteome. "
            f"section_32={rep['section_32']}"
        )

    write_outputs(reports, isolation)
    if not isolation:
        return incomplete("isolation", "protected paths are dirty")

    print("D1_MEASURED_DECLINE=REPORT")
    print(
        "Measured D_lambda on the frozen D0 trace. "
        "No band was chosen. Not Erickson. Not a dish PASS."
    )
    return 0


def measure(
    wid: str,
    rows: list[dict[str, str]],
    start_i: int,
    end_i: int,
    budget: float | None,
    label: str,
) -> dict[str, object]:
    t0 = float(rows[start_i]["t"])
    t1 = float(rows[end_i]["t"])
    c0 = float(rows[start_i]["C_mM"])
    c1 = float(rows[end_i]["C_mM"])
    lam0 = float(rows[start_i]["lambda_per_h"])
    n0 = int(rows[start_i]["N"])
    n1 = int(rows[end_i]["N"])
    v0 = float(rows[start_i]["sum_V_um3"])
    v1 = float(rows[end_i]["sum_V_um3"])
    integral = trapezoid_hours(rows, start_i, end_i)
    d_lambda = integral / math.log(2.0)
    count = math.log2(n1 / n0) if n0 > 0 and n1 > 0 else float("nan")
    vol = math.log2(v1 / v0) if v0 > 0 and v1 > 0 else float("nan")
    return {
        "id": wid,
        "t0_s": t0,
        "C_t0_mM": c0,
        "lambda_t0_per_h": lam0,
        "t1_s": t1,
        "C_t1_mM": c1,
        "D_lambda": d_lambda,
        "count_doublings": count,
        "volume_ratio_note": vol,
        "N_t0": n0,
        "N_t1": n1,
        "sum_V_t0": v0,
        "sum_V_t1": v1,
        "budget_estimate": budget,
        "section_32": label,
        "integral_hours": integral,
    }


def trapezoid_hours(rows: list[dict[str, str]], start_i: int, end_i: int) -> float:
    total = 0.0
    for i in range(start_i, end_i):
        t0 = float(rows[i]["t"])
        t1 = float(rows[i + 1]["t"])
        lam0 = float(rows[i]["lambda_per_h"])
        lam1 = float(rows[i + 1]["lambda_per_h"])
        total += 0.5 * (lam0 + lam1) * (t1 - t0) / 3600.0
    return total


def first_index(rows: list[dict[str, str]], pred) -> int | None:
    for i, row in enumerate(rows):
        if pred(row):
            return i
    return None


def write_outputs(reports: list[dict[str, object]], isolation: bool) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / "d1_windows.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "id",
                "t0_s",
                "C_t0_mM",
                "lambda_t0_per_h",
                "t1_s",
                "C_t1_mM",
                "D_lambda",
                "count_doublings",
                "volume_ratio_note",
                "N_t0",
                "N_t1",
                "sum_V_t0",
                "sum_V_t1",
                "budget_estimate",
                "section_32",
            ],
        )
        writer.writeheader()
        for rep in reports:
            row = {k: rep[k] for k in writer.fieldnames}
            writer.writerow(row)
    payload = {
        "status_label": "D1_MEASURED_DECLINE",
        "object": "LANE_D_LIFE_CYCLE",
        "device": "LANE_D_CLOSED_BATH",
        "report_not_gate": True,
        "death": False,
        "erickson": False,
        "new_java": False,
        "new_trace": False,
        "k_s_mM": K_S,
        "trace": "examples/LaneD_LifeCycle/results/d0_grow.csv",
        "windows": reports,
        "isolation": isolation,
        "hygiene": "ok" if isolation else "incomplete",
    }
    (RESULTS / "d1_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )


def isolation_ok() -> bool:
    existing = [p for p in ISOLATION_PATHS if (ROOT / p).exists()]
    if not existing:
        return False
    try:
        proc = subprocess.run(
            ["git", "diff", "HEAD", "--exit-code", "--", *existing],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return False
    if proc.returncode == 0:
        return True
    sys.stderr.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    return False


def incomplete(killer: str, detail: str) -> int:
    print(f"D1_MEASURED_DECLINE=INCOMPLETE {killer}")
    print(detail)
    return 1


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
