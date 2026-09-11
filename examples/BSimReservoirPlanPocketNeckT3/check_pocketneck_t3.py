#!/usr/bin/env python3
"""PocketNeck-T3 occupancy checker.

Cannot see a NARMA target. No ridge, no Mackey-Glass, no waveform AUC.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
CONFIG = HERE / "config"
FORBIDDEN = ("narma", "mackey-glass", "mackey_glass", "waveform_auc", "ridge_lambda")
ARM_ORDER = (
    "FIELD_W100_FLUSH",
    "FIELD_W20_L20",
    "LIVE_W100_FLUSH",
    "LIVE_W20_L20",
)
T0_MEAN_R = 0.0318
T0_REL_TOL = 0.15
T2_W20_MEAN_R = 0.1253


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    for token in FORBIDDEN:
        if token in blob:
            raise SystemExit("PocketNeck-T3 checker cannot see a NARMA target.")
    for path in HERE.rglob("*"):
        name = path.name.lower()
        if path.is_file() and any(token in name for token in ("narma", "mackey")):
            raise SystemExit(f"PocketNeck-T3 checker refused forbidden file name: {path}")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def occupancy_flag(mean_r: float) -> str:
    if not math.isfinite(mean_r):
        return "NA"
    if mean_r > 0.95:
        return "SATURATED"
    if mean_r >= 0.05:
        return "ALIVE"
    return "DEAD"


def t0_flush_match(mean_r: float) -> bool:
    if not math.isfinite(mean_r) or T0_MEAN_R <= 0:
        return False
    return abs(mean_r - T0_MEAN_R) / T0_MEAN_R <= T0_REL_TOL


def print_theory() -> None:
    j_max = 636000.0
    v = 100.0 * 100.0 * 10.0
    k = 0.0033
    conv = 602.2
    c1 = j_max / (k * v * conv)
    c05 = 0.5 * c1
    h = (c05 ** 2) / (1.6 ** 2 + c05 ** 2)
    dt = 0.02
    d = 159.0
    dx = 5.0
    k_ftcs = d * dt / (dx * dx)
    print("PocketNeck-T3 THEORY (well-mixed closed algebra + FTCS gate)")
    print(f"  J_max={j_max:.3e} molecules/s  V={v:.0f} um^3")
    print(
        f"  C_ss(u=1)={c1:.6f} uM (target 3.2)  "
        f"C_ss(u=0.5)={c05:.6f} uM  H={h:.4f}  flag={occupancy_flag(h)}"
    )
    print(
        f"  dt={dt:.4f} s  D dt/dx^2={k_ftcs:.4f}  6D dt/dx^2={6.0 * k_ftcs:.4f} "
        "(ENGINEERING, not a clock retune)"
    )
    print("  Hypothesis (predeclared): FIELD_W100_FLUSH DEAD T0-class mean_R~0.0318;")
    print("  FIELD_W20_L20 ALIVE T2-class mean_R~0.125; living occupancy follows field.")
    print("  If FIELD_W100_FLUSH is ALIVE: stop, fix mask/bus. Do not retune K.")
    print("  If FIELD_W20_L20 is DEAD with a clean mass budget: write DEAD, do not retune.")
    print("  Occupancy voxels = pocket only. Cells read AHL; they do not write AHL.")
    print("  LIVING_DECISION=STOP_AFTER_T3_OCCUPANCY")


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def fmt(value: object) -> str:
    if value is None or value == "":
        return "NA"
    if isinstance(value, str):
        if value.lower() in {"nan", "na", "none"}:
            return "NA"
        try:
            number = float(value)
        except ValueError:
            return value
        value = number
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.6g}"
    return str(value)


def collect_rows() -> list[dict[str, str]]:
    screen = RESULTS / "occupancy_screen.csv"
    if screen.exists():
        return load_csv(screen)
    rows = []
    runs = RESULTS / "runs"
    if not runs.exists():
        return rows
    for arm in ARM_ORDER:
        path = runs / arm / "summary.csv"
        if path.exists():
            rows.extend(load_csv(path))
    return rows


def last_cell_sample(arm: str) -> dict[str, str]:
    path = RESULTS / "runs" / arm / "timeseries.csv"
    if not path.exists():
        return {}
    last = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                n = int(float(row.get("N", "0")))
            except (TypeError, ValueError):
                n = 0
            if n > 0:
                last = row
    return last


def occupancy_rows(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    out = {}
    for row in rows:
        if row.get("Label", "OCCUPANCY") != "OCCUPANCY":
            continue
        if str(row.get("Smoke", "false")).lower() in {"true", "1"}:
            continue
        out[row["Arm"]] = row
    return out


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "nan"))
    except (TypeError, ValueError):
        return float("nan")


def living_text(step: dict[str, dict[str, str]]) -> tuple[str, str]:
    decision = "STOP_AFTER_T3_OCCUPANCY"
    if "FIELD_W100_FLUSH" not in step:
        return decision, "missing FIELD_W100_FLUSH"
    flush = step["FIELD_W100_FLUSH"]
    occ = {arm: step[arm]["Occupancy"] for arm in ARM_ORDER if arm in step}
    mean_r = as_float(flush, "PocketMeanR")
    if occ.get("FIELD_W100_FLUSH") in {"ALIVE", "SATURATED"}:
        return (
            decision,
            "FIELD_W100_FLUSH is not DEAD. The wall-mask / bus advection is wrong. "
            "Stop. Fix the ticker wrap before interpreting W20 or living arms. "
            "Do not retune K.",
        )
    if not t0_flush_match(mean_r):
        return (
            decision,
            f"FIELD_W100_FLUSH mean_R={mean_r:.6g} is not T0/T2 class "
            f"(0.0318 ±15%). Fix the mask/bus. Do not celebrate a neck.",
        )
    neck = step.get("FIELD_W20_L20")
    if neck is None:
        return decision, "missing FIELD_W20_L20"
    neck_flag = neck["Occupancy"]
    if neck_flag == "SATURATED":
        return (
            decision,
            "FIELD_W20_L20 is SATURATED. Unexpected on an open neck. "
            "Fix conversion/source, not K.",
        )
    if neck_flag == "DEAD":
        return (
            decision,
            "FIELD_W20_L20 is DEAD on the BSim explicit stencil. "
            "Do not retune K, k, or J_max. Compare to T2 BDF (mean_R~0.125, dt, "
            "explicit vs implicit). Occupancy class was not recovered in Java.",
        )
    bits = [
        f"FIELD_W100_FLUSH C={as_float(flush, 'PocketMeanAHL_uM'):.4g} µM, "
        f"mean_R={mean_r:.4g} (T0 class, DEAD)",
        f"FIELD_W20_L20 C={as_float(neck, 'PocketMeanAHL_uM'):.4g} µM, "
        f"mean_R={as_float(neck, 'PocketMeanR'):.4g} (T2 class target ~{T2_W20_MEAN_R}, {neck_flag})",
    ]
    for live, field in (
        ("LIVE_W100_FLUSH", "FIELD_W100_FLUSH"),
        ("LIVE_W20_L20", "FIELD_W20_L20"),
    ):
        if live not in step:
            continue
        live_flag = occ[live]
        field_flag = occ[field]
        n_end = step[live].get("N_end", "NA")
        spill = step[live].get("Spillover", "NA")
        last = last_cell_sample(live)
        last_t = last.get("t_s", "NA")
        last_n = last.get("N", "NA")
        last_cr = last.get("cell_mean_R", "NA")
        bits.append(
            f"{live} field-gate {live_flag} (N_end={n_end}, spillover={spill}; "
            f"last N>0 at t={last_t}s N={last_n} cell_R={last_cr})"
        )
        if live_flag != field_flag and {live_flag, field_flag} <= {"ALIVE", "DEAD"}:
            return (
                decision,
                "Living arm occupancy disagrees with the matching field arm. "
                "Cells are leaking through walls or sitting in the bus — fix BC, "
                "do not retune K. " + "; ".join(bits),
            )
    return (
        decision,
        "BSim field port recovered the T2 leak story: flush door DEAD, W20 neck "
        "ALIVE. Living receivers (growth off, no AHL write) follow the field class. "
        + "; ".join(bits)
        + ". Stop after T3 occupancy. Do not start T4, NARMA, or LuxI.",
    )


def transport_status(step: dict[str, dict[str, str]]) -> str:
    field_ok = True
    for arm in ("FIELD_W100_FLUSH", "FIELD_W20_L20"):
        if arm not in step:
            return "NA"
        row = step[arm]
        vs = as_float(row, "MassResidualVsDominant")
        cpos = str(row.get("C_nonnegative", "")).lower() in {"true", "1"}
        if not (math.isfinite(vs) and vs <= 0.01 and cpos):
            field_ok = False
    live_note = True
    for arm in ("LIVE_W100_FLUSH", "LIVE_W20_L20"):
        if arm not in step:
            continue
        wrote = str(step[arm].get("CellsWriteAHL", "false")).lower() in {"true", "1"}
        if wrote:
            live_note = False
    if field_ok and live_note:
        return "VALIDATED_FOR_SCREENING"
    return "MASS_BUDGET_FAIL"


def write_screen_csv(rows: list[dict[str, str]]) -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / "occupancy_screen.csv"
    if not rows:
        raise SystemExit("no occupancy rows to write")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            if str(row.get("Smoke", "false")).lower() in {"true", "1"}:
                continue
            writer.writerow(row)
    return path


def write_hashes() -> list[tuple[str, str]]:
    hashes = []
    files = [
        CONFIG / "sim_config.properties",
        CONFIG / "field_w100.properties",
        CONFIG / "field_w20.properties",
        CONFIG / "live_w100.properties",
        CONFIG / "live_w20.properties",
        HERE / "bsim" / "BSimPocketNeckT3.java",
        HERE / "check_pocketneck_t3.py",
        HERE / "run_occupancy_screen.py",
        HERE / "compile_and_run.cmd",
    ]
    for path in files:
        hashes.append((path.relative_to(HERE).as_posix(), sha256_file(path)))
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketNeck-T3 SHA-256\n\n| File | SHA-256 |\n|---|---|\n"
        + "\n".join(f"| `{name}` | `{digest}` |" for name, digest in hashes)
        + "\n",
        encoding="utf-8",
    )
    return hashes


def write_report(rows: list[dict[str, str]]) -> Path:
    step = occupancy_rows(rows)
    status = transport_status(step)
    decision, text = living_text(step)
    hashes = write_hashes()
    lines = [
        "# PocketNeck-T3 occupancy screen",
        "",
        "BSim Java chip (wall-mask + bus advection) + living receivers,",
        "growth off. Not a task scout. No NARMA, no ridge, no Overall",
        "PASS/FAIL on a task.",
        "",
        "Geometry freeze: `examples/PocketDish/GEOMETRY_NECK.md`.",
        "Chip plan: `examples/PocketDish/CHIP_PLAN.md`.",
        "Clone source: `examples/HybridDish/bsim/BSimHybridDish.java` (not edited).",
        "This job: `examples/BSimReservoirPlanPocketNeckT3/PROTOCOL.md`.",
        "",
        f"`TRANSPORT_MODEL_STATUS={status}`",
        "",
    ]
    for arm in ARM_ORDER:
        flag = step[arm]["Occupancy"] if arm in step else "NA"
        lines.append(f"`OCCUPANCY_{arm}={flag}`")
    lines += [
        "",
        "`LIVING_DECISION=STOP_AFTER_T3_OCCUPANCY`",
        "",
        text,
        "",
        "## Result (occupancy — not a task score)",
        "",
        "T0/T2 reduced Python: flush 100 µm door DEAD (`mean_R=0.0318`);",
        "Kim-class `W20_L20` ALIVE (`mean_R=0.1253`). This screen asks whether",
        "that leak story survives the real BSim field and bacteria.",
        "",
    ]
    if "FIELD_W100_FLUSH" in step:
        flush = step["FIELD_W100_FLUSH"]
        lines.append(
            f"- `FIELD_W100_FLUSH`: mean AHL = {fmt(flush.get('PocketMeanAHL_uM'))} µM, "
            f"mean R = {fmt(flush.get('PocketMeanR'))}. **{flush['Occupancy']}.** "
            "Must match T0/T2 class (0.0318 ±15%). If this arm is ALIVE, the mask is wrong."
        )
    if "FIELD_W20_L20" in step:
        row = step["FIELD_W20_L20"]
        lines.append(
            f"- `FIELD_W20_L20`: mean AHL = {fmt(row.get('PocketMeanAHL_uM'))} µM, "
            f"mean R = {fmt(row.get('PocketMeanR'))}. **{row['Occupancy']}.** "
            f"Neck-mean AHL = {fmt(row.get('NeckMeanAHL_uM'))} µM (not in occupancy). "
            "T2 class target ~0.125; explicit stencil need not match BDF digit for digit."
        )
    for arm in ("LIVE_W100_FLUSH", "LIVE_W20_L20"):
        if arm not in step:
            continue
        row = step[arm]
        last = last_cell_sample(arm)
        lines.append(
            f"- `{arm}`: pocket-field mean R = {fmt(row.get('PocketMeanR'))}. "
            f"**{row['Occupancy']}.** N {fmt(row.get('N_start'))}→{fmt(row.get('N_end'))}, "
            f"spillover={fmt(row.get('Spillover'))}. "
            f"Last N>0 at t={fmt(last.get('t_s'))} s, N={fmt(last.get('N'))}, "
            f"cell mean R={fmt(last.get('cell_mean_R'))}, cell mean L={fmt(last.get('cell_mean_L'))} "
            "(diagnostic; gate is the field). Cells did not write AHL. "
            "Motile run-tumble + an open door + growth off empties the garage by "
            "declared spillover on the 9000 s clock; that is T4's N(t) question, "
            "not a K retune."
        )
    lines += [
        "",
        "Occupancy uses **pocket boxes only**. Growth is off. Do not retune K.",
        "Do not start NARMA or T4.",
        "",
        "## STEP_ON occupancy (u=0.5, T=9000 s, dt=0.02 s, pocket means)",
        "",
        "| Arm | cells | W (µm) | mean AHL (µM) | mean R | occupancy | plateau R | neck AHL | bus AHL | N start/end | spillover | residual vs dominant |",
        "|---|---|---:|---:|---:|---|---|---:|---:|---|---:|---:|",
    ]
    for arm in ARM_ORDER:
        if arm not in step:
            continue
        row = step[arm]
        cells = "50" if arm.startswith("LIVE") else "none"
        lines.append(
            f"| {arm} | {cells} | {fmt(row.get('W_um'))} | "
            f"{fmt(row.get('PocketMeanAHL_uM'))} | {fmt(row.get('PocketMeanR'))} | "
            f"{row.get('Occupancy')} | {row.get('PlateauR')} | "
            f"{fmt(row.get('NeckMeanAHL_uM'))} | {fmt(row.get('BusMeanAHL_uM'))} | "
            f"{fmt(row.get('N_start'))}/{fmt(row.get('N_end'))} | "
            f"{fmt(row.get('Spillover'))} | {fmt(row.get('MassResidualVsDominant'))} |"
        )
    lines += [
        "",
        "## Mass budget (STEP_ON)",
        "",
        "Production is the T0 spray can. Residual ≤1% of the dominant flux on field arms.",
        "Living arms: same budget plus cells did not write AHL.",
        "",
        "| Arm | production | remaining | decay | bus outlet | interface | wall sink | residual vs dominant | C ≥ 0 | walls ~0 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for arm in ARM_ORDER:
        if arm not in step:
            continue
        row = step[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('InjectedMass_molecules'))} | "
            f"{fmt(row.get('RemainingMass_molecules'))} | {fmt(row.get('DecayLoss_molecules'))} | "
            f"{fmt(row.get('BoundaryLoss_molecules'))} | {fmt(row.get('InterfaceLeak_molecules'))} | "
            f"{fmt(row.get('WallSink_molecules'))} | {fmt(row.get('MassResidualVsDominant'))} | "
            f"{row.get('C_nonnegative')} | {row.get('WallsZero')} |"
        )
    lines += [
        "",
        "## Preview",
        "",
        "`compile_and_run.cmd config\\live_w20.properties preview`",
        "",
        "Camera pulls back by bound Y (200 µm), not thickness 10. Close the window;",
        "headless occupancy is the score. `BSim.preview()` does not stop at 9000 s.",
        "",
        "## SHA-256 (configs and main Java)",
        "",
        "| File | SHA-256 |",
        "|---|---|",
    ]
    for name, digest in hashes:
        lines.append(f"| `{name}` | `{digest}` |")
    lines += [
        "",
        "## Stop",
        "",
        "T3 is living occupancy + preview only. Do not start T4 growth-on, NARMA,",
        "or LuxI unless a new frozen prompt says so. Do not retune HybridDish K",
        "or T0 J_max. Do not edit HybridDish Java or any GATE_EVIDENCE.md.",
        "",
        f"`TRANSPORT_MODEL_STATUS={status}`",
    ]
    for arm in ARM_ORDER:
        flag = step[arm]["Occupancy"] if arm in step else "NA"
        lines.append(f"`OCCUPANCY_{arm}={flag}`")
    lines += [
        "`LIVING_DECISION=STOP_AFTER_T3_OCCUPANCY`",
        "",
    ]
    path = RESULTS / "OCCUPANCY_SCREEN.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> int:
    refuse_narma()
    parser = argparse.ArgumentParser()
    parser.add_argument("--theory", action="store_true")
    args = parser.parse_args()
    if args.theory:
        print_theory()
        return 0
    rows = collect_rows()
    if not rows:
        raise SystemExit("missing occupancy rows; run run_occupancy_screen.py first")
    if any("narma" in key.lower() for row in rows for key in row):
        raise SystemExit("checker refused a NARMA column")
    write_screen_csv(rows)
    path = write_report(rows)
    print(f"wrote {path}")
    step = occupancy_rows(rows)
    status = transport_status(step)
    print(f"TRANSPORT_MODEL_STATUS={status}")
    for arm in ARM_ORDER:
        if arm in step:
            print(f"OCCUPANCY_{arm}={step[arm]['Occupancy']}")
    print("LIVING_DECISION=STOP_AFTER_T3_OCCUPANCY")
    _decision, text = living_text(step)
    print(text)
    if "FIELD_W100_FLUSH" in step:
        flush = step["FIELD_W100_FLUSH"]
        if flush["Occupancy"] != "DEAD" or not t0_flush_match(as_float(flush, "PocketMeanR")):
            return 3
    if "FIELD_W20_L20" in step and step["FIELD_W20_L20"]["Occupancy"] == "DEAD":
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
