#!/usr/bin/env python3
"""PocketNeck-T2 occupancy checker.

Cannot see a NARMA target. No ridge, no Mackey-Glass, no waveform AUC.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import transport_model as tm  # noqa: E402

RESULTS = HERE / "results"
CONFIGS = HERE / "configs"
FORBIDDEN = ("narma", "mackey-glass", "mackey_glass", "waveform_auc", "ridge_lambda")
ARM_ORDER = ("W100_FLUSH", "W50_L20", "W20_L20", "W10_L20")


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    for token in FORBIDDEN:
        if token in blob:
            raise SystemExit("PocketNeck-T2 checker cannot see a NARMA target.")
    for path in HERE.rglob("*"):
        name = path.name.lower()
        if path.is_file() and any(token in name for token in ("narma", "mackey")):
            raise SystemExit(f"PocketNeck-T2 checker refused forbidden file name: {path}")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def theory_table() -> dict[str, object]:
    c1 = tm.closed_form_css(1.0)
    c05 = tm.closed_form_css(0.5)
    h05 = (c05 ** 2) / (tm.K_HILL ** 2 + c05 ** 2)
    rows = {arm: tm.transport_scales(arm) for arm in ARM_ORDER}
    return {
        "tau_AHL_s": 1.0 / tm.K_DECAY,
        "tau_R_s": tm.TAU_R,
        "tau_L_s": tm.TAU_L,
        "J_max": tm.J_MAX,
        "V_um3": tm.POCKET_V,
        "C_ss_u1_uM": c1,
        "C_ss_u05_uM": c05,
        "H_wellmixed_u05": h05,
        "wellmixed_occupancy": tm.occupancy_flag(h05),
        "arms": rows,
        "hypothesis": {
            "W100_FLUSH": "DEAD, T0 OPEN_BUS_3 class mean_R ≈ 0.0318",
            "W50_L20": "unknown; leak∝W hypothesis, near occupancy gate",
            "W20_L20": "unknown; literature constriction class",
            "W10_L20": "unknown; aggressive, not 3 µm hex",
        },
    }


def print_theory() -> None:
    data = theory_table()
    print("PocketNeck-T2 THEORY (well-mixed closed algebra + transport scales)")
    print(f"  J_max={data['J_max']:.3e} molecules/s  V={data['V_um3']:.0f} um^3")
    print(
        f"  C_ss(u=1)={data['C_ss_u1_uM']:.6f} uM (target 3.2)  "
        f"C_ss(u=0.5)={data['C_ss_u05_uM']:.6f} uM  H={data['H_wellmixed_u05']:.4f}  "
        f"flag={data['wellmixed_occupancy']}"
    )
    print(
        f"  tau_AHL={data['tau_AHL_s']:.4f} s  tau_R={data['tau_R_s']:.1f} s  "
        f"tau_L={data['tau_L_s']:.1f} s"
    )
    print(
        "  L^2/D pocket 100 um = {:.3f} s  Da = {:.4f}  bus Pe = {:.4f}".format(
            data["arms"]["W100_FLUSH"]["L2_over_D_s"],
            data["arms"]["W100_FLUSH"]["Damkohler"],
            data["arms"]["W100_FLUSH"]["Pe"],
        )
    )
    print("  Hypothesis (predeclared): W100_FLUSH DEAD T0-class; narrower W may ALIVE.")
    print("  If all neck arms DEAD: stop. Do not retune K, k, J_max, or L_n.")
    print("  Occupancy voxels = pocket only. Neck is not Danino full-edge feed.")


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def fmt(value: object) -> str:
    if value is None or value == "":
        return "NA"
    if isinstance(value, str):
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


def occupancy_rows(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    out = {}
    for row in rows:
        if row.get("Label", "OCCUPANCY") != "OCCUPANCY":
            continue
        if row.get("Drive", "STEP_ON") != "STEP_ON":
            continue
        out[row["Arm"]] = row
    return out


def living_text(step: dict[str, dict[str, str]]) -> tuple[str, str]:
    flush = step["W100_FLUSH"]
    occ = {arm: step[arm]["Occupancy"] for arm in ARM_ORDER if arm in step}
    mean_r = float(flush["PocketMeanR"])
    if occ.get("W100_FLUSH") == "ALIVE" or occ.get("W100_FLUSH") == "SATURATED":
        return (
            "STOP_AFTER_NECK_SCREEN",
            "W100_FLUSH is not DEAD. The flush clone is wrong. Do not celebrate a neck. "
            "Fix geometry/stencil before interpreting W50/W20/W10.",
        )
    if not tm.t0_flush_match(mean_r):
        return (
            "STOP_AFTER_NECK_SCREEN",
            f"W100_FLUSH mean_R={mean_r:.6g} is not T0 OPEN_BUS_3 class (0.0318 ±10%). "
            "Fix the flush clone before interpreting necks.",
        )
    neck_flags = [occ[arm] for arm in ("W50_L20", "W20_L20", "W10_L20") if arm in occ]
    if any(flag == "SATURATED" for flag in neck_flags):
        return (
            "STOP_AFTER_NECK_SCREEN",
            "A neck arm is SATURATED. Unexpected on an open neck. Fix conversion/source, not K.",
        )
    if any(flag == "ALIVE" for flag in neck_flags):
        alive = [arm for arm in ("W50_L20", "W20_L20", "W10_L20") if occ.get(arm) == "ALIVE"]
        bits = []
        for arm in ("W50_L20", "W20_L20", "W10_L20"):
            if arm not in step:
                continue
            bits.append(
                f"{arm} mean_R={float(step[arm]['PocketMeanR']):.4g} "
                f"C={float(step[arm]['PocketMeanAHL_uM']):.4g} µM"
            )
        flush_c = float(flush["PocketMeanAHL_uM"])
        flush_iface = float(flush["InterfaceLeak_molecules"])
        flush_decay = float(flush["DecayLoss_molecules"])
        return (
            "STOP_AFTER_NECK_SCREEN",
            "Leak story: a Kim-class neck occupies PocketHill on HybridDish clocks while "
            "the flush 100 µm door does not. "
            f"W100_FLUSH C={flush_c:.4g} µM, mean_R={mean_r:.4g} (T0 OPEN_BUS_3 class), "
            f"interface {flush_iface:.3g} vs decay {flush_decay:.3g}. "
            + "; ".join(bits)
            + ". Narrower W cuts leak (interface falls, decay takes a larger share). "
            "That is the door story, not a K retune. "
            f"ALIVE arm(s): {', '.join(alive)}. "
            "Do not start NARMA, LuxI, or a living pack. Do not add more widths.",
        )
    if neck_flags and all(flag == "DEAD" for flag in neck_flags):
        return (
            "STOP_AFTER_NECK_SCREEN",
            "All three neck arms are DEAD. Neck is not enough on HybridDish clocks. "
            "Do not retune K, k, J_max, or L_n. Do not add W=5 µm after looking. "
            "A later named prompt may ask hours + seed (PocketFill), not a quieter K.",
        )
    return (
        "STOP_AFTER_NECK_SCREEN",
        "Occupancy table incomplete. Run the screen first.",
    )


def write_report(rows: list[dict[str, str]]) -> Path:
    step = occupancy_rows(rows)
    conv_rows = [row for row in rows if row.get("Label") == "CONVERGENCE"]
    status = rows[0].get("TRANSPORT_MODEL_STATUS", "NA") if rows else "NA"
    smoke = rows[0].get("Smoke", "False") if rows else "NA"
    decision, text = living_text(step) if "W100_FLUSH" in step else ("STOP_AFTER_NECK_SCREEN", "missing W100_FLUSH")
    hashes = []
    for path in [
        CONFIGS / "dish.json",
        CONFIGS / "arms.json",
        CONFIGS / "drives.json",
        HERE / "check_pocketneck_t2.py",
        HERE / "run_occupancy_screen.py",
        HERE / "transport_model.py",
        HERE / "test_mass_budget.py",
    ]:
        hashes.append((path.relative_to(HERE).as_posix(), sha256_file(path)))
    conv_path = RESULTS / "convergence.json"
    conv = json.loads(conv_path.read_text(encoding="utf-8")) if conv_path.exists() else {}
    form_path = RESULTS / "closed_form.json"
    form = json.loads(form_path.read_text(encoding="utf-8")) if form_path.exists() else {}

    lines = [
        "# PocketNeck-T2 occupancy screen",
        "",
        "Reduced 2-D field + Hill/L surrogate. Not a task scout. No NARMA,",
        "no ridge, no Overall PASS/FAIL on a task.",
        "",
        "Geometry freeze: `examples/PocketDish/GEOMETRY_NECK.md`.",
        "T0 job protocol: `examples/BSimReservoirPlanPocketDishT0/PROTOCOL.md`.",
        "This job: `examples/BSimReservoirPlanPocketNeckT2/PROTOCOL.md`.",
        "",
        f"`TRANSPORT_MODEL_STATUS={status}`",
        "",
    ]
    for arm in ARM_ORDER:
        flag = step[arm]["Occupancy"] if arm in step else "NA"
        lines.append(f"`OCCUPANCY_{arm}={flag}`")
    lines += [
        "",
        "`LIVING_DECISION=STOP_AFTER_NECK_SCREEN`",
        "",
        text,
        "",
        "## Result (occupancy — not a task score)",
        "",
        "On HybridDish clocks, T0 opened the entire 100 µm wall onto a 3 µm/s bus",
        "and PocketHill was DEAD (`OPEN_BUS_3` mean_R = 0.0318). This screen asks",
        "whether a Kim 2016 crevice-class neck occupies the same Hill with the",
        "same T0 spray can. The neck is inserted between pocket and bus; pocket",
        "volume is not eaten.",
        "",
    ]
    if "W100_FLUSH" in step:
        flush = step["W100_FLUSH"]
        lines.append(
            f"- `W100_FLUSH`: mean AHL = {fmt(flush.get('PocketMeanAHL_uM'))} µM, "
            f"mean R = {fmt(flush.get('PocketMeanR'))}. **{flush['Occupancy']}.** "
            "Must match T0 `OPEN_BUS_3` class (0.0318 ±10%). If this arm is ALIVE, "
            "the clone is wrong."
        )
    for arm in ("W50_L20", "W20_L20", "W10_L20"):
        if arm not in step:
            continue
        row = step[arm]
        lines.append(
            f"- `{arm}`: mean AHL = {fmt(row.get('PocketMeanAHL_uM'))} µM, "
            f"mean R = {fmt(row.get('PocketMeanR'))}. **{row['Occupancy']}.** "
            f"Neck-mean AHL = {fmt(row.get('NeckMeanAHL_uM'))} µM "
            f"(not in occupancy)."
        )
    lines += [
        "",
        "Occupancy uses **pocket voxels only**. Neck and bus AHL are diagnostics.",
        "Do not retune K, k, J_max, or L_n. Do not start NARMA.",
        "",
        "## Closed-form J_max (units, not a retune)",
        "",
        f"- well-mixed C(u=1) = {fmt(form.get('C_ss_u1_uM'))} uM (target 3.2)",
        f"- well-mixed C(u=0.5) = {fmt(form.get('C_ss_u05_uM'))} uM",
        f"- well-mixed H(u=0.5) = {fmt(form.get('H_wellmixed_u05'))} -> `{form.get('wellmixed_flag', 'NA')}`",
        "",
        "## STEP_ON occupancy (u=0.5, T=9000 s, pocket means)",
        "",
        "| Arm | W (µm) | L_n (µm) | mean AHL (µM) | mean R | mean L | occupancy | plateau R | neck AHL | bus AHL | corner/interior (closed-end) | Pe | Da | mass residual vs dominant |",
        "|---|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for arm in ARM_ORDER:
        if arm not in step:
            continue
        row = step[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('W_um'))} | {fmt(row.get('L_n_um'))} | "
            f"{fmt(row.get('PocketMeanAHL_uM'))} | {fmt(row.get('PocketMeanR'))} | "
            f"{fmt(row.get('PocketMeanL'))} | {row.get('Occupancy')} | {row.get('PlateauR')} | "
            f"{fmt(row.get('NeckMeanAHL_uM'))} | {fmt(row.get('BusMeanAHL_uM'))} | "
            f"{fmt(row.get('CornerInteriorRatio_closed_end'))} | {fmt(row.get('Pe'))} | "
            f"{fmt(row.get('Damkohler'))} | {fmt(row.get('MassResidualVsDominant'))} |"
        )
    lines += [
        "",
        "## Mass budget (STEP_ON)",
        "",
        "T0 lesson: `W100_FLUSH` interface should dominate decay (~2.3×10⁹ vs ~1.1×10⁹ class).",
        "Interface is pocket→bus on the flush arm and pocket→neck on necked arms.",
        "Boundary loss is the bus outlet (inlet C=0).",
        "",
        "| Arm | production | remaining | decay | bus outlet | interface | residual fraction | residual vs dominant | C ≥ 0 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for arm in ARM_ORDER:
        if arm not in step:
            continue
        row = step[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('InjectedMass_molecules'))} | "
            f"{fmt(row.get('RemainingMass_molecules'))} | {fmt(row.get('DecayLoss_molecules'))} | "
            f"{fmt(row.get('BoundaryLoss_molecules'))} | {fmt(row.get('InterfaceLeak_molecules'))} | "
            f"{fmt(row.get('MassResidualFraction'))} | {fmt(row.get('MassResidualVsDominant'))} | "
            f"{row.get('C_nonnegative')} |"
        )
    lines += [
        "",
        "## Corners (artifacts, not a K retune)",
        "",
        "T0 closed-end ratio on the flush door was ~0.98 (slightly colder, not a hotspot).",
        "",
        "| Arm | SW | SE | NW | NE | interior | closed-end ratio | all-four ratio | wall/interior |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for arm in ARM_ORDER:
        if arm not in step:
            continue
        row = step[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('C_SW'))} | {fmt(row.get('C_SE'))} | {fmt(row.get('C_NW'))} | "
            f"{fmt(row.get('C_NE'))} | {fmt(row.get('C_interior'))} | "
            f"{fmt(row.get('CornerInteriorRatio_closed_end'))} | "
            f"{fmt(row.get('CornerInteriorRatio_all_four'))} | {fmt(row.get('WallInteriorRatio'))} |"
        )
    lines += [
        "",
        "## Numerical status",
        "",
        f"- smoke: {smoke}",
        f"- optional CONVERGENCE (W10_L20, dx=2.5 µm): `{json.dumps(conv)}`",
        f"- N(t): NA on all reduced rows (no cells)",
        "",
    ]
    if conv_rows:
        crow = conv_rows[0]
        lines += [
            "CONVERGENCE row is **not** occupancy. Do not replace a coarse DEAD with a finer ALIVE.",
            "",
            f"- `W10_L20` dx=2.5: mean R = {fmt(crow.get('PocketMeanR'))}, "
            f"flag {crow.get('Occupancy')} (label `{crow.get('Label')}`)",
            "",
        ]
    elif conv.get("skipped"):
        lines += [
            "Optional dx=2.5 `W10_L20` CONVERGENCE was **skipped**.",
            "",
        ]
    lines += [
        "## Maps",
        "",
        "PNG maps of pocket (+ neck + bus) C and R are under `results/maps/`.",
        "Four-arm STEP_ON R strip: `results/maps/STEP_ON_R_four_arms.png`.",
        "",
        "## SHA-256 (configs and checker)",
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
        "T2 is occupancy / leak vs door width only. Do not start NARMA. Do not",
        "start week-2 living Java, circle twin, or QS seed unless a new frozen",
        "prompt says so. Do not retune HybridDish K, n, tau, or T0 J_max.",
        "",
        f"`TRANSPORT_MODEL_STATUS={status}`",
    ]
    for arm in ARM_ORDER:
        flag = step[arm]["Occupancy"] if arm in step else "NA"
        lines.append(f"`OCCUPANCY_{arm}={flag}`")
    lines += [
        "`LIVING_DECISION=STOP_AFTER_NECK_SCREEN`",
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
    csv_path = RESULTS / "occupancy_screen.csv"
    if not csv_path.exists():
        raise SystemExit("missing results/occupancy_screen.csv; run run_occupancy_screen.py first")
    rows = load_csv(csv_path)
    if any("narma" in key.lower() for row in rows for key in row):
        raise SystemExit("checker refused a NARMA column")
    path = write_report(rows)
    print(f"wrote {path}")
    step = occupancy_rows(rows)
    print(f"TRANSPORT_MODEL_STATUS={rows[0].get('TRANSPORT_MODEL_STATUS')}")
    for arm in ARM_ORDER:
        if arm in step:
            print(f"OCCUPANCY_{arm}={step[arm]['Occupancy']}")
    print("LIVING_DECISION=STOP_AFTER_NECK_SCREEN")
    _decision, text = living_text(step)
    print(text)
    if "W100_FLUSH" in step:
        flush = step["W100_FLUSH"]
        if flush["Occupancy"] != "DEAD" or not tm.t0_flush_match(float(flush["PocketMeanR"])):
            return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
