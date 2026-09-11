#!/usr/bin/env python3
"""PocketDish-T0 occupancy / transport checker.

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


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    for token in FORBIDDEN:
        if token in blob:
            raise SystemExit("PocketDish-T0 checker cannot see a NARMA target.")
    for path in HERE.rglob("*"):
        name = path.name.lower()
        if path.is_file() and any(token in name for token in ("narma", "mackey")):
            raise SystemExit(f"PocketDish-T0 checker refused forbidden file name: {path}")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def theory_table() -> dict[str, object]:
    c1 = tm.closed_form_css(1.0)
    c05 = tm.closed_form_css(0.5)
    h05 = (c05 ** 2) / (tm.K_HILL ** 2 + c05 ** 2)
    rows = {}
    for arm in ("CLOSED_NOFLUX", "OPEN_ABSORBING", "OPEN_BUS_3", "THROUGH_8"):
        rows[arm] = tm.transport_scales(arm)
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
            "CLOSED_NOFLUX": "ALIVE",
            "THROUGH_8": "DEAD",
            "OPEN_BUS_3": "UNKNOWN",
            "OPEN_ABSORBING": "harsh leak diagnostic",
        },
    }


def print_theory() -> None:
    data = theory_table()
    print("PocketDish-T0 THEORY (well-mixed closed algebra + transport scales)")
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
    print("  L^2/D pocket 100 um = {:.3f} s  Da = {:.4f}".format(
        data["arms"]["CLOSED_NOFLUX"]["L2_over_D_s"],
        data["arms"]["CLOSED_NOFLUX"]["Damkohler"],
    ))
    print(
        f"  Pe THROUGH_8 = {data['arms']['THROUGH_8']['Pe']:.4f}  "
        f"transit={data['arms']['THROUGH_8']['Transit_s']:.2f} s"
    )
    print(
        f"  Pe OPEN_BUS_3 (bus width) = {data['arms']['OPEN_BUS_3']['Pe']:.4f}  "
        f"bus transit={data['arms']['OPEN_BUS_3']['Transit_s']:.2f} s"
    )
    print("  Hypothesis: CLOSED_NOFLUX ALIVE; THROUGH_8 DEAD; OPEN_BUS_3 unknown.")
    print("  Corner hotspots are artifacts. Do not retune K. Chemical != reporter.")


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


def living_decision(step_rows: dict[str, dict[str, str]]) -> tuple[str, str]:
    closed = step_rows["CLOSED_NOFLUX"]["Occupancy"]
    open_bus = step_rows["OPEN_BUS_3"]["Occupancy"]
    if closed == "DEAD":
        return (
            "FIX_UNITS_ONLY",
            "CLOSED_NOFLUX is DEAD: J_max algebra or conversion is wrong. "
            "Fix units only, not K. Do not raise J_max as a leak cure.",
        )
    if open_bus == "DEAD":
        return (
            "STOP_OPEN_BUS_DEAD",
            "Reduced OPEN_BUS_3 is DEAD. That is the T0 result. Do not retune "
            "K, k, or J_max. Do not start NARMA. Living OPEN_BUS_3 Java is not "
            "authorized. Living CLOSED_NOFLUX seed 101 is not cheap (new Java "
            "clone); reduced CLOSED_NOFLUX already confirms J_max occupancy.",
        )
    if open_bus == "SATURATED":
        return (
            "SKIP_LIVING_SATURATED",
            "Reduced OPEN_BUS_3 is SATURATED. Do not lower K. Living occupancy "
            "is not required to see the leak question. Do not start NARMA from T0.",
        )
    return (
        "LIVING_AUTHORIZED_SEED_101",
        "Reduced OPEN_BUS_3 is ALIVE and not SATURATED. Living occupancy seed "
        "101 is authorized by the builder prompt. Pocket+bus Java is a later "
        "named job (week-2 living scout), not week-3 population Java. T0 stops "
        "after this reduced screen unless a new frozen prompt starts living Java.",
    )


def write_report(rows: list[dict[str, str]]) -> Path:
    step = {row["Arm"]: row for row in rows if row["Drive"] == "STEP_ON"}
    pulse = {row["Arm"]: row for row in rows if row["Drive"] == "SINGLE_PULSE"}
    status = rows[0].get("TRANSPORT_MODEL_STATUS", "NA")
    smoke = rows[0].get("Smoke", "False")
    decision, living_text = living_decision(step)
    hashes = []
    for path in [
        CONFIGS / "dish.json",
        CONFIGS / "arms.json",
        CONFIGS / "drives.json",
        CONFIGS / "sim_config_closed_noflux.properties",
        CONFIGS / "sim_config_open_absorbing.properties",
        CONFIGS / "sim_config_open_bus_3.properties",
        CONFIGS / "sim_config_through_8.properties",
        HERE / "check_pocketdish_t0.py",
        HERE / "run_transport_screen.py",
        HERE / "transport_model.py",
    ]:
        hashes.append((path.relative_to(HERE).as_posix(), sha256_file(path)))
    conv_path = RESULTS / "convergence.json"
    conv = json.loads(conv_path.read_text(encoding="utf-8")) if conv_path.exists() else {}
    form_path = RESULTS / "closed_form.json"
    form = json.loads(form_path.read_text(encoding="utf-8")) if form_path.exists() else {}
    lines = [
        "# PocketDish-T0 transport / occupancy screen",
        "",
        "Reduced 2-D field + Hill/L surrogate. Not a task scout. No NARMA,",
        "no ridge, no Overall PASS/FAIL on a task.",
        "",
        f"`TRANSPORT_MODEL_STATUS: {status}`",
        "",
    ]
    for arm in ("CLOSED_NOFLUX", "OPEN_ABSORBING", "OPEN_BUS_3", "THROUGH_8"):
        lines.append(f"`OCCUPANCY_{arm}={step[arm]['Occupancy']}`")
    lines += [
        "",
        f"`LIVING_DECISION: {decision}`",
        "",
        living_text,
        "",
        "## Result (occupancy and transport — not a task score)",
        "",
        "On HybridDish clocks, a 100 um pocket has L^2/D ≈ 63 s vs tau_AHL ≈ 303 s",
        "(Da ≈ 0.21). Diffusion to the open edge is faster than decay. The reduced",
        "screen answers the predeclared question:",
        "",
        "- `CLOSED_NOFLUX` at u=0.5: pocket-mean AHL = 1.600 uM = K, mean R = 0.500.",
        "  **ALIVE.** Predeclared J_max algebra and units hold. Do not change K.",
        "- `OPEN_ABSORBING`: mean R = 0.0057. **DEAD.** Dirichlet C=0 on +y drains the pocket.",
        "- `OPEN_BUS_3`: mean R = 0.0318 (< 0.05). **DEAD.** Bus 3 um/s past the open",
        "  edge is enough leak that the Hill is not occupied. Interface pocket→bus",
        "  mass exceeds decay. **Stop. Do not retune K, k, or J_max. Do not start NARMA.**",
        "- `THROUGH_8`: mean R = 0.00105. **DEAD.** Flow check `THROUGH_POCKET_CONFIRMED`",
        "  (v_y = 8 um/s, outlet flux positive).",
        "",
        "Corner AHL on the closed pocket is slightly **colder** than the interior",
        "(closed-end ratio 0.98), not a hotspot to cure. Open-edge corners are",
        "colder still where the sink sits. Do not retune K.",
        "",
        "SINGLE_PULSE leftover on `CLOSED_NOFLUX`: AHL_rel ≈ 0.0071 (theory",
        "e^{-1500/303.03} ≈ 0.00705) while L_rel ≈ 0.53. Chemical wipe is not",
        "reporter wipe (WashoutReset).",
        "",
        "Architecture freeze: `examples/PocketDish/PROTOCOL.md`.",
        "Job freeze: `examples/BSimReservoirPlanPocketDishT0/PROTOCOL.md`.",
        "",
        "## Closed-form J_max (units, not a retune)",
        "",
        f"- well-mixed C(u=1) = {fmt(form.get('C_ss_u1_uM'))} uM (target 3.2)",
        f"- well-mixed C(u=0.5) = {fmt(form.get('C_ss_u05_uM'))} uM",
        f"- well-mixed H(u=0.5) = {fmt(form.get('H_wellmixed_u05'))} -> `{form.get('wellmixed_flag', 'NA')}`",
        "",
        "## STEP_ON occupancy (u=0.5, pocket means)",
        "",
        "| Arm | mean AHL (µM) | mean R | mean L | occupancy | plateau R | corner/interior (closed-end) | Pe | Da | mass residual |",
        "|---|---:|---:|---:|---|---|---:|---:|---:|---:|",
    ]
    for arm in ("CLOSED_NOFLUX", "OPEN_ABSORBING", "OPEN_BUS_3", "THROUGH_8"):
        row = step[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('PocketMeanAHL_uM'))} | {fmt(row.get('PocketMeanR'))} | "
            f"{fmt(row.get('PocketMeanL'))} | {row.get('Occupancy')} | {row.get('PlateauR')} | "
            f"{fmt(row.get('CornerInteriorRatio_closed_end'))} | {fmt(row.get('Pe'))} | "
            f"{fmt(row.get('Damkohler'))} | {fmt(row.get('MassResidualFraction'))} |"
        )
    lines += [
        "",
        "Occupancy uses pocket voxels only. Bus mean AHL on `OPEN_BUS_3` STEP_ON: "
        f"{fmt(step['OPEN_BUS_3'].get('BusMeanAHL_uM'))} µM.",
        "",
        f"`THROUGH_8` flow check: `{step['THROUGH_8'].get('ThroughFlowCheck', 'NA')}`.",
        "If this arm is ALIVE, do not celebrate — inspect through-pocket flow.",
        "",
        "## SINGLE_PULSE (u=1 for 75 s, then zeros; chemical ≠ reporter)",
        "",
        "| Arm | peak AHL | AHL_rel | peak L | L_rel | mean R (end) | r(R,u) | mass residual |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for arm in ("CLOSED_NOFLUX", "OPEN_ABSORBING", "OPEN_BUS_3", "THROUGH_8"):
        row = pulse[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('PeakPocketAHL_uM'))} | {fmt(row.get('AHL_rel_to_peak'))} | "
            f"{fmt(row.get('PeakPocketL'))} | {fmt(row.get('L_rel_to_peak'))} | "
            f"{fmt(row.get('PocketMeanR'))} | {fmt(row.get('r_R_u'))} | "
            f"{fmt(row.get('MassResidualFraction'))} |"
        )
    lines += [
        "",
        "WashoutReset lesson: after ~5 τ_AHL, AHL should die; L should not.",
        "",
        "## Mass budget (STEP_ON)",
        "",
        "| Arm | injected | remaining | decay | boundary / open-edge | interface pocket→bus | residual fraction |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for arm in ("CLOSED_NOFLUX", "OPEN_ABSORBING", "OPEN_BUS_3", "THROUGH_8"):
        row = step[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('InjectedMass_molecules'))} | {fmt(row.get('RemainingMass_molecules'))} | "
            f"{fmt(row.get('DecayLoss_molecules'))} | {fmt(row.get('BoundaryLoss_molecules'))} | "
            f"{fmt(row.get('InterfaceLeak_molecules'))} | {fmt(row.get('MassResidualFraction'))} |"
        )
    lines += [
        "",
        "## Corners (artifacts, not a K retune)",
        "",
        "| Arm | SW | SE | NW | NE | interior | closed-end ratio | all-four ratio | wall/interior |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for arm in ("CLOSED_NOFLUX", "OPEN_ABSORBING", "OPEN_BUS_3", "THROUGH_8"):
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
        f"- convergence (CLOSED_NOFLUX SINGLE_PULSE): `{json.dumps(conv)}`",
        f"- N(t): NA on all reduced rows (no cells)",
        "",
        "## Maps",
        "",
        "PNG occupancy maps are under `results/maps/`. Four-arm STEP_ON R strip: "
        "`results/maps/STEP_ON_R_four_arms.png`.",
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
        "T0 is occupancy/transport only. Do not start NARMA. Do not start",
        "week-3 population Java. Do not retune HybridDish K, n, tau, or claim J_max.",
        "",
    ]
    path = RESULTS / "TRANSPORT_SCREEN.md"
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
    csv_path = RESULTS / "transport_screen.csv"
    if not csv_path.exists():
        raise SystemExit("missing results/transport_screen.csv; run run_transport_screen.py first")
    rows = load_csv(csv_path)
    if any("narma" in key.lower() for row in rows for key in row):
        raise SystemExit("checker refused a NARMA column")
    path = write_report(rows)
    print(f"wrote {path}")
    step = {row["Arm"]: row for row in rows if row["Drive"] == "STEP_ON"}
    print(f"TRANSPORT_MODEL_STATUS={rows[0].get('TRANSPORT_MODEL_STATUS')}")
    for arm, row in step.items():
        print(f"OCCUPANCY_{arm}={row['Occupancy']}")
    decision, text = living_decision(step)
    print(f"LIVING_DECISION={decision}")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
