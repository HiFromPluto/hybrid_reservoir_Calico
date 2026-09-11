#!/usr/bin/env python3
"""PocketNeck-T4 population checker.

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
ARM_ORDER = ("LIVE_W20_GROWTH", "LIVE_W100_GROWTH")
POP_FLAGS = {"EMPTY", "CLAMPED", "EQUILIBRIUM", "TRANSIENT"}


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    for token in FORBIDDEN:
        if token in blob:
            raise SystemExit("PocketNeck-T4 checker cannot see a NARMA target.")
    for path in HERE.rglob("*"):
        name = path.name.lower()
        if path.is_file() and any(token in name for token in ("narma", "mackey")):
            raise SystemExit(f"PocketNeck-T4 checker refused forbidden file name: {path}")


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
    growth = 4.0 * math.pi / 1800.0
    print("PocketNeck-T4 THEORY (growth-on N(t) / spillover; T3 chip)")
    print(f"  J_max={j_max:.3e} molecules/s  V={v:.0f} um^3")
    print(
        f"  C_ss(u=1)={c1:.6f} uM (target 3.2)  "
        f"C_ss(u=0.5)={c05:.6f} uM  H={h:.4f}  flag={occupancy_flag(h)}"
    )
    print(
        f"  dt={dt:.4f} s  D dt/dx^2={k_ftcs:.4f}  6D dt/dx^2={6.0 * k_ftcs:.4f} "
        "(ENGINEERING, not a clock retune)"
    )
    print(
        f"  GROWTH_RATE={growth:.6g}  T_gen=1800 s  K_clamp=4000  "
        "p_removal HybridDish, acid off"
    )
    print("  Hypothesis (predeclared): LIVE_W20_GROWTH occupancy ALIVE.")
    print("  Population unknown: EMPTY | CLAMPED | EQUILIBRIUM | TRANSIENT.")
    print("  Do not retune K, J_max, or growth after seeing N.")
    print("  Chemostat sentence only if EQUILIBRIUM.")
    print("  LIVING_DECISION=STOP_AFTER_T4_POPULATION")


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
    path = RESULTS / "runs" / arm / "n_trace.csv"
    if not path.exists():
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
        label = row.get("Label", "POPULATION")
        if label not in {"POPULATION", "OCCUPANCY"}:
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


def as_int(row: dict[str, str], key: str) -> int:
    try:
        return int(float(row.get(key, "nan")))
    except (TypeError, ValueError):
        return -1


def chemostat_sentence(pop: str) -> str:
    if pop == "EQUILIBRIUM":
        return (
            "N(t) equilibrated off the clamp on the 18000 s clock. "
            "Chemostat sentence is allowed for this arm."
        )
    if pop == "CLAMPED":
        return (
            "N reached the safety cap. The cap is still the story. "
            "Do not call this a physical chemostat."
        )
    if pop == "EMPTY":
        return (
            "Growth lost to the neck: same physics as T3 plus birth, "
            "not a K retune. Do not call this a chemostat."
        )
    return (
        "N(t) is still rising or falling at 18000 s (TRANSIENT). "
        "Write that. Do not extend the clock. Do not call this a chemostat."
    )


def living_text(step: dict[str, dict[str, str]]) -> tuple[str, str]:
    decision = "STOP_AFTER_T4_POPULATION"
    if "LIVE_W20_GROWTH" not in step:
        return decision, "missing LIVE_W20_GROWTH"
    w20 = step["LIVE_W20_GROWTH"]
    occ = w20["Occupancy"]
    pop = w20.get("Population", "NA")
    mean_r = as_float(w20, "PocketMeanR")
    bits = [
        f"LIVE_W20_GROWTH occupancy {occ} mean_R={mean_r:.4g}, "
        f"population {pop}, N {fmt(w20.get('N_start'))}->{fmt(w20.get('N_end'))} "
        f"Nmax={fmt(w20.get('N_max'))} spill={fmt(w20.get('Spillover'))} "
        f"births={fmt(w20.get('Births'))} clamp={fmt(w20.get('ClampRemovals'))}"
    ]
    if occ == "DEAD":
        return (
            decision,
            "LIVE_W20_GROWTH occupancy is DEAD. Stop. Do not raise J_max. "
            + "; ".join(bits),
        )
    if occ == "SATURATED":
        return (
            decision,
            "LIVE_W20_GROWTH occupancy is SATURATED. Unexpected on an open neck. "
            "Fix conversion/source, not K. " + "; ".join(bits),
        )
    if "LIVE_W100_GROWTH" in step:
        w100 = step["LIVE_W100_GROWTH"]
        bits.append(
            f"LIVE_W100_GROWTH occupancy {w100['Occupancy']} "
            f"population {w100.get('Population', 'NA')} "
            f"N {fmt(w100.get('N_start'))}->{fmt(w100.get('N_end'))} "
            f"Nmax={fmt(w100.get('N_max'))}"
        )
    return decision, " ".join(bits) + ". " + chemostat_sentence(pop)


def transport_status(step: dict[str, dict[str, str]]) -> str:
    ok = True
    for arm, row in step.items():
        vs = as_float(row, "MassResidualVsDominant")
        cpos = str(row.get("C_nonnegative", "")).lower() in {"true", "1"}
        wrote = str(row.get("CellsWriteAHL", "false")).lower() in {"true", "1"}
        if not (math.isfinite(vs) and vs <= 0.01 and cpos and not wrote):
            ok = False
    if not step:
        return "NA"
    return "VALIDATED_FOR_SCREENING" if ok else "MASS_BUDGET_FAIL"


def write_screen_csv(rows: list[dict[str, str]]) -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / "population_screen.csv"
    if not rows:
        raise SystemExit("no population rows to write")
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
        CONFIG / "live_w20_growth.properties",
        CONFIG / "live_w100_growth.properties",
        HERE / "bsim" / "BSimPocketNeckT4.java",
        HERE / "check_pocketneck_t4.py",
        HERE / "run_population_screen.py",
        HERE / "compile_and_run.cmd",
    ]
    for path in files:
        hashes.append((path.relative_to(HERE).as_posix(), sha256_file(path)))
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketNeck-T4 SHA-256\n\n| File | SHA-256 |\n|---|---|\n"
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
        "# PocketNeck-T4 population screen",
        "",
        "BSim Java chip (T3 wall-mask + bus advection) + living receivers,",
        "growth on, HybridDish p_removal clamp at K=4000, acid off.",
        "Not a task scout. No NARMA, no ridge, no Overall PASS/FAIL on a task.",
        "",
        "Geometry freeze: `examples/PocketDish/GEOMETRY_NECK.md`.",
        "Chip plan: `examples/PocketDish/CHIP_PLAN.md`.",
        "T3 standing (cite only): `examples/PocketDish/T3_STANDING.md`.",
        "This job: `examples/BSimReservoirPlanPocketNeckT4/PROTOCOL.md`.",
        "",
        f"`TRANSPORT_MODEL_STATUS={status}`",
        "",
    ]
    for arm in ARM_ORDER:
        occ = step[arm]["Occupancy"] if arm in step else "NA"
        pop = step[arm].get("Population", "NA") if arm in step else "NA"
        lines.append(f"`OCCUPANCY_{arm}={occ}`")
        lines.append(f"`POPULATION_{arm}={pop}`")
    lines += [
        "",
        "`LIVING_DECISION=STOP_AFTER_T4_POPULATION`",
        "",
        text,
        "",
        "## Result (population — not a task score)",
        "",
        "T3 growth-off emptied the garage (N=50->0 by spillover). T4 asks",
        "whether growth-on on the same W20_L20 chip keeps a population,",
        "hits the safety clamp, or still goes empty.",
        "",
    ]
    for arm in ARM_ORDER:
        if arm not in step:
            continue
        row = step[arm]
        pop = row.get("Population", "NA")
        last = last_cell_sample(arm)
        lines.append(
            f"- `{arm}`: pocket-field mean R = {fmt(row.get('PocketMeanR'))}. "
            f"**{row['Occupancy']}.** Population **{pop}.** "
            f"N(0)={fmt(row.get('N_start'))}, N_max={fmt(row.get('N_max'))}, "
            f"N_end={fmt(row.get('N_end'))}, spillover={fmt(row.get('Spillover'))}, "
            f"births={fmt(row.get('Births'))}, clamp-removals={fmt(row.get('ClampRemovals'))}. "
            f"Last N>0 at t={fmt(last.get('t_s'))} s, N={fmt(last.get('N'))}. "
            f"{chemostat_sentence(pop)} Cells did not write AHL."
        )
    lines += [
        "",
        "Occupancy uses **pocket boxes only**. Do not retune K or growth.",
        "Do not start NARMA, LuxI, or T5.",
        "",
        "## STEP_ON (u=0.5, T=18000 s, dt=0.02 s, pocket means)",
        "",
        "| Arm | W (um) | mean AHL (uM) | mean R | occupancy | population | N start/max/end | spillover | births | clamp | residual vs dominant |",
        "|---|---:|---:|---:|---|---|---|---:|---:|---:|---:|",
    ]
    for arm in ARM_ORDER:
        if arm not in step:
            continue
        row = step[arm]
        lines.append(
            f"| {arm} | {fmt(row.get('W_um'))} | "
            f"{fmt(row.get('PocketMeanAHL_uM'))} | {fmt(row.get('PocketMeanR'))} | "
            f"{row.get('Occupancy')} | {row.get('Population')} | "
            f"{fmt(row.get('N_start'))}/{fmt(row.get('N_max'))}/{fmt(row.get('N_end'))} | "
            f"{fmt(row.get('Spillover'))} | {fmt(row.get('Births'))} | "
            f"{fmt(row.get('ClampRemovals'))} | {fmt(row.get('MassResidualVsDominant'))} |"
        )
    lines += [
        "",
        "## Mass budget (STEP_ON)",
        "",
        "Production is the T0 spray can. Residual <=1% of the dominant flux.",
        "Cells did not write AHL.",
        "",
        "| Arm | production | remaining | decay | bus outlet | interface | wall sink | residual vs dominant | C >= 0 | walls ~0 |",
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
        "## N(t) traces",
        "",
        "Per-arm `results/runs/<ARM>/n_trace.csv` and `timeseries.csv`.",
        "",
        "## Preview",
        "",
        "`compile_and_run.cmd config\\live_w20_growth.properties preview`",
        "",
        "Cyan AC, green cells, yellow AHL, garage outline. Close the window;",
        "headless 18000 s is the score. `BSim.preview()` does not stop at 18000 s.",
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
        "T4 is growth-on N(t) / spillover only. Do not start T5 optical,",
        "NARMA, or LuxI unless a new frozen prompt says so. Do not retune",
        "HybridDish K, T0 J_max, or growth. Do not edit HybridDish Java",
        "or any GATE_EVIDENCE.md.",
        "",
        f"`TRANSPORT_MODEL_STATUS={status}`",
    ]
    for arm in ARM_ORDER:
        occ = step[arm]["Occupancy"] if arm in step else "NA"
        pop = step[arm].get("Population", "NA") if arm in step else "NA"
        lines.append(f"`OCCUPANCY_{arm}={occ}`")
        lines.append(f"`POPULATION_{arm}={pop}`")
    lines += [
        "`LIVING_DECISION=STOP_AFTER_T4_POPULATION`",
        "",
    ]
    path = RESULTS / "POPULATION_SCREEN.md"
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
        raise SystemExit("missing population rows; run run_population_screen.py first")
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
            print(f"POPULATION_{arm}={step[arm].get('Population', 'NA')}")
    print("LIVING_DECISION=STOP_AFTER_T4_POPULATION")
    _decision, text = living_text(step)
    print(text)
    if "LIVE_W20_GROWTH" in step:
        w20 = step["LIVE_W20_GROWTH"]
        if w20["Occupancy"] == "DEAD":
            return 4
        pop = w20.get("Population", "")
        if pop not in POP_FLAGS:
            return 5
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
