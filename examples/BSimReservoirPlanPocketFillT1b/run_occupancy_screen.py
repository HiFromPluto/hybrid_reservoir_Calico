#!/usr/bin/env python3
"""PocketFill-T1b seed occupancy screen. No NARMA."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import fill_model as fm  # noqa: E402

RESULTS = HERE / "results"
MAPS = RESULTS / "maps"
ARM_ORDER = ("VOLUMETRIC_ZERO", "SEED_PULSE", "SEED_BATH")


def refuse_narma() -> None:
    if "narma" in " ".join(sys.argv).lower():
        raise SystemExit("PocketFill-T1b cannot see a NARMA target.")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_mean(t: np.ndarray, series: np.ndarray, window: float) -> float:
    mask = t >= (t[-1] - window)
    return float(np.mean(series[mask]))


def run_arm(arm_id: str) -> dict:
    spec = fm.ARMS[arm_id]
    print(f"{arm_id} integrating {spec['t_s']} s seed={spec['seed']} ...", flush=True)
    geom = fm.make_bus_geometry(fm.K_DANINO, 1.0)
    out = fm.integrate_volumetric(
        geom,
        float(spec["t_s"]),
        seed=str(spec["seed"]),
        bath_um=fm.BATH_UM,
    )
    window = float(spec["occupancy_window_s"])
    c_mean_t = np.array([fm.pocket_mean(geom, row) for row in out["C"]])
    la_t = out["Y"][:, 3]
    lux_t = out["Y"][:, 0]
    ahl_in_t = out["Y"][:, 1]
    h_la_t = np.array([float(fm.hill(x, fm.QS_KMLA)) for x in la_t])
    h_c_t = np.array([float(fm.hill(x, fm.QS_KMLA)) for x in c_mean_t])
    mean_c = last_mean(out["t"], c_mean_t, window)
    mean_la = last_mean(out["t"], la_t, window)
    mean_lux = last_mean(out["t"], lux_t, window)
    mean_h_la = last_mean(out["t"], h_la_t, window)
    mean_h_c = last_mean(out["t"], h_c_t, window)
    t_win = out["t"][-1] - window
    i_win = int(np.argmin(np.abs(out["t"] - t_win)))
    rising = bool(la_t[-1] > la_t[i_win])
    flag = fm.occupancy_flag(mean_h_la)
    budget = fm.field_budget(geom, out["sol"], str(spec["seed"]), float(spec["t_s"]))
    ts_path = RESULTS / f"timeseries_{arm_id}.csv"
    with ts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t_s", "C_pocket_mean_uM", "AHL_in_uM", "LuxI", "LA", "H_C", "H_LA"])
        for i, t in enumerate(out["t"]):
            writer.writerow(
                [t, c_mean_t[i], ahl_in_t[i], lux_t[i], la_t[i], h_c_t[i], h_la_t[i]]
            )
    print(
        f"    {arm_id} H(LA)={mean_h_la:.6g} C={mean_c:.6g} LA={mean_la:.6g} "
        f"flag={flag} residual={budget['residual_frac']:.3g}",
        flush=True,
    )
    return {
        "arm": arm_id,
        "kind": spec["kind"],
        "seed": spec["seed"],
        "mean_C": mean_c,
        "mean_LA": mean_la,
        "mean_LuxI": mean_lux,
        "H_LA": mean_h_la,
        "H_C": mean_h_c,
        "occupancy": flag,
        "LA_rising": rising,
        "t_end": float(out["t"][-1]),
        "C_map": fm.unpack(geom, out["C"][-1]),
        "mass_residual_frac": budget["residual_frac"],
        "mass_class": budget["mass_class"],
        "C_nonnegative": budget["C_nonnegative"],
        "initial_molecules": budget["initial_molecules"],
        "membrane_molecules": budget["membrane_molecules"],
        "bath_molecules": budget["bath_molecules"],
        "remaining_molecules": budget["remaining_molecules"],
        "decay_molecules": budget["decay_molecules"],
        "leak_molecules": budget["leak_molecules"],
        "residual_molecules": budget["residual_molecules"],
    }


def save_map(path: Path, field: np.ndarray, title: str) -> None:
    import matplotlib.pyplot as plt

    MAPS.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 3.2))
    im = ax.imshow(np.nan_to_num(field.T, nan=0.0), origin="lower", aspect="auto")
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def fmt(value: float, digits: str = ".6g") -> str:
    if value is None or (isinstance(value, float) and value != value):
        return "NA"
    return format(value, digits)


def main() -> None:
    refuse_narma()
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    for arm_id in ARM_ORDER:
        rec = run_arm(arm_id)
        rows.append(rec)
        save_map(MAPS / f"{rec['arm']}_C.png", rec["C_map"], f"{rec['arm']} C (µM)")

    csv_path = RESULTS / "occupancy_screen.csv"
    keys = [
        "arm",
        "kind",
        "seed",
        "mean_C",
        "mean_LA",
        "mean_LuxI",
        "H_LA",
        "H_C",
        "occupancy",
        "LA_rising",
        "t_end",
        "mass_residual_frac",
        "mass_class",
        "C_nonnegative",
        "initial_molecules",
        "membrane_molecules",
        "bath_molecules",
        "remaining_molecules",
        "decay_molecules",
        "leak_molecules",
        "residual_molecules",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for rec in rows:
            writer.writerow({k: rec.get(k, "") for k in keys})

    hashes = {
        "configs/dish.json": sha256_file(HERE / "configs" / "dish.json"),
        "configs/arms.json": sha256_file(HERE / "configs" / "arms.json"),
        "fill_model.py": sha256_file(HERE / "fill_model.py"),
        "run_occupancy_screen.py": sha256_file(HERE / "run_occupancy_screen.py"),
        "check_pocketfill_t1b.py": sha256_file(HERE / "check_pocketfill_t1b.py"),
    }
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketFill-T1b SHA-256\n\n"
        + "\n".join(f"| `{k}` | `{v}` |" for k, v in hashes.items())
        + "\n",
        encoding="utf-8",
    )

    by = {rec["arm"]: rec for rec in rows}
    zero, pulse, bath = by["VOLUMETRIC_ZERO"], by["SEED_PULSE"], by["SEED_BATH"]
    mass_ok = all(rec["mass_class"] == "LE1PCT" and rec["C_nonnegative"] for rec in rows)
    status = "VALIDATED_FOR_SCREENING" if mass_ok else "NOT_VALIDATED"

    if pulse["occupancy"] == "DEAD" and bath["occupancy"] in ("ALIVE", "SATURATED"):
        reading = (
            "nucleation needs a maintained bath; leak still wins. "
            "That is not a PocketFill occupancy win."
        )
    elif pulse["occupancy"] in ("ALIVE", "SATURATED"):
        reading = "Track C occupancy exists (SEED_PULSE). Still no NARMA from this job."
    elif pulse["occupancy"] == "DEAD" and bath["occupancy"] == "DEAD":
        reading = "pulse and bath both DEAD. Do not retune QS_KMLA or the seed."
    else:
        reading = "see flags."

    md = RESULTS / "OCCUPANCY_SCREEN.md"
    lines = [
        "# PocketFill-T1b occupancy screen",
        "",
        "Reduced OPEN_BUS_3. No NARMA. T1 volumetric recipe + predeclared 0.05 µM seed.",
        "Do not retune QS_KMLA, N_pack, or the seed. Do not start NARMA.",
        "",
        "| Arm | mean C (µM) | mean LA | mean LuxI | H(C) | H(LA) | LA rising | mass residual | Flag |",
        "|---|---:|---:|---:|---:|---:|---|---:|---|",
    ]
    for rec in rows:
        lines.append(
            f"| `{rec['arm']}` | {fmt(rec['mean_C'])} | {fmt(rec['mean_LA'])} | "
            f"{fmt(rec['mean_LuxI'])} | {fmt(rec['H_C'])} | {fmt(rec['H_LA'])} | "
            f"{rec['LA_rising']} | {fmt(rec['mass_residual_frac'], '.3g')} | **{rec['occupancy']}** |"
        )
    lines += [
        "",
        f"`OCCUPANCY_VOLUMETRIC_ZERO={zero['occupancy']}`",
        f"`OCCUPANCY_SEED_PULSE={pulse['occupancy']}`",
        f"`OCCUPANCY_SEED_BATH={bath['occupancy']}`",
        f"`TRANSPORT_MODEL_STATUS={status}`",
        "`LIVING_DECISION=STOP_AFTER_T1B_SEED`",
        "",
        "VOLUMETRIC_ZERO must be T1-class DEAD. SEED_PULSE is the claim arm. "
        "SEED_BATH is ENGINEERING / D1g-class diagnostic.",
        "",
        reading,
        "",
        "If SEED_PULSE is DEAD: stop. Do not retune. Do not start NARMA.",
        "",
    ]
    md.write_text("\n".join(lines), encoding="utf-8")
    payload = {rec["arm"]: rec["occupancy"] for rec in rows}
    payload["TRANSPORT_MODEL_STATUS"] = status
    payload["LIVING_DECISION"] = "STOP_AFTER_T1B_SEED"
    (RESULTS / "occupancy_flags.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(md.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
