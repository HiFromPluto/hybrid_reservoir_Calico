#!/usr/bin/env python3
"""PocketFill-T1c door-vs-circuit occupancy. No NARMA."""

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
ARM_ORDER = ("VOLUMETRIC_FLUSH", "VOLUMETRIC_CLOSED", "VOLUMETRIC_W20")


def refuse_narma() -> None:
    if "narma" in " ".join(sys.argv).lower():
        raise SystemExit("PocketFill-T1c cannot see a NARMA target.")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_mean(t: np.ndarray, series: np.ndarray, window: float) -> float:
    mask = t >= (t[-1] - window)
    return float(np.mean(series[mask]))


def run_arm(arm_id: str) -> dict:
    spec = fm.ARMS[arm_id]
    print(f"{arm_id} integrating {spec['t_s']} s geom={spec['kind']} ...", flush=True)
    geom = fm.make_geometry(str(spec["kind"]))
    out = fm.integrate_volumetric(geom, float(spec["t_s"]))
    window = float(spec["occupancy_window_s"])
    c_mean_t = np.array([fm.pocket_mean(geom, row) for row in out["C"]])
    la_t = out["Y"][:, 3]
    lux_t = out["Y"][:, 0]
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
    budget = fm.field_budget(geom, out["sol"], float(spec["t_s"]))
    ts_path = RESULTS / f"timeseries_{arm_id}.csv"
    with ts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t_s", "C_pocket_mean_uM", "LuxI", "LA", "H_C", "H_LA"])
        for i, t in enumerate(out["t"]):
            writer.writerow([t, c_mean_t[i], lux_t[i], la_t[i], h_c_t[i], h_la_t[i]])
    print(
        f"    {arm_id} H(LA)={mean_h_la:.6g} C={mean_c:.6g} flag={flag} "
        f"residual={budget['residual_frac']:.3g} leak={budget['leak_molecules']:.3g}",
        flush=True,
    )
    return {
        "arm": arm_id,
        "kind": spec["kind"],
        "mean_C": mean_c,
        "mean_LA": mean_la,
        "mean_LuxI": mean_lux,
        "H_LA": mean_h_la,
        "H_C": mean_h_c,
        "occupancy": flag,
        "LA_rising": rising,
        "mass_residual_frac": budget["residual_frac"],
        "mass_class": budget["mass_class"],
        "C_nonnegative": budget["C_nonnegative"],
        "leak_molecules": budget["leak_molecules"],
        "decay_molecules": budget["decay_molecules"],
        "membrane_molecules": budget["membrane_molecules"],
    }


def fmt(value: float) -> str:
    return format(value, ".6g")


def main() -> None:
    refuse_narma()
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = [run_arm(arm_id) for arm_id in ARM_ORDER]
    keys = [
        "arm", "kind", "mean_C", "mean_LA", "mean_LuxI", "H_LA", "H_C",
        "occupancy", "LA_rising", "mass_residual_frac", "mass_class",
        "C_nonnegative", "leak_molecules", "decay_molecules", "membrane_molecules",
    ]
    with (RESULTS / "occupancy_screen.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for rec in rows:
            writer.writerow({k: rec.get(k, "") for k in keys})

    hashes = {
        "configs/dish.json": sha256_file(HERE / "configs" / "dish.json"),
        "configs/arms.json": sha256_file(HERE / "configs" / "arms.json"),
        "fill_model.py": sha256_file(HERE / "fill_model.py"),
        "run_occupancy_screen.py": sha256_file(HERE / "run_occupancy_screen.py"),
    }
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketFill-T1c SHA-256\n\n"
        + "\n".join(f"| `{k}` | `{v}` |" for k, v in hashes.items())
        + "\n",
        encoding="utf-8",
    )
    by = {rec["arm"]: rec for rec in rows}
    flush, closed, w20 = by["VOLUMETRIC_FLUSH"], by["VOLUMETRIC_CLOSED"], by["VOLUMETRIC_W20"]
    mass_ok = all(rec["mass_class"] == "LE1PCT" and rec["C_nonnegative"] for rec in rows)
    status = "VALIDATED_FOR_SCREENING" if mass_ok else "NOT_VALIDATED"
    closed_flag = closed["occupancy"]
    if flush["occupancy"] in ("ALIVE", "SATURATED"):
        decision = "STOP_PORT_BROKE"
        reading = "FLUSH ALIVE: T1 port broke. Do not celebrate."
    elif closed_flag == "DEAD":
        decision = "STOP_AFTER_T1C_CLOSED_DEAD"
        reading = (
            "CLOSED is DEAD. Cited pack does not occupy even without a street. "
            "Stop Track C. Do not ease QS_KMLA."
        )
    elif w20["occupancy"] == "DEAD":
        decision = "STOP_AFTER_T1C_LEAK_BUDGET"
        reading = (
            "CLOSED occupies and W20 does not. Leak budget is the story. "
            "Do not add extra W overnight."
        )
    else:
        decision = "STOP_AFTER_T1C_W20_ALIVE"
        reading = (
            "CLOSED and W20 occupy. Track C occupancy exists on the living neck. "
            "Still no NARMA from this job."
        )

    lines = [
        "# PocketFill-T1c occupancy screen",
        "",
        "Door vs circuit. Same QS_* as T1/T1b. No NARMA. Do not retune QS_KMLA.",
        "",
        "| Arm | geom | mean C (µM) | mean LA | H(LA) | leak mol | Flag |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for rec in rows:
        lines.append(
            f"| `{rec['arm']}` | {rec['kind']} | {fmt(rec['mean_C'])} | "
            f"{fmt(rec['mean_LA'])} | {fmt(rec['H_LA'])} | "
            f"{fmt(rec['leak_molecules'])} | **{rec['occupancy']}** |"
        )
    lines += [
        "",
        f"`OCCUPANCY_VOLUMETRIC_FLUSH={flush['occupancy']}`",
        f"`OCCUPANCY_VOLUMETRIC_CLOSED={closed['occupancy']}`",
        f"`OCCUPANCY_VOLUMETRIC_W20={w20['occupancy']}`",
        f"`TRANSPORT_MODEL_STATUS={status}`",
        f"`LIVING_DECISION={decision}`",
        "",
        reading,
        "",
    ]
    (RESULTS / "OCCUPANCY_SCREEN.md").write_text("\n".join(lines), encoding="utf-8")
    payload = {rec["arm"]: rec["occupancy"] for rec in rows}
    payload["TRANSPORT_MODEL_STATUS"] = status
    payload["LIVING_DECISION"] = decision
    (RESULTS / "occupancy_flags.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print((RESULTS / "OCCUPANCY_SCREEN.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
