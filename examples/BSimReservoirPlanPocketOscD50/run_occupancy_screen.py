#!/usr/bin/env python3
"""PocketOsc-D50 closed occupancy. No NARMA."""

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

import osc_model as om  # noqa: E402

RESULTS = HERE / "results"
ARM_ORDER = ("D05_H165_CLOSED", "D005_H10_REPLAY")


def refuse_narma() -> None:
    if "narma" in " ".join(sys.argv).lower():
        raise SystemExit("PocketOsc-D50 cannot see a NARMA target.")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_mean(t: np.ndarray, series: np.ndarray, window: float) -> float:
    mask = t >= (t[-1] - window)
    return float(np.mean(series[mask]))


def run_arm(arm_id: str) -> dict:
    spec = om.ARMS[arm_id]
    d = float(spec["d"])
    height = float(spec["height_um"])
    duration = float(spec["t_s"])
    window = float(spec["occupancy_window_s"])
    print(f"{arm_id} integrating {duration} s d={d} h={height} ...", flush=True)
    out = om.integrate_closed(d, height, duration)
    pack = out["pack"]
    la_t = out["Y"][:, 3]
    lux_t = out["Y"][:, 0]
    c_t = np.asarray(out["C"], dtype=float)
    h_la_t = np.array([float(om.hill_qs(x)) for x in la_t])
    h_c_t = np.array([float(om.hill_qs(x)) for x in c_t])
    mean_c = last_mean(out["t"], c_t, window)
    mean_la = last_mean(out["t"], la_t, window)
    mean_lux = last_mean(out["t"], lux_t, window)
    mean_h_la = last_mean(out["t"], h_la_t, window)
    mean_h_c = last_mean(out["t"], h_c_t, window)
    t_win = out["t"][-1] - window
    i_win = int(np.argmin(np.abs(out["t"] - t_win)))
    rising = bool(la_t[-1] > la_t[i_win])
    flag = om.occupancy_flag(mean_h_la)
    budget = om.field_budget(out, duration)
    ts_path = RESULTS / f"timeseries_{arm_id}.csv"
    with ts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t_s", "C_uM", "LuxI", "AHL_in", "AiiA", "LA", "H_C", "H_LA"])
        for i, t in enumerate(out["t"]):
            writer.writerow(
                [
                    t,
                    c_t[i],
                    lux_t[i],
                    out["Y"][i, 1],
                    out["Y"][i, 2],
                    la_t[i],
                    h_c_t[i],
                    h_la_t[i],
                ]
            )
    print(
        f"    {arm_id} H(LA)={mean_h_la:.6g} C={mean_c:.6g} flag={flag} "
        f"residual={budget['residual_frac']:.3g} leak={budget['leak_molecules']:.3g}",
        flush=True,
    )
    return {
        "arm": arm_id,
        "role": spec["role"],
        "d": pack["d"],
        "d_over_one_minus_d": pack["d_over_one_minus_d"],
        "height_um": pack["height_um"],
        "V_trap_um3": pack["V_trap_um3"],
        "N_implied_engineering": pack["N_implied_engineering"],
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
        "arm", "role", "d", "d_over_one_minus_d", "height_um", "V_trap_um3",
        "N_implied_engineering", "mean_C", "mean_LA", "mean_LuxI", "H_LA", "H_C",
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
        "osc_model.py": sha256_file(HERE / "osc_model.py"),
        "run_occupancy_screen.py": sha256_file(HERE / "run_occupancy_screen.py"),
    }
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketOsc-D50 SHA-256\n\n"
        + "\n".join(f"| `{k}` | `{v}` |" for k, v in hashes.items())
        + "\n",
        encoding="utf-8",
    )
    by = {rec["arm"]: rec for rec in rows}
    claim, replay = by["D05_H165_CLOSED"], by["D005_H10_REPLAY"]
    mass_ok = all(rec["mass_class"] == "LE1PCT" and rec["C_nonnegative"] for rec in rows)
    status = "VALIDATED_FOR_SCREENING" if mass_ok else "NOT_VALIDATED"
    if replay["occupancy"] in ("ALIVE", "SATURATED"):
        decision = "STOP_PORT_BROKE"
        reading = "REPLAY ALIVE: T1c-class port broke. Do not celebrate."
    elif claim["occupancy"] == "DEAD":
        decision = "STOP_AFTER_D50"
        reading = (
            "CLAIM d=0.5 at 1.65 um is DEAD. Cited monolayer does not occupy "
            "optical LA. Stop. Do not ease QS_KMLA. Do not raise d. No NARMA."
        )
    elif claim["occupancy"] == "SATURATED":
        decision = "STOP_AFTER_D50_SATURATED"
        reading = (
            "CLAIM SATURATED. Occupancy exists; not a useful graded clock from "
            "this job. Still no period-vs-flow. Do not lower QS_KMLA."
        )
    else:
        decision = "STOP_AFTER_D50_ALIVE"
        reading = (
            "CLAIM ALIVE at cited d=0.5, 1.65 um, closed, no bath, no AC. "
            "Occupancy exists. Still no period-vs-flow from this chat. No NARMA."
        )

    lines = [
        "# PocketOsc-D50 occupancy screen",
        "",
        "Cited monolayer vs T1c-class pack. Same QS_*. Closed. No NARMA.",
        "Do not retune QS_KMLA. Do not raise d. Do not add 0.95 um.",
        "",
        "| Arm | d | h (µm) | d/(1-d) | V_trap (µm³) | mean C (µM) | mean LA | H(LA) | Flag |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for rec in rows:
        lines.append(
            f"| `{rec['arm']}` | {fmt(rec['d'])} | {fmt(rec['height_um'])} | "
            f"{fmt(rec['d_over_one_minus_d'])} | {fmt(rec['V_trap_um3'])} | "
            f"{fmt(rec['mean_C'])} | {fmt(rec['mean_LA'])} | {fmt(rec['H_LA'])} | "
            f"**{rec['occupancy']}** |"
        )
    lines += [
        "",
        f"`OCCUPANCY_D05_H165_CLOSED={claim['occupancy']}`",
        f"`OCCUPANCY_D005_H10_REPLAY={replay['occupancy']}`",
        f"`TRANSPORT_MODEL_STATUS={status}`",
        f"`LIVING_DECISION={decision}`",
        "",
        f"Claim: d={fmt(claim['d'])}, d/(1-d)={fmt(claim['d_over_one_minus_d'])}, "
        f"V_trap={fmt(claim['V_trap_um3'])} µm³, "
        f"N_implied={fmt(claim['N_implied_engineering'])} (ENGINEERING).",
        f"Replay leak={fmt(replay['leak_molecules'])}; claim leak={fmt(claim['leak_molecules'])}.",
        f"Claim mass residual={fmt(claim['mass_residual_frac'])}; "
        f"LA_rising={claim['LA_rising']}.",
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
