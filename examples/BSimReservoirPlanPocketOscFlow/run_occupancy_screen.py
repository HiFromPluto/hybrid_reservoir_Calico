#!/usr/bin/env python3
"""PocketOsc-Flow occupancy + period. No NARMA."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import osc_model as om  # noqa: E402

RESULTS = HERE / "results"
ARM_ORDER = ("CLOSED_D50", "FLOW_180", "FLOW_296")


def refuse_narma() -> None:
    if "narma" in " ".join(sys.argv).lower():
        raise SystemExit("PocketOsc-Flow cannot see a NARMA target.")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_mean(t: np.ndarray, series: np.ndarray, t_end: float, window: float) -> float:
    mask = (t >= (t_end - window)) & (t <= t_end + 1e-9)
    return float(np.mean(series[mask]))


def run_arm(arm_id: str) -> dict:
    spec = om.ARMS[arm_id]
    d = float(spec["d"])
    height = float(spec["height_um"])
    v = float(spec["v_um_min"]) / 60.0
    duration = om.T_END
    print(f"{arm_id} integrating {duration} s d={d} h={height} v={v:.4g} um/s ...", flush=True)
    out = om.integrate_arm(d, height, v, duration)
    pack = out["pack"]
    la_t = out["Y"][:, 3]
    lux_t = out["Y"][:, 0]
    c_t = np.asarray(out["C"], dtype=float)
    h_la_t = np.array([float(om.hill_qs(x)) for x in la_t])
    peaks = om.la_peaks(out["t"], la_t)
    mean_c = last_mean(out["t"], c_t, duration, om.OCC_WINDOW)
    mean_la = last_mean(out["t"], la_t, duration, om.OCC_WINDOW)
    mean_lux = last_mean(out["t"], lux_t, duration, om.OCC_WINDOW)
    mean_h_la = last_mean(out["t"], h_la_t, duration, om.OCC_WINDOW)
    d50_h = last_mean(out["t"], h_la_t, om.D50_DURATION, om.OCC_WINDOW)
    flag = om.occupancy_flag(mean_h_la)
    d50_flag = om.occupancy_flag(d50_h)
    budget = om.field_budget(out, duration)
    ts_path = RESULTS / f"timeseries_{arm_id}.csv"
    with ts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t_s", "C_uM", "LuxI", "AHL_in", "AiiA", "LA", "H_LA"])
        for i, t in enumerate(out["t"]):
            writer.writerow([t, c_t[i], lux_t[i], out["Y"][i, 1], out["Y"][i, 2], la_t[i], h_la_t[i]])
    period = peaks["period_s"]
    period_min = period / 60.0 if math.isfinite(period) else float("nan")
    print(
        f"    {arm_id} H(LA)={mean_h_la:.6g} flag={flag} mu={out['mu_bus']:.6g} "
        f"n_peaks={peaks['n_peaks']} T_min={period_min:.6g} residual={budget['residual_frac']:.3g}",
        flush=True,
    )
    return {
        "arm": arm_id,
        "role": spec["role"],
        "d": pack["d"],
        "height_um": pack["height_um"],
        "v_um_s": v,
        "v_um_min": float(spec["v_um_min"]),
        "mu_bus_s": out["mu_bus"],
        "mean_C": mean_c,
        "mean_LA": mean_la,
        "mean_LuxI": mean_lux,
        "H_LA": mean_h_la,
        "H_LA_d50_window": d50_h,
        "occupancy": flag,
        "occupancy_d50_window": d50_flag,
        "n_peaks": peaks["n_peaks"],
        "period_s": period,
        "period_min": period_min,
        "mass_residual_frac": budget["residual_frac"],
        "mass_class": budget["mass_class"],
        "C_nonnegative": budget["C_nonnegative"],
        "leak_molecules": budget["leak_molecules"],
        "decay_molecules": budget["decay_molecules"],
        "membrane_molecules": budget["membrane_molecules"],
    }


def fmt(value: float) -> str:
    if isinstance(value, float) and not math.isfinite(value):
        return "NA"
    return format(value, ".6g")


def main() -> None:
    refuse_narma()
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = [run_arm(arm_id) for arm_id in ARM_ORDER]
    keys = [
        "arm", "role", "d", "height_um", "v_um_s", "v_um_min", "mu_bus_s",
        "mean_C", "mean_LA", "mean_LuxI", "H_LA", "H_LA_d50_window",
        "occupancy", "occupancy_d50_window", "n_peaks", "period_s", "period_min",
        "mass_residual_frac", "mass_class", "C_nonnegative",
        "leak_molecules", "decay_molecules", "membrane_molecules",
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
        "PROTOCOL.md": sha256_file(HERE / "PROTOCOL.md"),
    }
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketOsc-Flow SHA-256\n\n"
        + "\n".join(f"| `{k}` | `{v}` |" for k, v in hashes.items())
        + "\n",
        encoding="utf-8",
    )
    by = {rec["arm"]: rec for rec in rows}
    closed, lo, hi = by["CLOSED_D50"], by["FLOW_180"], by["FLOW_296"]
    mass_ok = all(rec["mass_class"] == "LE1PCT" and rec["C_nonnegative"] for rec in rows)
    status = "VALIDATED_FOR_SCREENING" if mass_ok else "NOT_VALIDATED"

    def has_period(rec: dict) -> bool:
        return rec["n_peaks"] >= 3 and math.isfinite(float(rec["period_s"]))

    if closed["occupancy_d50_window"] == "DEAD":
        decision = "STOP_PORT_BROKE"
        reading = "CLOSED D50 window DEAD: D50 port broke. Do not celebrate flow."
    elif lo["occupancy"] == "DEAD" or hi["occupancy"] == "DEAD":
        decision = "STOP_AFTER_FLOW_DEAD"
        reading = (
            "An open-flow arm is DEAD. Cited monolayer does not stay occupied "
            "with the frozen bus leak. Stop. Do not ease QS_KMLA. Do not raise d. "
            "Do not close the door and call that identity."
        )
    elif not has_period(lo) or not has_period(hi):
        decision = "STOP_AFTER_FLOW_NO_IDENTITY"
        reading = (
            "Open arms occupy but periods do not exist on both flows. "
            "Stop. Do not add delay. Do not ease QS_*."
        )
    elif float(hi["period_s"]) > float(lo["period_s"]):
        decision = "STOP_AFTER_FLOW_IDENTITY"
        reading = (
            "T(296)>T(180) in the paper's direction. Occupancy exists with a bus. "
            "Matching 55/90 min is not a retune ticket. Still no GFP maps. No NARMA."
        )
    else:
        decision = "STOP_AFTER_FLOW_NO_IDENTITY"
        reading = (
            "Periods exist but do not order as the paper (T should lengthen with flow). "
            "Stop. Do not ease QS_*. Do not retune mu_bus."
        )

    lines = [
        "# PocketOsc-Flow occupancy screen",
        "",
        "Period vs channel flow at cited d=0.5, 1.65 um. Same QS_*. No NARMA.",
        "mu_bus is ENGINEERING two-resistance; not fitted to 55/90 min.",
        "",
        "| Arm | v (µm/min) | µ_bus (1/s) | mean C (µM) | H(LA) | n_peaks | T (min) | Flag |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for rec in rows:
        lines.append(
            f"| `{rec['arm']}` | {fmt(rec['v_um_min'])} | {fmt(rec['mu_bus_s'])} | "
            f"{fmt(rec['mean_C'])} | {fmt(rec['H_LA'])} | {rec['n_peaks']} | "
            f"{fmt(rec['period_min'])} | **{rec['occupancy']}** |"
        )
    lines += [
        "",
        f"`OCCUPANCY_CLOSED_D50={closed['occupancy']}`",
        f"`OCCUPANCY_CLOSED_D50_WINDOW={closed['occupancy_d50_window']}`",
        f"`OCCUPANCY_FLOW_180={lo['occupancy']}`",
        f"`OCCUPANCY_FLOW_296={hi['occupancy']}`",
        f"`TRANSPORT_MODEL_STATUS={status}`",
        f"`LIVING_DECISION={decision}`",
        "",
        f"Closed D50-window H(LA)={fmt(closed['H_LA_d50_window'])} "
        f"(expect ALIVE). Closed leak={fmt(closed['leak_molecules'])}.",
        f"FLOW_180 T={fmt(lo['period_min'])} min; FLOW_296 T={fmt(hi['period_min'])} min.",
        "",
        reading,
        "",
    ]
    (RESULTS / "OCCUPANCY_SCREEN.md").write_text("\n".join(lines), encoding="utf-8")
    payload = {rec["arm"]: rec["occupancy"] for rec in rows}
    payload["CLOSED_D50_WINDOW"] = closed["occupancy_d50_window"]
    payload["TRANSPORT_MODEL_STATUS"] = status
    payload["LIVING_DECISION"] = decision
    (RESULTS / "occupancy_flags.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print((RESULTS / "OCCUPANCY_SCREEN.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
