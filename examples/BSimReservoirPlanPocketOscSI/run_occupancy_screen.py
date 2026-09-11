#!/usr/bin/env python3
"""PocketOsc-SI bulk delay-DDE screen. Period vs mu. No NARMA."""

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

import si_dde as dde  # noqa: E402

RESULTS = HERE / "results"
ARM_ORDER = ("MU_0", "MU_025", "MU_05", "MU_1", "MU_2", "MU_CSTR_MIN")


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("PocketOsc-SI cannot see a NARMA target.")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_arm(arm_id: str) -> dict:
    spec = dde.ARMS[arm_id]
    d = float(spec["d"])
    mu = float(spec["mu"])
    print(f"{arm_id} integrating T={dde.T_END} d={d} mu={mu} D1=0 ...", flush=True)
    out = dde.integrate_bulk(d, mu)
    y = out["Y"]
    peaks = dde.i_peaks(out["t"], y[:, 1])
    means = dde.post_means(out["t"], y)
    period = peaks["period"]
    osc = bool(peaks["n_peaks"] >= dde.MIN_PEAKS and math.isfinite(period))
    ts_path = RESULTS / f"timeseries_{arm_id}.csv"
    with ts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t", "A", "I", "H_i", "H_e"])
        for i, t in enumerate(out["t"]):
            writer.writerow([t, y[i, 0], y[i, 1], y[i, 2], y[i, 3]])
    print(
        f"    {arm_id} mean_I={means['mean_I']:.6g} mean_A={means['mean_A']:.6g} "
        f"mean_Hi={means['mean_Hi']:.6g} mean_He={means['mean_He']:.6g} "
        f"n_peaks={peaks['n_peaks']} T={period if math.isfinite(period) else 'NO_PERIOD'} "
        f"OSC={osc}",
        flush=True,
    )
    return {
        "arm": arm_id,
        "role": spec["role"],
        "d": d,
        "mu": mu,
        "D1": 0.0,
        "mean_A": means["mean_A"],
        "mean_I": means["mean_I"],
        "mean_Hi": means["mean_Hi"],
        "mean_He": means["mean_He"],
        "n_peaks": peaks["n_peaks"],
        "period": period,
        "osc": osc,
        "flag": "OSC" if osc else "NO_PERIOD",
    }


def fmt(value: float) -> str:
    if isinstance(value, float) and not math.isfinite(value):
        return "NA"
    return format(value, ".6g")


def living_decision(cover: list[dict]) -> tuple[str, str]:
    osc_arms = [rec for rec in cover if rec["osc"]]
    if not osc_arms:
        return (
            "STOP_AFTER_SI_BULK_DEAD",
            "No COVER arm oscillates. Cited SI bulk delay-DDE is dead on this "
            "predeclared mu line. Stop. Do not retune alpha, tau, or gamma_*. "
            "Do not revive mu_bus. No GFP maps. No NARMA.",
        )
    if len(osc_arms) < 2:
        return (
            "STOP_AFTER_SI_NO_IDENTITY",
            "Fewer than two COVER arms oscillate, so T(mu) cannot be scored. "
            "Stop. Do not add mu points. Do not fit 55/90 min.",
        )
    periods = [float(rec["period"]) for rec in osc_arms]
    increasing = all(periods[i] < periods[i + 1] for i in range(len(periods) - 1))
    if increasing:
        labels = ", ".join(f"{rec['arm']} T={fmt(rec['period'])}" for rec in osc_arms)
        return (
            "STOP_AFTER_SI_BULK_ALIVE",
            "Identity holds: at least two COVER arms OSC and T strictly increases "
            f"with mu ({labels}). Model twin ALIVE. Still not Fig. 2c. Still no "
            "GFP maps. No NARMA. Stop. No spatial add-on.",
        )
    return (
        "STOP_AFTER_SI_NO_IDENTITY",
        "COVER arms oscillate but T does not strictly lengthen with mu. "
        "Stop. Do not retune. Do not fit 55/90 min. No NARMA.",
    )


def main() -> None:
    refuse_narma()
    RESULTS.mkdir(parents=True, exist_ok=True)
    dde.print_ics()
    rows = [run_arm(arm_id) for arm_id in ARM_ORDER]
    keys = [
        "arm", "role", "d", "mu", "D1", "mean_A", "mean_I", "mean_Hi", "mean_He",
        "n_peaks", "period", "osc", "flag",
    ]
    with (RESULTS / "occupancy_screen.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for rec in rows:
            writer.writerow({k: rec.get(k, "") for k in keys})

    hashes = {
        "PROTOCOL.md": sha256_file(HERE / "PROTOCOL.md"),
        "configs/params.json": sha256_file(HERE / "configs" / "params.json"),
        "configs/arms.json": sha256_file(HERE / "configs" / "arms.json"),
        "si_dde.py": sha256_file(HERE / "si_dde.py"),
        "run_occupancy_screen.py": sha256_file(HERE / "run_occupancy_screen.py"),
    }
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketOsc-SI SHA-256\n\n"
        + "\n".join(f"| `{k}` | `{v}` |" for k, v in hashes.items())
        + "\n",
        encoding="utf-8",
    )
    cover = [rec for rec in rows if rec["arm"] in dde.COVER_IDS]
    decision, reading = living_decision(cover)

    lines = [
        "# PocketOsc-SI occupancy screen",
        "",
        "SI bulk delay-DDE at d=0.5. Period vs mu. TAKEN table. No NARMA.",
        "Time is SI-scaled. Readout is I (LuxI). Do not fit 55/90 min.",
        "ICs ENGINEERING: A=I=0, Hi=He=1, Hi(t<0)=1.",
        "",
        "| Arm | role | mu | mean I | mean A | mean Hi | mean He | n_peaks | T | Flag |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for rec in rows:
        lines.append(
            f"| `{rec['arm']}` | {rec['role']} | {fmt(rec['mu'])} | "
            f"{fmt(rec['mean_I'])} | {fmt(rec['mean_A'])} | {fmt(rec['mean_Hi'])} | "
            f"{fmt(rec['mean_He'])} | {rec['n_peaks']} | {fmt(rec['period'])} | "
            f"**{rec['flag']}** |"
        )
    osc_cover = [rec["arm"] for rec in cover if rec["osc"]]
    lines += [
        "",
        f"`OSC_COVER={','.join(osc_cover) if osc_cover else 'NONE'}`",
        f"`LIVING_DECISION={decision}`",
        "",
        reading,
        "",
    ]
    (RESULTS / "OCCUPANCY_SCREEN.md").write_text("\n".join(lines), encoding="utf-8")
    payload = {rec["arm"]: rec["flag"] for rec in rows}
    payload["OSC_COVER"] = osc_cover
    payload["LIVING_DECISION"] = decision
    (RESULTS / "occupancy_flags.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print((RESULTS / "OCCUPANCY_SCREEN.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
