#!/usr/bin/env python3
"""PocketFill-T1 occupancy screen. No NARMA."""

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

import fill_model as fm  # noqa: E402

RESULTS = HERE / "results"
MAPS = RESULTS / "maps"


def refuse_narma() -> None:
    if "narma" in " ".join(sys.argv).lower():
        raise SystemExit("PocketFill-T1 cannot see a NARMA target.")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_mean(t: np.ndarray, series: np.ndarray, window: float) -> float:
    mask = t >= (t[-1] - window)
    return float(np.mean(series[mask]))


def run_point(arm_id: str) -> dict:
    spec = fm.ARMS[arm_id]
    geom = fm.make_bus_geometry(fm.K_HYBRID, float(spec["D_pocket_factor"]))
    out = fm.integrate_point(geom, float(spec["u"]), float(spec["t_s"]))
    c_end = out["C"][-1]
    r_end = out["R"][-1]
    mean_c = fm.pocket_mean(geom, c_end)
    mean_r = fm.pocket_mean(geom, r_end)
    flag = fm.occupancy_flag(mean_r)
    field = fm.unpack(geom, c_end)
    return {
        "arm": arm_id,
        "kind": spec["kind"],
        "mean_C": mean_c,
        "mean_R": mean_r,
        "mean_LA": float("nan"),
        "H_LA": float("nan"),
        "H_C": float(fm.hill(mean_c, fm.K_HILL)),
        "occupancy": flag,
        "t_end": float(out["t"][-1]),
        "C_map": field,
        "geom": geom,
    }


def run_volumetric() -> dict:
    spec = fm.ARMS["VOLUMETRIC_LUXI"]
    geom = fm.make_bus_geometry(fm.K_DANINO, 1.0)
    out = fm.integrate_volumetric(geom, float(spec["t_s"]))
    window = float(spec["occupancy_window_s"])
    c_mean_t = np.array([fm.pocket_mean(geom, row) for row in out["C"]])
    la_t = out["Y"][:, 3]
    lux_t = out["Y"][:, 0]
    h_la_t = np.array([float(fm.hill(x, fm.QS_KMLA)) for x in la_t])
    h_c_t = np.array([float(fm.hill(x, fm.QS_KMLA)) for x in c_mean_t])
    mean_c = last_mean(out["t"], c_mean_t, window)
    mean_la = last_mean(out["t"], la_t, window)
    mean_h_la = last_mean(out["t"], h_la_t, window)
    mean_h_c = last_mean(out["t"], h_c_t, window)
    rising = bool(la_t[-1] > la_t[int(np.argmin(np.abs(out["t"] - (out["t"][-1] - window))))])
    flag = fm.occupancy_flag(mean_h_la)
    return {
        "arm": "VOLUMETRIC_LUXI",
        "kind": spec["kind"],
        "mean_C": mean_c,
        "mean_R": mean_h_c,
        "mean_LA": mean_la,
        "mean_LuxI": last_mean(out["t"], lux_t, window),
        "H_LA": mean_h_la,
        "H_C": mean_h_c,
        "occupancy": flag,
        "LA_rising": rising,
        "t_end": float(out["t"][-1]),
        "C_map": fm.unpack(geom, out["C"][-1]),
        "t": out["t"],
        "C_mean_t": c_mean_t,
        "LA_t": la_t,
        "geom": geom,
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


def main() -> None:
    refuse_narma()
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    empty = run_point("EMPTY_POINT")
    plug = run_point("PACKED_PLUG")
    vol = run_volumetric()
    for rec in (empty, plug, vol):
        rows.append(rec)
        save_map(MAPS / f"{rec['arm']}_C.png", rec["C_map"], f"{rec['arm']} C (µM)")

    csv_path = RESULTS / "occupancy_screen.csv"
    keys = ["arm", "kind", "mean_C", "mean_R", "mean_LA", "H_LA", "H_C", "occupancy", "t_end"]
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
    }
    (RESULTS / "INPUT_SHA256.md").write_text(
        "# PocketFill-T1 SHA-256\n\n"
        + "\n".join(f"| `{k}` | `{v}` |" for k, v in hashes.items())
        + "\n",
        encoding="utf-8",
    )

    md = RESULTS / "OCCUPANCY_SCREEN.md"
    lines = [
        "# PocketFill-T1 occupancy screen",
        "",
        "Reduced OPEN_BUS_3. No NARMA. D1g QS_* write the field.",
        "Do not retune QS_KMLA, N_pack, or T0 Jmax.",
        "",
        "| Arm | mean C (µM) | HybridDish R or H(C) | mean LA | H(LA) | Flag |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for rec in rows:
        la = rec.get("mean_LA")
        hla = rec.get("H_LA")
        la_s = "NA" if la is None or (isinstance(la, float) and la != la) else f"{la:.6g}"
        hla_s = "NA" if hla is None or (isinstance(hla, float) and hla != hla) else f"{hla:.6g}"
        lines.append(
            f"| `{rec['arm']}` | {rec['mean_C']:.6g} | {rec['mean_R']:.6g} | {la_s} | {hla_s} | **{rec['occupancy']}** |"
        )
    lines += [
        "",
        f"`OCCUPANCY_EMPTY_POINT={empty['occupancy']}`",
        f"`OCCUPANCY_PACKED_PLUG={plug['occupancy']}`",
        f"`OCCUPANCY_VOLUMETRIC_LUXI={vol['occupancy']}`",
        "",
        "EMPTY_POINT must be DEAD (T0 class). VOLUMETRIC primary gate is H(LA) over last 7200 s.",
        f"LA still rising: {vol.get('LA_rising')}.",
        "",
        "If VOLUMETRIC_LUXI is DEAD: stop. Do not retune. Do not start NARMA.",
        "",
    ]
    md.write_text("\n".join(lines), encoding="utf-8")
    print(md.read_text(encoding="utf-8"))
    payload = {rec["arm"]: rec["occupancy"] for rec in rows}
    (RESULTS / "occupancy_flags.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
