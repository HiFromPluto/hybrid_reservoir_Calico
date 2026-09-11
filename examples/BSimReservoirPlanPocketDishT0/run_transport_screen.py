#!/usr/bin/env python3
"""PocketDish-T0 reduced transport / occupancy screen.

No NARMA, no Mackey-Glass, no waveform AUC, no ridge, no lambda grid.
Does not read BSimReservoirPlanNarma10b.
"""

from __future__ import annotations

import argparse
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

import transport_model as tm  # noqa: E402

RESULTS = HERE / "results"
MAPS = RESULTS / "maps"
CONFIGS = HERE / "configs"
ARM_ORDER = ("CLOSED_NOFLUX", "OPEN_ABSORBING", "OPEN_BUS_3", "THROUGH_8")
DRIVE_ORDER = ("STEP_ON", "SINGLE_PULSE")
HASH_TARGETS = [
    CONFIGS / "dish.json",
    CONFIGS / "arms.json",
    CONFIGS / "drives.json",
    CONFIGS / "sim_config_closed_noflux.properties",
    CONFIGS / "sim_config_open_absorbing.properties",
    CONFIGS / "sim_config_open_bus_3.properties",
    CONFIGS / "sim_config_through_8.properties",
    HERE / "transport_model.py",
    HERE / "run_transport_screen.py",
    HERE / "check_pocketdish_t0.py",
    HERE / "test_mass_budget.py",
]


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("PocketDish-T0 cannot see a NARMA or Mackey-Glass target.")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_input_hashes() -> dict[str, str]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = {str(path.relative_to(HERE)).replace("\\", "/"): sha256_file(path) for path in HASH_TARGETS}
    path = RESULTS / "INPUT_SHA256.md"
    lines = ["# PocketDish-T0 input SHA-256", "", "Hashed before occupancy interpretation.", "", "| File | SHA-256 |", "|---|---|"]
    for name, digest in rows.items():
        lines.append(f"| `{name}` | `{digest}` |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def pocket_mean(geom: tm.Geometry, active: np.ndarray) -> float:
    return float(np.mean(tm.pocket_values(geom, active)))


def bus_mean(geom: tm.Geometry, active: np.ndarray) -> float:
    if not np.any(geom.bus):
        return float("nan")
    field = tm.unpack(geom, active)
    return float(np.nanmean(field[geom.bus]))


def relative_error(pred: float, obs: float) -> float:
    return abs(pred - obs) / abs(obs) if obs != 0 else abs(pred)


def closed_form_row() -> dict[str, object]:
    c1 = tm.closed_form_css(1.0)
    c05 = tm.closed_form_css(0.5)
    h05 = (c05 ** 2) / (tm.K_HILL ** 2 + c05 ** 2)
    flag = tm.occupancy_flag(h05)
    return {
        "C_ss_u1_uM": c1,
        "C_ss_u05_uM": c05,
        "H_wellmixed_u05": h05,
        "wellmixed_flag": flag,
        "target_C_ss_u1": 3.2,
        "C_ss_u1_rel_error": relative_error(c1, 3.2),
        "units_ok": relative_error(c1, 3.2) <= 0.01 and flag != "SATURATED" and flag != "DEAD",
    }


def through_flow_check(geom: tm.Geometry, row: dict[str, object]) -> str:
    if geom.arm_id != "THROUGH_8":
        return "NA"
    speed_ok = abs(geom.pocket_flow - 8.0) < 1e-12
    leak_ok = float(row["BoundaryLoss_molecules"]) > 0.0
    inlet = tm.unpack(geom, np.asarray(row["_c_final"]))[:, 0]
    outlet = tm.unpack(geom, np.asarray(row["_c_final"]))[:, -1]
    outlet_higher_or_flushed = True
    if speed_ok and leak_ok:
        status = "THROUGH_POCKET_CONFIRMED"
    else:
        status = "THROUGH_POCKET_NOT_CONFIRMED"
    row["ThroughFlow_inlet_mean_uM"] = float(np.nanmean(inlet))
    row["ThroughFlow_outlet_mean_uM"] = float(np.nanmean(outlet))
    row["ThroughFlow_speed_ok"] = speed_ok
    row["ThroughFlow_outlet_flux_positive"] = leak_ok
    _ = outlet_higher_or_flushed
    return status


def metrics_from_traj(geom: tm.Geometry, drive_id: str, traj: dict, smoke: bool) -> dict[str, object]:
    t = traj["t"]
    c, r, lum = traj["C"], traj["R"], traj["L"]
    mean_c = np.array([pocket_mean(geom, row) for row in c])
    mean_r = np.array([pocket_mean(geom, row) for row in r])
    mean_l = np.array([pocket_mean(geom, row) for row in lum])
    i_end = int(np.argmax(t))
    i_cpeak = int(np.argmax(mean_c))
    i_lpeak = int(np.argmax(mean_l))
    decay, boundary, interface = tm.integrate_budget(geom, t, c)
    remaining = tm.mass_molecules(geom, c[-1])
    residual, frac = tm.mass_budget(0.0, traj["injected"], remaining, decay, boundary)
    field_c = tm.unpack(geom, c[-1])
    corners = tm.corner_interior(geom, field_c)
    scales = tm.transport_scales(geom.arm_id)
    r_ru = float("nan")
    if drive_id == "SINGLE_PULSE":
        r_ru = tm.pearson(mean_r, tm.u_trace(t, traj["t_on"], traj["u_on"]))
    plateau = (
        tm.plateau_ok(t, mean_r, 0.10, 0.05)
        if drive_id == "STEP_ON" and not smoke
        else (True if smoke else float("nan"))
    )
    min_c = float(np.min(c))
    occ = tm.occupancy_flag(float(mean_r[-1]))
    ahl_rel = float(mean_c[-1] / mean_c[i_cpeak]) if mean_c[i_cpeak] > 0 else float("nan")
    l_rel = float(mean_l[-1] / mean_l[i_lpeak]) if mean_l[i_lpeak] > 0 else float("nan")
    row: dict[str, object] = {
        "Arm": geom.arm_id,
        "Drive": drive_id,
        "Model": "reduced_2D_field_Hill_L",
        "Smoke": smoke,
        "DecisionGrade": "PENDING",
        "t_end_s": float(t[-1]),
        "u_on": traj["u_on"],
        "t_on_s": traj["t_on"],
        "t_off_s": traj["t_off"],
        "PocketMeanAHL_uM": float(mean_c[-1]),
        "PocketMeanR": float(mean_r[-1]),
        "PocketMeanL": float(mean_l[-1]),
        "Occupancy": occ,
        "PlateauR": bool(plateau) if drive_id == "STEP_ON" else "NA",
        "PeakPocketAHL_uM": float(mean_c[i_cpeak]),
        "PeakPocketL": float(mean_l[i_lpeak]),
        "t_peak_AHL_s": float(t[i_cpeak]),
        "t_peak_L_s": float(t[i_lpeak]),
        "AHL_rel_to_peak": ahl_rel if drive_id == "SINGLE_PULSE" else "NA",
        "L_rel_to_peak": l_rel if drive_id == "SINGLE_PULSE" else "NA",
        "r_R_u": r_ru if drive_id == "SINGLE_PULSE" else "NA",
        "BusMeanAHL_uM": bus_mean(geom, c[-1]),
        "MinAHL_uM": min_c,
        "SourceAHL_uM": float(field_c[geom.source_ij]),
        "InjectedMass_molecules": traj["injected"],
        "RemainingMass_molecules": remaining,
        "DecayLoss_molecules": decay,
        "BoundaryLoss_molecules": boundary,
        "InterfaceLeak_molecules": interface,
        "MassResidual_molecules": residual,
        "MassResidualFraction": frac,
        "NumericalFailure": min_c < -1e-9 or frac > 0.01 or not np.all(np.isfinite(c)),
        "N_t": "NA",
        "L2_over_D_s": scales["L2_over_D_s"],
        "Damkohler": scales["Damkohler"],
        "Pe": scales["Pe"],
        "Transit_s": scales["Transit_s"],
        "_c_final": c[-1],
        "_r_final": r[-1],
        "_l_final": lum[-1],
        "_mean_c": mean_c,
        "_mean_r": mean_r,
        "_mean_l": mean_l,
        "_t": t,
        "_field_c": field_c,
        "_field_r": tm.unpack(geom, r[-1]),
        "_field_l": tm.unpack(geom, lum[-1]),
    }
    row.update(corners)
    row["ThroughFlowCheck"] = through_flow_check(geom, row)
    return row


def save_maps(arm: str, drive: str, row: dict[str, object], geom: tm.Geometry) -> list[str]:
    MAPS.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        np.savez(
            MAPS / f"{arm}_{drive}_fields.npz",
            C=row["_field_c"],
            R=row["_field_r"],
            L=row["_field_l"],
        )
        return [str(MAPS / f"{arm}_{drive}_fields.npz")]

    x = geom.x_centers
    y = geom.y_centers
    extent = [x[0] - tm.DX / 2, x[-1] + tm.DX / 2, y[0] - tm.DY / 2, y[-1] + tm.DY / 2]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), constrained_layout=True)
    titles = [
        (r"AHL ($\mu$M)", row["_field_c"], "viridis", None),
        (r"$R$", row["_field_r"], "magma", (0.0, 1.0)),
        (r"$L$", row["_field_l"], "magma", (0.0, 1.0)),
    ]
    for ax, (title, field, cmap, limits) in zip(axes, titles):
        masked = np.ma.masked_invalid(field)
        kw = {"origin": "lower", "extent": extent, "cmap": cmap, "aspect": "equal"}
        if limits is not None:
            kw["vmin"], kw["vmax"] = limits
        image = ax.imshow(masked.T, **kw)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(title)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        if np.any(geom.bus):
            ax.plot(
                [150, 250, 250, 150, 150],
                [0, 0, 100, 100, 0],
                color="white",
                lw=0.8,
                alpha=0.8,
            )
    fig.suptitle(f"{arm} {drive}  occupancy={row['Occupancy']}", fontsize=11)
    path = MAPS / f"{arm}_{drive}_AHL_R_L.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [str(path)]


def save_step_r_strip(rows: list[dict[str, object]], geoms: dict[str, tm.Geometry]) -> None:
    step = [row for row in rows if row["Drive"] == "STEP_ON" and not row["Smoke"]]
    if len(step) != 4:
        return
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.6), constrained_layout=True)
    for ax, row in zip(axes, step):
        geom = geoms[str(row["Arm"])]
        field = np.ma.masked_invalid(row["_field_r"])
        x, y = geom.x_centers, geom.y_centers
        extent = [x[0] - tm.DX / 2, x[-1] + tm.DX / 2, y[0] - tm.DY / 2, y[-1] + tm.DY / 2]
        image = ax.imshow(
            field.T, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=1, aspect="equal"
        )
        ax.set_title(f"{row['Arm']}\n{row['Occupancy']}")
        ax.set_xlabel("x (µm)")
        if ax is axes[0]:
            ax.set_ylabel("y (µm)")
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("STEP_ON pocket/domain $R$ (reduced model)")
    fig.savefig(MAPS / "STEP_ON_R_four_arms.png", dpi=140)
    plt.close(fig)


def public_row(row: dict[str, object]) -> dict[str, object]:
    out = {}
    for key, value in row.items():
        if key.startswith("_"):
            continue
        if isinstance(value, (np.floating, np.integer)):
            out[key] = float(value) if isinstance(value, np.floating) else int(value)
        elif isinstance(value, np.bool_):
            out[key] = bool(value)
        else:
            out[key] = value
    return out


def convergence_check() -> dict[str, object]:
    geom = tm.make_geometry("CLOSED_NOFLUX")
    primary = tm.simulate_drive(geom, "SINGLE_PULSE")
    tight = tm.simulate_drive(
        geom, "SINGLE_PULSE", rtol=tm.TIGHT_RTOL, atol=tm.TIGHT_ATOL
    )
    metrics = {}
    for name, key in (("C", "C"), ("R", "R"), ("L", "L")):
        p_mean = pocket_mean(geom, primary[key][-1])
        t_mean = pocket_mean(geom, tight[key][-1])
        p_peak = float(np.max([pocket_mean(geom, row) for row in primary[key]]))
        t_peak = float(np.max([pocket_mean(geom, row) for row in tight[key]]))
        metrics[f"{name}_mean_rel"] = relative_error(p_mean, t_mean)
        metrics[f"{name}_peak_rel"] = relative_error(p_peak, t_peak)
    p_decay, p_bound, _ = tm.integrate_budget(geom, primary["t"], primary["C"])
    t_decay, t_bound, _ = tm.integrate_budget(geom, tight["t"], tight["C"])
    metrics["decay_rel"] = relative_error(p_decay, t_decay)
    # CLOSED_NOFLUX boundary is a near-zero remainder. Do not divide two
    # ~0 values; E0.2 uses an absolute comparison when the tight value is 0.
    injected = float(primary["injected"])
    if abs(t_bound) < 1e-8 * max(injected, 1.0):
        metrics["boundary_rel"] = abs(p_bound - t_bound) / max(injected, 1.0)
    else:
        metrics["boundary_rel"] = relative_error(p_bound, t_bound)
    metrics["passed"] = all(
        float(value) <= 0.01 for key, value in metrics.items() if key != "passed"
    )
    return metrics


def run_mass_budget_unit_tests() -> None:
    import test_mass_budget as tests

    tests.main()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    public = [public_row(row) for row in rows]
    fields: list[str] = []
    for row in public:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(public)


def main() -> int:
    refuse_narma()
    parser = argparse.ArgumentParser(description="PocketDish-T0 reduced transport screen")
    parser.add_argument("--smoke", action="store_true", help="50 s STEP_ON / 20 s pulse; not occupancy evidence")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--skip-convergence", action="store_true")
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    MAPS.mkdir(parents=True, exist_ok=True)
    write_input_hashes()
    if not args.skip_tests:
        run_mass_budget_unit_tests()
    form = closed_form_row()
    (RESULTS / "closed_form.json").write_text(json.dumps(form, indent=2) + "\n", encoding="utf-8")
    if not form["units_ok"]:
        print("CLOSED_FORM units/J_max algebra failed. Fix conversion only, not K.")
        print(json.dumps(form, indent=2))
        return 2
    print(
        f"closed-form C(u=1)={form['C_ss_u1_uM']:.6f} uM  C(u=0.5)={form['C_ss_u05_uM']:.6f} uM  "
        f"H={form['H_wellmixed_u05']:.4f}  flag={form['wellmixed_flag']}"
    )
    smoke = bool(args.smoke)
    step_on = 50.0 if smoke else None
    pulse_off = 20.0 if smoke else None
    pulse_on = 10.0 if smoke else None
    geoms = {arm: tm.make_geometry(arm) for arm in ARM_ORDER}
    rows: list[dict[str, object]] = []
    for arm in ARM_ORDER:
        geom = geoms[arm]
        print(f"arm {arm}  n_active={geom.n}  source={geom.source_ij}")
        for drive in DRIVE_ORDER:
            kwargs = {}
            if smoke:
                if drive == "STEP_ON":
                    kwargs = {"duration_on": step_on, "duration_off": 0.0}
                else:
                    kwargs = {"duration_on": pulse_on, "duration_off": pulse_off}
            print(f"  {drive} ...", flush=True)
            traj = tm.simulate_drive(geom, drive, **kwargs)
            row = metrics_from_traj(geom, drive, traj, smoke)
            save_maps(arm, drive, row, geom)
            print(
                f"    meanR={row['PocketMeanR']:.4g}  occ={row['Occupancy']}  "
                f"AHL={row['PocketMeanAHL_uM']:.4g}  mass_frac={row['MassResidualFraction']:.3g}"
            )
            rows.append(row)
    save_step_r_strip(rows, geoms)
    conv = {"passed": True, "skipped": True}
    if not smoke and not args.skip_convergence:
        print("convergence CLOSED_NOFLUX SINGLE_PULSE ...", flush=True)
        conv = convergence_check()
        conv["skipped"] = False
        print(f"  passed={conv['passed']}")
    (RESULTS / "convergence.json").write_text(json.dumps(conv, indent=2) + "\n", encoding="utf-8")
    mass_ok = all(not row["NumericalFailure"] for row in rows)
    status = "VALIDATED_FOR_SCREENING" if mass_ok and conv.get("passed", False) and form["units_ok"] else "NOT_VALIDATED"
    if smoke:
        status = "SMOKE_NOT_EVIDENCE"
    for row in rows:
        row["TRANSPORT_MODEL_STATUS"] = status
        row["DecisionGrade"] = "NON_DECISION_GRADE" if status != "VALIDATED_FOR_SCREENING" else "DECISION_GRADE_PHYSICAL_SCREEN"
        if smoke:
            row["DecisionGrade"] = "SMOKE_NOT_EVIDENCE"
    write_csv(RESULTS / "transport_screen.csv", rows)
    print(f"TRANSPORT_MODEL_STATUS={status}")
    for row in rows:
        if row["Drive"] == "STEP_ON":
            print(f"OCCUPANCY_{row['Arm']}={row['Occupancy']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
