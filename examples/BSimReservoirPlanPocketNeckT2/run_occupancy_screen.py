#!/usr/bin/env python3
"""PocketNeck-T2 reduced occupancy screen.

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
ARM_ORDER = ("W100_FLUSH", "W50_L20", "W20_L20", "W10_L20")
HASH_TARGETS = [
    CONFIGS / "dish.json",
    CONFIGS / "arms.json",
    CONFIGS / "drives.json",
    HERE / "transport_model.py",
    HERE / "run_occupancy_screen.py",
    HERE / "check_pocketneck_t2.py",
    HERE / "test_mass_budget.py",
]


def refuse_narma() -> None:
    blob = " ".join(sys.argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("PocketNeck-T2 cannot see a NARMA or Mackey-Glass target.")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_input_hashes() -> dict[str, str]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = {str(path.relative_to(HERE)).replace("\\", "/"): sha256_file(path) for path in HASH_TARGETS}
    path = RESULTS / "INPUT_SHA256.md"
    lines = ["# PocketNeck-T2 input SHA-256", "", "Hashed before occupancy interpretation.", "", "| File | SHA-256 |", "|---|---|"]
    for name, digest in rows.items():
        lines.append(f"| `{name}` | `{digest}` |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def pocket_mean(geom: tm.Geometry, active: np.ndarray) -> float:
    return float(np.mean(tm.pocket_values(geom, active)))


def mask_mean(geom: tm.Geometry, active: np.ndarray, mask: np.ndarray) -> float:
    if not np.any(mask):
        return float("nan")
    field = tm.unpack(geom, active)
    return float(np.nanmean(field[mask]))


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


def metrics_from_traj(geom: tm.Geometry, drive_id: str, traj: dict, smoke: bool) -> dict[str, object]:
    t = traj["t"]
    c, r, lum = traj["C"], traj["R"], traj["L"]
    mean_c = np.array([pocket_mean(geom, row) for row in c])
    mean_r = np.array([pocket_mean(geom, row) for row in r])
    mean_l = np.array([pocket_mean(geom, row) for row in lum])
    decay, boundary, interface = tm.integrate_budget(geom, t, c)
    remaining = tm.mass_molecules(geom, c[-1])
    residual, frac_inj = tm.mass_budget(0.0, traj["injected"], remaining, decay, boundary)
    dominant, frac_dom = tm.residual_vs_dominant(residual, traj["injected"], decay, boundary, interface)
    field_c = tm.unpack(geom, c[-1])
    corners = tm.corner_interior(geom, field_c)
    scales = tm.transport_scales(geom.arm_id)
    plateau = tm.plateau_ok(t, mean_r, 0.10, 0.05) if not smoke else True
    min_c = float(np.min(c))
    occ = tm.occupancy_flag(float(mean_r[-1]))
    numerical = min_c < -1e-9 or frac_dom > 0.01 or frac_inj > 0.01 or not np.all(np.isfinite(c))
    row: dict[str, object] = {
        "Arm": geom.arm_id,
        "Drive": drive_id,
        "Model": "reduced_2D_field_Hill_L",
        "Label": "OCCUPANCY",
        "Smoke": smoke,
        "dx_um": geom.dx,
        "W_um": float(tm.ARMS[geom.arm_id]["W_um"]),
        "L_n_um": float(tm.ARMS[geom.arm_id]["L_n_um"]),
        "DecisionGrade": "PENDING",
        "t_end_s": float(t[-1]),
        "u_on": traj["u_on"],
        "t_on_s": traj["t_on"],
        "PocketMeanAHL_uM": float(mean_c[-1]),
        "PocketMeanR": float(mean_r[-1]),
        "PocketMeanL": float(mean_l[-1]),
        "Occupancy": occ,
        "PlateauR": bool(plateau),
        "NeckMeanAHL_uM": mask_mean(geom, c[-1], geom.neck),
        "BusMeanAHL_uM": mask_mean(geom, c[-1], geom.bus),
        "MinAHL_uM": min_c,
        "SourceAHL_uM": float(field_c[geom.source_ij]),
        "InjectedMass_molecules": traj["injected"],
        "RemainingMass_molecules": remaining,
        "DecayLoss_molecules": decay,
        "BoundaryLoss_molecules": boundary,
        "InterfaceLeak_molecules": interface,
        "MassResidual_molecules": residual,
        "MassResidualFraction": frac_inj,
        "MassResidualVsDominant": frac_dom,
        "DominantFlux_molecules": dominant,
        "C_nonnegative": min_c >= -1e-9,
        "NumericalFailure": numerical,
        "N_t": "NA",
        "L2_over_D_s": scales["L2_over_D_s"],
        "Damkohler": scales["Damkohler"],
        "Pe": scales["Pe"],
        "Transit_s": scales["Transit_s"],
        "FlushT0Class": tm.t0_flush_match(float(mean_r[-1])) if geom.arm_id == "W100_FLUSH" and not smoke else "NA",
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
    return row


def _outline_pocket(ax) -> None:
    ax.plot([150, 250, 250, 150, 150], [0, 0, 100, 100, 0], color="white", lw=0.8, alpha=0.85)


def _outline_neck(ax, w: float) -> None:
    x0 = 200.0 - 0.5 * w
    ax.plot([x0, x0 + w, x0 + w, x0, x0], [100, 100, 120, 120, 100], color="white", lw=0.8, alpha=0.85)


def save_maps(arm: str, row: dict[str, object], geom: tm.Geometry) -> list[str]:
    MAPS.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        np.savez(
            MAPS / f"{arm}_STEP_ON_fields.npz",
            C=row["_field_c"],
            R=row["_field_r"],
            L=row["_field_l"],
        )
        return [str(MAPS / f"{arm}_STEP_ON_fields.npz")]

    x = geom.x_centers
    y = geom.y_centers
    extent = [x[0] - geom.dx / 2, x[-1] + geom.dx / 2, y[0] - geom.dy / 2, y[-1] + geom.dy / 2]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), constrained_layout=True)
    titles = [
        (r"AHL ($\mu$M)", row["_field_c"], "viridis", None),
        (r"$R$", row["_field_r"], "magma", (0.0, 1.0)),
    ]
    w = float(row["W_um"])
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
        _outline_pocket(ax)
        if np.any(geom.neck):
            _outline_neck(ax, w)
    fig.suptitle(f"{arm} STEP_ON  occupancy={row['Occupancy']}", fontsize=11)
    path = MAPS / f"{arm}_STEP_ON_C_R.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [str(path)]


def save_r_strip(rows: list[dict[str, object]], geoms: dict[str, tm.Geometry]) -> None:
    occ = [row for row in rows if row["Label"] == "OCCUPANCY" and not row["Smoke"]]
    if len(occ) != 4:
        return
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8), constrained_layout=True)
    for ax, row in zip(axes, occ):
        geom = geoms[str(row["Arm"])]
        field = np.ma.masked_invalid(row["_field_r"])
        x, y = geom.x_centers, geom.y_centers
        extent = [x[0] - geom.dx / 2, x[-1] + geom.dx / 2, y[0] - geom.dy / 2, y[-1] + geom.dy / 2]
        image = ax.imshow(field.T, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=1, aspect="equal")
        ax.set_title(f"{row['Arm']}\n{row['Occupancy']}")
        ax.set_xlabel("x (µm)")
        if ax is axes[0]:
            ax.set_ylabel("y (µm)")
        _outline_pocket(ax)
        if np.any(geom.neck):
            _outline_neck(ax, float(row["W_um"]))
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("STEP_ON $R$ (reduced model; occupancy = pocket voxels)")
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


def run_mass_budget_unit_tests() -> None:
    import test_mass_budget as tests

    tests.main()


def run_arm(arm: str, *, smoke: bool, sample_dt: float) -> tuple[tm.Geometry, dict[str, object]]:
    geom = tm.make_geometry(arm)
    tm.print_voxel_check(geom)
    print(f"arm {arm}  n_active={geom.n}  source={geom.source_ij}  n_pocket={int(geom.pocket.sum())}", flush=True)
    kwargs: dict = {}
    if smoke:
        kwargs = {"duration_on": 50.0, "duration_off": 0.0}
    print("  STEP_ON ...", flush=True)
    traj = tm.simulate_drive(geom, "STEP_ON", sample_dt=sample_dt, **kwargs)
    row = metrics_from_traj(geom, "STEP_ON", traj, smoke)
    save_maps(arm, row, geom)
    neck_s = row["NeckMeanAHL_uM"]
    neck_txt = "NA" if not math.isfinite(float(neck_s)) else f"{float(neck_s):.4g}"
    print(
        f"    meanR={row['PocketMeanR']:.4g}  occ={row['Occupancy']}  "
        f"AHL={row['PocketMeanAHL_uM']:.4g}  neckC={neck_txt}  "
        f"mass_dom={row['MassResidualVsDominant']:.3g}",
        flush=True,
    )
    return geom, row


def run_convergence_w10(sample_dt: float) -> dict[str, object]:
    print("CONVERGENCE W10_L20 dx=2.5 um (not occupancy) ...", flush=True)
    tm.apply_grid(2.5)
    try:
        geom = tm.make_geometry("W10_L20")
        tm.print_voxel_check(geom)
        traj = tm.simulate_drive(geom, "STEP_ON", sample_dt=sample_dt)
        row = metrics_from_traj(geom, "STEP_ON", traj, smoke=False)
        row["Label"] = "CONVERGENCE"
        row["Arm"] = "W10_L20_dx2.5"
        save_maps("W10_L20_dx2.5", row, geom)
        print(
            f"    CONVERGENCE meanR={row['PocketMeanR']:.4g}  occ={row['Occupancy']} "
            f"(do not replace the dx=5 table)",
            flush=True,
        )
        return row
    finally:
        tm.restore_primary_grid()


def main() -> int:
    refuse_narma()
    parser = argparse.ArgumentParser(description="PocketNeck-T2 reduced occupancy screen")
    parser.add_argument("--smoke", action="store_true", help="50 s STEP_ON; not occupancy evidence")
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument("--skip-convergence", action="store_true", help="skip optional W10_L20 dx=2.5 CONVERGENCE")
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
    sample_dt = 20.0
    rows: list[dict[str, object]] = []
    geoms: dict[str, tm.Geometry] = {}

    geom, flush = run_arm("W100_FLUSH", smoke=smoke, sample_dt=sample_dt)
    geoms["W100_FLUSH"] = geom
    rows.append(flush)
    if not smoke:
        if flush["Occupancy"] == "ALIVE" or flush["Occupancy"] == "SATURATED" or not flush["FlushT0Class"]:
            write_csv(RESULTS / "occupancy_screen.csv", rows)
            print("STOP. W100_FLUSH is not T0 OPEN_BUS_3 class. Fix the flush clone before interpreting necks.")
            print(f"OCCUPANCY_W100_FLUSH={flush['Occupancy']}")
            print(f"PocketMeanR={flush['PocketMeanR']}")
            return 3

    for arm in ARM_ORDER[1:]:
        geom, row = run_arm(arm, smoke=smoke, sample_dt=sample_dt)
        geoms[arm] = geom
        rows.append(row)

    save_r_strip(rows, geoms)

    conv_meta: dict[str, object] = {"ran": False, "skipped": True}
    if not smoke and not args.skip_convergence:
        conv_row = run_convergence_w10(sample_dt)
        rows.append(conv_row)
        conv_meta = {
            "ran": True,
            "skipped": False,
            "dx_um": 2.5,
            "PocketMeanR": conv_row["PocketMeanR"],
            "Occupancy": conv_row["Occupancy"],
            "note": "CONVERGENCE only. Do not replace the dx=5 occupancy table.",
        }
    elif args.skip_convergence:
        conv_meta = {"ran": False, "skipped": True, "reason": "--skip-convergence"}
    (RESULTS / "convergence.json").write_text(json.dumps(conv_meta, indent=2) + "\n", encoding="utf-8")

    occ_rows = [row for row in rows if row["Label"] == "OCCUPANCY"]
    mass_ok = all(not row["NumericalFailure"] for row in occ_rows)
    c_ok = all(bool(row["C_nonnegative"]) for row in occ_rows)
    plateau_ok = all(bool(row["PlateauR"]) for row in occ_rows) if not smoke else True
    status = "VALIDATED_FOR_SCREENING" if mass_ok and c_ok and plateau_ok and form["units_ok"] else "NOT_VALIDATED"
    if smoke:
        status = "SMOKE_NOT_EVIDENCE"
    for row in rows:
        row["TRANSPORT_MODEL_STATUS"] = status
        if smoke:
            row["DecisionGrade"] = "SMOKE_NOT_EVIDENCE"
        elif row["Label"] == "CONVERGENCE":
            row["DecisionGrade"] = "CONVERGENCE_NOT_OCCUPANCY"
        elif status != "VALIDATED_FOR_SCREENING":
            row["DecisionGrade"] = "NON_DECISION_GRADE"
        else:
            row["DecisionGrade"] = "DECISION_GRADE_PHYSICAL_SCREEN"
    write_csv(RESULTS / "occupancy_screen.csv", rows)
    print(f"TRANSPORT_MODEL_STATUS={status}")
    for row in occ_rows:
        print(f"OCCUPANCY_{row['Arm']}={row['Occupancy']}")
    print("LIVING_DECISION=STOP_AFTER_NECK_SCREEN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
