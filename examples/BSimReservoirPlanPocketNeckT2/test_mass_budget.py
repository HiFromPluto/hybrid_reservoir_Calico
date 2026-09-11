#!/usr/bin/env python3
"""Mass / voxel / flush-clone tests for PocketNeck-T2. No NARMA."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import transport_model as tm  # noqa: E402

ARM_ORDER = ("W100_FLUSH", "W50_L20", "W20_L20", "W10_L20")


def _assert_rel(name: str, value: float, ceiling: float) -> None:
    if not math.isfinite(value) or value > ceiling:
        raise RuntimeError(f"mass-budget unit test {name} failed: {value} > {ceiling}")


def _assert(name: str, cond: bool, detail: str = "") -> None:
    if not cond:
        raise RuntimeError(f"mass-budget unit test {name} failed {detail}".rstrip())


def main() -> None:
    # E: well-mixed closed C_ss(u=1) algebra (T0; pocket V unchanged).
    _assert_rel("E_css", abs(tm.closed_form_css(1.0) - 3.2) / 3.2, 0.01)
    _assert_rel("E_css_u05", abs(tm.closed_form_css(0.5) - 1.6) / 1.6, 0.01)

    # Voxel integers and pocket volume not eaten.
    expected_width_seen = []
    for arm in ARM_ORDER:
        geom = tm.make_geometry(arm)
        row = tm.print_voxel_check(geom)
        _assert("pocket_V", bool(row["pocket_volume_ok"]), arm)
        if tm.ARMS[arm]["neck"]:
            _assert("W_dx_int", bool(row["W_over_dx_integer"]), arm)
            _assert("Ln_dx_int", bool(row["Ln_over_dx_integer"]), arm)
            _assert("Ln_over_dx_4", abs(float(row["Ln_over_dx"]) - 4.0) < 1e-9, arm)
            expected_width_seen.append(int(round(float(row["W_over_dx"]))))
            _assert("flush_has_no_neck", True)
        else:
            _assert("flush_n_neck_0", int(row["n_neck"]) == 0, arm)
            _assert("flush_iface_100um", int(row["n_interface"]) == 20, f"{row['n_interface']}")
    _assert("neck_widths_10_4_2", expected_width_seen == [10, 4, 2], str(expected_width_seen))

    # Short W100_FLUSH STEP_ON residual (T0 identity). Not occupancy evidence.
    geom = tm.make_geometry("W100_FLUSH")
    traj = tm.simulate_drive(geom, "STEP_ON", duration_on=200.0, duration_off=0.0, sample_dt=20.0)
    decay, boundary, iface = tm.integrate_budget(geom, traj["t"], traj["C"])
    remaining = tm.mass_molecules(geom, traj["C"][-1])
    residual, frac = tm.mass_budget(0.0, traj["injected"], remaining, decay, boundary)
    _dom, frac_dom = tm.residual_vs_dominant(residual, traj["injected"], decay, boundary, iface)
    _assert_rel("flush_short_vs_injected", frac, 0.01)
    _assert_rel("flush_short_vs_dominant", frac_dom, 0.01)
    _assert("C_nonneg_short", float(np.min(traj["C"])) >= -1e-9)

    # Mixing a second horizon into the identity must not close (T0 test D).
    traj2 = tm.simulate_drive(geom, "STEP_ON", duration_on=50.0, duration_off=0.0, sample_dt=10.0)
    mixed = traj["injected"] + traj2["injected"]
    _, wrong = tm.mass_budget(0.0, mixed, remaining, decay, boundary)
    if wrong <= 0.01:
        raise RuntimeError("mass-budget unit test D failed: mixed-horizon identity unexpectedly closed")

    print("PocketNeck-T2 mass-budget / voxel tests passed.")


if __name__ == "__main__":
    main()
