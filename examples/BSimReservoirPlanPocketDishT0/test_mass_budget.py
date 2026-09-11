#!/usr/bin/env python3
"""Mass-budget unit tests A-E for PocketDish-T0. No NARMA."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import transport_model as tm  # noqa: E402


def _assert_rel(name: str, value: float, ceiling: float) -> None:
    if not math.isfinite(value) or value > ceiling:
        raise RuntimeError(f"mass-budget unit test {name} failed: {value} > {ceiling}")


def main() -> None:
    k = tm.K_DECAY

    # A: no decay, no source — mass conserved.
    tm.K_DECAY = 0.0  # type: ignore[misc]
    try:
        geom_a = tm.make_geometry("CLOSED_NOFLUX")
        rng = np.random.default_rng(0)
        blob = rng.random(geom_a.n) * 0.4 + 0.1
        mass0 = tm.mass_molecules(geom_a, blob)
        out = tm.integrate_concentration(geom_a, lambda _t: 0.0, 400.0, c0=blob, sample_dt=20.0)
        mass1 = tm.mass_molecules(geom_a, out["C"][-1])
        _assert_rel("A", abs(mass1 - mass0) / mass0, 1e-6)
    finally:
        tm.K_DECAY = k  # type: ignore[misc]

    # B: uniform C, decay only — exponential.
    geom = tm.make_geometry("CLOSED_NOFLUX")
    c_uniform = np.full(geom.n, 2.0)
    duration = 400.0
    mass0 = tm.mass_molecules(geom, c_uniform)
    out = tm.integrate_concentration(geom, lambda _t: 0.0, duration, c0=c_uniform, sample_dt=10.0)
    expected = mass0 * math.exp(-k * duration)
    got = tm.mass_molecules(geom, out["C"][-1])
    _assert_rel("B", abs(got - expected) / expected, 1e-5)

    # C: closed SINGLE_PULSE budget residual; boundary ~ 0.
    geom = tm.make_geometry("CLOSED_NOFLUX")
    traj = tm.simulate_drive(geom, "SINGLE_PULSE")
    decay, boundary, _iface = tm.integrate_budget(geom, traj["t"], traj["C"])
    remaining = tm.mass_molecules(geom, traj["C"][-1])
    _residual, frac = tm.mass_budget(0.0, traj["injected"], remaining, decay, boundary)
    _assert_rel("C", frac, 0.01)
    _assert_rel("C_boundary_closed", abs(boundary) / max(traj["injected"], 1.0), 0.01)

    # D: mixing STEP_ON injection into the pulse identity must not close.
    step = tm.simulate_drive(geom, "STEP_ON", duration_on=200.0, duration_off=0.0)
    pulse = traj
    decay_s, bound_s, _ = tm.integrate_budget(geom, step["t"], step["C"])
    rem_s = tm.mass_molecules(geom, step["C"][-1])
    _, frac_s = tm.mass_budget(0.0, step["injected"], rem_s, decay_s, bound_s)
    _assert_rel("D_step", frac_s, 0.01)
    mixed_injected = step["injected"] + pulse["injected"]
    decay_p, bound_p, _ = tm.integrate_budget(geom, pulse["t"], pulse["C"])
    rem_p = tm.mass_molecules(geom, pulse["C"][-1])
    _, wrong = tm.mass_budget(0.0, mixed_injected, rem_p, decay_p, bound_p)
    if wrong <= 0.01:
        raise RuntimeError("mass-budget unit test D failed: mixed-horizon identity unexpectedly closed")

    # E: well-mixed closed C_ss(u=1) algebra.
    _assert_rel("E_css", abs(tm.closed_form_css(1.0) - 3.2) / 3.2, 0.01)
    print("PocketDish-T0 mass-budget unit tests A-E passed.")


if __name__ == "__main__":
    main()
