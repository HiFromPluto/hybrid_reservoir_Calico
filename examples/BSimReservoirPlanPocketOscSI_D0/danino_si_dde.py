#!/usr/bin/env python3
"""Danino 2010 SI bulk delay-DDE oracle (D0).

Transcribed from 41586_2010_BFnature08753_MOESM346_ESM.pdf Modeling,
not from DANINO_SI_MODEL.md or PocketOsc-SI.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
CONFIGS = HERE / "configs"

PRODUCTION_SI_HILL = "si_hill"
PRODUCTION_POCKETOSCSI_SPLIT = "pocketoscsi_split"


def load_json(name: str) -> dict:
    with (CONFIGS / name).open(encoding="utf-8") as handle:
        return json.load(handle)


def production(h_tau: float, delta: float, alpha: float, k1: float, form: str) -> float:
    h2 = h_tau * h_tau
    denom = 1.0 + k1 * h2
    if form == PRODUCTION_SI_HILL:
        return (delta + alpha * h2) / denom
    if form == PRODUCTION_POCKETOSCSI_SPLIT:
        return delta + alpha * h2 / denom
    raise ValueError(f"unknown production form {form!r}")


class DelayTape:
    """Piecewise dense interpolant for Hi(t). Raises on a history gap.

    Hi(t<0) is the frozen constant history. The t=0 state is y0, which
    may differ from that history.
    """

    def __init__(self, y0: np.ndarray, hi_history: float):
        self.y0 = np.asarray(y0, dtype=float)
        self.hi_history = float(hi_history)
        self.segments: list[tuple[float, float, object]] = []

    def add(self, t_left: float, t_right: float, sol: object) -> None:
        if t_right < t_left:
            raise RuntimeError(f"empty segment [{t_left}, {t_right}]")
        if self.segments:
            prev_end = self.segments[-1][1]
            if t_left < prev_end - 1e-12:
                raise RuntimeError(
                    f"overlapping delay segments ending {prev_end} starting {t_left}"
                )
            if t_left > prev_end + 1e-8:
                raise RuntimeError(
                    f"delay interpolant gap ({prev_end}, {t_left}); history was not skipped silently"
                )
        self.segments.append((float(t_left), float(t_right), sol))

    def state(self, t: float) -> np.ndarray:
        if t <= 1e-15:
            return self.y0.copy()
        for t0, t1, sol in reversed(self.segments):
            if t0 - 1e-14 <= t <= t1 + 1e-12:
                tt = min(max(t, t0), t1)
                return np.asarray(sol.sol(tt), dtype=float)
        raise RuntimeError(
            f"delay interpolant: no recorded history at t={t}; lookup cannot skip"
        )

    def hi(self, t: float) -> float:
        if t < 0.0:
            return self.hi_history
        return float(self.state(t)[2])


class DaninoSIDDE:
    def __init__(
        self,
        params: dict | None = None,
        production_form: str | None = None,
        rtol: float | None = None,
        atol: float | None = None,
        max_step: float | None = None,
    ) -> None:
        self.params = params if params is not None else load_json("params.json")
        p = self.params
        self.C_A = float(p["C_A"])
        self.C_I = float(p["C_I"])
        self.delta = float(p["delta"])
        self.alpha = float(p["alpha"])
        self.tau = float(p["tau"])
        self.k = float(p["k"])
        self.k1 = float(p["k1"])
        self.b = float(p["b"])
        self.gamma_A = float(p["gamma_A"])
        self.gamma_I = float(p["gamma_I"])
        self.gamma_H = float(p["gamma_H"])
        self.f = float(p["f"])
        self.g = float(p["g"])
        self.d0 = float(p["d0"])
        self.D_mem = float(p["D_membrane"])
        integ = p["integrator"]
        self.production_form = production_form or str(p["production_form"])
        self.rtol = float(integ["rtol"] if rtol is None else rtol)
        self.atol = float(integ["atol"] if atol is None else atol)
        self.max_step = float(integ["max_step"] if max_step is None else max_step)

    def p_of(self, h_tau: float) -> float:
        return production(h_tau, self.delta, self.alpha, self.k1, self.production_form)

    def rhs(self, y: np.ndarray, h_tau: float, d: float, mu: float) -> np.ndarray:
        a, i, h_i, h_e = y
        p = self.p_of(h_tau)
        dens = 1.0 - (d / self.d0) ** 4
        denom_f = 1.0 + self.f * (a + i)
        da = self.C_A * dens * p - self.gamma_A * a / denom_f
        di = self.C_I * dens * p - self.gamma_I * i / denom_f
        dhi = (
            self.b * i / (1.0 + self.k * i)
            - (self.gamma_H * a * h_i) / (1.0 + self.g * a)
            + self.D_mem * (h_e - h_i)
        )
        dhe = -(d / (1.0 - d)) * self.D_mem * (h_e - h_i) - mu * h_e
        return np.array([da, di, dhi, dhe], dtype=float)

    def integrate(
        self,
        d: float,
        mu: float,
        y0: np.ndarray,
        hi_history: float,
        t_end: float,
        sample_dt: float,
    ) -> dict:
        y0 = np.asarray(y0, dtype=float)
        if y0.shape != (4,):
            raise ValueError("y0 must be (A, I, Hi, He)")
        tape = DelayTape(y0, float(hi_history))
        tau = self.tau

        def fun(t: float, y: np.ndarray) -> np.ndarray:
            return self.rhs(y, tape.hi(t - tau), d, mu)

        t_left = 0.0
        y_left = y0.copy()
        while t_left < t_end - 1e-15:
            t_right = min(t_left + tau, t_end)
            sol = solve_ivp(
                fun,
                (t_left, t_right),
                y_left,
                method="BDF",
                rtol=self.rtol,
                atol=self.atol,
                dense_output=True,
                max_step=self.max_step,
            )
            if not sol.success:
                raise RuntimeError(sol.message)
            tape.add(t_left, t_right, sol)
            y_left = sol.y[:, -1]
            t_left = t_right

        n = int(round(t_end / sample_dt))
        t_eval = np.linspace(0.0, t_end, n + 1)
        y_out = np.vstack([tape.state(t) for t in t_eval])
        return {
            "t": t_eval,
            "Y": y_out,
            "d": float(d),
            "mu": float(mu),
            "D1": 0.0,
            "production_form": self.production_form,
            "y0": y0.copy(),
            "H_i_history": float(hi_history),
            "rtol": self.rtol,
            "atol": self.atol,
            "negative_state": bool(np.any(y_out < -1e-12)),
        }
