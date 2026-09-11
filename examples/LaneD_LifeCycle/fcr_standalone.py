"""Standalone Erickson 2017 FCR (Eq. 5 / S47). Not the D0 bath."""

from __future__ import annotations

import math
from dataclasses import dataclass


PHI_RB0 = 0.049
GAMMA = 11.02
LAMBDA_C = 1.17
LAMBDA_F = 0.45
LAMBDA_I = 0.90


def phi_rb_star(lam: float) -> float:
    return PHI_RB0 + lam / GAMMA


def sigma_star(lam: float) -> float:
    return lam / phi_rb_star(lam)


def mu_of(lam: float) -> float:
    return lam / (1.0 - lam / LAMBDA_C)


def chi_rb(sigma: float) -> float:
    return PHI_RB0 / (1.0 - sigma / GAMMA)


def chi_cat(sigma: float) -> float:
    return 1.0 - (sigma / LAMBDA_C) * chi_rb(sigma)


def dsigma_eq5(sigma: float, mu_f: float) -> float:
    return sigma * (mu_f * chi_cat(sigma) - sigma * chi_rb(sigma))


def dsigma_s47(sigma: float, mu_f: float, sigma_f: float) -> float:
    return mu_f * sigma * (1.0 - sigma / sigma_f) / (1.0 - sigma / GAMMA)


def sigma0_s57(sigma_i: float, lam_i: float, lam_f: float) -> float:
    """Drop form of SI S57 (paper: flux and growth rate drop)."""
    return sigma_i * (LAMBDA_C - lam_i) / (LAMBDA_C - lam_f)


def s49_lhs(sigma: float, sigma0: float, sigma_f: float, mu_f: float) -> float:
    """Implicit S49: mu_f * t = ln[(σ/σ0) * (1-σ0/σf) / (1-σ/σf)].

    SI S48–S49 after the γ factors cancel in the plan/S47 form.
    """
    a = math.log(sigma / sigma0)
    b = math.log((1.0 - sigma0 / sigma_f) / (1.0 - sigma / sigma_f))
    return (a + b) / mu_f


def rk4_sigma(sigma0: float, mu_f: float, sigma_f: float, t_end: float, dt: float):
    n = int(round(t_end / dt))
    t = 0.0
    s = sigma0
    ts = [0.0]
    ss = [s]
    for _ in range(n):
        k1 = dsigma_eq5(s, mu_f)
        k2 = dsigma_eq5(s + 0.5 * dt * k1, mu_f)
        k3 = dsigma_eq5(s + 0.5 * dt * k2, mu_f)
        k4 = dsigma_eq5(s + dt * k3, mu_f)
        s = s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t += dt
        ts.append(t)
        ss.append(s)
    return ts, ss


def _proteome_derivs(sigma: float, phi_rb: float, phi_cat: float, mu_f: float):
    """σ from Eq. 5; fractions from dilution (Note 4.4f / 5.2)."""
    lam = sigma * phi_rb
    return (
        dsigma_eq5(sigma, mu_f),
        lam * (chi_rb(sigma) - phi_rb),
        lam * (chi_cat(sigma) - phi_cat),
    )


def rk4_proteome(sigma0: float, phi_rb0: float, phi_cat0: float, mu_f: float,
                 t_end: float, dt: float):
    n = int(round(t_end / dt))
    t = 0.0
    s, rb, cat = sigma0, phi_rb0, phi_cat0
    ts = [0.0]
    ss = [s]
    rbs = [rb]
    cats = [cat]
    for _ in range(n):
        k1 = _proteome_derivs(s, rb, cat, mu_f)
        k2 = _proteome_derivs(
            s + 0.5 * dt * k1[0], rb + 0.5 * dt * k1[1],
            cat + 0.5 * dt * k1[2], mu_f,
        )
        k3 = _proteome_derivs(
            s + 0.5 * dt * k2[0], rb + 0.5 * dt * k2[1],
            cat + 0.5 * dt * k2[2], mu_f,
        )
        k4 = _proteome_derivs(
            s + dt * k3[0], rb + dt * k3[1], cat + dt * k3[2], mu_f,
        )
        s += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        rb += (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        cat += (dt / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2])
        t += dt
        ts.append(t)
        ss.append(s)
        rbs.append(rb)
        cats.append(cat)
    return ts, ss, rbs, cats


@dataclass(frozen=True)
class Downshift:
    lambda_i: float = LAMBDA_I
    lambda_f: float = LAMBDA_F
    phi_rb0: float = PHI_RB0
    gamma: float = GAMMA
    lambda_c: float = LAMBDA_C

    @property
    def sigma_i(self) -> float:
        return sigma_star(self.lambda_i)

    @property
    def sigma_f(self) -> float:
        return sigma_star(self.lambda_f)

    @property
    def mu_f(self) -> float:
        return mu_of(self.lambda_f)

    @property
    def sigma0(self) -> float:
        return sigma0_s57(self.sigma_i, self.lambda_i, self.lambda_f)

    @property
    def phi_rb_i(self) -> float:
        return phi_rb_star(self.lambda_i)

    @property
    def lambda0(self) -> float:
        return self.sigma0 * self.phi_rb_i
