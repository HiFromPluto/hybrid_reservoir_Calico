#!/usr/bin/env python3
"""D3_1_FCR_FIGURE. Standalone Eq. 5 / S47. Not the D0 dish."""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

from fcr_standalone import (
    GAMMA,
    LAMBDA_C,
    LAMBDA_F,
    LAMBDA_I,
    PHI_RB0,
    Downshift,
    dsigma_eq5,
    dsigma_s47,
    rk4_sigma,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_JSON = HERE / "configs" / "protocol_d31.json"
OUT = HERE / "results" / "d31_summary.json"

ISOLATION = [
    "examples/LaneD_LifeCycle/PROTOCOL_D0.md",
    "examples/LaneD_LifeCycle/results/d0_grow.csv",
    "src/bsim/laned/LaneDClosedBath.java",
    "src/bsim/laned/LaneDClosedBathCell.java",
    "examples/BacteriumFromScratch/ChassisParameters.java",
    "examples/BacteriumFromScratch/EcoliRodCell.java",
    "examples/BacteriumFromScratch/NutrientField.java",
]


def git_isolation() -> bool:
    r = subprocess.run(
        ["git", "diff", "HEAD", "--", *ISOLATION],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    return r.returncode == 0 and not r.stdout.strip()


def s49_time(sigma: float, sigma0: float, sigma_f: float, mu_f: float) -> float:
    """Integral of S47 (SI S48–S49; OCR of S49 is unreadable, so this is the antiderivative)."""
    return (
        math.log(sigma / sigma0)
        + (1.0 - sigma_f / GAMMA)
        * math.log((1.0 - sigma0 / sigma_f) / (1.0 - sigma / sigma_f))
    ) / mu_f


def main() -> int:
    spec = json.loads(PROTOCOL_JSON.read_text(encoding="utf-8"))
    assert spec["frozen_before_traces"] is True
    assert spec["death"] is False
    assert spec["dish"] is False
    assert spec["fit"] is False

    no_fit = (
        spec["phi_rb0"] == PHI_RB0
        and spec["gamma_per_h"] == GAMMA
        and spec["lambda_c_per_h"] == LAMBDA_C
        and spec["lambda_f_per_h"] == LAMBDA_F
        and spec["lambda_i_per_h"] == LAMBDA_I
    )

    d = Downshift()
    ts, ss = rk4_sigma(d.sigma0, d.mu_f, d.sigma_f, spec["t_end_h"], spec["dt_h"])
    s_end = ss[-1]
    lam0 = d.lambda0
    lam_end = s_end * (PHI_RB0 + LAMBDA_F / GAMMA)

    max_ode = 0.0
    max_s49 = 0.0
    for t, s in zip(ts[1:], ss[1:]):
        e5 = dsigma_eq5(s, d.mu_f)
        e47 = dsigma_s47(s, d.mu_f, d.sigma_f)
        max_ode = max(max_ode, abs(e5 - e47) / max(abs(e47), 1e-12))
        t_imp = s49_time(s, d.sigma0, d.sigma_f, d.mu_f)
        max_s49 = max(max_s49, abs(t_imp - t) / max(t, 1e-12))

    t_lf = abs(lam_end - LAMBDA_F) / LAMBDA_F <= spec["delta_lambda"]
    t_sigma = abs(s_end - d.sigma_f) / d.sigma_f <= spec["delta_sigma"]
    t_drop = lam0 < LAMBDA_I and lam_end >= lam0
    t_ode = max_ode <= spec["delta_ode"]
    t_s49 = max_s49 <= spec["delta_s49"]
    isolation = git_isolation()

    honesty = (
        spec["device"] == "D3_1_FCR_STANDALONE"
        and spec["death"] is False
        and spec["dish"] is False
        and spec["fit"] is False
    )

    holds = all(
        [honesty, no_fit, t_lf, t_sigma, t_drop, t_ode, t_s49, isolation]
    )

    summary = {
        "gate": "D3_1_FCR_FIGURE",
        "device": "D3_1_FCR_STANDALONE",
        "death": False,
        "dish": False,
        "fit": False,
        "phi_rb0": PHI_RB0,
        "gamma_per_h": GAMMA,
        "lambda_c_per_h": LAMBDA_C,
        "lambda_i_per_h": LAMBDA_I,
        "lambda_f_per_h": LAMBDA_F,
        "mu_f_per_h": d.mu_f,
        "sigma_i": d.sigma_i,
        "sigma0": d.sigma0,
        "sigma_f": d.sigma_f,
        "lambda0": lam0,
        "lambda_end": lam_end,
        "sigma_end": s_end,
        "max_eq5_s47_rel": max_ode,
        "max_s49_rel": max_s49,
        "T_LF": t_lf,
        "T_SIGMA": t_sigma,
        "T_DROP": t_drop,
        "T_EQ5_S47": t_ode,
        "T_S49": t_s49,
        "T_NOFIT": no_fit,
        "T_OVERLAY": "SCOPE_NOTE",
        "honesty": honesty,
        "isolation": isolation,
        "verdict": "PORT_HOLDS" if holds else "FAIL",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(
        "device=D3_1_FCR_STANDALONE death=OFF dish=OFF fit=OFF "
        f"lambda_i={LAMBDA_I} lambda_f={LAMBDA_F} "
        f"lambda0={lam0:.6f} lambda_end={lam_end:.6f} "
        f"T_OVERLAY=SCOPE_NOTE"
    )
    print(f"D3_1_FCR_FIGURE={summary['verdict']}")
    if not holds:
        print(json.dumps({k: summary[k] for k in (
            "T_LF", "T_SIGMA", "T_DROP", "T_EQ5_S47", "T_S49",
            "T_NOFIT", "honesty", "isolation", "max_eq5_s47_rel",
            "max_s49_rel",
        )}, indent=2))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
