#!/usr/bin/env python3
"""D3_2_CUTOFF_CONVERGENCE. Finite-time FCR handoff vs lambda_cut -> 0."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from fcr_standalone import Downshift, rk4_sigma

HERE = Path(__file__).resolve().parent
SPEC_PATH = HERE / "configs" / "protocol_d32.json"
OUT = HERE / "results" / "d32_summary.json"


def main() -> int:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
    assert spec["frozen_before_traces"] is True
    rows = []
    for lam_c in spec["lambda_cut_per_h"]:
        d = Downshift(lambda_f=lam_c)
        _, ss = rk4_sigma(d.sigma0, d.mu_f, d.sigma_f, spec["t_h"], spec["dt_h"])
        rows.append(
            {
                "lambda_cut": lam_c,
                "mu_f": d.mu_f,
                "sigma0": d.sigma0,
                "sigma_f_star": d.sigma_f,
                "sigma_T": ss[-1],
            }
        )

    last3 = rows[-3:]
    rels = []
    for a, b in zip(last3, last3[1:]):
        rels.append(abs(b["sigma_T"] - a["sigma_T"]) / max(abs(a["sigma_T"]), 1e-12))
    converged = all(r <= spec["delta_cv"] for r in rels)

    # Walk test: sigma_T should not scale like |log lambda_cut| across last three.
    span = abs(last3[-1]["sigma_T"] - last3[0]["sigma_T"]) / max(
        abs(last3[0]["sigma_T"]), 1e-12
    )
    converged = converged and span <= spec["delta_cv"]

    verdict = "CONVERGES" if converged else "KILL"
    summary = {
        "gate": "D3_2_CUTOFF_CONVERGENCE",
        "device": "D3_1_FCR_STANDALONE",
        "death": False,
        "dish": False,
        "rows": rows,
        "last3_rel": rels,
        "last3_span": span,
        "delta_cv": spec["delta_cv"],
        "verdict": verdict,
        "d4_still_forbidden": True,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        "device=D3_1_FCR_STANDALONE death=OFF dish=OFF "
        f"last3_rel={rels} span={span:.6f} d4_forbidden=true"
    )
    print(f"D3_2_CUTOFF_CONVERGENCE={verdict}")
    return 0 if converged else 1


if __name__ == "__main__":
    sys.exit(main())
