#!/usr/bin/env python3
"""D3_3_PROTEOME_TABLE. Fig. 4b,c vs Supplementary Table 2. Not the D0 dish."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from fcr_standalone import (
    GAMMA,
    LAMBDA_C,
    LAMBDA_F,
    LAMBDA_I,
    PHI_RB0,
    Downshift,
    phi_rb_star,
    rk4_proteome,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_JSON = HERE / "configs" / "protocol_d33.json"
OUT = HERE / "results" / "d33_summary.json"

ISOLATION = [
    "examples/LaneD_LifeCycle/PROTOCOL_D0.md",
    "examples/LaneD_LifeCycle/results/d0_grow.csv",
    "src/bsim/laned/LaneDClosedBath.java",
    "src/bsim/laned/LaneDClosedBathCell.java",
    "examples/BacteriumFromScratch/ChassisParameters.java",
    "examples/BacteriumFromScratch/EcoliRodCell.java",
    "examples/BacteriumFromScratch/NutrientField.java",
]

SECTOR_KEY = {"C up": "up", "C down": "down", "C flat": "flat"}


def git_isolation() -> bool:
    r = subprocess.run(
        ["git", "diff", "HEAD", "--", *ISOLATION],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    return r.returncode == 0 and not r.stdout.strip()


def parse_num(s: str) -> float:
    s = (s or "").strip()
    if not s:
        return 0.0
    return float(s.replace(",", "."))


def load_table(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter=";"))
    times = [parse_num(x) for x in rows[1][3:11]]
    sums = {k: [0.0] * 8 for k in ("up", "down", "flat")}
    n = 0
    for r in rows[2:]:
        if len(r) < 11:
            continue
        key = SECTOR_KEY[r[2]]
        n += 1
        for i, v in enumerate(r[3:11]):
            sums[key][i] += parse_num(v)
    return times, sums, n


def nearest(ts, ys, t):
    i = min(range(len(ts)), key=lambda j: abs(ts[j] - t))
    return ys[i]


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

    table_path = HERE / spec["table_file"]
    times, sums, n_genes = load_table(table_path)
    t_identity = (
        n_genes == spec["n_genes_expected"]
        and times == [-0.25, 0.0, 0.333, 1.0, 2.0, 3.0, 3.5, 4.0]
        and spec["table_identity"].startswith("Erickson 2017 Supplementary Table 2")
        and spec["biselli_phi_x"] == "phi_Rb"
        and abs(spec["biselli_phi_x_prime_per_h"] - 1.0 / GAMMA) < 1e-12
    )

    col_totals = [sums["up"][i] + sums["down"][i] + sums["flat"][i] for i in range(8)]
    t_sum = all(abs(tot - 100.0) <= spec["delta_sum_pp"] for tot in col_totals)

    pre_idx = [times.index(t) for t in spec["pre_times_h"]]
    c_pre = {
        k: sum(sums[k][i] for i in pre_idx) / len(pre_idx)
        for k in ("up", "down", "flat")
    }
    t_flat = all(
        abs(sums["flat"][i] - c_pre["flat"]) <= spec["delta_flat_pp"]
        for i in range(8)
    )

    post_i = times.index(spec["slope_post_time_h"])
    phi_rb_i = phi_rb_star(LAMBDA_I)
    phi_rb_f = phi_rb_star(LAMBDA_F)
    phi_cat_i = 1.0 - LAMBDA_I / LAMBDA_C
    phi_cat_f = 1.0 - LAMBDA_F / LAMBDA_C

    c_dn_prime = (sums["down"][post_i] - c_pre["down"]) / (phi_rb_f - phi_rb_i)
    c_dn0 = c_pre["down"] - c_dn_prime * phi_rb_i
    c_up_prime = (sums["up"][post_i] - c_pre["up"]) / (phi_cat_f - phi_cat_i)
    c_up0 = c_pre["up"] - c_up_prime * phi_cat_i

    d = Downshift()
    ts, ss, rbs, cats = rk4_proteome(
        d.sigma0, d.phi_rb_i, phi_cat_i, d.mu_f,
        spec["t_end_h"], spec["dt_h"],
    )

    rows = []
    max_dn = 0.0
    max_up = 0.0
    for t in spec["gate_times_h"]:
        i = times.index(t)
        phi_rb = nearest(ts, rbs, t)
        phi_cat = nearest(ts, cats, t)
        pred_dn = c_dn0 + c_dn_prime * phi_rb
        pred_up = c_up0 + c_up_prime * phi_cat
        err_dn = abs(pred_dn - sums["down"][i])
        err_up = abs(pred_up - sums["up"][i])
        max_dn = max(max_dn, err_dn)
        max_up = max(max_up, err_up)
        rows.append(
            {
                "t_h": t,
                "C_down": sums["down"][i],
                "C_up": sums["up"][i],
                "C_flat": sums["flat"][i],
                "phi_rb": phi_rb,
                "phi_cat": phi_cat,
                "pred_down": pred_dn,
                "pred_up": pred_up,
                "err_down_pp": err_dn,
                "err_up_pp": err_up,
            }
        )

    t_sector = max_dn <= spec["delta_pp"] and max_up <= spec["delta_pp"]
    isolation = git_isolation()
    honesty = (
        spec["device"] == "D3_1_FCR_STANDALONE"
        and spec["death"] is False
        and spec["dish"] is False
        and spec["fit"] is False
    )
    t_bisselli = spec["biselli_phi_x"] == "phi_Rb"

    holds = all(
        [
            honesty, no_fit, t_identity, t_sum, t_flat, t_sector,
            t_bisselli, isolation,
        ]
    )

    summary = {
        "gate": "D3_3_PROTEOME_TABLE",
        "device": "D3_1_FCR_STANDALONE",
        "death": False,
        "dish": False,
        "fit": False,
        "n_genes": n_genes,
        "col_totals": col_totals,
        "C_pre": c_pre,
        "C_down0": c_dn0,
        "C_down_prime": c_dn_prime,
        "C_up0": c_up0,
        "C_up_prime": c_up_prime,
        "phi_rb_i": phi_rb_i,
        "phi_rb_f": phi_rb_f,
        "phi_cat_i": phi_cat_i,
        "phi_cat_f": phi_cat_f,
        "max_err_down_pp": max_dn,
        "max_err_up_pp": max_up,
        "delta_pp": spec["delta_pp"],
        "rows": rows,
        "T_IDENTITY": t_identity,
        "T_SUM": t_sum,
        "T_FLAT": t_flat,
        "T_NOFIT": no_fit,
        "T_SECTOR": t_sector,
        "T_BISSELLI": t_bisselli,
        "honesty": honesty,
        "isolation": isolation,
        "verdict": "PORT_HOLDS" if holds else "FAIL",
        "lifts_t_overlay": bool(holds),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(
        "device=D3_1_FCR_STANDALONE death=OFF dish=OFF fit=OFF "
        f"n_genes={n_genes} max_err_down={max_dn:.3f} "
        f"max_err_up={max_up:.3f} pp"
    )
    print(f"D3_3_PROTEOME_TABLE={summary['verdict']}")
    if not holds:
        print(json.dumps({k: summary[k] for k in (
            "T_IDENTITY", "T_SUM", "T_FLAT", "T_NOFIT", "T_SECTOR",
            "T_BISSELLI", "honesty", "isolation",
            "max_err_down_pp", "max_err_up_pp",
        )}, indent=2))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
