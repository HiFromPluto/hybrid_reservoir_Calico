#!/usr/bin/env python3
"""Wringe pack 1 batch helper: NARMA-10 comparators on frozen Narma10b voxels.

Does not run BSim. Closed val-slice (lambda on rows 88..109 only).
MATLAB peers can reproduce one driven folder with report_narma_pack.m.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent.parent
NARMA_DIR = EXAMPLES / "BSimReservoirPlanNarma10b"
STAGE6_DIR = EXAMPLES / "BSimReservoirPlanStage6"
sys.path.insert(0, str(NARMA_DIR))
import check_narma10b as c  # noqa: E402

SEEDS = (111, 222, 333)
STAGE6_SEEDS = (101, 202, 303)
HORIZONS = (1, 5, 10)
GATE = {
    "driven": (0.9311744825832927, 0.9227617723822614, 0.9309934976088081),
    "brownian": (1.1624567155075367, 1.1622904218745869, 1.1623717712221238),
    "silent": (1.1621883174398704, 1.1621883174398704, 1.1621883174398704),
    "field": (1.0288648225259103, 1.0288648225259103, 1.0288648225259103),
}
MATCH_TOL = 1e-3
OUT_DIR = NARMA_DIR / "results"


def evaluate_var(X_all, y_all):
    """Closed-slice ridge; test length is inferred (horizons drop tail windows)."""
    n_test = len(y_all) - c.WASHOUT - c.TRAIN
    if n_test < 1:
        raise ValueError("not enough rows for train+test")
    X = X_all[c.WASHOUT :]
    y = y_all[c.WASHOUT :]
    lam, grid = c.select_lambda(X, y)
    metrics = c.fit_eval(X, y, lam)
    metrics["lambda_grid"] = grid
    metrics["n_test"] = n_test
    return metrics


def horizon_target(y_next, h):
    n_valid = len(y_next) - h + 1
    y_h = y_next[h - 1 : h - 1 + n_valid]
    return n_valid, y_h


def intercept_persistence(y_next):
    train = y_next[c.WASHOUT : c.WASHOUT + c.TRAIN]
    test = y_next[c.WASHOUT + c.TRAIN :]
    intercept = np.full_like(test, float(np.mean(train)))
    persist = y_next[c.WASHOUT + c.TRAIN - 1 : c.WASHOUT + c.TRAIN - 1 + c.TEST]
    return {
        "train_mean": float(np.mean(train)),
        "test_mean": float(np.mean(test)),
        "test_std": float(np.std(test)),
        "intercept_nrmse": c.nrmse(test, intercept),
        "persistence_nrmse": c.nrmse(test, persist),
    }


def load_xy(folder, arm, u, y):
    if arm == "brownian":
        matrix = c.read_window_matrix(
            Path(folder) / "voxels.csv", c.BROWNIAN_MEAN_PREFIXES, c.BROWNIAN_LAST_PREFIXES
        )
    else:
        matrix = c.read_window_matrix(
            Path(folder) / "voxels.csv", c.BIOLOGY_MEAN_PREFIXES, c.BIOLOGY_LAST_PREFIXES
        )
    field = None
    if arm == "driven":
        field = c.read_window_matrix(
            Path(folder) / "voxels.csv", c.FIELD_MEAN_PREFIXES, ()
        )
    return matrix, field


def analyse_arm(folder, arm, u, y, do_mc=True, do_horizons=True):
    matrix, field = load_xy(folder, arm, u, y)
    X = matrix["X"]
    h1 = c.evaluate_readout(X, y, matrix["windows"])
    result = {
        "folder": str(folder),
        "arm": arm,
        "n_features": len(matrix["feature_names"]),
        "h1_nrmse": h1["test_nrmse"],
        "h1_r2": h1["test_r2"],
        "h1_lambda": h1["lambda"],
        "field_h1_nrmse": None,
        "mc": None,
        "horizons": [],
        "field_horizons": [],
    }
    if field is not None:
        fh1 = c.evaluate_readout(field["X"], y, field["windows"])
        result["field_h1_nrmse"] = fh1["test_nrmse"]
        result["field_h1_r2"] = fh1["test_r2"]
        result["field_h1_lambda"] = fh1["lambda"]
    if do_mc:
        mc = c.memory_capacity(X, u, matrix["windows"])
        result["mc"] = mc["mc"]
        result["mc_r2"] = [d["test_r2"] for d in mc["delays"]]
        result["k0_r2"] = mc["k0_diagnostic"]["test_r2"]
    if do_horizons:
        for h in HORIZONS:
            n_valid, y_h = horizon_target(y, h)
            fit = evaluate_var(X[:n_valid], y_h)
            result["horizons"].append({
                "h": h,
                "n_test": fit["n_test"],
                "test_nrmse": fit["test_nrmse"],
                "test_r2": fit["test_r2"],
                "lambda": fit["lambda"],
            })
            if field is not None:
                ffit = evaluate_var(field["X"][:n_valid], y_h)
                result["field_horizons"].append({
                    "h": h,
                    "n_test": ffit["n_test"],
                    "test_nrmse": ffit["test_nrmse"],
                    "test_r2": ffit["test_r2"],
                    "lambda": ffit["lambda"],
                })
    return result


def mean_se(vals):
    return c.mean_se(vals)


def check_gate(computed):
    mismatches = []
    for arm, expected in GATE.items():
        got = computed[arm]
        for i, (g, e) in enumerate(zip(got, expected)):
            if abs(g - e) > MATCH_TOL:
                mismatches.append(f"{arm} seed{SEEDS[i]}: got {g:.6f} expected {e:.6f}")
    return mismatches


def fmt_row(values):
    return " | ".join(f"{v:.4f}" for v in values)


def write_figures(pack, base):
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    labels = ["Driven", "Field-only", "Brownian", "Silent", "Intercept", "Persistence"]
    vals = [
        pack["mean"]["driven"],
        pack["mean"]["field"],
        pack["mean"]["brownian"],
        pack["mean"]["silent"],
        pack["baselines"]["intercept_nrmse"],
        pack["baselines"]["persistence_nrmse"],
    ]
    colors = ["#2a6f4e", "#6b7280", "#6b7280", "#6b7280", "#b45309", "#b45309"]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    ax.bar(labels, vals, color=colors)
    ax.axhline(1.0, color="#444", linestyle=":", linewidth=1, label="NRMSE = 1")
    ax.set_ylabel("Test NRMSE")
    ax.set_title("NARMA-10 h=1 NRMSE vs Wringe comparators")
    ax.set_ylim(0, max(vals) * 1.15)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_comparators.png", dpi=160)
    plt.close()

    k = list(range(1, 21))
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    for arm, tone, label in (
        ("driven", "#2a6f4e", "Driven"),
        ("brownian", "#6b7280", "Brownian"),
        ("silent", "#b45309", "Silent"),
    ):
        curves = pack["mc_r2"][arm]
        mean_curve = np.mean(curves, axis=0)
        ax.plot(k, mean_curve, "o-", color=tone, linewidth=1.4, markersize=4, label=label)
    ax.axhline(0, color="#444", linewidth=1)
    ax.set_xlabel("Delay k (windows)")
    ax.set_ylabel("Test R² for u[n−k]")
    ax.set_title("Linear memory capacity (mean across 3 seeds)")
    ax.set_xticks(k)
    ax.set_ylim(-0.35, 1.05)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_mc.png", dpi=160)
    plt.close()

    hs = list(HORIZONS)
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    driven_h = [pack["horizon_mean"]["driven"][h] for h in hs]
    field_h = [pack["horizon_mean"]["field"][h] for h in hs]
    intercept = pack["baselines"]["intercept_nrmse"]
    ax.plot(hs, driven_h, "o-", color="#2a6f4e", linewidth=1.6, markersize=6, label="Driven biology")
    ax.plot(hs, field_h, "s--", color="#6b7280", linewidth=1.4, markersize=6, label="Field-only AHL")
    ax.axhline(intercept, color="#b45309", linestyle=":", linewidth=1.2, label="Intercept baseline")
    ax.set_xlabel("Teacher-forced horizon h (windows)")
    ax.set_ylabel("Test NRMSE")
    ax.set_title("Teacher-forced NARMA-10 horizon (not free-run)")
    ax.set_xticks(hs)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_horizon.png", dpi=160)
    plt.close()


def write_markdown(pack, path):
    b = pack["baselines"]
    m = pack["mean"]
    se = pack["se"]
    lines = [
        "# Wringe pack 1 — NARMA-10 reporting layer",
        "",
        "Analysis only. Frozen Narma10b voxels, closed val-slice. No BSim rerun.",
        "Stage 6 NARMA-10 stays **PASS**. Narma10b stays **PASS**. Neither",
        "`GATE_EVIDENCE.md` Overall line was edited.",
        "",
        "Wringe et al. 2024 ([arXiv:2405.06561](https://arxiv.org/abs/2405.06561))",
        "comparators on this dish: intercept, persistence, linear MC, teacher-forced",
        "horizons 1 / 5 / 10. Nonlinear IPC / kernel rank / generalisation rank are",
        "not in this pack.",
        "",
        "## Flagship claim (unchanged)",
        "",
        "Biology beat Brownian, silent, and field-only on **h=1** NARMA-10 NRMSE.",
        "Driven **0.928 ± 0.003** vs Brownian **1.162** vs silent **1.162** vs",
        "field-only **1.029**. Recomputed here; matches Narma10b `GATE_EVIDENCE.md`",
        "within 1e-3.",
        "",
        "0.93 is a **weak predictor** (test R² ≈ 0.15). Literature persistence on",
        "NARMA-10 is ~0.83 (Wringe table 9, 1000-point runs). Driven 0.93 does",
        "**not** beat that published persistence baseline. That is not a rewrite",
        "of PASS.",
        "",
        "## Six-comparator table (h=1 test NRMSE)",
        "",
        "| Comparator | 111 | 222 | 333 | mean ± s.e. |",
        "|---|---|---|---|---|",
    ]
    for name, key in (
        ("Driven biology", "driven"),
        ("Field-only AHL", "field"),
        ("Brownian", "brownian"),
        ("Silent", "silent"),
    ):
        vals = pack["h1"][key]
        lines.append(
            f"| {name} | {vals[0]:.4f} | {vals[1]:.4f} | {vals[2]:.4f} | "
            f"{m[key]:.4f} ± {se[key]:.4f} |"
        )
    lines.extend([
        f"| Intercept (train-mean constant) | {b['intercept_nrmse']:.4f} | {b['intercept_nrmse']:.4f} | {b['intercept_nrmse']:.4f} | {b['intercept_nrmse']:.4f} |",
        f"| Persistence (teacher-forced y[n]→y[n+1]) | {b['persistence_nrmse']:.4f} | {b['persistence_nrmse']:.4f} | {b['persistence_nrmse']:.4f} | {b['persistence_nrmse']:.4f} |",
        "",
        f"Intercept uses train y_next mean {b['train_mean']:.4f} on test (test mean "
        f"{b['test_mean']:.4f}). That is Wringe’s NRMSE=1 case, evaluated honestly",
        "on this split — train mean ≠ test mean, so the constant scores ~1.16, not 1.00.",
        "Persistence and intercept do not use voxels; they are identical across seeds.",
        "",
        "Literature persistence on NARMA-10 is ~0.83 (Wringe table 9, 1000-point runs).",
        f"Our 50-point test persistence is **{b['persistence_nrmse']:.4f}** — weaker than",
        "the literature number (short split, not a 1000-step trajectory). Driven 0.93",
        f"is below *this* persistence ({b['persistence_nrmse']:.2f}) but does **not**",
        "beat the published ~0.83. Do not copy 0.826 from Wringe table 9 onto this dish.",
        "",
        "![h=1 NRMSE vs comparators](fig_wringe_comparators.png)",
        "",
    ])
    if pack.get("stage6_h1"):
        s6 = pack["stage6_h1"]
        lines.extend([
            "## Stage 6, not this run (closed slice on Stage 6 voxels)",
            "",
            "Optional second column. Same closed ridge, seeds 101/202/303. Not a",
            "re-simulation. Published Stage 6 GATE_EVIDENCE used a leaky val slice",
            "for field-only; numbers here are recomputed closed.",
            "",
            "| Comparator | 101 | 202 | 303 | mean |",
            "|---|---|---|---|---|",
        ])
        for name, key in (
            ("Driven biology", "driven"),
            ("Field-only AHL", "field"),
            ("Brownian", "brownian"),
            ("Silent", "silent"),
        ):
            vals = s6[key]
            mean = float(np.mean(vals))
            lines.append(
                f"| {name} | {vals[0]:.4f} | {vals[1]:.4f} | {vals[2]:.4f} | {mean:.4f} |"
            )
        lines.append("")
    lines.extend([
        "## Linear memory capacity",
        "",
        "Same ridge, reconstruct `u[n−k]` for k=1..20. MC = Σ max(0, test R²_k).",
        "k=0 (current pulse) is printed, not added. Independent lambda per delay.",
        "",
        "| Arm | MC k=1..20 | k=0 test R² (not added) |",
        "|---|---|---|",
    ])
    for arm in ("driven", "brownian", "silent"):
        mc = pack["mc"][arm]
        k0 = pack["k0"][arm]
        mc_m, mc_se = mean_se(mc)
        lines.append(
            f"| {arm} | {mc[0]:.3f} / {mc[1]:.3f} / {mc[2]:.3f} "
            f"(mean {mc_m:.3f} ± {mc_se:.3f}) | "
            f"{k0[0]:.3f} / {k0[1]:.3f} / {k0[2]:.3f} |"
        )
    lines.extend([
        "",
        "Driven MC ≈ 1.2 with useful lags ~5. Brownian and silent MC ≈ 0.",
        "",
        "![Linear MC R² vs delay](fig_wringe_mc.png)",
        "",
        "## Teacher-forced horizons",
        "",
        "Window n predicts `y[n+h] = y_next[n+h-1]`. Not closed-loop free-run;",
        "u still drove the dish. Test count 50 / 46 / 41 (windows 150..(199−h+1)).",
        "Independent lambda.",
        "Diagnostics. If h=5 and h=10 sit near intercept, the dish is not carrying",
        "a long NARMA state — only a short fading memory, which matches the MC curve.",
        "",
        "| h | n_test | Driven NRMSE (111/222/333) | mean | Field-only mean |",
        "|---|---|---|---|---|",
    ])
    for h in HORIZONS:
        driven_vals = pack["horizon"]["driven"][h]
        field_vals = pack["horizon"]["field"][h]
        dmean = float(np.mean(driven_vals))
        fmean = float(np.mean(field_vals))
        n_test = {1: 50, 5: 46, 10: 41}[h]
        lines.append(
            f"| {h} | {n_test} | {driven_vals[0]:.4f} / {driven_vals[1]:.4f} / "
            f"{driven_vals[2]:.4f} | {dmean:.4f} | {fmean:.4f} |"
        )
    intercept = b["intercept_nrmse"]
    h5 = float(np.mean(pack["horizon"]["driven"][5]))
    h10 = float(np.mean(pack["horizon"]["driven"][10]))
    f5 = float(np.mean(pack["horizon"]["field"][5]))
    f10 = float(np.mean(pack["horizon"]["field"][10]))
    horizon_note = (
        f"h=5 driven {h5:.3f} / field {f5:.3f} sit at or above intercept "
        f"({intercept:.3f}): no horizon skill. h=10 driven {h10:.3f} / field "
        f"{f10:.3f} drop together; biology is within 0.01 of the delay line. "
        "NRMSE is scaled by each horizon’s own test-target std and n_test "
        "(46 / 41), so h=10 is not a better NARMA engine than h=1. Short fading "
        "memory, consistent with MC ≈ 5 useful lags."
    )
    lines.extend([
        "",
        horizon_note,
        "",
        "![Teacher-forced horizon NRMSE](fig_wringe_horizon.png)",
        "",
        "## MATLAB reproduction (one driven seed)",
        "",
        "From `examples/HybridDish`:",
        "",
        "```matlab",
        "addpath('matlab')",
        "root = fullfile('..', 'BSimReservoirPlanNarma10b');",
        "T = load_voxels(fullfile(root, 'results', 'narma10b_driven_seed222'));",
        "[X, names] = window_features(T, 'biology');",
        "target = readtable(fullfile(root, 'narma10_target.csv'), ...",
        "    'Delimiter', ';', 'FileType', 'text');",
        "S = ridge_closed(X, target.y_next);",
        "B = baselines(target.y_next);",
        "MC = memory_capacity(X, target.u);",
        "H5 = horizon_narma(X, target.y_next, 5);",
        "```",
        "",
        "Or `R = report_narma_pack('Only', \"driven_seed222\")`.",
        "",
        "Ridge: washout 40 / train 110 / test 50; lambda on post-washout rows",
        "88..109 only; bias unregularised; train pop-std; ties take larger lambda.",
        "Biology ridge is 408 channels. Voxel AHL is field-only, not in biology.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    u, y, digest = c.load_target(NARMA_DIR / "narma10_target.csv", NARMA_DIR / "input_ahl_narma200.txt")
    print(f"u_sha256={digest}")
    baselines = intercept_persistence(y)
    print(f"intercept_nrmse={baselines['intercept_nrmse']:.6f}")
    print(f"persistence_nrmse={baselines['persistence_nrmse']:.6f}")

    computed = {arm: [] for arm in ("driven", "brownian", "silent", "field")}
    mc = {arm: [] for arm in ("driven", "brownian", "silent")}
    k0 = {arm: [] for arm in ("driven", "brownian", "silent")}
    mc_r2 = {arm: [] for arm in ("driven", "brownian", "silent")}
    horizon = {"driven": {h: [] for h in HORIZONS}, "field": {h: [] for h in HORIZONS}}
    runs = []

    for arm in ("driven", "brownian", "silent"):
        for seed in SEEDS:
            folder = OUT_DIR / f"narma10b_{arm}_seed{seed}"
            print(f"=== {arm} seed{seed} ===")
            result = analyse_arm(folder, arm, u, y)
            runs.append(result)
            computed[arm].append(result["h1_nrmse"])
            if result["field_h1_nrmse"] is not None:
                computed["field"].append(result["field_h1_nrmse"])
            mc[arm].append(result["mc"])
            k0[arm].append(result["k0_r2"])
            mc_r2[arm].append(result["mc_r2"])
            if arm == "driven":
                for item in result["horizons"]:
                    horizon["driven"][item["h"]].append(item["test_nrmse"])
                for item in result["field_horizons"]:
                    horizon["field"][item["h"]].append(item["test_nrmse"])
            print(f"  h1={result['h1_nrmse']:.6f} mc={result['mc']:.3f}")

    mismatches = check_gate(computed)
    if mismatches:
        print("MISMATCH vs Narma10b GATE_EVIDENCE (tol 1e-3):")
        for line in mismatches:
            print(" ", line)
        raise SystemExit("h=1 numbers do not match GATE_EVIDENCE; stopping.")
    print("h=1 matches Narma10b GATE_EVIDENCE within 1e-3")

    stage6_h1 = None
    if all((STAGE6_DIR / "results" / f"stage6_{arm}_seed{seed}" / "voxels.csv").exists()
           for arm in ("driven", "brownian", "silent") for seed in STAGE6_SEEDS):
        stage6_h1 = {arm: [] for arm in ("driven", "brownian", "silent", "field")}
        for arm in ("driven", "brownian", "silent"):
            for seed in STAGE6_SEEDS:
                folder = STAGE6_DIR / "results" / f"stage6_{arm}_seed{seed}"
                print(f"=== Stage 6 {arm} seed{seed} (h=1 only) ===")
                result = analyse_arm(folder, arm, u, y, do_mc=False, do_horizons=False)
                stage6_h1[arm].append(result["h1_nrmse"])
                if result["field_h1_nrmse"] is not None:
                    stage6_h1["field"].append(result["field_h1_nrmse"])

    mean = {k: mean_se(v)[0] for k, v in computed.items()}
    se = {k: mean_se(v)[1] for k, v in computed.items()}
    horizon_mean = {
        "driven": {h: float(np.mean(horizon["driven"][h])) for h in HORIZONS},
        "field": {h: float(np.mean(horizon["field"][h])) for h in HORIZONS},
    }
    pack = {
        "baselines": baselines,
        "h1": computed,
        "mean": mean,
        "se": se,
        "mc": mc,
        "k0": k0,
        "mc_r2": mc_r2,
        "horizon": horizon,
        "horizon_mean": horizon_mean,
        "stage6_h1": stage6_h1,
        "u_sha256": digest,
    }

    def json_default(obj):
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(type(obj))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "wringe_narma_pack.json"
    json_path.write_text(json.dumps({
        **{k: v for k, v in pack.items() if k != "mc_r2"},
        "mc_r2": {arm: np.asarray(curves).tolist() for arm, curves in mc_r2.items()},
        "runs": [{k: v for k, v in run.items() if k not in ("mc_r2",)} | {
            "mc_r2": run.get("mc_r2"),
        } for run in runs],
    }, indent=2, default=json_default) + "\n", encoding="utf-8")

    write_figures(pack, OUT_DIR)
    md_path = OUT_DIR / "WRINGE_NARMA_PACK.md"
    write_markdown(pack, md_path)
    print(f"markdown={md_path}")
    print(f"json={json_path}")
    print(
        f"driven={mean['driven']:.4f} field={mean['field']:.4f} "
        f"brownian={mean['brownian']:.4f} silent={mean['silent']:.4f} "
        f"intercept={baselines['intercept_nrmse']:.4f} "
        f"persistence={baselines['persistence_nrmse']:.4f}"
    )


if __name__ == "__main__":
    main()
