#!/usr/bin/env python3
"""Wringe pack 3: Mackey–Glass one-pager on frozen BenchA voxels.

Does not run BSim. Closed val-slice (lambda on rows 88..109 only).
MATLAB peers can reproduce one driven folder with report_mg_pack.m.
Does not edit any GATE_EVIDENCE.md Overall line.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent.parent
BENCHA = EXAMPLES / "BSimReservoirPlanBenchA"
sys.path.insert(0, str(BENCHA))
import check_bencha as c  # noqa: E402

SEEDS = (101, 202, 303)
HORIZONS = (1, 5, 10)
MATCH_TOL = 1e-3
MG_SHA = "e06810ea70f61bbe4baf61939f93b1e1cae9d92d9eee8d757865c4935dbac780"
OUT_DIR = BENCHA / "results"

GATE_MG = {
    "driven": (0.3306, 0.2152, 0.2347),
    "brownian": (1.0415, 1.0410, 1.0431),
    "silent": (1.0416, 1.0416, 1.0416),
    "field": (0.0229, 0.0229, 0.0229),
}
GATE_MG_MEAN = {
    "driven": 0.2602,
    "brownian": 1.0418,
    "silent": 1.0416,
    "field": 0.0229,
}
GATE_LORENZ = {
    "driven": (0.9902, 0.9965, 0.9819),
    "brownian": (0.9161, 1.1799, 1.2016),
    "field": (0.0951, 0.0951, 0.0951),
}


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


def horizon_target(y_next, h):
    n_valid = len(y_next) - h + 1
    y_h = y_next[h - 1 : h - 1 + n_valid]
    return n_valid, y_h


def evaluate_var(X_all, y_all):
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


def load_xy(folder, arm):
    voxels = Path(folder) / "voxels.csv"
    if arm == "brownian":
        matrix = c.read_window_matrix(
            voxels, c.BROWNIAN_MEAN_PREFIXES, c.BROWNIAN_LAST_PREFIXES
        )
    else:
        matrix = c.read_window_matrix(
            voxels, c.BIOLOGY_MEAN_PREFIXES, c.BIOLOGY_LAST_PREFIXES
        )
    field = None
    if arm == "driven":
        field = c.read_window_matrix(voxels, c.FIELD_MEAN_PREFIXES, ())
    return matrix, field


def analyse_arm(folder, arm, y, do_horizons=True):
    matrix, field = load_xy(folder, arm)
    h1 = c.evaluate_readout(matrix["X"], y, matrix["windows"])
    result = {
        "folder": str(folder),
        "arm": arm,
        "n_features": len(matrix["feature_names"]),
        "h1_nrmse": h1["test_nrmse"],
        "h1_r2": h1["test_r2"],
        "h1_lambda": h1["lambda"],
        "field_h1_nrmse": None,
        "horizons": [],
        "field_horizons": [],
    }
    if field is not None:
        fh1 = c.evaluate_readout(field["X"], y, field["windows"])
        result["field_h1_nrmse"] = fh1["test_nrmse"]
        result["field_h1_r2"] = fh1["test_r2"]
        result["field_h1_lambda"] = fh1["lambda"]
    if do_horizons:
        for h in HORIZONS:
            n_valid, y_h = horizon_target(y, h)
            fit = evaluate_var(matrix["X"][:n_valid], y_h)
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


def check_gate(computed, expected, seeds=SEEDS):
    mismatches = []
    for arm, exp in expected.items():
        got = computed[arm]
        for i, (g, e) in enumerate(zip(got, exp)):
            if abs(g - e) > MATCH_TOL:
                mismatches.append(
                    f"{arm} seed{seeds[i]}: got {g:.6f} expected {e:.6f}"
                )
    return mismatches


def write_figure(pack, base):
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
    colors = ["#2a6f4e", "#6b7280", "#9ca3af", "#9ca3af", "#b45309", "#b45309"]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    bars = ax.bar(labels, vals, color=colors, edgecolor="#222", linewidth=0.4)
    driven_se = pack["se"]["driven"]
    ax.errorbar(0, vals[0], yerr=driven_se, color="#222", capsize=3.5,
                elinewidth=0.9, lw=0, zorder=3)
    ax.axhline(1.0, color="#444", linestyle=":", linewidth=1, label="NRMSE = 1")
    ax.set_ylabel("Test NRMSE")
    ax.set_title("Mackey–Glass one-step NRMSE vs Wringe comparators")
    ax.set_ylim(0, max(vals) * 1.18)
    for bar, val in zip(bars, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2, val + 0.03,
            f"{val:.3f}", ha="center", va="bottom", fontsize=8,
        )
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_mg_comparators.png", dpi=160)
    plt.close()


def write_markdown(pack, path):
    b = pack["baselines"]
    m = pack["mean"]
    se = pack["se"]
    lz = pack["lorenz"]
    lines = [
        "# Wringe pack 3 — Mackey–Glass one-pager",
        "",
        "Analysis only. Frozen BenchA Mackey–Glass voxels, closed val-slice.",
        "No BSim rerun. No kinetics retune. Packs 1–2 stay closed. Nonlinear",
        "KR is not this pack. No `GATE_EVIDENCE.md` Overall line was edited:",
        "BenchA Mackey–Glass stays **PASS**, BenchA Lorenz stays **FAIL**,",
        "waveform stays **FAIL**, Stage 6 / Narma10b NARMA-10 stay **PASS**.",
        "",
        "Wringe et al. 2024 ([arXiv:2405.06561](https://arxiv.org/abs/2405.06561))",
        "§5.1 *prediction of known dynamics* on this dish. Teacher-forced",
        "**one-step** readout: window *n* predicts raw `x[n+1]`. The dish is",
        "still driven by affine `u[n] ∈ [0, 0.5]`. This is **not**",
        "Jaeger-style autonomous (closed-loop) Mackey–Glass.",
        "",
        "## Claim (do not exaggerate)",
        "",
        f"Driven biology **{m['driven']:.2f} ± {se['driven']:.3f}** vs Brownian",
        f"**{m['brownian']:.2f}** vs silent **{m['silent']:.2f}** is a real PASS",
        "against the nulls. Field-only **0.023** is much stronger. Same pattern",
        "as waveform: biology tracks the drive; the plume already *is* `x`.",
        "",
        "Why field wins: `u` is affine `x`, the target is `x[n+1]`. A linear",
        "voxel-AHL map is almost a delay line on the MG state. That is Wringe",
        "prediction of known dynamics, not a mystery reservoir.",
        "",
        f"Persistence on this file is **{b['persistence_nrmse']:.3f}**: copying",
        f"`x[n]` beats the 408-feature hybrid readout ({m['driven']:.2f}). Driven",
        "does not beat persistence or the delay line. Silent equals intercept",
        f"({b['intercept_nrmse']:.3f}); sources off, the ridge is a constant.",
        "",
        f"Driven {m['driven']:.2f} is still a decent **one-step** score versus",
        "Brownian. It is not an ESN free-run (those quote ~0.1–0.4 on",
        "closed-loop MG with trained recurrent weights). Do not paste a",
        "literature MG number onto this 50-window test.",
        "",
        "NARMA-10 remains the **only** task where biology beat field-only.",
        "",
        "## Six-comparator table (h=1 test NRMSE)",
        "",
        "Target: `mg_target.csv` column `x_next`. SHA-256 of the 12-decimal `u`",
        f"sequence: `{pack['u_sha256']}`. Discrete map dt=1, τ=17.",
        "",
        "| Comparator | 101 | 202 | 303 | mean ± s.e. |",
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
        f"| Intercept (train-mean constant) | {b['intercept_nrmse']:.4f} | "
        f"{b['intercept_nrmse']:.4f} | {b['intercept_nrmse']:.4f} | "
        f"{b['intercept_nrmse']:.4f} |",
        f"| Persistence (`x[n]` → `x[n+1]`) | {b['persistence_nrmse']:.4f} | "
        f"{b['persistence_nrmse']:.4f} | {b['persistence_nrmse']:.4f} | "
        f"{b['persistence_nrmse']:.4f} |",
        "",
        f"Intercept uses train `x_next` mean {b['train_mean']:.4f} on test "
        f"(test mean {b['test_mean']:.4f}, test pop-std {b['test_std']:.4f}).",
        "Persistence and intercept do not use voxels; they are identical",
        "across seeds. Computed on this file; not a literature MG number.",
        "",
        "Recomputed h=1 matches BenchA `GATE_EVIDENCE.md` Mackey–Glass table",
        "within 1e-3. Silent dirs are Stage 6 copies; they were not rerun.",
        "",
        "![h=1 NRMSE vs comparators](fig_wringe_mg_comparators.png)",
        "",
        "## Teacher-forced horizons (optional, not free-run)",
        "",
        "Window *n* predicts `x[n+h]`. `u` still drove the dish. Test count",
        "50 / 46 / 41. Independent lambda. Field is expected to dominate:",
        "the delay line is the MG state.",
        "",
        "| h | n_test | Driven NRMSE (101/202/303) | mean | Field-only mean |",
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
    f5 = float(np.mean(pack["horizon"]["field"][5]))
    f10 = float(np.mean(pack["horizon"]["field"][10]))
    d5 = float(np.mean(pack["horizon"]["driven"][5]))
    d10 = float(np.mean(pack["horizon"]["driven"][10]))
    lines.extend([
        "",
        f"h=5 / h=10: field {f5:.3f} / {f10:.3f} vs driven {d5:.3f} / {d10:.3f}.",
        "Teacher-forced only. Linear MC of this drive is **not** reported as",
        "comparable to NARMA MC (~1.2): MG `u` is a smooth chaotic series, not",
        "Uniform[0, 0.5]. If reconstructed, that diagnostic is “reconstruct",
        "delayed drive, not white-input memory.”",
        "",
        "## Lorenz’63 — FAIL box (not this figure)",
        "",
        "Same closed ridge. No SKIP=50. Not retuned. Not a new plot.",
        "",
        f"Driven mean **{lz['mean']['driven']:.4f}** vs Brownian",
        f"**{lz['mean']['brownian']:.4f}**: gate 1 **FAIL** (not every seed:",
        f"seed 101 driven {lz['h1']['driven'][0]:.4f} does not beat Brownian",
        f"{lz['h1']['brownian'][0]:.4f}). Field-only **{lz['mean']['field']:.3f}**:",
        "Δt=0.02 is easy for a delay line. Silent is not the gate.",
        "",
        "| Arm | 101 | 202 | 303 | mean |",
        "|---|---|---|---|---|",
        f"| Driven | {lz['h1']['driven'][0]:.4f} | {lz['h1']['driven'][1]:.4f} | "
        f"{lz['h1']['driven'][2]:.4f} | {lz['mean']['driven']:.4f} |",
        f"| Brownian | {lz['h1']['brownian'][0]:.4f} | {lz['h1']['brownian'][1]:.4f} | "
        f"{lz['h1']['brownian'][2]:.4f} | {lz['mean']['brownian']:.4f} |",
        f"| Field-only | {lz['h1']['field'][0]:.4f} | {lz['h1']['field'][1]:.4f} | "
        f"{lz['h1']['field'][2]:.4f} | {lz['mean']['field']:.4f} |",
        "",
        "Matches BenchA Lorenz GATE within 1e-3. Overall remains FAIL.",
        "",
        "## MATLAB reproduction (one driven seed)",
        "",
        "From `examples/HybridDish`:",
        "",
        "```matlab",
        "addpath('matlab')",
        "root = fullfile('..', 'BSimReservoirPlanBenchA');",
        "T = load_voxels(fullfile(root, 'results', 'mg_driven_seed101'));",
        "[X, names] = window_features(T, 'biology');",
        "target = readtable(fullfile(root, 'mg_target.csv'), ...",
        "    'Delimiter', ';', 'FileType', 'text');",
        "S = ridge_closed(X, target.x_next);",
        "B = baselines(target.x_next, 'Label', 'Mackey-Glass');",
        "[Xf, ~] = window_features(T, 'field');",
        "Sf = ridge_closed(Xf, target.x_next);",
        "```",
        "",
        "Or `R = report_mg_pack('Only', \"driven_seed101\")`.",
        "",
        "Expected seed 101: driven NRMSE 0.3306, field-only 0.0229.",
        "Biology ridge is 408 channels (`Receiver_R_*`, `Lum_Mean_*`,",
        "`Input_Driven_Death_*`). Voxel AHL is field-only, not in biology.",
        "Lambda on post-washout rows 88..109 only; then refit on 110 train.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    u, y = c.load_target(BENCHA / "mg_target.csv", "mg")
    digest = c.sha256_u(u)
    if digest != MG_SHA:
        raise SystemExit(f"MG u SHA-256 drifted: {digest}")
    print(f"u_sha256={digest}")
    baselines = intercept_persistence(y)
    print(f"intercept_nrmse={baselines['intercept_nrmse']:.6f}")
    print(f"persistence_nrmse={baselines['persistence_nrmse']:.6f}")

    computed = {arm: [] for arm in ("driven", "brownian", "silent", "field")}
    horizon = {"driven": {h: [] for h in HORIZONS}, "field": {h: [] for h in HORIZONS}}
    runs = []

    for arm in ("driven", "brownian", "silent"):
        for seed in SEEDS:
            folder = OUT_DIR / f"mg_{arm}_seed{seed}"
            print(f"=== MG {arm} seed{seed} ===")
            result = analyse_arm(folder, arm, y, do_horizons=(arm == "driven"))
            runs.append(result)
            computed[arm].append(result["h1_nrmse"])
            if result["field_h1_nrmse"] is not None:
                computed["field"].append(result["field_h1_nrmse"])
            if arm == "driven":
                for item in result["horizons"]:
                    horizon["driven"][item["h"]].append(item["test_nrmse"])
                for item in result["field_horizons"]:
                    horizon["field"][item["h"]].append(item["test_nrmse"])
            extra = ""
            if result["field_h1_nrmse"] is not None:
                extra = f" field={result['field_h1_nrmse']:.6f}"
            print(f"  h1={result['h1_nrmse']:.6f} n_features={result['n_features']}{extra}")

    mismatches = check_gate(computed, GATE_MG)
    if mismatches:
        print("MISMATCH vs BenchA MG GATE_EVIDENCE (tol 1e-3):")
        for line in mismatches:
            print(" ", line)
        raise SystemExit("MG h=1 numbers do not match GATE_EVIDENCE; stopping.")
    print("MG h=1 matches GATE_EVIDENCE within 1e-3")

    print("=== Lorenz FAIL box (h=1 only) ===")
    lorenz_u, lorenz_y = c.load_target(BENCHA / "lorenz_target.csv", "lorenz")
    del lorenz_u
    lorenz = {arm: [] for arm in ("driven", "brownian", "field")}
    for arm in ("driven", "brownian"):
        for seed in SEEDS:
            folder = OUT_DIR / f"lorenz_{arm}_seed{seed}"
            result = analyse_arm(folder, arm, lorenz_y, do_horizons=False)
            lorenz[arm].append(result["h1_nrmse"])
            if result["field_h1_nrmse"] is not None:
                lorenz["field"].append(result["field_h1_nrmse"])
            print(f"  Lorenz {arm} seed{seed} h1={result['h1_nrmse']:.6f}")
    lz_mismatch = check_gate(lorenz, GATE_LORENZ)
    if lz_mismatch:
        print("MISMATCH vs BenchA Lorenz GATE_EVIDENCE (tol 1e-3):")
        for line in lz_mismatch:
            print(" ", line)
        raise SystemExit("Lorenz h=1 numbers do not match GATE_EVIDENCE; stopping.")
    print("Lorenz h=1 matches GATE_EVIDENCE within 1e-3")

    mean = {k: c.mean_se(v)[0] for k, v in computed.items()}
    se = {k: c.mean_se(v)[1] for k, v in computed.items()}
    for key, expected in GATE_MG_MEAN.items():
        if abs(mean[key] - expected) > MATCH_TOL:
            raise SystemExit(
                f"MG mean {key} {mean[key]:.6f} != GATE {expected:.4f}"
            )

    pack = {
        "baselines": baselines,
        "h1": computed,
        "mean": mean,
        "se": se,
        "horizon": horizon,
        "horizon_mean": {
            "driven": {h: float(np.mean(horizon["driven"][h])) for h in HORIZONS},
            "field": {h: float(np.mean(horizon["field"][h])) for h in HORIZONS},
        },
        "u_sha256": digest,
        "lorenz": {
            "h1": lorenz,
            "mean": {k: float(np.mean(v)) for k, v in lorenz.items()},
        },
        "teacher_forced": True,
        "not_autonomous_mg": True,
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
    json_path = OUT_DIR / "wringe_mg_pack.json"
    json_path.write_text(json.dumps({
        **pack,
        "runs": runs,
    }, indent=2, default=json_default) + "\n", encoding="utf-8")

    write_figure(pack, OUT_DIR)
    md_path = OUT_DIR / "WRINGE_MG_PACK.md"
    write_markdown(pack, md_path)
    print(f"markdown={md_path}")
    print(f"json={json_path}")
    print(f"figure={OUT_DIR / 'fig_wringe_mg_comparators.png'}")
    print(
        f"driven={mean['driven']:.4f} field={mean['field']:.4f} "
        f"brownian={mean['brownian']:.4f} silent={mean['silent']:.4f} "
        f"intercept={baselines['intercept_nrmse']:.4f} "
        f"persistence={baselines['persistence_nrmse']:.4f}"
    )
    print(
        f"lorenz driven={pack['lorenz']['mean']['driven']:.4f} "
        f"brownian={pack['lorenz']['mean']['brownian']:.4f} "
        f"field={pack['lorenz']['mean']['field']:.4f}"
    )
    print("No GATE_EVIDENCE.md Overall line was edited.")


if __name__ == "__main__":
    main()
