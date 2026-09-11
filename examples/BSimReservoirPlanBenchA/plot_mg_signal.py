#!/usr/bin/env python3
"""Mackey-Glass report figure from frozen BenchA CSVs.

Uses check_bencha.py ridge helpers. Does not rerun BSim or retune lambda.
Stops if recomputed test NRMSE drifts from GATE_EVIDENCE by more than 1e-4.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from check_bencha import (
    BIOLOGY_LAST_PREFIXES,
    BIOLOGY_MEAN_PREFIXES,
    BROWNIAN_LAST_PREFIXES,
    BROWNIAN_MEAN_PREFIXES,
    FIELD_MEAN_PREFIXES,
    HERE,
    TRAIN,
    WASHOUT,
    load_target,
    mean_se,
    nrmse,
    read_window_matrix,
    ridge_fit,
    ridge_predict,
    select_lambda,
    standardize,
)

SEEDS = (101, 202, 303)
FROZEN = {
    "brownian": (1.0415, 1.0410, 1.0431),
    "silent": (1.0416, 1.0416, 1.0416),
    "driven": (0.3306, 0.2152, 0.2347),
    "field-only": (0.0229, 0.0229, 0.0229),
}
FROZEN_MEAN = {
    "brownian": 1.0418,
    "silent": 1.0416,
    "driven": 0.2602,
    "field-only": 0.0229,
}
FROZEN_SE = {
    "brownian": 0.0006,
    "silent": 0.0,
    "driven": 0.0357,
    "field-only": 0.0,
}

CAPTION = (
    "Mackey–Glass τ=17 on the frozen Stage 6 dish (one sample per 300 s window, "
    "u affine-mapped into [0,0.5] on the AHL AC). Ridge on 408 biology features "
    "(R, L, input-driven deaths). Test windows 150–199. Driven mean test NRMSE "
    "0.260 ± 0.036 vs Brownian 1.042 vs silent 1.042. Field-only voxel AHL is "
    "0.023: the chemical delay line still beats the hybrid readout; the PASS "
    "is vs Brownian, not vs the field. Panel A shows the test target with all "
    "three driven-seed predictions as thin navy lines of the same hue (seed 101 "
    "labelled; 202 and 303 have test NRMSE 0.215 and 0.235) and the AHL "
    "field-only readout as a dashed gray line."
)


def run_dir(arm: str, seed: int) -> Path:
    return HERE / "results" / f"mg_{arm}_seed{seed}"


def voxels_path(arm: str, seed: int) -> Path:
    source_arm = "driven" if arm == "field" else arm
    path = run_dir(source_arm, seed) / "voxels.csv"
    if path.exists():
        return path
    if arm == "silent":
        # Checker default is results/mg_silent_seed*; Stage 6 copies are the fallback.
        for candidate in (
            HERE / "results" / f"silent_seed{seed}" / "voxels.csv",
            HERE.parent / "BSimReservoirPlanStage6" / "results" / f"stage6_silent_seed{seed}" / "voxels.csv",
        ):
            if candidate.exists():
                return candidate
    raise FileNotFoundError(path)


def load_features(arm: str, seed: int):
    path = voxels_path(arm, seed)
    if arm == "brownian":
        return read_window_matrix(path, BROWNIAN_MEAN_PREFIXES, BROWNIAN_LAST_PREFIXES)
    if arm == "field":
        return read_window_matrix(path, FIELD_MEAN_PREFIXES, ())
    return read_window_matrix(path, BIOLOGY_MEAN_PREFIXES, BIOLOGY_LAST_PREFIXES)


def fit_test_predictions(X_all, y_all):
    X = X_all[WASHOUT:]
    y = y_all[WASHOUT:]
    lam, _ = select_lambda(X, y)
    X_train, y_train = X[:TRAIN], y[:TRAIN]
    X_test, y_test = X[TRAIN:], y[TRAIN:]
    (X_train_z, X_test_z), _, _ = standardize(X_train, X_test)
    weights = ridge_fit(X_train_z, y_train, lam)
    pred_test = ridge_predict(X_test_z, weights)
    return {
        "lambda": lam,
        "test_nrmse": nrmse(y_test, pred_test),
        "pred_test": np.asarray(pred_test, dtype=float),
        "y_test": np.asarray(y_test, dtype=float),
    }


def check_frozen(label: str, values: list[float]) -> None:
    frozen = FROZEN[label]
    mismatches = []
    for seed, got, expected in zip(SEEDS, values, frozen):
        if abs(got - expected) > 1e-4:
            mismatches.append(
                f"{label} seed {seed}: recomputed {got:.6f} vs frozen {expected:.4f}"
            )
    if mismatches:
        print("NRMSE mismatch > 1e-4; not drawing.", file=sys.stderr)
        for line in mismatches:
            print(line, file=sys.stderr)
        raise SystemExit(1)


def style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


def draw(n_test, y_test, driven_preds, driven_nrmse, field_pred, field_nrmse, means, ses):
    driven_color = "#1f4e79"
    field_color = "#5a5a5a"
    bar_colors = {
        "brownian": "#c8c8c8",
        "silent": "#9a9a9a",
        "driven": driven_color,
        "field-only": "#7a7a7a",
    }

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(7.0, 5.35),
        gridspec_kw={"height_ratios": [1.75, 1.0], "hspace": 0.40},
    )

    # Thin same-hue traces first; they are not legend entries.
    for seed in (202, 303):
        ax_a.plot(
            n_test, driven_preds[seed], color=driven_color, lw=0.75,
            alpha=0.45, zorder=2,
        )
    ax_a.plot(
        n_test, driven_preds[101], color=driven_color, lw=1.35, zorder=3,
        label=f"driven seed 101 (NRMSE {driven_nrmse[101]:.3f})",
    )
    ax_a.plot(
        n_test, y_test, color="black", lw=1.8, label="target $x[n+1]$", zorder=4,
    )
    ax_a.plot(
        n_test, field_pred, color=field_color, lw=1.35, ls=(0, (3.2, 1.6)),
        zorder=5, label=f"AHL field-only (NRMSE {field_nrmse:.3f})",
    )
    ax_a.set_xlim(150, 199)
    ax_a.set_xlabel("window $n$")
    ax_a.set_ylabel("$x[n+1]$")
    ax_a.set_ylim(0.32, 1.28)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    # Upper left is empty: the test window starts at the trough.
    ax_a.legend(frameon=False, loc="upper left", handlelength=2.4, borderaxespad=0.15)
    ax_a.text(-0.09, 1.06, "A", transform=ax_a.transAxes, fontsize=12, fontweight="bold", va="top")

    labels = ("brownian", "silent", "driven", "field-only")
    xs = np.arange(len(labels))
    heights = [means[name] for name in labels]
    yerr = [ses[name] for name in labels]
    bars = ax_b.bar(
        xs, heights, color=[bar_colors[name] for name in labels],
        width=0.68, edgecolor="black", linewidth=0.5, zorder=2,
    )
    for x, height, err in zip(xs, heights, yerr):
        if err >= 0.005:
            ax_b.errorbar(
                x, height, yerr=err, color="black", capsize=3.5,
                elinewidth=0.9, lw=0, zorder=3,
            )
    ax_b.axhline(1.0, color="black", ls=":", lw=0.9, zorder=1)
    ax_b.text(3.48, 1.0, "mean predictor", ha="right", va="bottom", fontsize=7.5, color="#333333")
    ax_b.set_xticks(xs, labels)
    ax_b.set_ylabel("test NRMSE")
    ax_b.set_ylim(0, 1.28)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    for bar, height, err in zip(bars, heights, yerr):
        offset = 0.035 + (err if err >= 0.005 else 0.0)
        ax_b.text(
            bar.get_x() + bar.get_width() / 2, height + offset,
            f"{height:.3f}", ha="center", va="bottom", fontsize=8,
        )
    ax_b.text(-0.09, 1.06, "B", transform=ax_b.transAxes, fontsize=12, fontweight="bold", va="top")

    return fig


def main() -> None:
    _, y = load_target(HERE / "mg_target.csv", "mg")
    n_test = np.arange(150, 200)

    nrmse_by_arm = {"brownian": [], "silent": [], "driven": [], "field-only": []}
    driven_preds = {}
    field_preds = []
    y_test = None

    for seed in SEEDS:
        driven = load_features("driven", seed)
        brownian = load_features("brownian", seed)
        silent = load_features("silent", seed)
        field = load_features("field", seed)
        d = fit_test_predictions(driven["X"], y)
        b = fit_test_predictions(brownian["X"], y)
        s = fit_test_predictions(silent["X"], y)
        f = fit_test_predictions(field["X"], y)
        nrmse_by_arm["driven"].append(d["test_nrmse"])
        nrmse_by_arm["brownian"].append(b["test_nrmse"])
        nrmse_by_arm["silent"].append(s["test_nrmse"])
        nrmse_by_arm["field-only"].append(f["test_nrmse"])
        driven_preds[seed] = d["pred_test"]
        field_preds.append(f["pred_test"])
        y_test = d["y_test"]

    for label, values in nrmse_by_arm.items():
        check_frozen(label, values)

    field_pred = field_preds[0]
    if not all(np.allclose(pred, field_pred) for pred in field_preds):
        print("field-only predictions differ across driven seeds; using seed 101", file=sys.stderr)

    means = {}
    for label, values in nrmse_by_arm.items():
        mean, _ = mean_se(values)
        means[label] = mean
    # Error bars from the frozen evidence table (driven 0.0357; others ~0).
    ses = dict(FROZEN_SE)

    out_dir = HERE / "results" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_path = out_dir / "mg_test_predictions.csv"
    with pred_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "n", "x_next", "pred_driven_101", "pred_driven_202",
            "pred_driven_303", "pred_field_only",
        ])
        for i, window in enumerate(n_test):
            writer.writerow([
                int(window),
                f"{y_test[i]:.12f}",
                f"{driven_preds[101][i]:.12f}",
                f"{driven_preds[202][i]:.12f}",
                f"{driven_preds[303][i]:.12f}",
                f"{field_pred[i]:.12f}",
            ])

    style()
    fig = draw(
        n_test, y_test, driven_preds,
        {seed: nrmse_by_arm["driven"][i] for i, seed in enumerate(SEEDS)},
        field_pred, nrmse_by_arm["field-only"][0],
        means, ses,
    )
    pdf_path = out_dir / "mg_signal_nrmse.pdf"
    png_path = out_dir / "mg_signal_nrmse.png"
    caption_path = out_dir / "mg_signal_nrmse_caption.txt"
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    caption_path.write_text(CAPTION + "\n", encoding="utf-8")
    plt.close(fig)

    print("recomputed test NRMSE")
    print(f"{'arm':<12} {'101':>8} {'202':>8} {'303':>8} {'mean':>8}   frozen")
    for label in ("brownian", "silent", "driven", "field-only"):
        vals = nrmse_by_arm[label]
        frozen = FROZEN[label]
        print(
            f"{label:<12} {vals[0]:8.4f} {vals[1]:8.4f} {vals[2]:8.4f} "
            f"{means[label]:8.4f}   {frozen[0]:.4f} / {frozen[1]:.4f} / {frozen[2]:.4f}  "
            f"mean {FROZEN_MEAN[label]:.4f}"
        )
    print(f"pdf={pdf_path}")
    print(f"png={png_path}")
    print(f"caption={caption_path}")
    print(f"predictions={pred_path}")
    print("CAPTION")
    print(CAPTION.encode("ascii", "replace").decode("ascii"))


if __name__ == "__main__":
    main()
