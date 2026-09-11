#!/usr/bin/env python3
"""Waveform2c publication figures. Analysis and plotting only.

Reuses frozen check_waveform2c.py. Does not run Java or BSim.
Does not edit Waveform GATE_EVIDENCE.md. Field >= driven is not rewritten.
RAW_U_8 = 1.0 is the answer-key ceiling, not a biology failure.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_waveform2c as W  # noqa: E402

OUT = HERE / "results" / "figures"
WASHOUT, TRAIN, TEST = W.WASHOUT, W.TRAIN, W.TEST
NUM_WINDOWS = W.NUM_WINDOWS
N_TEST_BLOCKS = W.TEST_BLOCKS
CLASS_NAMES = W.CLASS_NAMES

EXPECTED = {
    "driven_mean": 0.8352,
    "driven_101": 0.8389,
    "field": 0.8864,
    "brownian_mean": 0.5369,
    "silent": 0.5000,
    "moments": 0.5000,
    "raw_u_8": 1.0000,
    "mean_R_101": 0.1906,
    "r_Ru_101": 0.8889,
}
SANITY_TOL = 1e-3

C_UNI = "#0072B2"
C_DOWN = "#D55E00"
C_UP = "#009E73"
C_U = "#000000"
C_R = "#0072B2"
C_L = "#882255"
C_AHL = "#56B4E9"
C_FIELD = "#009E73"
C_F408 = "#D55E00"
C_BROWN = "#E69F00"
C_SILENT = "#999999"
C_MOM = "#CC79A7"
C_RAW = "#000000"
CLASS_COLOR = {0: C_UNI, 1: C_DOWN, 2: C_UP}


def style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
    })


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    pdf = OUT / f"{stem}.pdf"
    png = OUT / f"{stem}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(f"wrote {pdf.name} {png.name}")
    return pdf, png


def assert_close(name, got, expected):
    if abs(got - expected) > SANITY_TOL:
        raise SystemExit(
            f"STOP sanity {name}: {got:.6f} vs {expected:.4f} "
            f"(|Δ|={abs(got - expected):.4g} > {SANITY_TOL})"
        )


def split_bands(ax, n=NUM_WINDOWS):
    ax.axvspan(-0.5, W.WASHOUT - 0.5, color="#F0F0F0", zorder=0, lw=0)
    ax.axvspan(W.WASHOUT - 0.5, W.WASHOUT + W.TRAIN - 0.5, color="#E8F4EA", zorder=0, lw=0)
    ax.axvspan(W.WASHOUT + W.TRAIN - 0.5, n - 0.5, color="#FFF4E0", zorder=0, lw=0)


def test_block_scores(X_all, y_window):
    """Closed ridge on washout-stripped windows; return test-block pooled scores."""
    y_all = np.asarray(y_window, dtype=int)
    X = np.asarray(X_all, dtype=float)[W.WASHOUT:]
    y = y_all[W.WASHOUT:]
    lam, _ = W.select_lambda_windows(X, y)
    X_train, y_train = X[:W.TRAIN], y[:W.TRAIN]
    X_test, y_test = X[W.TRAIN:], y[W.TRAIN:]
    (X_train_z, X_test_z), _, _ = W.standardize(X_train, X_test)
    test_scores, _ = W.ovr_scores(X_train_z, y_train, X_test_z, lam)
    block = W.pool_block_scores(y_test, test_scores, W.TEST)
    return {
        "lambda": lam,
        "test_scores": test_scores,
        "y_test": y_test,
        "block": block,
        "block_macro": block["macro_ovr_auc"],
        "block_pred": np.array(block["pred"], dtype=int),
        "block_y": W.block_labels(y_all)[W.WASHOUT_BLOCKS + W.TRAIN_BLOCKS:],
    }


def figure1():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1.4, 1]})
    t = np.arange(8)
    names = CLASS_NAMES
    colors = (C_UNI, C_DOWN, C_UP)
    for seq, name, color in zip(W.TEMPLATES, names, colors):
        axes[0].plot(t, seq, "-o", color=color, ms=4, lw=1.4, label=name)
    axes[0].set_xticks(t)
    axes[0].set_xlabel("Sample in 8-window block")
    axes[0].set_ylabel("Command u")
    axes[0].set_ylim(-0.04, 0.54)
    axes[0].set_title("A  Canonical order templates (phase 0)")
    axes[0].legend(frameon=False, loc="upper right")
    shared = np.array(W.SHARED_MULTISET, dtype=float)
    axes[1].hist(shared, bins=np.linspace(-0.025, 0.525, 12), color="#CCCCCC",
                 edgecolor="black", lw=0.6)
    axes[1].axvline(0.25, color="black", ls="--", lw=0.8, label="mean 0.25")
    axes[1].set_xlabel("u value")
    axes[1].set_ylabel("Count in one block")
    axes[1].set_title("B  Shared histogram (moments matched)")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save(fig, "fig1_templates")


def figure2(u, y_block):
    fig, ax = plt.subplots(figsize=(7.2, 2.5))
    split_bands(ax)
    n = np.arange(NUM_WINDOWS)
    for block in range(W.NUM_BLOCKS):
        sl = slice(block * 8, (block + 1) * 8)
        ax.plot(n[sl], u[sl], color=CLASS_COLOR[int(y_block[block])], lw=1.0)
    ax.set_xlabel("Window")
    ax.set_ylabel("Command u")
    ax.set_ylim(-0.04, 0.54)
    ax.set_title("Injected AHL command (50 blocks, phase identically 0)")
    handles = [
        plt.Line2D([0], [0], color=CLASS_COLOR[k], lw=1.5, label=CLASS_NAMES[k])
        for k in range(3)
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right", ncol=3)
    fig.tight_layout()
    save(fig, "fig2_injected_signal")


def figure3(driven_aucs, field_auc, brownian_aucs):
    labels = [
        "RAW_U_8\n(answer key)",
        "Field AHL",
        "Driven F408",
        "Brownian",
        "MOMENTS",
        "Silent",
    ]
    means = [
        EXPECTED["raw_u_8"],
        field_auc,
        float(np.mean(driven_aucs)),
        float(np.mean(brownian_aucs)),
        EXPECTED["moments"],
        EXPECTED["silent"],
    ]
    colors = [C_RAW, C_FIELD, C_F408, C_BROWN, C_MOM, C_SILENT]
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    x = np.arange(len(labels))
    ax.bar(x, means, color=colors, width=0.72, edgecolor="none")
    ax.scatter([2, 2, 2], driven_aucs, color="black", s=14, zorder=3)
    ax.scatter([3, 3, 3], brownian_aucs, color="black", s=14, zorder=3)
    ax.axhline(0.5, color="#666666", ls=":", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Test-block macro OVR AUC")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Waveform2c test-block ranking (seeds 101/202/303 as dots)")
    fig.tight_layout()
    save(fig, "fig3_block_auc")


def figure4(block_y, driven_pred, field_pred, driven_macro, field_macro):
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 3.8), sharex=True)
    idx = np.arange(N_TEST_BLOCKS)
    for ax, pred, title, color, auc in (
        (axes[0], driven_pred, f"Driven F408  AUC {driven_macro:.3f}", C_F408, driven_macro),
        (axes[1], field_pred, f"Field AHL  AUC {field_macro:.3f}", C_FIELD, field_macro),
    ):
        ax.scatter(idx, block_y, color="black", s=28, marker="o", label="true class", zorder=3)
        ax.scatter(idx, pred, color=color, s=36, marker="x", linewidths=1.4, label="argmax", zorder=4)
        ax.set_yticks([0, 1, 2], CLASS_NAMES)
        ax.set_ylim(-0.5, 2.5)
        ax.set_title(title)
        ax.legend(frameon=False, loc="upper right", ncol=2)
    axes[1].set_xlabel("Test block (16 blocks)")
    fig.tight_layout()
    save(fig, "fig4_test_blocks")


def figure5(u, mean_r, mean_l, mean_ahl, occ):
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    split_bands(ax)
    n = np.arange(NUM_WINDOWS)
    ax.plot(n, u, color=C_U, lw=0.8, label="u")
    ax.plot(n, mean_ahl / (np.max(mean_ahl) + 1e-15) * 0.5, color=C_AHL, lw=0.9,
            label="mean AHL (scaled to 0.5)")
    ax.set_ylabel("u / scaled AHL")
    ax.set_ylim(-0.04, 0.56)
    ax2 = ax.twinx()
    ax2.plot(n, mean_r, color=C_R, lw=1.1, label="mean R")
    ax2.plot(n, mean_l, color=C_L, lw=1.1, label="mean L")
    ax2.set_ylabel("Hill R, luminescence L")
    ax.set_xlabel("Window")
    ax.set_title(
        f"Living state, seed 101  mean R={occ['mean_R']:.3f}  "
        f"r(R,u)={occ['r_meanR_u']:.2f}  {occ['occupancy']}"
    )
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, loc="upper right", ncol=2)
    fig.tight_layout()
    save(fig, "fig5_living_state")


def figure6(per_class):
    names = list(CLASS_NAMES)
    driven = [per_class["driven"][n] for n in names]
    field = [per_class["field"][n] for n in names]
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(5.4, 2.8))
    ax.bar(x - 0.18, driven, 0.36, color=C_F408, label="Driven F408 mean")
    ax.bar(x + 0.18, field, 0.36, color=C_FIELD, label="Field AHL")
    ax.axhline(0.5, color="#666666", ls=":", lw=0.8)
    ax.set_xticks(x, names)
    ax.set_ylabel("Test-block OVR AUC")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Per-class test-block AUC (driven = mean of 3 seeds)")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "fig6_per_class")


def main():
    style()
    frozen = W.verify_frozen_inputs()
    u = np.asarray(frozen["u"], dtype=float)
    y_window = np.asarray(frozen["y_window"], dtype=int)
    y_block = W.block_labels(y_window)

    driven_dir = W.available_run("driven", 101)
    if driven_dir is None:
        raise SystemExit("STOP missing driven seed 101 voxels")

    driven_aucs = []
    brownian_aucs = []
    field_auc = None
    per_class_driven = {n: [] for n in CLASS_NAMES}
    driven_101 = None
    field_101 = None

    for seed in W.SEEDS:
        ddir = W.available_run("driven", seed)
        bdir = W.available_run("brownian", seed)
        if ddir is None or bdir is None:
            raise SystemExit(f"STOP missing voxels for seed {seed}")
        driven_run = W.read_run(ddir, "driven")
        brown_run = W.read_run(bdir, "brownian")
        d_eval = test_block_scores(driven_run["matrix"]["X"], y_window)
        b_eval = test_block_scores(brown_run["matrix"]["X"], y_window)
        f_eval = test_block_scores(driven_run["field"]["X"], y_window)
        driven_aucs.append(d_eval["block_macro"])
        brownian_aucs.append(b_eval["block_macro"])
        field_auc = f_eval["block_macro"]
        for name in CLASS_NAMES:
            per_class_driven[name].append(d_eval["block"]["per_class_ovr_auc"][name])
        if seed == 101:
            driven_101 = d_eval
            field_101 = f_eval

    assert_close("driven_101", driven_aucs[0], EXPECTED["driven_101"])
    assert_close("driven_mean", float(np.mean(driven_aucs)), EXPECTED["driven_mean"])
    assert_close("field", field_auc, EXPECTED["field"])
    assert_close("brownian_mean", float(np.mean(brownian_aucs)), EXPECTED["brownian_mean"])

    voxels = W.load_voxel_arrays(driven_dir / "voxels.csv")
    occ = W.occupancy(voxels, u)
    assert_close("mean_R_101", occ["mean_R"], EXPECTED["mean_R_101"])
    assert_close("r_Ru_101", occ["r_meanR_u"], EXPECTED["r_Ru_101"])
    mean_r = np.mean(W.window_mean_matrix(voxels["receiver"], voxels["window"]), axis=1)
    mean_l = np.mean(W.window_mean_matrix(voxels["lum"], voxels["window"]), axis=1)
    mean_ahl = np.mean(W.window_mean_matrix(voxels["ahl"], voxels["window"]), axis=1)

    figure1()
    figure2(u, y_block)
    figure3(driven_aucs, field_auc, brownian_aucs)
    figure4(
        driven_101["block_y"],
        driven_101["block_pred"],
        field_101["block_pred"],
        driven_101["block_macro"],
        field_101["block_macro"],
    )
    figure5(u, mean_r, mean_l, mean_ahl, occ)
    figure6({
        "driven": {n: float(np.mean(per_class_driven[n])) for n in CLASS_NAMES},
        "field": field_101["block"]["per_class_ovr_auc"],
    })

    csv_path = OUT / "w2c_test_block_predictions.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["test_block", "true", "driven101_pred", "field_pred"])
        for i, (yt, yd, yf) in enumerate(
            zip(driven_101["block_y"], driven_101["block_pred"], field_101["block_pred"])
        ):
            writer.writerow([i, int(yt), int(yd), int(yf)])
    print(f"wrote {csv_path.name}")
    print("Waveform2c figures OK")


if __name__ == "__main__":
    main()
