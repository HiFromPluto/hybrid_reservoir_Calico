#!/usr/bin/env python3
"""Lorenz L3 publication figures. Analysis and plotting only.

Reuses the frozen closed ridge in check_lorenz_l3.py / screen_lorenz_l3.py
/ check_bencha.py. Does not run Java or BSim. Does not change k, target,
kinetics, or GATE_EVIDENCE.md. Unmasked RL is descriptive, not biology.
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
sys.path.insert(0, str(HERE.parent / "BSimReservoirPlanBenchA"))

import check_bencha as BC  # noqa: E402
import check_lorenz_l3 as L3  # noqa: E402
from screen_lorenz_l3 import delay_matrix  # noqa: E402

OUT = HERE / "results" / "figures"
WASHOUT, TRAIN, TEST = L3.WASHOUT, L3.TRAIN, L3.TEST
NUM_WINDOWS = L3.NUM_WINDOWS
N_TEST = np.arange(WASHOUT + TRAIN, NUM_WINDOWS)

EXPECTED = {
    "ar": 0.8233,
    "f408": 0.8256,
    "field": 0.7973,
    "brownian": 1.0076,
    "silent": 1.0076,
    "masked": 0.7967,
    "unmasked": 0.3477,
}
SANITY_TOL = 1e-3

# Wong colour-blind palette.
C_TRUE = "#000000"
C_PERSIST = "#E69F00"
C_AR = "#0072B2"
C_FIELD = "#009E73"
C_F408 = "#D55E00"
C_UNMASKED = "#CC79A7"
C_AHL = "#56B4E9"
C_L = "#882255"
C_INTERCEPT = "#999999"


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


def fit_closed_predictions(X_all, y_all):
    """Train-fit closed ridge; test predictions. Same lambda rule as checkers."""
    X_all = np.asarray(X_all, dtype=float)
    y_all = np.asarray(y_all, dtype=float)
    X = X_all[WASHOUT:]
    y = y_all[WASHOUT:]
    lam, _ = BC.select_lambda(X, y)
    X_train, y_train = X[:TRAIN], y[:TRAIN]
    X_test, y_test = X[TRAIN:], y[TRAIN:]
    (X_train_z, X_test_z), _, _ = BC.standardize(X_train, X_test)
    weights = BC.ridge_fit(X_train_z, y_train, lam)
    pred_test = BC.ridge_predict(X_test_z, weights)
    return {
        "lambda": lam,
        "test_nrmse": BC.nrmse(y_test, pred_test),
        "pred_test": np.asarray(pred_test, dtype=float),
        "y_test": np.asarray(y_test, dtype=float),
    }


def assert_close(name, got, expected):
    if abs(got - expected) > SANITY_TOL:
        raise SystemExit(
            f"STOP sanity {name}: {got:.6f} vs {expected:.4f} "
            f"(|Δ|={abs(got - expected):.4g} > {SANITY_TOL})"
        )


def split_bands(ax, full=True):
    if full:
        ax.axvspan(-0.5, WASHOUT - 0.5, color="#F0F0F0", zorder=0, lw=0)
        ax.axvspan(WASHOUT - 0.5, WASHOUT + TRAIN - 0.5, color="#E8F4EA", zorder=0, lw=0)
        ax.axvspan(WASHOUT + TRAIN - 0.5, NUM_WINDOWS - 0.5, color="#FFF4E0", zorder=0, lw=0)
    else:
        ax.axvspan(WASHOUT + TRAIN - 0.5, NUM_WINDOWS - 0.5, color="#FFF4E0", zorder=0, lw=0)


def save(fig, stem):
    pdf = OUT / f"{stem}.pdf"
    png = OUT / f"{stem}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(f"wrote {pdf.name} {png.name}")
    return pdf, png


def load_clock_csv():
    path = HERE / "results" / "lorenz_l3_clock_screen.csv"
    if not path.exists():
        raise SystemExit("STOP missing clock table; run screen_lorenz_l3.py (no BSim)")
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return rows


def clock_xy(rows, task, model):
    ks, dts, scores = [], [], []
    for row in rows:
        if row["task"] != task or row["model"] != model:
            continue
        ks.append(int(row["k"]))
        dts.append(float(row["dt_sample"]))
        scores.append(float(row["test_nrmse"]))
    order = np.argsort(ks)
    return np.array(ks)[order], np.array(dts)[order], np.array(scores)[order]


def figure1(rows):
    auto = {
        "intercept": clock_xy(rows, "AUTO_X", "TRAIN_INTERCEPT"),
        "persist": clock_xy(rows, "AUTO_X", "PERSISTENCE_X_N"),
        "ar": clock_xy(rows, "AUTO_X", "LINEAR_AR_VALIDATION_SELECTED"),
    }
    cross = {
        "intercept": clock_xy(rows, "CROSS_Y", "TRAIN_INTERCEPT"),
        "persist": clock_xy(rows, "CROSS_Y", "PERSISTENCE_Y_N"),
        "delay": clock_xy(rows, "CROSS_Y", "LINEAR_X_DELAY_VALIDATION_SELECTED"),
    }
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 5.8), sharex=True)
    fig.suptitle("Intermediate clock, not a living win", fontsize=11, y=0.995)

    def draw(ax, series, delay_key, delay_label, panel, note, annotate_k):
        _, dt, y = series["intercept"]
        ax.plot(dt, y, color=C_INTERCEPT, marker="s", ms=4.5, lw=1.1,
                label="train-intercept")
        _, dt, y = series["persist"]
        ax.plot(dt, y, color=C_PERSIST, marker="o", ms=5, lw=1.2,
                label="persistence")
        _, dt, y = series[delay_key]
        ax.plot(dt, y, color=C_AR, marker="D", ms=5, lw=1.3,
                label=delay_label)
        ax.axhline(0.30, color="#0072B2", ls="--", lw=0.9, alpha=0.85)
        ax.axhline(0.95, color="#D55E00", ls="--", lw=0.9, alpha=0.85)
        ax.text(0.023, 0.12, "TRIVIAL 0.30", fontsize=7, color="#0072B2")
        ax.text(0.023, 1.00, "DESTROYED 0.95", fontsize=7, color="#D55E00")
        for xx in (0.02, 0.20, 1.00):
            ax.axvline(xx, color="#444444", ls=":", lw=0.8, zorder=1)
        if annotate_k:
            ax.text(0.02, 1.72, "k=1", ha="center", fontsize=7, color="#333333")
            ax.text(0.20, 1.72, "k=10 selected", ha="center", fontsize=7, color="#333333")
            ax.text(1.00, 1.72, "k=50", ha="center", fontsize=7, color="#333333")
        ax.set_xscale("log")
        ax.set_ylim(0.0, 1.82)
        ax.set_ylabel("test NRMSE")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.12, 1.08, panel, transform=ax.transAxes, fontsize=12, fontweight="bold")
        ax.set_title(note, loc="left", fontsize=9)
        ax.legend(frameon=False, loc="lower right", handlelength=1.8)

    draw(axes[0], auto, "ar", "selected AR of $x$", "A",
         "AUTO_X  $x[n]\\to x[n+1]$", True)
    draw(axes[1], cross, "delay", "selected $x$-delay", "B",
         "CROSS_Y  $x[n]\\to y[n+1]$  (Module-1 only; not BSim'd)", False)
    axes[1].set_xlabel("Lorenz $\\Delta t$ per 300 s dish window")
    axes[1].set_xticks([0.02, 0.10, 0.20, 0.40, 0.80, 1.00])
    axes[1].set_xticklabels(["0.02", "0.10", "0.20", "0.40", "0.80", "1.00"])
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save(fig, "fig1_clock_map")


def figure2(n, x, x_next, u):
    fig, (ax_x, ax_u, ax_s) = plt.subplots(
        3, 1, figsize=(6.4, 6.4),
        gridspec_kw={"height_ratios": [1.05, 0.72, 1.15], "hspace": 0.42},
    )
    split_bands(ax_x)
    split_bands(ax_u)
    ax_x.plot(n, x, color=C_TRUE, lw=1.15)
    ax_u.plot(n, u, color=C_AHL, lw=1.1)
    ax_x.set_xlim(-1, 200)
    ax_u.set_xlim(-1, 200)
    ax_x.set_ylabel("$x[n]$")
    ax_u.set_ylabel("affine $u$")
    ax_u.set_ylim(-0.02, 0.55)
    ax_x.set_xlabel("")
    ax_u.set_xlabel("window $n$")
    ax_x.spines["top"].set_visible(False)
    ax_x.spines["right"].set_visible(False)
    ax_u.spines["top"].set_visible(False)
    ax_u.spines["right"].set_visible(False)
    for ax in (ax_x, ax_u):
        ax.text(20, 0.92, "washout", fontsize=7, color="#666666", ha="center",
                transform=ax.get_xaxis_transform())
        ax.text(95, 0.92, "train", fontsize=7, color="#3d6b45", ha="center",
                transform=ax.get_xaxis_transform())
        ax.text(175, 0.92, "test", fontsize=7, color="#8a5a00", ha="center",
                transform=ax.get_xaxis_transform())
    ax_x.text(-0.12, 1.10, "A", transform=ax_x.transAxes, fontsize=12, fontweight="bold")
    ax_x.set_title("Injected scalar drive (not a reservoir score)", loc="left", fontsize=9)

    train = slice(WASHOUT, WASHOUT + TRAIN)
    test = slice(WASHOUT + TRAIN, NUM_WINDOWS)
    ax_s.scatter(x[train], x_next[train], s=14, c="#0072B2", alpha=0.7,
                 edgecolors="none", label="train")
    ax_s.scatter(x[test], x_next[test], s=18, c="#D55E00", marker="s", alpha=0.85,
                 edgecolors="none", label="test")
    lo = min(float(np.min(x)), float(np.min(x_next)))
    hi = max(float(np.max(x)), float(np.max(x_next)))
    pad = 0.06 * (hi - lo)
    ax_s.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#888888",
              ls=":", lw=0.9, zorder=0)
    ax_s.set_xlim(lo - pad, hi + pad)
    ax_s.set_ylim(lo - pad, hi + pad)
    ax_s.set_aspect("equal", adjustable="box")
    ax_s.set_xlabel("$x[n]$")
    ax_s.set_ylabel("$x[n+1]$")
    ax_s.spines["top"].set_visible(False)
    ax_s.spines["right"].set_visible(False)
    ax_s.legend(frameon=False, loc="lower right")
    ax_s.text(-0.12, 1.08, "B", transform=ax_s.transAxes, fontsize=12, fontweight="bold")
    ax_s.set_title("AUTO_X map at $k=10$", loc="left", fontsize=9)
    return save(fig, "fig2_injected_signal")


def figure3(n_test, y_test, preds, scores, mean_r_test, appendix=False):
    fig, ax = plt.subplots(figsize=(6.8, 3.55))
    ax.axvspan(n_test[0] - 0.5, n_test[-1] + 0.5, color="#FFF4E0", zorder=0, lw=0)
    ax.plot(n_test, y_test, color=C_TRUE, lw=1.85, zorder=5, label="true $x[n+1]$")
    ax.plot(n_test, preds["persist"], color=C_PERSIST, lw=1.05, zorder=3,
            label=f"persistence $x[n]$  {scores['persist']:.4f}")
    ax.plot(n_test, preds["ar"], color=C_AR, lw=1.35, zorder=4,
            label=f"AR $m=10$  {scores['ar']:.4f}")
    ax.plot(n_test, preds["field"], color=C_FIELD, lw=1.25, ls=(0, (3.2, 1.4)),
            zorder=4, label=f"field-only AHL  {scores['field']:.4f}")
    ax.plot(n_test, preds["f408"], color=C_F408, lw=1.35, zorder=4,
            label=f"driven F408  {scores['f408']:.4f}")
    if appendix:
        ax.plot(n_test, preds["unmasked"], color=C_UNMASKED, lw=1.15,
                ls=(0, (1.4, 1.0)), zorder=3,
                label=f"unmasked RL  {scores['unmasked']:.4f}")
    ax2 = ax.twinx()
    ax2.plot(n_test, mean_r_test, color=C_AHL, lw=0.7, alpha=0.7, zorder=2)
    ax2.set_ylabel("mean $R$ (state, not a predictor)", color=C_AHL)
    ax2.tick_params(axis="y", labelcolor=C_AHL)
    ax2.set_ylim(-0.05, 1.05)
    ax.set_xlim(149.5, 199.5)
    ax.set_xlabel("test window $n$")
    ax.set_ylabel("$x[n+1]$")
    ax.spines["top"].set_visible(False)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0),
              handlelength=2.3, borderaxespad=0.15, labelspacing=0.25)
    title = "Test follow at $k=10$ AUTO_X"
    if appendix:
        title = "Appendix: occupancy-free kinetic surrogate, not a living-cell score"
    ax.set_title(title, loc="left", fontsize=9)
    stem = "figA_unmasked_overlay" if appendix else "fig3_test_follow"
    return save(fig, stem)


def figure4(y_test, preds, scores):
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), sharex=True, sharey=True)
    panels = (
        ("AR $m=10$", preds["ar"], scores["ar"], C_AR),
        ("field-only AHL", preds["field"], scores["field"], C_FIELD),
        ("driven F408", preds["f408"], scores["f408"], C_F408),
    )
    lo = min(float(np.min(y_test)), min(float(np.min(preds[k])) for k in ("ar", "field", "f408")))
    hi = max(float(np.max(y_test)), max(float(np.max(preds[k])) for k in ("ar", "field", "f408")))
    pad = 0.08 * (hi - lo)
    lim = (lo - pad, hi + pad)
    for ax, (title, pred, nrmse, color) in zip(axes, panels):
        ax.plot(lim, lim, color="#888888", ls=":", lw=0.9, zorder=0)
        ax.scatter(y_test, pred, s=16, c=color, alpha=0.85, edgecolors="none", zorder=2)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{title}\nNRMSE {nrmse:.4f}", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("true $x[n+1]$")
    axes[0].set_ylabel("predicted $x[n+1]$")
    fig.suptitle("Test residuals at identical limits", fontsize=10, y=1.04)
    fig.tight_layout()
    return save(fig, "fig4_residual_scatter")


def figure5(n, u, x, mean_ahl, mean_r, mean_l, r_ru):
    fig, (ax_t, ax_s) = plt.subplots(
        2, 1, figsize=(6.4, 5.4),
        gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.38},
    )
    split_bands(ax_t)
    ax_t.plot(n, mean_r, color=C_F408, lw=1.2, label="mean $R$")
    ax_t.plot(n, mean_l, color=C_L, lw=1.2, label="mean $L$")
    ax_a = ax_t.twinx()
    ax_a.plot(n, mean_ahl, color=C_AHL, lw=1.0, alpha=0.9, label="mean AHL")
    ax_t.set_xlim(-1, 200)
    ax_t.set_xlabel("window $n$")
    ax_t.set_ylabel("receiver / reporter")
    ax_a.set_ylabel("mean AHL ($\\mu$M)")
    ax_t.spines["top"].set_visible(False)
    h1, l1 = ax_t.get_legend_handles_labels()
    h2, l2 = ax_a.get_legend_handles_labels()
    ax_t.legend(h1 + h2, l1 + l2, frameon=False, loc="upper right")
    ax_t.text(20, 0.92, "washout", fontsize=7, color="#666666", ha="center",
              transform=ax_t.get_xaxis_transform())
    ax_t.text(95, 0.92, "train", fontsize=7, color="#3d6b45", ha="center",
              transform=ax_t.get_xaxis_transform())
    ax_t.text(175, 0.92, "test", fontsize=7, color="#8a5a00", ha="center",
              transform=ax_t.get_xaxis_transform())
    ax_t.text(-0.12, 1.08, "A", transform=ax_t.transAxes, fontsize=12, fontweight="bold")
    ax_t.set_title("Living state follows the drive (occupancy, not a score)", loc="left", fontsize=9)

    ax_s.scatter(u, mean_r, s=12, c=C_F408, alpha=0.7, edgecolors="none", label="mean $R$")
    ax_s.scatter(u, mean_l, s=12, c=C_L, alpha=0.55, edgecolors="none", label="mean $L$")
    ax_s.set_xlabel("affine drive $u$")
    ax_s.set_ylabel("window-mean state")
    ax_s.spines["top"].set_visible(False)
    ax_s.spines["right"].set_visible(False)
    ax_s.legend(frameon=False, loc="upper left")
    ax_s.text(0.04, 0.92, f"$r(\\mathrm{{mean\\ }}R, u)={r_ru:.2f}$",
              transform=ax_s.transAxes, fontsize=8)
    ax_s.text(-0.12, 1.08, "B", transform=ax_s.transAxes, fontsize=12, fontweight="bold")
    return save(fig, "fig5_living_state")


def write_pred_csv(path, n_test, x, x_next, u, preds, mean_r, mean_l, mean_ahl):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "n", "x_n", "x_next", "u", "pred_persist", "pred_ar", "pred_field",
            "pred_f408", "pred_masked_rl", "pred_unmasked_rl", "mean_R",
            "mean_L", "mean_AHL",
        ])
        for i, window in enumerate(n_test):
            w = int(window)
            writer.writerow([
                w,
                f"{x[w]:.12f}",
                f"{x_next[w]:.12f}",
                f"{u[w]:.12f}",
                f"{preds['persist'][i]:.12f}",
                f"{preds['ar'][i]:.12f}",
                f"{preds['field'][i]:.12f}",
                f"{preds['f408'][i]:.12f}",
                f"{preds['masked'][i]:.12f}",
                f"{preds['unmasked'][i]:.12f}",
                f"{mean_r[w]:.12f}",
                f"{mean_l[w]:.12f}",
                f"{mean_ahl[w]:.12f}",
            ])


def write_report(scores, r_ru, mean_r):
    path = HERE / "results" / "LORENZ_L3_FIGURES.md"
    lines = [
        "# Lorenz L3 figures",
        "",
        "Analysis and plotting only. Frozen $k=10$ AUTO_X scout, seed 101.",
        "Closed ridge unchanged. Living-layer FAIL is not rewritten.",
        "Unmasked RL is occupancy-free and is not a living-cell score.",
        "The F408 vs AR $0.002$ tick is kept at four decimals.",
        "",
        f"u SHA-256 `{L3.SELECTED_SHA}`.",
        "",
        "Reproduce from the repo root:",
        "",
        "```",
        "python examples/BSimReservoirPlanLorenzL3/plot_lorenz_l3.py",
        "```",
        "",
        "## Sanity (within 1e-3 of the scout)",
        "",
        "| Quantity | Figure pack | Scout |",
        "|---|---|---|",
        f"| AR $m=10$ | {scores['ar']:.4f} | 0.8233 |",
        f"| driven F408 | {scores['f408']:.4f} | 0.8256 |",
        f"| field | {scores['field']:.4f} | 0.7973 |",
        f"| Brownian | {scores['brownian']:.4f} | 1.0076 |",
        f"| silent | {scores['silent']:.4f} | 1.0076 |",
        f"| occupancy-masked RL | {scores['masked']:.4f} | 0.7967 |",
        f"| unmasked RL | {scores['unmasked']:.4f} | 0.3477 |",
        "",
        "Test-window series:",
        "`results/figures/l3_k10_test_predictions.csv`.",
        "",
        "## Figure 1 — clock map",
        "",
        "- Files: `results/figures/fig1_clock_map.pdf`, `.png`",
        "- Claim: an intermediate Lorenz sample clock exists between",
        "  trivial $\\Delta t=0.02$ and SKIP=50 destruction; $k=10$ AUTO_X",
        "  is that clock.",
        "- Not: a living-layer win, a BSim score, or a CROSS_Y scout.",
        "  Panel B is Module-1 only and was not sent to the dish.",
        "",
        "## Figure 2 — what was injected",
        "",
        "- Files: `results/figures/fig2_injected_signal.pdf`, `.png`",
        "- Claim: the dish saw one affine-mapped Lorenz $x$ sequence at",
        "  $k=10$, and AUTO_X is the $(x[n],x[n+1])$ map of that sequence.",
        "- Not: a reservoir score or an attractor embedding figure.",
        "",
        "## Figure 3 — test follow",
        "",
        "- Files: `results/figures/fig3_test_follow.pdf`, `.png`",
        "- Claim: on test windows, driven F408 tracks $x[n+1]$ about as",
        f"  well as AR $m=10$ ({scores['f408']:.4f} vs {scores['ar']:.4f}),",
        f"  not better; field-only is at least as tight ({scores['field']:.4f}).",
        "- Not: a living-layer pass. Mean $R$ on the twin axis is dish",
        "  state, not a predictor. The $0.002$ tick is visible in the legend.",
        "",
        "## Figure 4 — residual comparison",
        "",
        "- Files: `results/figures/fig4_residual_scatter.pdf`, `.png`",
        "- Claim: AR, field, and F408 occupy the same residual scale;",
        "  F408 is not a tighter map once axes are shared.",
        "- Not: a cropped zoom that hides the AR–F408 gap.",
        "",
        "## Figure 5 — living state follows the drive",
        "",
        "- Files: `results/figures/fig5_living_state.pdf`, `.png`",
        f"- Claim: occupancy is ALIVE ($r(\\mathrm{{mean\\ }}R,u)={r_ru:.2f}$,",
        f"  mean $R={mean_r:.3f}$); $L$ lags $R$.",
        "- Not: a performance claim, and not evidence that biology beats AR.",
        "",
        "## Appendix — unmasked kinetic surrogate",
        "",
        "- Files: `results/figures/figA_unmasked_overlay.pdf`, `.png`",
        f"- Claim: occupancy-free R/L kinetics on the AHL field score",
        f"  {scores['unmasked']:.4f} on this target.",
        "- Not: a living-cell result. Empty voxels are included. Do not",
        "  present $0.35$ as biology beating AR.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {path}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    u, x, x_next, y_now, y_next = L3.load_selected_target()
    del y_now, y_next
    digest = L3.sha256_u(u)
    if digest != L3.SELECTED_SHA:
        raise SystemExit(f"STOP u hash {digest}")
    n = np.arange(NUM_WINDOWS)
    target = x_next

    baselines = L3.legal_baselines(x, x_next, np.zeros_like(x))
    ar = fit_closed_predictions(delay_matrix(x, 10), target)
    persist_test = x[WASHOUT + TRAIN:]
    persist_nrmse = BC.nrmse(target[WASHOUT + TRAIN:], persist_test)
    assert_close("AR m=10", ar["test_nrmse"], EXPECTED["ar"])
    assert_close("AR from legal_baselines", baselines["best_test"], EXPECTED["ar"])

    driven_dir = HERE / "results" / "l3_k10_driven_seed101"
    brown_dir = HERE / "results" / "l3_k10_brownian_seed101"
    silent_dir = L3.SILENT_DIRS[101]
    if not L3.csv_ok(driven_dir) or not L3.csv_ok(brown_dir) or not L3.csv_ok(silent_dir):
        raise SystemExit("STOP CSV completeness")

    f408_mat = BC.read_window_matrix(
        driven_dir / "voxels.csv", BC.BIOLOGY_MEAN_PREFIXES, BC.BIOLOGY_LAST_PREFIXES
    )
    field_mat = BC.read_window_matrix(driven_dir / "voxels.csv", BC.FIELD_MEAN_PREFIXES, ())
    brown_mat = BC.read_window_matrix(
        brown_dir / "voxels.csv", BC.BROWNIAN_MEAN_PREFIXES, BC.BROWNIAN_LAST_PREFIXES
    )
    silent_mat = BC.read_window_matrix(
        silent_dir / "voxels.csv", BC.BIOLOGY_MEAN_PREFIXES, BC.BIOLOGY_LAST_PREFIXES
    )
    f408 = fit_closed_predictions(f408_mat["X"], target)
    field = fit_closed_predictions(field_mat["X"], target)
    brown = fit_closed_predictions(brown_mat["X"], target)
    silent = fit_closed_predictions(silent_mat["X"], target)

    voxels = L3.load_voxel_arrays(driven_dir / "voxels.csv")
    occ = L3.occupancy(voxels, u)
    sr, sl = L3.integrate_surrogate(voxels["absolute_time"], voxels["ahl"])
    mask = (voxels["den"] != 0).astype(float)
    unmasked_X = np.column_stack([
        L3.window_mean_matrix(sr, voxels["window"]),
        L3.window_mean_matrix(sl, voxels["window"]),
    ])
    masked_X = np.column_stack([
        L3.window_mean_matrix(sr * mask, voxels["window"]),
        L3.window_mean_matrix(sl * mask, voxels["window"]),
    ])
    unmasked = fit_closed_predictions(unmasked_X, target)
    masked = fit_closed_predictions(masked_X, target)

    assert_close("driven F408", f408["test_nrmse"], EXPECTED["f408"])
    assert_close("field", field["test_nrmse"], EXPECTED["field"])
    assert_close("Brownian", brown["test_nrmse"], EXPECTED["brownian"])
    assert_close("silent", silent["test_nrmse"], EXPECTED["silent"])
    assert_close("occupancy-masked RL", masked["test_nrmse"], EXPECTED["masked"])
    assert_close("unmasked RL", unmasked["test_nrmse"], EXPECTED["unmasked"])
    print("SANITY PASS")
    print(f"  AR={ar['test_nrmse']:.4f} F408={f408['test_nrmse']:.4f} "
          f"field={field['test_nrmse']:.4f} persist={persist_nrmse:.4f}")
    print(f"  brown={brown['test_nrmse']:.4f} silent={silent['test_nrmse']:.4f} "
          f"masked={masked['test_nrmse']:.4f} unmasked={unmasked['test_nrmse']:.4f}")

    mean_r = np.mean(L3.window_mean_matrix(voxels["receiver"], voxels["window"]), axis=1)
    mean_l = np.mean(L3.window_mean_matrix(voxels["lum"], voxels["window"]), axis=1)
    mean_ahl = np.mean(L3.window_mean_matrix(voxels["ahl"], voxels["window"]), axis=1)
    scores = {
        "ar": ar["test_nrmse"],
        "f408": f408["test_nrmse"],
        "field": field["test_nrmse"],
        "persist": persist_nrmse,
        "brownian": brown["test_nrmse"],
        "silent": silent["test_nrmse"],
        "masked": masked["test_nrmse"],
        "unmasked": unmasked["test_nrmse"],
    }
    preds = {
        "persist": persist_test,
        "ar": ar["pred_test"],
        "field": field["pred_test"],
        "f408": f408["pred_test"],
        "masked": masked["pred_test"],
        "unmasked": unmasked["pred_test"],
    }
    y_test = target[WASHOUT + TRAIN:]

    csv_path = OUT / "l3_k10_test_predictions.csv"
    write_pred_csv(csv_path, N_TEST, x, x_next, u, preds, mean_r, mean_l, mean_ahl)
    print(f"wrote {csv_path.name}")

    clock_rows = load_clock_csv()
    figure1(clock_rows)
    figure2(n, x, x_next, u)
    figure3(N_TEST, y_test, preds, scores, mean_r[WASHOUT + TRAIN:], appendix=False)
    figure3(N_TEST, y_test, preds, scores, mean_r[WASHOUT + TRAIN:], appendix=True)
    figure4(y_test, preds, scores)
    figure5(n, u, x, mean_ahl, mean_r, mean_l, occ["r_meanR_u"])
    write_report(scores, occ["r_meanR_u"], occ["mean_R"])
    print("FIGURE PACK COMPLETE")
    print("living-layer FAIL unchanged; unmasked is not biology")


if __name__ == "__main__":
    main()
