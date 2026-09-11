#!/usr/bin/env python3
"""Lane A closed-benchmark figures from frozen standings. No Java.

Numbers are the standing three-decimal (or exact AUC) values. Does not
rewrite occupancy, carrier FAIL, Narma10b 0.928, or HybridDish GATE_EVIDENCE.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "figures"
# Committed gallery for the PocketDish report (results/ is gitignored).
OUT_COMMITTED = HERE.parent / "PocketDish" / "figures" / "lane_a"

C_RL = "#D55E00"
C_FIELD = "#009E73"
C_SILENT = "#999999"
C_CEIL = "#0072B2"
C_OCC = "#882255"
C_R = "#56B4E9"


def style() -> None:
    plt.rcParams.update(
        {
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
        }
    )


def save(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    OUT_COMMITTED.mkdir(parents=True, exist_ok=True)
    for dest in (OUT, OUT_COMMITTED):
        fig.savefig(dest / f"{stem}.pdf")
        fig.savefig(dest / f"{stem}.png", dpi=300)
    plt.close(fig)
    print(f"wrote {stem}.pdf {stem}.png")


def grouped(ax, labels, series, ylabel, title, y0=0.0, y1=None, hline=None, hlabel=None):
    x = np.arange(len(labels))
    n = len(series)
    width = min(0.22, 0.8 / max(n, 1))
    offset = (n - 1) * width / 2
    for i, (name, vals, color) in enumerate(series):
        ax.bar(x - offset + i * width, vals, width, label=name, color=color, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(y0, y1)
    if hline is not None:
        ax.axhline(hline, color="#666666", ls="--", lw=0.8, zorder=1)
        if hlabel:
            ax.text(len(labels) - 0.45, hline + 0.02 * (y1 or 1), hlabel, fontsize=7, color="#666666")
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fig1_nrmse_tasks() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    grouped(
        ax,
        ["NARMA-10", "MG one-step", "Lorenz k=10"],
        [
            ("RL (living)", [0.895, 0.031, 0.439], C_RL),
            ("Field AHL", [0.945, 0.012, 0.797], C_FIELD),
            ("Silent", [1.066, 1.042, 1.008], C_SILENT),
            ("Ceiling", [0.840, 0.003, 0.823], C_CEIL),
        ],
        "Test NRMSE (lower is better)",
        "Lane A system extras — NRMSE",
        y0=0,
        y1=1.2,
        hline=1.0,
        hlabel="mean predictor",
    )
    fig.text(
        0.01,
        -0.02,
        "Ceiling: NARMA = 10-tap of u; MG = 10-tap of affine u (0.003); Lorenz = AR m=10 of x (0.823). "
        "Source: LANE_A_*_STANDING.md, 2026-08-28. Not Narma10b 0.928.",
        fontsize=7,
        color="#444444",
    )
    save(fig, "fig1_nrmse_tasks")


def fig2_waveform_auc() -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    labels = ["RL", "R", "L", "Field", "Silent", "MOMENTS", "RAW_U_8"]
    vals = [1.000, 1.000, 1.000, 0.820, 0.500, 0.500, 1.000]
    colors = [C_RL, C_R, C_OCC, C_FIELD, C_SILENT, "#CC79A7", "#000000"]
    ax.bar(labels, vals, color=colors, zorder=2)
    ax.axhline(0.5, color="#666666", ls="--", lw=0.8, zorder=1)
    ax.text(6.35, 0.52, "chance", fontsize=7, color="#666666")
    ax.set_ylabel("Test-block macro OVR AUC (higher is better)")
    ax.set_title("Lane A Waveform2c — 16 test blocks")
    ax.set_ylim(0, 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        -0.02,
        "u-only PASS: MOMENTS=0.500, RAW_U_8=1.000 (answer key). 16 blocks can saturate AUC at 1.0. "
        "Not Waveform2c 0.835. Source: LANE_A_WAVEFORM_STANDING.md.",
        fontsize=7,
        color="#444444",
    )
    save(fig, "fig2_waveform_auc")


def fig3_carrier() -> None:
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    labels = ["Field", "R", "L", "Silent R", "10-tap u"]
    vals = [0.444, 0.476, 0.571, 1.000, 0.544]
    colors = [C_FIELD, C_R, C_OCC, C_SILENT, C_CEIL]
    ax.bar(labels, vals, color=colors, zorder=2)
    ax.axhline(1.0, color="#666666", ls="--", lw=0.8, zorder=1)
    ax.set_ylabel("Test NRMSE reconstructing u(t)")
    ax.set_title("Lane A carrier — occupancy pulse train (FAIL vs field)")
    ax.set_ylim(0, 1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        -0.02,
        "Sample-level reconstruction of U_PULSE_TRAIN. Window-mean u VOID. Occupancy still PASS. "
        "Source: LANE_A_CARRIER_STANDING.md.",
        fontsize=7,
        color="#444444",
    )
    save(fig, "fig3_carrier_nrmse")


def fig4_occupancy() -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    labels = ["Pulse train", "NARMA", "MG", "Lorenz", "Waveform"]
    vals = [0.661, 0.359, 0.447, 0.401, 0.352]
    ax.bar(labels, vals, color=C_OCC, zorder=2)
    ax.axhline(0.05, color="#D55E00", ls="--", lw=0.9, zorder=1)
    ax.text(4.05, 0.07, "ALIVE 0.05", fontsize=7, color="#D55E00")
    ax.set_ylabel(r"Test-band mean $R$")
    ax.set_title("Lane A occupancy — all scored extras ALIVE")
    ax.set_ylim(0, 0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        -0.02,
        "Pulse-train occupancy gate used windows 4–7; NARMA/MG/Lorenz windows 150–199; "
        "Waveform windows 272–399. J_max never raised. Source: occupancy and extra standings.",
        fontsize=7,
        color="#444444",
    )
    save(fig, "fig4_occupancy")


def fig5_who_wins() -> None:
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    tasks = ["Pulse u", "NARMA-10", "MG", "Lorenz", "Waveform2c"]
    # 1 = living beats field, 0 = field wins or ties (carrier FAIL / MG)
    living = [0, 1, 0, 1, 1]
    field = [1, 0, 1, 0, 0]
    x = np.arange(len(tasks))
    ax.bar(x - 0.18, living, 0.36, label="RL beats field", color=C_RL, zorder=2)
    ax.bar(x + 0.18, field, 0.36, label="Field wins or ties", color=C_FIELD, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(tasks)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["no", "yes"])
    ax.set_ylabel("Named comparison")
    ax.set_title("Living-layer vs field — task-specific, not general")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        -0.02,
        "Pulse reconstruction and MG one-step live in the plume. NARMA, Lorenz k=10 AUTO_X, and "
        "Waveform2c are named notes on those u only. Not a general “cells beat the plume” claim.",
        fontsize=7,
        color="#444444",
    )
    save(fig, "fig5_living_vs_field")


def main() -> None:
    style()
    fig1_nrmse_tasks()
    fig2_waveform_auc()
    fig3_carrier()
    fig4_occupancy()
    fig5_who_wins()
    print(f"figures in {OUT} and {OUT_COMMITTED}")


if __name__ == "__main__":
    main()
