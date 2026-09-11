#!/usr/bin/env python3
"""Figures for the IEEE manuscript. Standing numbers only. No Java."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import csv

HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
LANE_A = HERE.parent.parent / "LaneA_OccupiedMillimetre"
LANE_A_CSV = LANE_A  # unused
CHEY = HERE.parent.parent / "LaneB_ChemotacticSpatial" / "results"
MC_JSON = LANE_A / "results" / "lane_a_memory_capacity_summary.json"

C_RL = "#D55E00"
C_FIELD = "#009E73"
C_SILENT = "#999999"
C_CEIL = "#0072B2"


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
        }
    )


def save(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    plt.close(fig)
    print("wrote", stem)


def fig_nrmse() -> None:
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    labs = ["NARMA-10", "MG", "Lorenz"]
    x = np.arange(len(labs))
    w = 0.2
    ax.bar(x - 1.5 * w, [0.895, 0.031, 0.439], w, label=r"$R\|L$", color=C_RL)
    ax.bar(x - 0.5 * w, [0.945, 0.012, 0.797], w, label="Field", color=C_FIELD)
    ax.bar(x + 0.5 * w, [1.066, 1.042, 1.008], w, label="Silent", color=C_SILENT)
    ax.bar(x + 1.5 * w, [0.840, 0.003, 0.823], w, label="Ceiling", color=C_CEIL)
    ax.axhline(1.0, color="#666", ls="--", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylabel("Test NRMSE")
    ax.set_ylim(0, 1.2)
    ax.legend(frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_nrmse")


def fig_waveform() -> None:
    fig, ax = plt.subplots(figsize=(3.5, 2.3))
    labs = [r"$R\|L$", "Field", "Silent", "MOMENTS", "RAW_U"]
    vals = [1.0, 0.820, 0.5, 0.5, 1.0]
    cols = [C_RL, C_FIELD, C_SILENT, "#CC79A7", "#000"]
    ax.bar(labs, vals, color=cols)
    ax.axhline(0.5, color="#666", ls="--", lw=0.6)
    ax.set_ylabel("Block macro AUC")
    ax.set_ylim(0, 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_waveform")


def fig_carrier() -> None:
    fig, ax = plt.subplots(figsize=(3.4, 2.3))
    labs = ["Field", r"$R$", r"$L$", "Silent", "10-tap"]
    vals = [0.444, 0.476, 0.571, 1.0, 0.544]
    ax.bar(labs, vals, color=[C_FIELD, C_RL, "#E69F00", C_SILENT, C_CEIL])
    ax.set_ylabel("Test NRMSE reconstructing $u$")
    ax.set_ylim(0, 1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_carrier")


def fig_occupancy() -> None:
    fig, ax = plt.subplots(figsize=(3.4, 2.2))
    labs = ["Pulse", "NARMA", "MG", "Lorenz", "Wave"]
    ax.bar(labs, [0.661, 0.359, 0.447, 0.401, 0.352], color="#882255")
    ax.axhline(0.05, color=C_RL, ls="--", lw=0.8)
    ax.set_ylabel(r"Test-band mean $R$")
    ax.set_ylim(0, 0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_occupancy")


def fig_axis() -> None:
    fig, ax = plt.subplots(figsize=(7.0, 2.15))
    names = ["Pulse $u$", "MG one-step", "NARMA-10", "Lorenz $k{=}10$", "Waveform2c"]
    living = [0, 0, 1, 1, 1]
    colors = [C_FIELD if v == 0 else C_RL for v in living]
    ax.barh(names[::-1], [1] * 5, color=colors[::-1], height=0.55)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_xlabel("Winner of living vs field  (green = field, orange = $R\\|L$)")
    for i, (n, v) in enumerate(zip(names[::-1], living[::-1])):
        ax.text(0.5, i, "field" if v == 0 else r"living $R\|L$", ha="center", va="center",
                color="white", fontsize=8, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    save(fig, "fig_axis")


def fig_hybriddish() -> None:
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    x = np.arange(4)
    w = 0.35
    hd = [0.928, 0.260, 0.826, 0.835]
    la = [0.895, 0.031, 0.439, 1.000]
    ax.bar(x - w / 2, hd, w, label="HybridDish living", color="#88CCEE")
    ax.bar(x + w / 2, la, w, label=r"Lane A $R\|L$", color=C_RL)
    ax.set_xticks(x)
    ax.set_xticklabels(["NARMA", "MG", "Lorenz", "WF AUC"])
    ax.set_ylabel("Living score (NRMSE / AUC)")
    ax.set_title("Predecessor vs N0 rebuild (different $u$)")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_hybriddish")


def fig_laneb_kappa() -> None:
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    labs = ["B0 22 s", "B1 60 s", "B2 1.1 s", "B3 25 s"]
    mot = [0.03, 0.010, 0.052, 0.053]
    fld = [0.931, 0.833, 0.409, 0.080]
    x = np.arange(len(labs))
    w = 0.35
    ax.bar(x - w / 2, mot, w, label=r"$|\kappa|$ motile", color=C_RL)
    ax.bar(x + w / 2, fld, w, label=r"$|\kappa|$ field", color=C_FIELD)
    ax.axhline(0.12, color="#666", ls="--", lw=0.7, label="occupation gate")
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylabel(r"Scene contrast $|\kappa|$")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_laneb_kappa")


def fig_peclet() -> None:
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    chi = np.array([0.9, 1.8, 3.5])
    kappa_ss = np.array([0.082, 0.164, 0.319])
    drift = np.array([30, 60, 121])
    ax.plot(chi, kappa_ss, "o-", color=C_FIELD, label=r"steady $\kappa$ (est.)")
    ax.axhline(0.12, color="#666", ls="--", lw=0.7, label="B3 gate 0.12")
    ax.axhline(0.053, color=C_RL, ls=":", lw=1.0, label=r"B3 motile $|\kappa|$ at 25 s")
    ax.set_xlabel(r"$\chi_0/\mu$")
    ax.set_ylabel(r"Predicted stripe $\kappa$")
    ax.set_xlim(0.6, 3.8)
    ax.set_ylim(0, 0.4)
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_peclet")


def fig_chey() -> None:
    fig, axes = plt.subplots(3, 1, figsize=(3.4, 4.15), sharex=True)
    files = [
        (CHEY / "b4_hold.csv", "Hold $L=0$"),
        (CHEY / "b4_step_up.csv", r"Step $0\to 0.3$"),
        (CHEY / "b4_step_down.csv", r"Step $0.3\to 0$"),
    ]
    for ax, (path, title) in zip(axes, files):
        t, y = [], []
        if path.exists():
            with path.open(encoding="utf-8") as h:
                for row in csv.DictReader(h):
                    t.append(float(row["t"]))
                    y.append(float(row["cheyP"]))
            ax.plot(t, y, color=C_RL, lw=1.0)
        ax.axvline(5.0, color="#666", ls="--", lw=0.6)
        ax.axhline(1.0, color="#ccc", lw=0.5)
        ax.set_title(title, loc="left", fontsize=8)
        ax.set_xlim(0, 15)
        ax.set_ylabel("CheY-P $y$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1].set_xlabel("$t$ (s)")
    fig.tight_layout(h_pad=0.35)
    save(fig, "fig_chey")


def fig_schematic() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.16, 2.35))
    # Lane A
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.add_patch(plt.Rectangle((0.4, 0.5), 9.2, 4.0, fill=False, lw=1.2))
    ax.add_patch(plt.Rectangle((3.3, 1.7), 3.4, 1.6, facecolor="#D55E00", alpha=0.25, edgecolor="#D55E00"))
    ax.plot(5.0, 2.5, "o", color="#0072B2", ms=7)
    ax.annotate("CENTER AHL", xy=(5.0, 2.5), xytext=(6.6, 3.7),
                fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.text(5.0, 2.15, "seed $400\\times200$", ha="center", fontsize=7, color="#D55E00")
    ax.set_title("Lane A  (motility OFF)")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r"$1000\times500\times10$ µm, FLOW$=0$, N0")
    for sp in ax.spines.values():
        sp.set_visible(False)
    # Lane B
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.add_patch(plt.Rectangle((0.4, 0.5), 9.2, 4.0, fill=False, lw=1.2))
    ax.add_patch(plt.Rectangle((0.4, 0.5), 2.0, 4.0, facecolor="#009E73", alpha=0.22, edgecolor="#009E73"))
    ax.text(1.4, 2.5, "slab\natt/rep", ha="center", va="center", fontsize=7)
    ax.annotate("", xy=(4.2, 2.5), xytext=(2.5, 2.5),
                arrowprops=dict(arrowstyle="->", color="#009E73", lw=1.2))
    ax.text(5.5, 2.5, r"$\nabla\ln L$  (log-sensing)", ha="center", fontsize=7)
    ax.set_title("Lane B  (AHL/Hill OFF)")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r"same box; $D=800$ µm$^2$/s; Barkai--Leibler")
    for sp in ax.spines.values():
        sp.set_visible(False)
    save(fig, "fig_schematic")


def fig_memory() -> None:
    if not MC_JSON.is_file():
        raise SystemExit(f"missing {MC_JSON}; run check_lane_a_memory_capacity.py first")
    blob = json.loads(MC_JSON.read_text(encoding="utf-8"))
    k = list(range(21))
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    ax.plot(k, blob["test_r2"]["RL"], "o-", color=C_RL, ms=3.5, lw=1.1, label=r"$R\|L$")
    ax.plot(k, blob["test_r2"]["FIELD"], "s-", color=C_FIELD, ms=3.2, lw=1.1, label="Field")
    ax.plot(k, blob["test_r2"]["SILENT_RL"], "x-", color=C_SILENT, ms=3.5, lw=0.9, label="Silent")
    ax.plot(k, blob["test_r2"]["DELAY_U_10"], "^-", color=C_CEIL, ms=3.2, lw=1.0, label="10-tap")
    ax.axhline(0.5, color="#666", ls="--", lw=0.6)
    ax.axhline(0.0, color="#ccc", lw=0.5)
    ax.set_xlabel(r"lag $k$")
    ax.set_ylabel(r"test $R^2(k)$")
    ax.set_xlim(-0.4, 20.4)
    ax.set_ylim(-0.55, 1.12)
    ax.set_xticks([0, 5, 10, 15, 20])
    ax.legend(frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_memory")


def fig_hd_narma() -> None:
    fig, ax = plt.subplots(figsize=(3.4, 2.3))
    labs = ["F408", "Field", "10-tap", "Brown.", "Silent"]
    ax.bar(labs, [0.928, 1.029, 0.683, 1.162, 1.162],
           color=[C_RL, C_FIELD, C_CEIL, C_SILENT, C_SILENT])
    ax.axhline(1.0, color="#666", ls="--", lw=0.6)
    ax.set_ylabel("Test NRMSE")
    ax.set_title("HybridDish Narma10b")
    ax.set_ylim(0, 1.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "fig_hd_narma")


def main() -> None:
    style()
    fig_nrmse()
    fig_waveform()
    fig_carrier()
    fig_occupancy()
    fig_axis()
    fig_hybriddish()
    fig_laneb_kappa()
    fig_peclet()
    fig_chey()
    fig_hd_narma()
    fig_memory()
    fig_schematic()
    print("out", OUT)


if __name__ == "__main__":
    main()
