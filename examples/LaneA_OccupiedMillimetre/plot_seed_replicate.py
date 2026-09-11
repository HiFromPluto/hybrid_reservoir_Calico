#!/usr/bin/env python3
"""IEEE seed-replicate figure. Standing numbers only. No Java.

Colour house style (global with plot_manuscript_figures.py):
  orange = living / placement seeds
  green  = field
  grey   = silent / target-adjacent
  blue   = ceiling only (not used here)
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "seed_replicate"
MS = HERE.parent / "PocketDish" / "manuscript_ieee" / "figures"

C_TARGET = "#111111"
C_FIELD = "#009E73"
C_RL = "#D55E00"
C_SILENT = "#999999"
C_SEED = {101: "#8C3100", 202: "#D55E00", 303: "#E69F00"}
C_LIVING = "#D55E00"


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "legend.fontsize": 6.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
        }
    )


def load_trace(task: str, seed: int):
    path = OUT / f"trace_{task}_seed{seed}.csv"
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return rows


def draw_forecast(ax, task, ylabel, title, living_key="RL"):
    t101 = load_trace(task, 101)
    n = np.array([float(r["n"]) for r in t101])
    y = np.array([float(r["y"]) for r in t101])
    field = np.array([float(r["FIELD"]) for r in t101])
    ax.plot(n, y, color=C_TARGET, lw=1.35, zorder=3)
    ax.plot(n, field, color=C_FIELD, ls="--", lw=1.15, zorder=2)
    for seed in (101, 202, 303):
        rows = load_trace(task, seed)
        pred = np.array([float(r[living_key]) for r in rows])
        ax.plot(n, pred, color=C_SEED[seed], lw=0.95 if seed == 101 else 0.75, zorder=1)
    ax.set_xlim(n[0], n[-1])
    ax.set_xlabel("window $n$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return n


def annotate_panel(ax, letter: str) -> None:
    ax.text(
        0.0,
        1.14,
        letter,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def main() -> None:
    summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
    tasks = summary["tasks"]
    style()

    fig = plt.figure(figsize=(7.16, 5.85))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.08, 1.0], hspace=0.55, wspace=0.34)

    ax_n = fig.add_subplot(gs[0, 0])
    ax_m = fig.add_subplot(gs[0, 1])
    ax_l = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_b = fig.add_subplot(gs[1, 1])
    ax_w = fig.add_subplot(gs[1, 2])

    na = tasks["narma10"]
    mg = tasks["mg"]
    lo = tasks["lorenz"]
    draw_forecast(
        ax_n,
        "narma10",
        "$y[n]$",
        rf"NARMA-10   $R\|L$ {na['living_mean']:.3f} / field {na['field']:.3f}",
    )
    draw_forecast(
        ax_m,
        "mg",
        "$x[n+1]$",
        rf"Mackey--Glass   $R\|L$ {mg['living_mean']:.3f} / field {mg['field']:.3f}",
    )
    draw_forecast(
        ax_l,
        "lorenz",
        "$x[n+10]$",
        rf"Lorenz $k{{=}}10$   $R\|L$ {lo['living_mean']:.3f} / field {lo['field']:.3f}",
    )

    t101 = load_trace("carrier", 101)
    idx = np.arange(len(t101))
    ax_c.plot(idx, [float(r["y"]) for r in t101], color=C_TARGET, lw=1.35, zorder=3)
    ax_c.plot(idx, [float(r["FIELD"]) for r in t101], color=C_FIELD, ls="--", lw=1.15, zorder=2)
    for seed in (101, 202, 303):
        rows = load_trace("carrier", seed)
        ax_c.plot(idx, [float(r["R"]) for r in rows], color=C_SEED[seed], lw=0.95 if seed == 101 else 0.75)
    car = tasks["carrier"]
    ax_c.set_xlabel("test sample")
    ax_c.set_ylabel("$u(t)$")
    ax_c.set_title(rf"Carrier pulse   $R$ {car['living_mean']:.3f} / field {car['field']:.3f}", pad=4)
    ax_c.set_xlim(0, len(t101) - 1)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)

    names = ["Pulse", "MG", "NARMA", "Lorenz"]
    keys = ["carrier", "mg", "narma10", "lorenz"]
    living_name = ["R", "RL", "RL", "RL"]
    silent_vals = []
    living_mu = []
    living_sd = []
    field_vals = []
    for key, lname in zip(keys, living_name):
        per = tasks[key]["per_seed"]
        if key == "carrier":
            liv = np.array([per[str(s)]["R"] for s in (101, 202, 303)])
            silent_vals.append(1.0)
        else:
            liv = np.array([per[str(s)][lname] for s in (101, 202, 303)])
            silent_vals.append(float(np.mean([per[str(s)]["SILENT_RL"] for s in (101, 202, 303)])))
        living_mu.append(liv.mean())
        living_sd.append(liv.std(ddof=1))
        field_vals.append(tasks[key]["field"])

    x = np.arange(len(names))
    w = 0.26
    ax_b.bar(x - w, silent_vals, w, color=C_SILENT, edgecolor="#666", lw=0.4, label="silent")
    ax_b.bar(
        x,
        living_mu,
        w,
        yerr=living_sd,
        color=C_LIVING,
        edgecolor="#8C3100",
        lw=0.4,
        error_kw={"ecolor": "black", "lw": 0.6, "capsize": 1.5},
        label=r"$R\|L$",
    )
    ax_b.bar(x + w, field_vals, w, color=C_FIELD, edgecolor="#00664A", lw=0.4, label="field")
    ax_b.axhline(1.0, color="#666", ls=":", lw=0.6)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(names)
    ax_b.set_ylabel("test NRMSE")
    ax_b.set_ylim(0, 1.28)
    ax_b.set_title("Living vs field vs silent", pad=4)
    ax_b.legend(frameon=False, fontsize=6, loc="upper right")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    for i, (mu, fld) in enumerate(zip(living_mu, field_vals)):
        if mu < 0.12 and fld < 0.12:
            ax_b.annotate(
                f"{mu:.3f}",
                xy=(i, mu),
                xytext=(i - 0.08, 0.22),
                fontsize=5.5,
                color=C_LIVING,
                ha="right",
                arrowprops=dict(arrowstyle="-", color=C_LIVING, lw=0.4),
            )
            ax_b.annotate(
                f"{fld:.3f}",
                xy=(i + w, fld),
                xytext=(i + w + 0.08, 0.34),
                fontsize=5.5,
                color="#00664A",
                ha="left",
                arrowprops=dict(arrowstyle="-", color="#00664A", lw=0.4),
            )
        else:
            ax_b.text(i, mu + 0.045, f"{mu:.3f}", ha="center", va="bottom", fontsize=5.5, color=C_LIVING)
            ax_b.text(i + w, fld + 0.03, f"{fld:.3f}", ha="center", va="bottom", fontsize=5.5, color="#00664A")

    wt = tasks["waveform"]
    perw = wt["per_seed"]
    brier_rl = float(np.mean([perw[str(s)]["RL_brier"] for s in (101, 202, 303)]))
    win_rl = float(np.mean([perw[str(s)]["RL_window_auc"] for s in (101, 202, 303)]))
    groups = ["Brier", "win. AUC", "block AUC"]
    living_f = [brier_rl, win_rl, wt["legacy_auc_rl"]]
    field_f = [wt["field"], wt["window_auc_field"], wt["legacy_auc_field"]]
    xf = np.arange(len(groups))
    wf = 0.36
    bars_l = ax_w.bar(
        xf - wf / 2,
        living_f,
        wf,
        color=C_LIVING,
        edgecolor="#8C3100",
        lw=0.4,
        label=r"$R\|L$",
    )
    ax_w.bar(xf + wf / 2, field_f, wf, color=C_FIELD, edgecolor="#00664A", lw=0.4, label="field")
    bars_l[2].set_alpha(0.45)
    ax_w.set_xticks(xf)
    ax_w.set_xticklabels(groups)
    ax_w.set_ylim(0, 1.22)
    ax_w.set_ylabel("score")
    ax_w.set_title("Waveform2c (Brier is the axis)", pad=4)
    ax_w.axhline(0.5, color="#666", ls=":", lw=0.5)
    ax_w.legend(frameon=False, fontsize=6, loc="upper left")
    ax_w.spines["top"].set_visible(False)
    ax_w.spines["right"].set_visible(False)
    for i, (lv, fv) in enumerate(zip(living_f, field_f)):
        ax_w.text(i - wf / 2, lv + 0.03, f"{lv:.3f}", ha="center", fontsize=5.5, color=C_LIVING)
        ax_w.text(i + wf / 2, fv + 0.03, f"{fv:.3f}", ha="center", fontsize=5.5, color="#00664A")

    handles = [
        Line2D([0], [0], color=C_TARGET, lw=1.4, label="target"),
        Line2D([0], [0], color=C_FIELD, ls="--", lw=1.15, label="field"),
        Line2D([0], [0], color=C_SEED[101], lw=1.0, label="seed 101"),
        Line2D([0], [0], color=C_SEED[202], lw=1.0, label="seed 202"),
        Line2D([0], [0], color=C_SEED[303], lw=1.0, label="seed 303"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
        fontsize=7,
    )

    for ax, letter in zip((ax_n, ax_m, ax_l, ax_c, ax_b, ax_w), "ABCDEF"):
        annotate_panel(ax, letter)

    OUT.mkdir(parents=True, exist_ok=True)
    MS.mkdir(parents=True, exist_ok=True)
    for dest in (OUT, MS):
        fig.savefig(dest / "fig_seed_replicate.pdf")
        fig.savefig(dest / "fig_seed_replicate.png", dpi=300)
    plt.close(fig)
    print("wrote", OUT / "fig_seed_replicate.png")
    print("AXIS", summary["LANE_A_SEED_REPLICATE"], summary["n_sign_match"], "/5")


if __name__ == "__main__":
    main()
