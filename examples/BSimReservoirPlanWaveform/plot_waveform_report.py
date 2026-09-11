#!/usr/bin/env python3
"""Report figures for waveform classification. Does not rerun BSim or retune."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_waveform as cw  # noqa: E402

OUT = HERE / "results" / "figures"
FROZEN_AUC = {
    "brownian": (0.5744, 0.6540, 0.5317),
    "silent": (0.5000, 0.5000, 0.5000),
    "driven": (0.7969, 0.8287, 0.8411),
    "field": (0.8603, 0.8603, 0.8603),
}
FROZEN_ACC = {
    "brownian": (0.3200, 0.3000, 0.3200),
    "silent": (0.3000, 0.3000, 0.3000),
    "driven": (0.7600, 0.7600, 0.6800),
    "field": (0.8600, 0.8600, 0.8600),
}

SINE = "#0072B2"
SQUARE = "#D55E00"
TRIANGLE = "#009E73"
CLASS_COLOR = (SINE, SQUARE, TRIANGLE)
ARM_COLOR = {
    "brownian": "#BBBBBB",
    "silent": "#888888",
    "driven": "#0072B2",
    "field": "#4D4D4D",
}


def assert_close(name, got, expected, tol=1e-4):
    for a, b in zip(got, expected):
        if abs(a - b) > tol:
            raise SystemExit(f"{name} mismatch: {got} vs frozen {expected}")


def test_scores(X_all, y):
    X = X_all[cw.WASHOUT :]
    yt = y[cw.WASHOUT :]
    lam, _ = cw.select_lambda(X, yt)
    X_train, y_train = X[: cw.TRAIN], yt[: cw.TRAIN]
    X_test, y_test = X[cw.TRAIN :], yt[cw.TRAIN :]
    (X_train_z, X_test_z), _, _ = cw.standardize(X_train, X_test)
    scores, _ = cw.ovr_scores(X_train_z, y_train, X_test_z, lam)
    macro, per = cw.macro_ovr_auc(y_test, scores)
    acc, pred, confusion = cw.accuracy_argmax(y_test, scores)
    return {
        "lambda": lam,
        "macro": macro,
        "per": per,
        "acc": acc,
        "pred": np.asarray(pred),
        "confusion": np.asarray(confusion),
        "scores": scores,
        "y_test": y_test,
    }


def load_all(y):
    out = {}
    for arm, prefixes in (
        ("brownian", (cw.BROWNIAN_MEAN_PREFIXES, cw.BROWNIAN_LAST_PREFIXES)),
        ("silent", (cw.BIOLOGY_MEAN_PREFIXES, cw.BIOLOGY_LAST_PREFIXES)),
        ("driven", (cw.BIOLOGY_MEAN_PREFIXES, cw.BIOLOGY_LAST_PREFIXES)),
    ):
        out[arm] = []
        for seed in cw.SEEDS:
            path = HERE / "results" / f"waveform_{arm}_seed{seed}" / "voxels.csv"
            matrix = cw.read_window_matrix(path, prefixes[0], prefixes[1])
            metrics = test_scores(matrix["X"], y)
            if arm == "driven":
                field = cw.read_window_matrix(path, cw.FIELD_MEAN_PREFIXES, ())
                metrics["field"] = test_scores(field["X"], y)
            out[arm].append(metrics)
    return out


def style():
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
    })


def fig_input(y, u):
    fig, axes = plt.subplots(
        2, 1, figsize=(7.2, 4.4), gridspec_kw={"height_ratios": [1.15, 2.4]}
    )
    w = np.arange(5)
    ax = axes[0]
    ax.plot(w, [cw.template_u(0, i) for i in w], "o-", color=SINE, lw=1.6, ms=5, label="sine")
    ax.plot(w, [cw.template_u(1, i) for i in w], "s-", color=SQUARE, lw=1.6, ms=5, label="square")
    ax.plot(w, [cw.template_u(2, i) for i in w], "^-", color=TRIANGLE, lw=1.6, ms=5, label="triangle")
    ax.set_xticks(w)
    ax.set_ylim(0, 0.55)
    ax.set_ylabel(r"$u$")
    ax.set_xlabel("index in 5-window block")
    ax.set_title("Frozen templates")
    ax.legend(frameon=False, ncol=3, loc="upper right")

    ax = axes[1]
    n = np.arange(len(u))
    for cls, color in enumerate(CLASS_COLOR):
        mask = y == cls
        ax.plot(n[mask], u[mask], ".", color=color, ms=4.5, label=cw.CLASS_NAMES[cls])
        ax.vlines(n[mask], 0, u[mask], color=color, lw=0.6, alpha=0.45)
    ax.axvspan(-0.5, 39.5, color="#f0f0f0", zorder=0)
    ax.axvspan(149.5, 199.5, color="#eeeeee", zorder=0)
    ax.axvline(39.5, color="#666666", lw=0.8, ls="--")
    ax.axvline(149.5, color="#666666", lw=0.8, ls="--")
    ax.set_xlim(-1, 200)
    ax.set_ylim(0, 0.55)
    ax.set_xlabel("window $n$")
    ax.set_ylabel(r"AHL $u[n]$")
    ax.set_title("Drive sequence (washout / train / test)")
    ax.legend(frameon=False, ncol=3, loc="upper right")
    fig.tight_layout()
    return fig


def fig_controls():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))
    names = ["brownian", "silent", "driven", "field"]
    labels = ["Brownian", "silent", "driven\nbiology", "AHL\nfield-only"]
    x = np.arange(len(names))

    ax = axes[0]
    means, ses = [], []
    for name in names:
        vals = np.array(FROZEN_AUC[name])
        means.append(float(np.mean(vals)))
        ses.append(float(np.std(vals, ddof=1) / np.sqrt(3)))
    colors = [ARM_COLOR[n] for n in names]
    ax.bar(x, means, yerr=ses, color=colors, capsize=3, width=0.72, zorder=2)
    ax.axhline(0.5, color="#333333", ls=":", lw=1, zorder=1, label="chance AUC 0.5")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel("test macro one-vs-rest AUC")
    ax.set_title("Primary metric")
    ax.legend(frameon=False, loc="upper left")
    for i, m in enumerate(means):
        ax.text(i, m + ses[i] + 0.015, f"{m:.2f}", ha="center", va="bottom", fontsize=8)

    ax = axes[1]
    means, ses = [], []
    for name in names:
        vals = np.array(FROZEN_ACC[name])
        means.append(float(np.mean(vals)))
        ses.append(float(np.std(vals, ddof=1) / np.sqrt(3)))
    ax.bar(x, means, yerr=ses, color=colors, capsize=3, width=0.72, zorder=2)
    ax.axhline(1.0 / 3.0, color="#333333", ls=":", lw=1, zorder=1, label="chance acc. 1/3")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("test accuracy (argmax)")
    ax.set_title("Diagnostic, not a gate")
    ax.legend(frameon=False, loc="upper left")
    for i, m in enumerate(means):
        ax.text(i, m + ses[i] + 0.02, f"{m:.2f}", ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    return fig


def fig_confusion(driven, field):
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.2),
                             gridspec_kw={"width_ratios": [1, 1, 1.15]})
    names = list(cw.CLASS_NAMES)
    for ax, mat, title in (
        (axes[0], driven["confusion"], "Driven seed 101"),
        (axes[1], field["confusion"], "AHL field-only"),
    ):
        im = ax.imshow(mat, cmap="Greys", vmin=0, vmax=20)
        ax.set_xticks(range(3), names, rotation=25, ha="right")
        ax.set_yticks(range(3), names)
        ax.set_xlabel("predicted")
        ax.set_ylabel("true")
        ax.set_title(title)
        for i in range(3):
            for j in range(3):
                val = int(mat[i, j])
                ax.text(j, i, str(val), ha="center", va="center",
                        color="white" if val > 12 else "black", fontsize=9)
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)

    ax = axes[2]
    x = np.arange(3)
    w = 0.36
    d_per = [driven["per"][k] for k in range(3)]
    f_per = [field["per"][k] for k in range(3)]
    ax.bar(x - w / 2, d_per, w, color="#0072B2", label="driven 101")
    ax.bar(x + w / 2, f_per, w, color="#4D4D4D", label="field-only")
    ax.axhline(0.5, color="#333333", ls=":", lw=1)
    ax.set_xticks(x, names)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("one-vs-rest AUC")
    ax.set_title("Per-class (test)")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    return fig


def fig_test_trace(y, u, driven, field):
    n = np.arange(cw.WASHOUT + cw.TRAIN, cw.NUM_WINDOWS)
    y_test = y[cw.WASHOUT + cw.TRAIN :]
    u_test = u[cw.WASHOUT + cw.TRAIN :]
    fig, axes = plt.subplots(
        2, 1, figsize=(7.2, 4.0), sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0]},
    )
    ax = axes[0]
    for start in range(0, 50, 5):
        cls = int(y_test[start])
        ax.axvspan(n[start] - 0.5, n[start] + 4.5, color=CLASS_COLOR[cls], alpha=0.18, zorder=0)
    ax.plot(n, u_test, color="#222222", lw=1.3)
    ax.set_ylim(0, 0.55)
    ax.set_ylabel(r"AHL $u[n]$")
    ax.set_title("Test windows 150–199")
    ax.legend(
        handles=[Patch(facecolor=CLASS_COLOR[k], alpha=0.35, label=cw.CLASS_NAMES[k]) for k in range(3)],
        frameon=False, ncol=3, loc="upper right",
    )

    ax = axes[1]
    jitter = 0.08
    ax.plot(n, y_test, color="#111111", lw=1.4, drawstyle="steps-mid", label="true")
    ax.plot(n, driven["pred"] - jitter, "o", color="#0072B2", ms=4.2, label="driven")
    ax.plot(n, field["pred"] + jitter, "s", color="#4D4D4D", ms=3.8, label="field-only")
    ax.set_yticks([0, 1, 2], list(cw.CLASS_NAMES))
    ax.set_ylim(-0.45, 2.45)
    ax.set_xlim(149.5, 199.5)
    ax.set_xlabel("window $n$")
    ax.set_ylabel("class")
    ax.legend(frameon=False, ncol=3, loc="upper right")
    fig.tight_layout()
    return fig


def main():
    style()
    y, _block, u_csv, _rows = cw.load_labels(HERE / "waveform_labels.csv")
    u = np.array(cw.load_sequence(HERE / "input_ahl_waveform200.txt"), dtype=float)
    data = load_all(y)

    driven_auc = tuple(m["macro"] for m in data["driven"])
    field_auc = tuple(m["field"]["macro"] for m in data["driven"])
    brown_auc = tuple(m["macro"] for m in data["brownian"])
    silent_auc = tuple(m["macro"] for m in data["silent"])
    assert_close("driven AUC", driven_auc, FROZEN_AUC["driven"])
    assert_close("field AUC", field_auc, FROZEN_AUC["field"])
    assert_close("brownian AUC", brown_auc, FROZEN_AUC["brownian"])
    assert_close("silent AUC", silent_auc, FROZEN_AUC["silent"])
    assert_close("driven acc", tuple(m["acc"] for m in data["driven"]), FROZEN_ACC["driven"])
    assert_close("field acc", tuple(m["field"]["acc"] for m in data["driven"]), FROZEN_ACC["field"])

    OUT.mkdir(parents=True, exist_ok=True)
    figs = {
        "waveform_input": fig_input(y, u),
        "waveform_controls": fig_controls(),
        "waveform_confusion": fig_confusion(data["driven"][0], data["driven"][0]["field"]),
        "waveform_test_trace": fig_test_trace(y, u, data["driven"][0], data["driven"][0]["field"]),
    }
    for name, fig in figs.items():
        fig.savefig(OUT / f"{name}.pdf")
        fig.savefig(OUT / f"{name}.png")
        plt.close(fig)

    pred_path = OUT / "waveform_test_predictions.csv"
    n = np.arange(cw.WASHOUT + cw.TRAIN, cw.NUM_WINDOWS)
    y_test = y[cw.WASHOUT + cw.TRAIN :]
    with pred_path.open("w", encoding="utf-8") as handle:
        handle.write("n;y;u;pred_driven_101;pred_field_only\n")
        for i, window in enumerate(n):
            handle.write(
                f"{window};{int(y_test[i])};{u[window]:.12f};"
                f"{int(data['driven'][0]['pred'][i])};"
                f"{int(data['driven'][0]['field']['pred'][i])}\n"
            )

    caption = OUT / "waveform_captions.txt"
    caption.write_text(
        "Figure 1. Frozen sine / square / triangle templates (top) and the "
        "200-window AHL drive u[n] coloured by class (bottom). Shaded bands "
        "are washout (0-39) and test (150-199). u is clipped to [0, 0.5]; "
        "acid is held at 0.5. Stage 6 dish, one AHL AC.\n\n"
        "Figure 2. Test macro one-vs-rest AUC (left; chance 0.5) and argmax "
        "accuracy (right; chance 1/3). Bars are mean +/- s.e. over seeds "
        "101/202/303. Driven biology 0.822 +/- 0.013 vs Brownian 0.587 +/- "
        "0.036 vs silent 0.500 vs AHL field-only 0.860. Field-only was a "
        "predeclared gate and is higher; biology still sits well above "
        "Brownian and chance. Accuracy is diagnostic, not a gate.\n\n"
        "Figure 3. Test-window confusion (rows true, columns predicted) for "
        "driven seed 101 and the AHL field-only readout, plus per-class "
        "one-vs-rest AUC. Square and triangle are easy for both; sine is "
        "the overlapping class (driven sine AUC 0.48 on seed 101).\n\n"
        "Figure 4. Test windows 150-199. Background tint is the true class "
        "block. Black: u[n]. Blue circles: driven argmax. Grey squares: "
        "field-only argmax (offset for visibility).\n",
        encoding="utf-8",
    )
    print("Wrote", OUT)
    for name in figs:
        print(" ", name + ".pdf/.png")
    print(" AUC check matched GATE_EVIDENCE.md")


if __name__ == "__main__":
    main()
