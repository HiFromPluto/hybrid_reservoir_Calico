#!/usr/bin/env python3
"""Wringe pack 2 batch helper: waveform F1 / confusion on frozen voxels.

Does not run BSim. Closed val-slice; one lambda for three OVR heads,
picked by highest macro OVR validation AUC. MATLAB peers can reproduce
one driven folder with classify_waveform.m / report_waveform_pack.m.
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
WAVE_DIR = EXAMPLES / "BSimReservoirPlanWaveform"
sys.path.insert(0, str(WAVE_DIR))
import check_waveform as cw  # noqa: E402

SEEDS = (101, 202, 303)
CLASS_NAMES = ("sine", "square", "triangle")
MATCH_TOL = 1e-3
OUT_DIR = WAVE_DIR / "results"

GATE_AUC = {
    "driven": (0.7969047619047619, 0.8287301587301587, 0.8411111111111111),
    "field": (0.8603174603174604, 0.8603174603174604, 0.8603174603174604),
    "brownian": (0.5743650793650793, 0.653968253968254, 0.5317460317460317),
    "silent": (0.5, 0.5, 0.5),
}
GATE_ACC = {
    "driven": (0.76, 0.76, 0.68),
    "field": (0.86, 0.86, 0.86),
    "brownian": (0.32, 0.30, 0.32),
    "silent": (0.30, 0.30, 0.30),
}
GATE_CONFUSION = {
    "driven": (
        [[6, 1, 8], [1, 13, 1], [1, 0, 19]],
        [[5, 5, 5], [0, 13, 2], [0, 0, 20]],
        [[9, 1, 5], [5, 8, 2], [2, 1, 17]],
    ),
    "field": (
        [[8, 7, 0], [0, 15, 0], [0, 0, 20]],
        [[8, 7, 0], [0, 15, 0], [0, 0, 20]],
        [[8, 7, 0], [0, 15, 0], [0, 0, 20]],
    ),
}


def prf_from_confusion(C):
    C = np.asarray(C, dtype=float)
    n = C.shape[0]
    precision = np.zeros(n)
    recall = np.zeros(n)
    f1 = np.zeros(n)
    support = C.sum(axis=1)
    for k in range(n):
        tp = C[k, k]
        fp = C[:, k].sum() - tp
        fn = C[k, :].sum() - tp
        precision[k] = tp / (tp + fp) if tp + fp > 0 else 0.0
        recall[k] = tp / (tp + fn) if tp + fn > 0 else 0.0
        f1[k] = (
            2 * precision[k] * recall[k] / (precision[k] + recall[k])
            if precision[k] + recall[k] > 0
            else 0.0
        )
    return {
        "precision": precision.tolist(),
        "recall": recall.tolist(),
        "f1": f1.tolist(),
        "support": [int(v) for v in support],
        "macro_f1": float(np.mean(f1)),
    }


def evaluate_folder(folder, arm, y):
    run = cw.read_run(folder, arm if arm != "field" else "driven", cw.NUM_WINDOWS, cw.EXPECTED_AUX_ROWS)
    if arm == "field":
        ev = cw.evaluate_readout(run["field"]["X"], y, run["field"]["windows"])
        n_features = len(run["field"]["feature_names"])
    else:
        ev = cw.evaluate_arm_run(run, y)["classification"]
        n_features = len(run["matrix"]["feature_names"])
    prf = prf_from_confusion(ev["test_confusion"])
    return {
        "folder": str(folder),
        "arm": arm,
        "n_features": n_features,
        "lambda": ev["lambda"],
        "macro_auc": ev["test_macro_ovr_auc"],
        "per_class_auc": ev["test_per_class_ovr_auc"],
        "accuracy": ev["test_accuracy_argmax"],
        "confusion": ev["test_confusion"],
        **prf,
    }


def mean_se(values):
    array = np.array(values, dtype=float)
    mean = float(np.mean(array))
    se = float(np.std(array, ddof=1) / np.sqrt(len(array))) if len(array) > 1 else 0.0
    return mean, se


def check_gate(computed):
    mismatches = []
    for arm in ("driven", "field", "brownian", "silent"):
        for i, seed in enumerate(SEEDS):
            if abs(computed[arm][i]["macro_auc"] - GATE_AUC[arm][i]) > MATCH_TOL:
                mismatches.append(
                    f"{arm} seed{seed} AUC {computed[arm][i]['macro_auc']:.6f} "
                    f"!= {GATE_AUC[arm][i]:.6f}"
                )
            if abs(computed[arm][i]["accuracy"] - GATE_ACC[arm][i]) > MATCH_TOL:
                mismatches.append(
                    f"{arm} seed{seed} acc {computed[arm][i]['accuracy']:.6f} "
                    f"!= {GATE_ACC[arm][i]:.6f}"
                )
        if arm in GATE_CONFUSION:
            for i, seed in enumerate(SEEDS):
                got = computed[arm][i]["confusion"]
                exp = GATE_CONFUSION[arm][i]
                if got != exp:
                    mismatches.append(f"{arm} seed{seed} confusion {got} != {exp}")
    return mismatches


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
    arms = ("driven", "field", "brownian", "silent")
    labels = ("Driven", "Field-only", "Brownian", "Silent")
    colors = ("#2a6f4e", "#6b7280", "#9ca3af", "#b45309")

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    auc_means = [pack["auc_mean"][a] for a in arms]
    acc_means = [pack["acc_mean"][a] for a in arms]
    axes[0].bar(labels, auc_means, color=colors)
    axes[0].axhline(0.5, color="#444", linestyle=":", linewidth=1, label="chance 0.5")
    axes[0].set_ylabel("Macro OVR AUC")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Test macro OVR AUC")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].bar(labels, acc_means, color=colors)
    axes[1].axhline(1.0 / 3.0, color="#444", linestyle=":", linewidth=1, label="chance 1/3")
    axes[1].set_ylabel("Argmax accuracy")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Test accuracy (not the gate)")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_waveform_auc_acc.png", dpi=160)
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8))
    for ax, arm, title in (
        (axes[0], "driven", "Driven biology, seed 101"),
        (axes[1], "field", "Field-only AHL (every seed)"),
    ):
        C = np.asarray(pack["runs"][arm][0]["confusion"], dtype=float)
        im = ax.imshow(C, cmap="Blues", vmin=0, vmax=20)
        ax.set_xticks([0, 1, 2], CLASS_NAMES)
        ax.set_yticks([0, 1, 2], CLASS_NAMES)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_title(title)
        ax.grid(False)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{int(C[i, j])}", ha="center", va="center",
                        color="white" if C[i, j] > 10 else "#111")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_waveform_confusion.png", dpi=160)
    plt.close()

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    x = np.arange(3)
    width = 0.35
    driven_f1 = pack["f1_mean"]["driven"]
    field_f1 = pack["f1_mean"]["field"]
    ax.bar(x - width / 2, driven_f1, width, color="#2a6f4e", label="Driven")
    ax.bar(x + width / 2, field_f1, width, color="#6b7280", label="Field-only")
    ax.set_xticks(x, CLASS_NAMES)
    ax.set_ylabel("F1")
    ax.set_ylim(0, 1.05)
    ax.set_title("Per-class F1 (mean of 3 seeds). Sine is the hard class.")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_waveform_f1.png", dpi=160)
    plt.close()


def write_markdown(pack, path):
    auc_m = pack["auc_mean"]
    auc_se = pack["auc_se"]
    acc_m = pack["acc_mean"]
    lines = [
        "# Wringe pack 2 — waveform F1 / confusion",
        "",
        "Analysis only. Frozen waveform voxels, closed val-slice, one-vs-rest",
        "ridge with lambda picked by highest macro OVR validation AUC. No BSim",
        "rerun. Waveform `GATE_EVIDENCE.md` Overall stays **FAIL** (gate 2:",
        "field-only AUC > driven). That line was not edited. Stage 6 NARMA-10",
        "stays **PASS**. Narma10b stays **PASS**. Wringe pack 1 is closed.",
        "",
        "Wringe et al. 2024 ([arXiv:2405.06561](https://arxiv.org/abs/2405.06561))",
        "§9.3.2 metrics on this dish: accuracy vs 1/3, macro one-vs-rest AUC vs",
        "0.5, confusion, per-class precision / recall / F1. This is not Track B.",
        "Track B was binary and failed.",
        "",
        "## Overall remains FAIL",
        "",
        "Field-only still slightly better on the primary metric: macro OVR AUC",
        f"**{auc_m['field']:.4f}** vs driven **{auc_m['driven']:.4f} ± {auc_se['driven']:.4f}**.",
        "The AHL plume already ranks the waveform. Accuracy vs 1/3 is real and",
        "is printed here; it is **not** the gate and does not override FAIL.",
        "",
        "Allowed claim: one-way AHL → Hill R → L sits well above 1/3 and",
        "Brownian, a little below the plume. NARMA-10 remains the only task",
        "where biology beat field-only.",
        "",
        "## Test macro OVR AUC and accuracy (must match GATE_EVIDENCE)",
        "",
        "Chance macro OVR AUC = 0.5. Chance accuracy = 1/3. Not Track B’s 0.5",
        "accuracy floor.",
        "",
        "| Arm | AUC 101 | 202 | 303 | mean ± s.e. | Acc 101 / 202 / 303 |",
        "|---|---|---|---|---|---|",
    ]
    for name, key in (
        ("Driven biology", "driven"),
        ("Field-only AHL", "field"),
        ("Brownian", "brownian"),
        ("Silent", "silent"),
    ):
        runs = pack["runs"][key]
        aucs = [r["macro_auc"] for r in runs]
        accs = [r["accuracy"] for r in runs]
        lines.append(
            f"| {name} | {aucs[0]:.4f} | {aucs[1]:.4f} | {aucs[2]:.4f} | "
            f"{auc_m[key]:.4f} ± {auc_se[key]:.4f} | "
            f"{accs[0]:.2f} / {accs[1]:.2f} / {accs[2]:.2f} |"
        )
    driven_acc_mean = acc_m["driven"]
    lines.extend([
        "",
        f"Driven accuracy mean **{driven_acc_mean:.2f}** vs 1/3. Field-only",
        "accuracy is 0.86 every seed. Brownian and silent sit at chance",
        "accuracy (~0.30–0.32). The gate is still AUC, and field-only wins it.",
        "",
        "![AUC and accuracy vs chance](fig_wringe_waveform_auc_acc.png)",
        "",
        "## Confusion (50 test windows)",
        "",
        "Rows = true (sine, square, triangle). Columns = predicted. Driven 101",
        "and field (every seed) match `GATE_EVIDENCE.md`.",
        "",
        "| Arm | seed | confusion |",
        "|---|---|---|",
    ])
    for arm in ("driven", "field", "brownian", "silent"):
        for seed, run in zip(SEEDS, pack["runs"][arm]):
            lines.append(f"| {arm} | {seed} | `{run['confusion']}` |")
    lines.extend([
        "",
        "![Confusion driven vs field](fig_wringe_waveform_confusion.png)",
        "",
        "## Per-class precision / recall / F1",
        "",
        "Zero-division → 0. Macro-F1 is the unweighted mean of the three class",
        "F1s (do not weight by 15/15/20). Support on test: sine 15, square 15,",
        "triangle 20.",
        "",
        "Sine is the hard class. Driven sine OVR AUC 0.48 / 0.58 / 0.67; field",
        "sine 0.61. Square and triangle are easy. F1 shows the same split.",
        "",
        "### Driven",
        "",
        "| Seed | sine P/R/F1 | square P/R/F1 | triangle P/R/F1 | macro-F1 |",
        "|---|---|---|---|---|",
    ])

    def prf_cell(run, k):
        return (
            f"{run['precision'][k]:.3f}/{run['recall'][k]:.3f}/{run['f1'][k]:.3f}"
        )

    for seed, run in zip(SEEDS, pack["runs"]["driven"]):
        lines.append(
            f"| {seed} | {prf_cell(run, 0)} | {prf_cell(run, 1)} | "
            f"{prf_cell(run, 2)} | {run['macro_f1']:.3f} |"
        )
    dmean = pack["f1_mean"]["driven"]
    lines.extend([
        f"| mean | F1 {dmean[0]:.3f} | F1 {dmean[1]:.3f} | F1 {dmean[2]:.3f} | "
        f"{pack['macro_f1_mean']['driven']:.3f} |",
        "",
        "### Field-only",
        "",
        "| Seed | sine P/R/F1 | square P/R/F1 | triangle P/R/F1 | macro-F1 |",
        "|---|---|---|---|---|",
    ])
    for seed, run in zip(SEEDS, pack["runs"]["field"]):
        lines.append(
            f"| {seed} | {prf_cell(run, 0)} | {prf_cell(run, 1)} | "
            f"{prf_cell(run, 2)} | {run['macro_f1']:.3f} |"
        )
    fmean = pack["f1_mean"]["field"]
    lines.extend([
        f"| mean | F1 {fmean[0]:.3f} | F1 {fmean[1]:.3f} | F1 {fmean[2]:.3f} | "
        f"{pack['macro_f1_mean']['field']:.3f} |",
        "",
        "### Brownian",
        "",
        "| Seed | sine P/R/F1 | square P/R/F1 | triangle P/R/F1 | macro-F1 |",
        "|---|---|---|---|---|",
    ])
    for seed, run in zip(SEEDS, pack["runs"]["brownian"]):
        lines.append(
            f"| {seed} | {prf_cell(run, 0)} | {prf_cell(run, 1)} | "
            f"{prf_cell(run, 2)} | {run['macro_f1']:.3f} |"
        )
    bmean = pack["f1_mean"]["brownian"]
    lines.extend([
        f"| mean | F1 {bmean[0]:.3f} | F1 {bmean[1]:.3f} | F1 {bmean[2]:.3f} | "
        f"{pack['macro_f1_mean']['brownian']:.3f} |",
        "",
        "### Silent",
        "",
        "Sources off. The ridge collapses to an intercept and always predicts",
        "sine (class 0) on every test window. Accuracy 0.30 = 15/50.",
        "",
        "| Seed | sine P/R/F1 | square P/R/F1 | triangle P/R/F1 | macro-F1 |",
        "|---|---|---|---|---|",
    ])
    for seed, run in zip(SEEDS, pack["runs"]["silent"]):
        lines.append(
            f"| {seed} | {prf_cell(run, 0)} | {prf_cell(run, 1)} | "
            f"{prf_cell(run, 2)} | {run['macro_f1']:.3f} |"
        )
    smean = pack["f1_mean"]["silent"]
    lines.extend([
        f"| mean | F1 {smean[0]:.3f} | F1 {smean[1]:.3f} | F1 {smean[2]:.3f} | "
        f"{pack['macro_f1_mean']['silent']:.3f} |",
        "",
        "![Per-class F1 driven vs field](fig_wringe_waveform_f1.png)",
        "",
        "Sine F1 is the weak cell on both arms. Square and triangle F1 stay",
        "high. That is the classification story; it does not change gate 2.",
        "",
        "## Per-class OVR AUC (from GATE_EVIDENCE, recomputed)",
        "",
        "| Arm | seed | sine | square | triangle |",
        "|---|---|---|---|---|",
    ])
    for arm in ("driven", "field", "brownian", "silent"):
        for seed, run in zip(SEEDS, pack["runs"][arm]):
            per = run["per_class_auc"]
            lines.append(
                f"| {arm} | {seed} | {per['sine']:.4f} | {per['square']:.4f} | "
                f"{per['triangle']:.4f} |"
            )
    lines.extend([
        "",
        "## MATLAB reproduction (one driven seed)",
        "",
        "From `examples/HybridDish`:",
        "",
        "```matlab",
        "addpath('matlab')",
        "root = fullfile('..', 'BSimReservoirPlanWaveform');",
        "S = classify_waveform( ...",
        "    fullfile(root, 'results', 'waveform_driven_seed101'), ...",
        "    fullfile(root, 'waveform_labels.csv'), 'Family', 'biology');",
        "Sf = classify_waveform( ...",
        "    fullfile(root, 'results', 'waveform_driven_seed101'), ...",
        "    fullfile(root, 'waveform_labels.csv'), 'Family', 'field');",
        "```",
        "",
        "Or `R = report_waveform_pack('Only', \"driven_seed101\")`.",
        "",
        "Expected seed 101: driven accuracy 0.76, macro AUC 0.7969, confusion",
        "`[[6,1,8],[1,13,1],[1,0,19]]`. Field accuracy 0.86, macro AUC 0.8603,",
        "confusion `[[8,7,0],[0,15,0],[0,0,20]]`.",
        "",
        "Ridge: washout 40 / train 110 / test 50; lambda on post-washout rows",
        "88..109 only; **one** lambda for three OVR heads, highest val macro",
        "AUC, ties take larger lambda; bias unregularised; train pop-std.",
        "Biology ridge is 408 channels. Voxel AHL is field-only, not in biology.",
        "",
        "Silent dirs are Stage 6 copies (sources off). They were not rerun.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    y, _, _, _ = cw.load_labels(WAVE_DIR / "waveform_labels.csv")
    reconstructed = ",".join(str(int(y[p * cw.WINDOWS_PER_BLOCK])) for p in range(40))
    if reconstructed != cw.CLASS_VECTOR:
        raise SystemExit("class vector drifted from PROTOCOL.md")

    computed = {arm: [] for arm in ("driven", "brownian", "silent", "field")}
    for arm in ("driven", "brownian", "silent"):
        for seed in SEEDS:
            folder = OUT_DIR / f"waveform_{arm}_seed{seed}"
            print(f"=== {arm} seed{seed} ===")
            result = evaluate_folder(folder, arm, y)
            computed[arm].append(result)
            print(
                f"  auc={result['macro_auc']:.4f} acc={result['accuracy']:.2f} "
                f"macroF1={result['macro_f1']:.3f} C={result['confusion']}"
            )
            if arm == "driven":
                field = evaluate_folder(folder, "field", y)
                computed["field"].append(field)
                print(
                    f"  field auc={field['macro_auc']:.4f} acc={field['accuracy']:.2f} "
                    f"macroF1={field['macro_f1']:.3f} C={field['confusion']}"
                )

    mismatches = check_gate(computed)
    if mismatches:
        print("MISMATCH vs waveform GATE_EVIDENCE (tol 1e-3):")
        for line in mismatches:
            print(" ", line)
        raise SystemExit("AUC/accuracy/confusion do not match GATE_EVIDENCE; stopping.")
    print("AUC / accuracy / confusion match GATE_EVIDENCE within 1e-3")

    auc_mean, auc_se, acc_mean = {}, {}, {}
    f1_mean, macro_f1_mean = {}, {}
    for arm, runs in computed.items():
        aucs = [r["macro_auc"] for r in runs]
        accs = [r["accuracy"] for r in runs]
        auc_mean[arm], auc_se[arm] = mean_se(aucs)
        acc_mean[arm], _ = mean_se(accs)
        f1s = np.array([r["f1"] for r in runs], dtype=float)
        f1_mean[arm] = f1s.mean(axis=0).tolist()
        macro_f1_mean[arm] = float(np.mean([r["macro_f1"] for r in runs]))

    pack = {
        "runs": computed,
        "auc_mean": auc_mean,
        "auc_se": auc_se,
        "acc_mean": acc_mean,
        "f1_mean": f1_mean,
        "macro_f1_mean": macro_f1_mean,
        "class_vector": reconstructed,
        "class_vector_sha256": cw.CLASS_SHA,
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
    json_path = OUT_DIR / "wringe_waveform_pack.json"
    json_path.write_text(
        json.dumps(pack, indent=2, default=json_default) + "\n", encoding="utf-8"
    )
    write_figures(pack, OUT_DIR)
    md_path = OUT_DIR / "WRINGE_WAVEFORM_PACK.md"
    write_markdown(pack, md_path)
    print(f"markdown={md_path}")
    print(f"json={json_path}")
    print(
        f"driven_auc={auc_mean['driven']:.4f} field_auc={auc_mean['field']:.4f} "
        f"driven_acc={acc_mean['driven']:.2f} field_acc={acc_mean['field']:.2f} "
        f"driven_macroF1={macro_f1_mean['driven']:.3f} "
        f"field_macroF1={macro_f1_mean['field']:.3f}"
    )
    print("Overall line in GATE_EVIDENCE.md was not edited. Waveform stays FAIL.")


if __name__ == "__main__":
    main()
