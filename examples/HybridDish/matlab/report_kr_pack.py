#!/usr/bin/env python3
"""Wringe pack 4: kernel-rank diagnostic on frozen Narma10b voxels.

Does not run BSim. No lambda, no target y. All 200 windows.
Primary families are 20×10 maps with S = N = 200. Official 408-D is a
ceiling, not a KR claim. MATLAB peers reproduce driven seed 222 with
kernel_rank.m / report_kr_pack.m. Does not edit any GATE_EVIDENCE.md
Overall line. Packs 1–3 stay closed. Stage 7 SVD stays FAIL.
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
sys.path.insert(0, str(NARMA_DIR))
import check_narma10b as c  # noqa: E402

SEEDS = (111, 222, 333)
THETA = (0.01, 0.05)
OUT_DIR = NARMA_DIR / "results"
FAMILIES = {
    "lum": ("Lum_Mean_",),
    "receiver": ("Receiver_R_",),
    "field": ("AHL_uM_",),
    "den": ("Den_",),
}


def kernel_rank(X, name=""):
    """Roy–Vetterli R_eff + frozen 1%/5% integer ranks. X is S×N."""
    X = np.asarray(X, dtype=float)
    n_win, n_feat = X.shape
    std = X.std(axis=0)
    n_zero = int(np.sum(std < 1e-12))
    std = np.where(std < 1e-12, 1.0, std)
    Z = (X - X.mean(axis=0)) / std
    sigma = np.linalg.svd(Z, compute_uv=False)
    sigma = np.clip(sigma, 0.0, None)
    total = float(np.sum(sigma))
    if total > 0:
        pr = sigma / total
        entropy = float(np.sum(pr[pr > 0] * np.log(pr[pr > 0])))
        r_eff = float(math.exp(-entropy))
    else:
        pr = np.zeros_like(sigma)
        r_eff = 0.0
    sig_max = float(sigma[0]) if len(sigma) else 0.0
    rank_theta = {}
    for theta in THETA:
        rank_theta[theta] = int(np.sum(sigma > theta * sig_max)) if sig_max > 0 else 0
    result = {
        "name": name,
        "S": int(n_win),
        "N": int(n_feat),
        "ceiling": int(min(n_win, n_feat)),
        "kr_claim": bool(n_win >= n_feat),
        "r_eff": r_eff,
        "rank_1pct": rank_theta[0.01],
        "rank_5pct": rank_theta[0.05],
        "n_sigma": int(len(sigma)),
        "n_zero_var": n_zero,
        "sigma_head": [float(v) for v in sigma[:8]],
    }
    claim = (
        "KR diagnostic (S >= N)"
        if result["kr_claim"]
        else f"NOT a KR claim (S={n_win} < N={n_feat}; ceiling {result['ceiling']})"
    )
    print(
        f"  kernel_rank {name}  S={n_win} N={n_feat}  R_eff={r_eff:.3f}  "
        f"rank@1%={result['rank_1pct']}  rank@5%={result['rank_5pct']}"
    )
    print(f"    {claim}")
    return result


def load_family(folder, prefixes, last_prefixes=()):
    return c.read_window_matrix(Path(folder) / "voxels.csv", prefixes, last_prefixes)


def mean_se(values):
    array = np.array(values, dtype=float)
    mean = float(np.mean(array))
    se = float(np.std(array, ddof=1) / np.sqrt(len(array))) if len(array) > 1 else 0.0
    return mean, se


def analyse_folder(folder, arm):
    folder = Path(folder)
    print(f"=== KR {arm} {folder.name} ===")
    out = {"folder": str(folder), "arm": arm}
    if arm == "brownian":
        den = load_family(folder, FAMILIES["den"])
        out["den"] = kernel_rank(den["X"], "Den_*")
        return out
    lum = load_family(folder, FAMILIES["lum"])
    rec = load_family(folder, FAMILIES["receiver"])
    bio = load_family(folder, c.BIOLOGY_MEAN_PREFIXES, c.BIOLOGY_LAST_PREFIXES)
    out["lum"] = kernel_rank(lum["X"], "Lum_Mean_*")
    out["receiver"] = kernel_rank(rec["X"], "Receiver_R_*")
    out["biology408"] = kernel_rank(bio["X"], "official 408")
    if arm == "driven":
        field = load_family(folder, FAMILIES["field"])
        out["field"] = kernel_rank(field["X"], "AHL_uM_*")
    return out


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
    # 408-D is censored from this figure (S < N).
    series = [
        ("Lum driven", pack["r_eff"]["lum"]["driven"], "#2a6f4e"),
        ("Lum silent", pack["r_eff"]["lum"]["silent"], "#b45309"),
        ("R driven", pack["r_eff"]["receiver"]["driven"], "#3f6f5a"),
        ("R silent", pack["r_eff"]["receiver"]["silent"], "#d97706"),
        ("Field driven", pack["r_eff"]["field"]["driven"], "#6b7280"),
        ("Den Brownian", pack["r_eff"]["den"]["brownian"], "#9ca3af"),
    ]
    labels = [row[0] for row in series]
    means = [float(np.mean(row[1])) for row in series]
    ses = [mean_se(row[1])[1] for row in series]
    colors = [row[2] for row in series]
    fig, ax = plt.subplots(figsize=(7.6, 3.9))
    xs = np.arange(len(labels))
    bars = ax.bar(xs, means, color=colors, edgecolor="#222", linewidth=0.4, width=0.72)
    ax.errorbar(xs, means, yerr=ses, fmt="none", color="#222", capsize=3,
                elinewidth=0.9, lw=0, zorder=3)
    ax.set_xticks(xs, labels, rotation=20, ha="right")
    ax.set_ylabel("Effective rank $R_{\\mathrm{eff}}$")
    ax.set_title("Kernel-rank diagnostic (20×10 maps, S = N = 200)")
    ymax = max(means) * 1.22 if means else 1.0
    ax.set_ylim(0, ymax)
    for bar, val in zip(bars, means):
        ax.text(
            bar.get_x() + bar.get_width() / 2, val + 0.03 * ymax,
            f"{val:.2f}", ha="center", va="bottom", fontsize=8,
        )
    ax.text(
        0.99, 0.97,
        "408-D censored: S=200 < N=408, rank ≤ 200. Not a KR claim.",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#333",
    )
    fig.tight_layout()
    fig.savefig(base / "fig_wringe_kr_effective_rank.png", dpi=160)
    plt.close()


def write_markdown(pack, path):
    lum_d = pack["r_eff"]["lum"]["driven"]
    lum_s = pack["r_eff"]["lum"]["silent"]
    lum_dm, lum_dse = mean_se(lum_d)
    lum_sm, lum_sse = mean_se(lum_s)
    rec_dm, _ = mean_se(pack["r_eff"]["receiver"]["driven"])
    field_m, _ = mean_se(pack["r_eff"]["field"]["driven"])
    den_m, _ = mean_se(pack["r_eff"]["den"]["brownian"])
    silent_zero = all(v == 0.0 for v in lum_s)
    n_zero = pack["table"]["lum"]["silent"][0]["n_zero_var"]
    if silent_zero:
        lum_note = (
            f"Silent Lum is the **zero map**: {n_zero}/{pack['table']['lum']['silent'][0]['N']} "
            "columns have zero variance (sources off; L ≡ 0 on every window). "
            f"Driven Lum $R_\\mathrm{{eff}}$ is {lum_dm:.1f}, not ≈ silent. "
            "On this dish the NARMA drive is what puts spatial structure into L. "
            "That is **not** Stage 7’s finding (silent 404 had 39 significant "
            "components on a different four-source matrix). Stage 7 FAIL is not "
            "reopened. High $R_\\mathrm{eff}$ is also not a computing claim: "
            f"Brownian Den is {den_m:.0f} and still NRMSE ≈ 1.16 on NARMA-10."
        )
    else:
        overlap = not (
            lum_dm + lum_dse < lum_sm - lum_sse or lum_sm + lum_sse < lum_dm - lum_dse
        )
        rel = abs(lum_dm - lum_sm) / max(lum_dm, lum_sm, 1e-12)
        if overlap or rel < 0.15:
            lum_note = (
                f"Silent Lum $R_\\mathrm{{eff}}$ {lum_sm:.2f} ≈ driven {lum_dm:.2f}. "
                "Rank is intrinsic population/field structure, not “the NARMA drive "
                "unfolded a rich kernel.” Stage 7 already found silent biology spanning "
                "the SVD rule (silent 404 had 39 significant components on a different "
                "matrix). That FAIL is not reopened."
            )
        else:
            lum_note = (
                f"Silent Lum $R_\\mathrm{{eff}}$ {lum_sm:.2f} vs driven {lum_dm:.2f} "
                f"(relative gap {rel:.0%}). Not a CHARC score, and not a rewrite of "
                "NARMA-10 NRMSE vs field."
            )

    def row(label, key, arm):
        vals = pack["table"][key][arm]
        r = [v["r_eff"] for v in vals]
        r1 = [v["rank_1pct"] for v in vals]
        r5 = [v["rank_5pct"] for v in vals]
        m, se = mean_se(r)
        return (
            f"| {label} | {vals[0]['N']} | {vals[0]['S']} | "
            f"{r[0]:.2f} / {r[1]:.2f} / {r[2]:.2f} | {m:.2f} ± {se:.2f} | "
            f"{r1[0]}/{r1[1]}/{r1[2]} | {r5[0]}/{r5[1]}/{r5[2]} |"
        )

    lines = [
        "# Wringe pack 4 — kernel-rank diagnostic",
        "",
        "Analysis only. Frozen Narma10b voxels. No BSim rerun. No kinetics",
        "retune. Packs 1–3 stay closed. No `GATE_EVIDENCE.md` Overall line was",
        "edited: Stage 6 / Narma10b NARMA-10 stay **PASS**, waveform stays",
        "**FAIL**, BenchA Mackey–Glass stays **PASS**, Lorenz stays **FAIL**,",
        "Stage 7 four-source SVD stays **FAIL** and is not reopened.",
        "",
        "Wringe et al. 2024 ([arXiv:2405.06561](https://arxiv.org/abs/2405.06561))",
        "§8.2 rank diagnostic, Dale 2019 adaptation, algorithm 2. **Not CHARC.**",
        "Not a gate. Dale/Wringe require $S \\ge N$. Official 408 features with",
        "200 windows cannot support a kernel-rank claim ($\\mathrm{rank} \\le",
        "\\min(S,N)=200$). That row is a ceiling, not a headline.",
        "",
        "NARMA $u \\sim$ Uniform[0, 0.5] is the closest existing drive to Dale’s",
        "random input. MG and waveform CSVs are not in this table (smooth chaos",
        "/ three templates). GR is not run (needs $U[-0.1,0.1]$, which we do",
        "not have). No new BSim with $U[-1,1]$.",
        "",
        "## What this is (and is not)",
        "",
        "State matrix $X$ is $S \\times N$ (rows = all 200 windows, columns =",
        "window-mean channels). Columns z-scored with this matrix’s mean /",
        "population std; zero-variance → 0. SVD. **No lambda, no target $y$,**",
        "no ridge split. Washout is not dropped: $S = 200$.",
        "",
        "Primary number: Roy & Vetterli effective rank",
        "$p_i = \\sigma_i / \\sum\\sigma$, $R_\\mathrm{eff} = \\exp(-\\sum p_i \\ln p_i)$.",
        "No cutoff for $R_\\mathrm{eff}$. Integer ranks at frozen $\\theta = 1\\%$",
        "and $5\\%$ of $\\sigma_\\max$ ($\\sigma_i > \\theta\\,\\sigma_\\max$) are",
        "sensitivity, not a θ hunt.",
        "",
        "The computing claim remains pack 1: NARMA-10 driven NRMSE vs field.",
        "KR does not replace it.",
        "",
        "## Effective rank (quote this)",
        "",
        "| Family | N | S | $R_\\mathrm{eff}$ 111 / 222 / 333 | mean ± s.e. | rank 1% | rank 5% |",
        "|---|---|---|---|---|---|---|",
        row("Lum_Mean_* driven (primary)", "lum", "driven"),
        row("Lum_Mean_* silent", "lum", "silent"),
        row("Receiver_R_* driven", "receiver", "driven"),
        row("Receiver_R_* silent", "receiver", "silent"),
        row("AHL_uM_* field (driven)", "field", "driven"),
        row("Den_* Brownian", "den", "brownian"),
        "",
        lum_note,
        "",
        f"Field AHL $R_\\mathrm{{eff}}$ {field_m:.2f} (rank@1% = rank@5% = 2): the "
        "plume is almost a rank-1 delay line. Receiver_R "
        f"{rec_dm:.1f} sits between field and Lum {lum_dm:.1f}, as expected for "
        "a fast Hill on that plume. rank@5% collapses Lum/R/field to 2; quote "
        "$R_\\mathrm{eff}$, which keeps the tail. Integer ranks are sensitivity.",
        "",
        f"Brownian Den $R_\\mathrm{{eff}}$ {den_m:.0f} is the highest bar. Spatial "
        "density noise fills the 20×10 map. High rank ≠ NARMA skill (pack 1 "
        "Brownian test NRMSE ≈ 1.16).",
        "",
        "![Effective rank by family; 408-D censored](fig_wringe_kr_effective_rank.png)",
        "",
        "## Official 408 — ceiling, not a KR claim",
        "",
        (
            "Driven 408-D: $S=200$, $N=408$, $\\min(S,N)=200$. "
            "$R_\\mathrm{eff}$ "
            f"{pack['table']['biology408']['driven'][0]['r_eff']:.2f} / "
            f"{pack['table']['biology408']['driven'][1]['r_eff']:.2f} / "
            f"{pack['table']['biology408']['driven'][2]['r_eff']:.2f}; "
            "rank@1% "
            f"{pack['table']['biology408']['driven'][0]['rank_1pct']}/"
            f"{pack['table']['biology408']['driven'][1]['rank_1pct']}/"
            f"{pack['table']['biology408']['driven'][2]['rank_1pct']}. "
            "Do not headline this as kernel rank. Silent 408 is the same ceiling."
        ),
        "",
        (
            "Silent 408 $R_\\mathrm{eff}$ "
            f"{pack['table']['biology408']['silent'][0]['r_eff']:.2f} / "
            f"{pack['table']['biology408']['silent'][1]['r_eff']:.2f} / "
            f"{pack['table']['biology408']['silent'][2]['r_eff']:.2f}."
        ),
        "",
        "## MATLAB reproduction (driven seed 222)",
        "",
        "From `examples/HybridDish`:",
        "",
        "```matlab",
        "addpath('matlab')",
        "root = fullfile('..', 'BSimReservoirPlanNarma10b');",
        "T = load_voxels(fullfile(root, 'results', 'narma10b_driven_seed222'));",
        "[Xl, ~] = window_features(T, 'lum');",
        "Kl = kernel_rank(Xl, 'Name', 'Lum_Mean_*');",
        "[Xr, ~] = window_features(T, 'receiver');",
        "Kr = kernel_rank(Xr, 'Name', 'Receiver_R_*');",
        "[Xf, ~] = window_features(T, 'field');",
        "Kf = kernel_rank(Xf, 'Name', 'AHL_uM_*');",
        "[Xb, ~] = window_features(T, 'biology');",
        "Kb = kernel_rank(Xb, 'Name', 'official 408');  % ceiling, not a claim",
        "```",
        "",
        "Or `R = report_kr_pack('Only', \"driven_seed222\")`.",
        "",
        f"Expected seed 222: Lum $R_\\mathrm{{eff}}$ "
        f"{pack['table']['lum']['driven'][1]['r_eff']:.3f}, "
        f"rank@1% {pack['table']['lum']['driven'][1]['rank_1pct']}, "
        f"rank@5% {pack['table']['lum']['driven'][1]['rank_5pct']}. "
        f"Field $R_\\mathrm{{eff}}$ "
        f"{pack['table']['field']['driven'][1]['r_eff']:.3f}. "
        f"408-D prints `NOT a KR claim`.",
        "",
        "All 200 windows. Population std. Frozen θ = 0.01 and 0.05 of σ_max.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    table = {
        "lum": {"driven": [], "silent": []},
        "receiver": {"driven": [], "silent": []},
        "field": {"driven": []},
        "den": {"brownian": []},
        "biology408": {"driven": [], "silent": []},
    }
    runs = []
    for arm in ("driven", "brownian", "silent"):
        for seed in SEEDS:
            folder = OUT_DIR / f"narma10b_{arm}_seed{seed}"
            result = analyse_folder(folder, arm)
            runs.append(result)
            if arm == "brownian":
                table["den"]["brownian"].append(result["den"])
            else:
                table["lum"][arm].append(result["lum"])
                table["receiver"][arm].append(result["receiver"])
                table["biology408"][arm].append(result["biology408"])
                if arm == "driven":
                    table["field"]["driven"].append(result["field"])

    for key, arm, n_feat in (
        ("lum", "driven", 200),
        ("lum", "silent", 200),
        ("receiver", "driven", 200),
        ("receiver", "silent", 200),
        ("field", "driven", 200),
        ("den", "brownian", 200),
    ):
        for item in table[key][arm]:
            if item["S"] != 200 or item["N"] != n_feat:
                raise SystemExit(
                    f"{key}/{arm} has S={item['S']} N={item['N']}, expected 200×{n_feat}"
                )
            if not item["kr_claim"]:
                raise SystemExit(f"{key}/{arm} unexpectedly S < N")
    for item in table["biology408"]["driven"] + table["biology408"]["silent"]:
        if item["N"] != 408 or item["S"] != 200 or item["kr_claim"]:
            raise SystemExit(f"408 row should be S=200 N=408 kr_claim=false, got {item}")

    r_eff = {
        key: {arm: [v["r_eff"] for v in runs_arm] for arm, runs_arm in arms.items()}
        for key, arms in table.items()
    }
    pack = {"table": table, "r_eff": r_eff, "theta": list(THETA), "seeds": list(SEEDS)}

    def json_default(obj):
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(type(obj))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "wringe_kr_pack.json"
    json_path.write_text(
        json.dumps({"runs": runs, **pack}, indent=2, default=json_default) + "\n",
        encoding="utf-8",
    )
    write_figure(pack, OUT_DIR)
    md_path = OUT_DIR / "WRINGE_KR_PACK.md"
    write_markdown(pack, md_path)
    print(f"markdown={md_path}")
    print(f"json={json_path}")
    print(f"figure={OUT_DIR / 'fig_wringe_kr_effective_rank.png'}")
    for key, arms in r_eff.items():
        for arm, vals in arms.items():
            m, se = mean_se(vals)
            print(f"{key:12} {arm:10} {vals[0]:.3f} / {vals[1]:.3f} / {vals[2]:.3f}  mean {m:.3f} ± {se:.3f}")
    print("No GATE_EVIDENCE.md Overall line was edited.")


if __name__ == "__main__":
    main()
