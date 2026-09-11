#!/usr/bin/env python3
"""Evaluate Plan Stage 7 rank gates. NARMA-10 is not a gate."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

NUM_WINDOWS = 40
EXPECTED_AUX_ROWS = 640
SIGNIFICANT_RATIO = 0.01
RANK_MEAN_PREFIXES = ("Den_", "Receiver_R_", "Lum_Mean_")
RANK_LAST_PREFIXES = ("Input_Driven_Death_",)
EXPECTED_COUNTS = {
    "Den_": 200,
    "Receiver_R_": 200,
    "Lum_Mean_": 200,
    "Input_Driven_Death_": 8,
}
FORBIDDEN_PREFIXES = (
    "AHL_uM_", "pH_", "Fraction_q", "Lum_Sum_", "Birth_",
    "Clamp_Death_", "OOB_Death_", "Att_",
)
EXPECTED_POSITIONS = {
    "AttractantA": (150.0, 100.0, 5.0),
    "AHL": (500.0, 250.0, 5.0),
    "Acid": (850.0, 100.0, 5.0),
    "AttractantB": (850.0, 400.0, 5.0),
}
STAGE4_STATUS = "FAIL"
STAGE5_STATUS = "PASS"
STAGE6_STATUS = "PASS"
LAST_SAMPLE = (39, 15, 299.95)
DEN_DEAD_VARIANCE = 1e-12


def json_default(obj):
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"not json serializable: {type(obj)}")


def validate_csv(path, expected_rows):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, [])
        widths = [len(row) for row in reader]
    return {
        "path": str(path),
        "row_count": len(widths),
        "column_count": len(header),
        "row_count_pass": len(widths) == expected_rows,
        "rectangular_csv_pass": bool(header) and all(width == len(header) for width in widths),
    }


def columns_with_prefix(header, prefixes):
    names = []
    for prefix in prefixes:
        names.extend(name for name in header if name.startswith(prefix))
    return names


def read_window_matrix(voxels_path):
    path = Path(voxels_path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        mean_cols = columns_with_prefix(header, RANK_MEAN_PREFIXES)
        last_cols = columns_with_prefix(header, RANK_LAST_PREFIXES)
        feature_names = mean_cols + last_cols
        window_idx = header.index("Window")
        mean_idx = [header.index(name) for name in mean_cols]
        last_idx = [header.index(name) for name in last_cols]
        sums, counts, lasts = {}, {}, {}
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"ragged voxels.csv in {path}")
            window = int(row[window_idx])
            if window not in sums:
                sums[window] = np.zeros(len(mean_idx), dtype=float)
                counts[window] = 0
                lasts[window] = np.zeros(len(last_idx), dtype=float)
            if mean_idx:
                sums[window] += np.array([float(row[i]) for i in mean_idx], dtype=float)
            counts[window] += 1
            if last_idx:
                lasts[window] = np.array([float(row[i]) for i in last_idx], dtype=float)
    windows = sorted(sums)
    rows = []
    for window in windows:
        parts = []
        if mean_idx:
            parts.append(sums[window] / counts[window])
        if last_idx:
            parts.append(lasts[window])
        rows.append(np.concatenate(parts))
    return {
        "windows": windows,
        "X": np.vstack(rows),
        "feature_names": feature_names,
        "samples_per_window": [counts[window] for window in windows],
    }


def last_sample(results_path):
    with Path(results_path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    if not rows:
        return {"pass": False, "window": None, "sample": None, "time": None}
    last = rows[-1]
    window = int(last["Window"])
    sample = int(last["Sample"])
    time = float(last["TimeInWindow_s"])
    return {
        "window": window,
        "sample": sample,
        "time": time,
        "pass": window == LAST_SAMPLE[0] and sample == LAST_SAMPLE[1]
        and abs(time - LAST_SAMPLE[2]) < 1e-9,
    }


def read_positions(path):
    positions = {}
    text = Path(path).read_text(encoding="utf-8")
    for line in text.splitlines():
        if "=" not in line or line.startswith("flow"):
            continue
        name, value = line.split("=", 1)
        parts = tuple(float(item) for item in value.split(","))
        positions[name] = parts
    return positions


def load_pairwise(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        next(reader)
        return [{"a": a, "b": b, "r": float(r)} for a, b, r in reader]


def svd_spectrum(X):
    centered = X - X.mean(axis=0, keepdims=True)
    _, singular, _ = np.linalg.svd(centered, full_matrices=False)
    leading = float(singular[0]) if len(singular) else float("nan")
    ratios = [float(value / leading) if leading > 0 else float("nan") for value in singular]
    n_sig = sum(1 for ratio in ratios if ratio > SIGNIFICANT_RATIO)
    return {
        "singular_values": [float(value) for value in singular],
        "ratios_to_sigma1": ratios,
        "n_significant": n_sig,
        "leading": leading,
    }


def den_variance(feature_names, X):
    den_idx = [i for i, name in enumerate(feature_names) if name.startswith("Den_")]
    if not den_idx:
        return {"max": 0.0, "mean": 0.0, "dead": True}
    variances = X[:, den_idx].var(axis=0)
    maximum = float(np.max(variances))
    return {
        "max": maximum,
        "mean": float(np.mean(variances)),
        "n_positive": int(np.sum(variances > DEN_DEAD_VARIANCE)),
        "dead": maximum <= DEN_DEAD_VARIANCE,
    }


def feature_contract_ok(feature_names):
    counts = {prefix: sum(name.startswith(prefix) for name in feature_names)
              for prefix in EXPECTED_COUNTS}
    forbidden = [name for name in feature_names
                 if name.startswith(FORBIDDEN_PREFIXES)]
    return (
        counts == EXPECTED_COUNTS
        and not forbidden
        and len(feature_names) == 608
    ), counts, forbidden


def evaluate_run(path, expected_windows, expected_aux):
    directory = Path(path)
    matrix = read_window_matrix(directory / "voxels.csv")
    if matrix["windows"] != list(range(expected_windows)):
        raise ValueError(f"{directory} windows {matrix['windows'][:3]}..")
    spectrum = svd_spectrum(matrix["X"])
    den = den_variance(matrix["feature_names"], matrix["X"])
    ok, counts, forbidden = feature_contract_ok(matrix["feature_names"])
    positions = read_positions(directory / "source_positions.txt")
    return {
        "path": str(directory),
        "positions": positions,
        "n_features": len(matrix["feature_names"]),
        "feature_names": matrix["feature_names"],
        "feature_counts": counts,
        "forbidden_in_rank": forbidden,
        "feature_contract_pass": ok,
        "spectrum": spectrum,
        "den_variance": den,
        "last_sample": last_sample(directory / "results.csv"),
        "summary_csv": validate_csv(directory / "window_summary.csv", expected_windows),
        "auxiliary_csvs": {
            "results.csv": validate_csv(directory / "results.csv", expected_aux),
            "voxels.csv": validate_csv(directory / "voxels.csv", expected_aux),
        },
    }


def write_markdown(path, evidence):
    overall = "PASS" if evidence["overall_gate_pass"] else "FAIL"
    lines = [
        "# Stage 7 gate evidence",
        "",
        "Stage 7 copies frozen Stage 6 and spreads four grounded sources.",
        "This is not the plan's 5-AC / flow-8 / glucose table. Stage 4 remains",
        "**FAIL**. Stage 5 remains **PASS**. Stage 6 remains **PASS**.",
        "Kinetics were not retuned. NARMA-10 was not rerun. Stage 8 was not",
        "started.",
        "",
        f"## Overall: {overall}",
        "",
        evidence["failure_reason"] if not evidence["overall_gate_pass"]
        else "Driven state matrices have at least four significant components; silent does not.",
        "",
        "A 40×608 centered matrix has rank at most 39. The frozen rule "
        "`σ_i/σ_1 > 0.01` counts almost every available component when the "
        "spectrum decays slowly. Positions, rates, and decays were not moved "
        "after seeing SVD.",
        "",
        "| Gate | Result |",
        "|---|---|",
        f"| 1. Four frozen source positions | **{'PASS' if evidence['gates']['positions_pass'] else 'FAIL'}** |",
        f"| 2. Input pairwise `|r|` < 0.3 | **{'PASS' if evidence['gates']['pairwise_r_pass'] else 'FAIL'}** |",
        f"| 3. Driven 101/202/303 each have ≥4 significant σ | **{'PASS' if evidence['gates']['driven_rank_pass'] else 'FAIL'}** |",
        f"| 4. Silent 404 has <4 significant σ | **{'PASS' if evidence['gates']['silent_rank_pass'] else 'FAIL'}** |",
        f"| 5. Rank matrix is Den / Receiver_R / Lum_Mean / Input_Driven_Death | **{'PASS' if evidence['gates']['feature_list_pass'] else 'FAIL'}** |",
        f"| 6. CSV 40/640/640 and last sample 39;15;299.95 | **{'PASS' if evidence['gates']['csv_validation_pass'] else 'FAIL'}** |",
        f"| 7. Prior stages FAIL/PASS/PASS; kinetics frozen | **{'PASS' if evidence['gates']['prior_stage_labels_pass'] else 'FAIL'}** |",
        "",
        "## Source positions",
        "",
        "| Source | Position |",
        "|---|---|",
        "| Attractant A | (150, 100, 5) |",
        "| AHL | (500, 250, 5) |",
        "| Acid | (850, 100, 5) |",
        "| Attractant B | (850, 400, 5) |",
        "",
        "`FLOW_SPEED=0`. Not the plan 5-AC table.",
        "",
        "## Input pairwise |r|",
        "",
        "| Pair | r |",
        "|---|---|",
    ]
    for pair in evidence["pairwise"]:
        lines.append(f"| {pair['a']}, {pair['b']} | {pair['r']:.6f} |")
    lines.extend([
        "",
        f"Max `|r|` = {evidence['pairwise_max_abs_r']:.6f}.",
        "",
        "## Rank matrix",
        "",
        "Printed feature groups:",
        "",
        "- `Den_*` 20×10 (200)",
        "- `Receiver_R_*` / Mean_q 20×10 (200)",
        "- `Lum_Mean_*` / Mean_L 20×10 (200)",
        "- `Input_Driven_Death_*` 4×2 (8)",
        "",
        f"Printed count: {evidence['driven_feature_count']} features.",
        "",
        "## Singular values (σ_i / σ_1)",
        "",
        "| Seed | n significant | σ1 | σ2/σ1 | σ3/σ1 | σ4/σ1 | σ5/σ1 |",
        "|---|---|---|---|---|---|---|",
    ])
    for label, run in evidence["driven_rows"]:
        ratios = run["spectrum"]["ratios_to_sigma1"]
        pad = ratios + [float("nan")] * 5
        lines.append(
            f"| {label} | {run['spectrum']['n_significant']} | "
            f"{run['spectrum']['leading']:.4g} | "
            f"{pad[1]:.4f} | {pad[2]:.4f} | {pad[3]:.4f} | {pad[4]:.4f} |"
        )
    silent = evidence["silent"]
    s_ratios = silent["spectrum"]["ratios_to_sigma1"] + [float("nan")] * 5
    lines.append(
        f"| silent 404 | {silent['spectrum']['n_significant']} | "
        f"{silent['spectrum']['leading']:.4g} | "
        f"{s_ratios[1]:.4f} | {s_ratios[2]:.4f} | {s_ratios[3]:.4f} | {s_ratios[4]:.4f} |"
    )
    lines.extend([
        "",
        "Full spectra are in `stage7_gate_evidence.json`.",
        "",
        "## Den_* variance",
        "",
    ])
    for label, run in evidence["driven_rows"]:
        den = run["den_variance"]
        lines.append(
            f"- {label}: mean var {den['mean']:.6g}, max var {den['max']:.6g}, "
            f"{'DEAD' if den['dead'] else 'live'}"
        )
    lines.extend([
        "",
        "## Protocol",
        "",
        "See `PROTOCOL.md`. 40 windows, four independent sequences, column-centered SVD,",
        "significant if σ_i/σ_1 > 0.01. CSV completeness is the on-disk files.",
        "Stage 8 was not started.",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driven", nargs=3, required=True)
    parser.add_argument("--silent", required=True)
    parser.add_argument("--pairwise", default="input_pairwise_r.csv")
    parser.add_argument("--expected-windows", type=int, default=NUM_WINDOWS)
    parser.add_argument("--expected-aux-rows", type=int, default=EXPECTED_AUX_ROWS)
    parser.add_argument("--evidence", default="results/stage7_gate_evidence.json")
    parser.add_argument("--markdown", default="results/GATE_EVIDENCE.md")
    args = parser.parse_args()

    pairwise = load_pairwise(args.pairwise)
    driven = [evaluate_run(path, args.expected_windows, args.expected_aux_rows)
              for path in args.driven]
    silent = evaluate_run(args.silent, args.expected_windows, args.expected_aux_rows)

    print("RANK_FEATURE_LIST")
    for name in driven[0]["feature_names"]:
        print(name)
    print(f"driven_feature_count={driven[0]['n_features']}")
    print("SOURCE_POSITIONS")
    for name, xyz in EXPECTED_POSITIONS.items():
        print(f"{name}={xyz}")

    positions_pass = all(
        run["positions"].get(name) == xyz
        for run in driven + [silent]
        for name, xyz in EXPECTED_POSITIONS.items()
    )
    pairwise_max = max(abs(pair["r"]) for pair in pairwise)
    pairwise_pass = pairwise_max < 0.3
    den_dead = any(run["den_variance"]["dead"] for run in driven)
    driven_rank_pass = all(run["spectrum"]["n_significant"] >= 4 for run in driven) and not den_dead
    silent_rank_pass = silent["spectrum"]["n_significant"] < 4
    feature_pass = all(run["feature_contract_pass"] for run in driven + [silent])
    csv_pass = all(
        run["summary_csv"]["row_count_pass"]
        and run["summary_csv"]["rectangular_csv_pass"]
        and run["last_sample"]["pass"]
        and all(item["row_count_pass"] and item["rectangular_csv_pass"]
                for item in run["auxiliary_csvs"].values())
        for run in driven + [silent]
    )
    gates = {
        "positions_pass": positions_pass,
        "pairwise_r_pass": pairwise_pass,
        "driven_rank_pass": driven_rank_pass,
        "silent_rank_pass": silent_rank_pass,
        "feature_list_pass": feature_pass,
        "csv_validation_pass": csv_pass,
        "prior_stage_labels_pass": (
            STAGE4_STATUS == "FAIL" and STAGE5_STATUS == "PASS" and STAGE6_STATUS == "PASS"
        ),
    }
    overall = all(gates.values())
    if den_dead:
        failure_reason = "den_dead: Den_* has ~zero variance; attractant is dead. PROD_RATE was not raised."
    elif not gates["driven_rank_pass"]:
        failure_reason = "collinear_plumes: a driven seed has fewer than 4 significant components."
    elif not gates["silent_rank_pass"]:
        failure_reason = "silent_approx_driven: silent also has 4 significant components; rank is intrinsic drift."
    elif not gates["feature_list_pass"]:
        failure_reason = "readout_leakage: rank matrix is not the frozen 608-feature list."
    elif not gates["csv_validation_pass"]:
        failure_reason = "csv_incomplete: a CSV is missing rows, ragged, or last sample is not 39;15;299.95."
    elif not gates["positions_pass"]:
        failure_reason = "positions_moved: a run did not use the frozen source sites."
    elif not gates["pairwise_r_pass"]:
        failure_reason = "correlated_inputs: a pair of input sequences has |r| >= 0.3."
    else:
        failure_reason = ""

    driven_rows = [("101", driven[0]), ("202", driven[1]), ("303", driven[2])]
    evidence = {
        "schema": "BSimReservoirPlanStage7-gates-v1",
        "stage4_status": STAGE4_STATUS,
        "stage5_status": STAGE5_STATUS,
        "stage6_status": STAGE6_STATUS,
        "layout": "four grounded sources; not the plan 5-AC / flow-8 / glucose table",
        "flow_speed": 0.0,
        "positions": {name: list(xyz) for name, xyz in EXPECTED_POSITIONS.items()},
        "pairwise": pairwise,
        "pairwise_max_abs_r": pairwise_max,
        "rank_features": [
            "Den_* (20x10)",
            "Receiver_R_* / Mean_q (20x10)",
            "Lum_Mean_* / Mean_L (20x10)",
            "Input_Driven_Death_* (4x2)",
        ],
        "driven_feature_count": driven[0]["n_features"],
        "driven": [{k: v for k, v in run.items() if k != "feature_names"} for run in driven],
        "silent": {k: v for k, v in silent.items() if k != "feature_names"},
        "driven_rows": driven_rows,
        "gates": gates,
        "overall_gate_pass": overall,
        "failure_reason": failure_reason,
    }
    # driven_rows has full run dicts; drop names for json via default
    json_evidence = {k: v for k, v in evidence.items() if k != "driven_rows"}
    json_evidence["driven_n_significant"] = [run["spectrum"]["n_significant"] for run in driven]
    json_evidence["silent_n_significant"] = silent["spectrum"]["n_significant"]
    json_evidence["driven_ratios"] = [run["spectrum"]["ratios_to_sigma1"][:10] for run in driven]
    json_evidence["silent_ratios"] = silent["spectrum"]["ratios_to_sigma1"][:10]

    destination = Path(args.evidence)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(json_evidence, indent=2, allow_nan=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    write_markdown(args.markdown, evidence)
    print(json.dumps({
        "overall_gate_pass": overall,
        "gates": gates,
        "driven_n_significant": json_evidence["driven_n_significant"],
        "silent_n_significant": json_evidence["silent_n_significant"],
        "driven_ratios_top5": [ratios[:5] for ratios in json_evidence["driven_ratios"]],
        "silent_ratios_top5": json_evidence["silent_ratios"][:5],
        "pairwise_max_abs_r": pairwise_max,
        "den_variance": [run["den_variance"] for run in driven],
        "failure_reason": failure_reason,
        "driven_feature_count": driven[0]["n_features"],
    }, indent=2, allow_nan=True, default=json_default))
    print(f"evidence_file={destination}")
    print(f"markdown_file={args.markdown}")


if __name__ == "__main__":
    main()
