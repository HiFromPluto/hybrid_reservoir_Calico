#!/usr/bin/env python3
"""Evaluate Plan Stage 5 feature-set gates. Stage 4 remains FAIL."""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path

RECEIVER_K_UM = 1.6
RECEIVER_HILL_N = 2.0
RECEIVER_TAU_S = 15.0
RECEIVER_T95_S = -math.log(0.05) * RECEIVER_TAU_S
DELTA_LUX = 1.0 / 1500.0
LUM_T95_S = -math.log(0.05) / DELTA_LUX
WINDOW_SECONDS = 300.0

ANALYSIS_CHANNELS = [
    "Mean_q",
    "Input_Driven_Deaths",
    "Mean_L",
]
EXPORT_ONLY = [
    "Population",
    "Births",
    "Total_Deaths",
    "Clamp_Deaths",
    "OOB_Deaths",
    "Lum_Sum",
    "Extracellular_AHL_uM_Mean",
    "Fraction_q_gt_0_5",
    "pH_Mean",
]
REQUIRED_COLUMNS = set(ANALYSIS_CHANNELS + EXPORT_ONLY + [
    "AHL_Input", "Acid_Input", "L_P10", "L_P50", "L_P90",
])


def pearson(x, y):
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    mx, my = statistics.mean(x), statistics.mean(y)
    dx = [value - mx for value in x]
    dy = [value - my for value in y]
    denominator = math.sqrt(sum(value * value for value in dx)
                            * sum(value * value for value in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denominator if denominator else float("nan")


def variance(values):
    return statistics.pvariance(values) if len(values) >= 2 else float("nan")


def validate_csv(path, expected_rows):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, [])
        widths = [len(row) for row in reader]
    return {
        "path": str(path),
        "row_count": len(widths),
        "column_count": len(header),
        "row_count_pass": len(widths) == expected_rows,
        "rectangular_csv_pass": bool(header) and all(
            width == len(header) for width in widths
        ),
    }


def read_run(path, expected_windows, expected_aux_rows):
    directory = Path(path)
    summary_path = directory / "window_summary.csv"
    with summary_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        header = set(reader.fieldnames or [])
        rows = list(reader)
    missing = sorted(REQUIRED_COLUMNS - header)
    numeric = []
    if not missing:
        for row in rows:
            item = {
                "window": int(row["Window"]),
                "ahl_input": float(row["AHL_Input"]),
                "acid_input": float(row["Acid_Input"]),
            }
            for name in ANALYSIS_CHANNELS + EXPORT_ONLY:
                item[name] = float(row[name])
            numeric.append(item)
    channels = {
        name: [row[name] for row in numeric] for name in ANALYSIS_CHANNELS
    } if numeric else {}
    matrix = {
        a: {b: pearson(channels[a], channels[b]) for b in ANALYSIS_CHANNELS}
        for a in ANALYSIS_CHANNELS
    } if channels else {}
    collinear = []
    for i, a in enumerate(ANALYSIS_CHANNELS):
        for b in ANALYSIS_CHANNELS[i + 1:]:
            value = matrix.get(a, {}).get(b, float("nan"))
            if math.isfinite(value) and abs(value) > 0.9:
                collinear.append({"a": a, "b": b, "r": value})
    summary_validation = validate_csv(summary_path, expected_windows)
    auxiliary = ({
        name: validate_csv(directory / name, expected_aux_rows)
        for name in ("results.csv", "voxels.csv")
    } if expected_aux_rows is not None else {})
    windows = [row["window"] for row in numeric]
    return {
        "path": str(summary_path),
        "label": rows[0].get("Run_Label", "") if rows else "",
        "replicate": rows[0].get("Stochastic_Replicate", "") if rows else "",
        "rows": numeric,
        "missing_required_columns": missing,
        "required_columns_pass": not missing,
        "summary_csv": summary_validation,
        "auxiliary_csvs": auxiliary,
        "channel_variance": {name: variance(values) for name, values in channels.items()},
        "correlation_matrix": matrix,
        "collinear_pairs": collinear,
        "receiver_r": pearson(
            [row["ahl_input"] for row in numeric],
            [row["Mean_q"] for row in numeric],
        ) if numeric else float("nan"),
        "death_r": pearson(
            [row["acid_input"] for row in numeric],
            [row["Input_Driven_Deaths"] for row in numeric],
        ) if numeric else float("nan"),
        "lum_vs_q_r": pearson(
            [row["Mean_L"] for row in numeric],
            [row["Mean_q"] for row in numeric],
        ) if numeric else float("nan"),
        "window_vs_q_r": pearson(windows, [row["Mean_q"] for row in numeric]) if numeric else float("nan"),
        "window_vs_l_r": pearson(windows, [row["Mean_L"] for row in numeric]) if numeric else float("nan"),
        "lag1_q_r": pearson(
            [row["Mean_q"] for row in numeric[:-1]],
            [row["Mean_q"] for row in numeric[1:]],
        ) if len(numeric) > 2 else float("nan"),
        "lag1_l_r": pearson(
            [row["Mean_L"] for row in numeric[:-1]],
            [row["Mean_L"] for row in numeric[1:]],
        ) if len(numeric) > 2 else float("nan"),
        "mean_L": statistics.mean(row["Mean_L"] for row in numeric) if numeric else float("nan"),
        "mid_mean_L": statistics.mean(
            row["Mean_L"] for row in numeric if row["ahl_input"] == 0.5
        ) if numeric else float("nan"),
    }


def grouped_mean(runs, key, input_key, predicate):
    values = [
        row[key] for run in runs for row in run["rows"]
        if predicate(row[input_key])
    ]
    return statistics.mean(values) if values else float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", nargs="+", required=True)
    parser.add_argument("--silent", required=True)
    parser.add_argument("--expected-windows", type=int, default=40)
    parser.add_argument("--expected-aux-rows", type=int, default=640)
    parser.add_argument("--evidence", default="results/stage5_gate_evidence.json")
    args = parser.parse_args()

    runs = [
        read_run(path, args.expected_windows, args.expected_aux_rows)
        for path in args.replicates
    ]
    silent = read_run(args.silent, args.expected_windows, args.expected_aux_rows)

    variance_passes = [
        all(math.isfinite(value) and value > 0 for value in run["channel_variance"].values())
        for run in runs
    ]
    collinear_pass = all(not run["collinear_pairs"] for run in runs)
    receiver_rs = [run["receiver_r"] for run in runs]
    receiver_r_pass = [
        math.isfinite(value) and abs(value) > 0.5 for value in receiver_rs
    ]
    death_rs = [run["death_r"] for run in runs]
    death_r_pass = [math.isfinite(value) and abs(value) > 0.5 for value in death_rs]
    lum_vs_q = [run["lum_vs_q_r"] for run in runs]
    lum_not_duplicate = [
        math.isfinite(value) and abs(value) < 0.9 for value in lum_vs_q
    ]
    low_response = grouped_mean(runs, "Mean_q", "ahl_input", lambda value: value <= 0.25)
    mid_response = grouped_mean(runs, "Mean_q", "ahl_input", lambda value: value == 0.5)
    high_response = grouped_mean(runs, "Mean_q", "ahl_input", lambda value: value >= 0.75)
    mid_occupancy = grouped_mean(
        runs, "Fraction_q_gt_0_5", "ahl_input", lambda value: value == 0.5
    )
    mid_mean_L = grouped_mean(runs, "Mean_L", "ahl_input", lambda value: value == 0.5)

    summary_csv_pass = all(
        run["required_columns_pass"]
        and run["summary_csv"]["row_count_pass"]
        and run["summary_csv"]["rectangular_csv_pass"]
        for run in runs + [silent]
    )
    auxiliary_csv_pass = all(
        validation["row_count_pass"] and validation["rectangular_csv_pass"]
        for run in runs + [silent] for validation in run["auxiliary_csvs"].values()
    )
    lum_usable = math.isfinite(mid_mean_L) and 0.05 <= mid_mean_L <= 0.95

    gates = {
        "replicate_count_pass": len(runs) >= 3,
        "nonzero_variance_pass": len(variance_passes) >= 3 and all(variance_passes),
        "no_collinear_analysis_pair_pass": collinear_pass,
        "receiver_correlation_pass": len(receiver_r_pass) >= 3 and all(receiver_r_pass),
        "receiver_mid_occupancy_pass": 0.2 <= mid_occupancy <= 0.8,
        "receiver_ordered_response_pass": low_response < mid_response < high_response,
        "input_driven_death_correlation_pass": len(death_r_pass) >= 3 and all(death_r_pass),
        "lum_not_duplicate_of_R_pass": len(lum_not_duplicate) >= 3 and all(lum_not_duplicate),
        "lum_usable_scale_pass": lum_usable,
        "summary_csv_validation_pass": summary_csv_pass,
        "auxiliary_csv_validation_pass": auxiliary_csv_pass,
        "stage4_remains_fail": True,
    }
    evidence = {
        "schema": "BSimReservoirPlanStage5-gates-v1",
        "stage4_status": "FAIL",
        "receiver": {
            "k_um": RECEIVER_K_UM,
            "hill_n": RECEIVER_HILL_N,
            "tau_s": RECEIVER_TAU_S,
            "t95_s": RECEIVER_T95_S,
            "status": "frozen Stage 3B; not retuned",
        },
        "luminescence": {
            "equation": "L += (ALPHA_LUX * R - DELTA_LUX * L) * dt",
            "delta_lux": DELTA_LUX,
            "tau_s": 1500.0,
            "t95_s": LUM_T95_S,
            "t95_semantics": "analytic -ln(0.05)/DELTA_LUX, not measured",
            "initialization": "L=0; 18000 s warmup sets the state",
            "alpha_lux": 1.0 / 1500.0,
        },
        "analysis_channels": ANALYSIS_CHANNELS,
        "export_only_not_claimed": EXPORT_ONLY,
        "dropped_columns": {
            "Att_Rep_voxels": "silent AC0/AC2 architecture; zero-variance fields",
            "AHL_molecules_um3": "linear duplicate of Extracellular_AHL_uM_Mean",
            "Receiver_R_Mean": "identical to Mean_q",
            "Acid_mM": "affine duplicate of pH via pH=7.1-acid_mM/2",
            "total_Death_voxels": "replaced by clamp / input-driven / OOB split",
            "window_AHL_uM_Mean": "analysis drop: |r|~0.98 with Mean_q because tau_R=15s << 300s window; kept as diagnostic/voxel field",
            "Fraction_q_gt_0_5": "analysis drop: |r|~0.93-0.95 with Mean_q; kept as occupancy diagnostic for the Stage 3B gate",
            "pH_Mean": "analysis drop: |r|~-0.87 to -0.91 with Input_Driven_Deaths; kept as voxel field",
        },
        "gate_rules": {
            "variance": "every analysis channel has non-zero variance in every replicate",
            "collinearity": "no analysis pair has |r| > 0.9",
            "receiver": "Stage 3B receiver gates still pass",
            "deaths": "|r(Acid_Input, Input_Driven_Deaths)| > 0.5; Total_Deaths not a pass condition",
            "lum": "|r(Mean_L, Mean_q)| < 0.9",
            "csv": "40 summary rows and 640 rectangular sample/voxel rows",
            "stage4": "Stage 4 remains FAIL",
        },
        "files": runs,
        "silent": {
            "path": silent["path"],
            "death_r": silent["death_r"],
            "sum_input_driven_deaths": sum(
                row["Input_Driven_Deaths"] for row in silent["rows"]
            ) if silent["rows"] else 0,
        },
        "receiver_correlations": receiver_rs,
        "death_correlations": death_rs,
        "lum_vs_q_correlations": lum_vs_q,
        "window_vs_q": [run["window_vs_q_r"] for run in runs],
        "window_vs_l": [run["window_vs_l_r"] for run in runs],
        "lag1_q": [run["lag1_q_r"] for run in runs],
        "lag1_l": [run["lag1_l_r"] for run in runs],
        "low_mean_response": low_response,
        "mid_mean_response": mid_response,
        "high_mean_response": high_response,
        "mid_occupancy": mid_occupancy,
        "mid_mean_L": mid_mean_L,
        "gates": gates,
        "overall_gate_pass": all(gates.values()),
    }

    destination = Path(args.evidence)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(evidence, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "overall_gate_pass": evidence["overall_gate_pass"],
        "gates": gates,
        "receiver_correlations": receiver_rs,
        "death_correlations": death_rs,
        "lum_vs_q_correlations": lum_vs_q,
        "window_vs_q": evidence["window_vs_q"],
        "window_vs_l": evidence["window_vs_l"],
        "lag1_q": evidence["lag1_q"],
        "lag1_l": evidence["lag1_l"],
        "mid_occupancy": mid_occupancy,
        "mid_mean_L": mid_mean_L,
        "collinear_pairs": [run["collinear_pairs"] for run in runs],
        "channel_variance": [run["channel_variance"] for run in runs],
        "correlation_matrices": [run["correlation_matrix"] for run in runs],
    }, indent=2, allow_nan=True))
    print(f"evidence_file={destination}")


if __name__ == "__main__":
    main()
