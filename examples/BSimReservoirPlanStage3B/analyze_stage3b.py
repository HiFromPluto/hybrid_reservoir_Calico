#!/usr/bin/env python3
"""Evaluate the predeclared spatial Plan Stage 3B production gate."""

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
WINDOW_SECONDS = 300.0
REQUIRED_COLUMNS = {"AHL_Input", "Mean_q", "Fraction_q_gt_0_5"}


def pearson(x, y):
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    mx, my = statistics.mean(x), statistics.mean(y)
    dx = [value - mx for value in x]
    dy = [value - my for value in y]
    denominator = math.sqrt(sum(value * value for value in dx)
                            * sum(value * value for value in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denominator if denominator else float("nan")


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
    numeric = [{
        "window": int(row["Window"]),
        "input": float(row["AHL_Input"]),
        "response": float(row["Mean_q"]),
        "fraction": float(row["Fraction_q_gt_0_5"]),
    } for row in rows] if not missing else []
    summary_validation = validate_csv(summary_path, expected_windows)
    auxiliary = ({
        name: validate_csv(directory / name, expected_aux_rows)
        for name in ("results.csv", "voxels.csv")
    } if expected_aux_rows is not None else {})
    return {
        "path": str(summary_path),
        "label": rows[0].get("Run_Label", "") if rows else "",
        "replicate": rows[0].get("Stochastic_Replicate", "") if rows else "",
        "rows": numeric,
        "required_columns": sorted(REQUIRED_COLUMNS),
        "missing_required_columns": missing,
        "required_columns_pass": not missing,
        "summary_csv": summary_validation,
        "auxiliary_csvs": auxiliary,
    }


def grouped_mean(runs, key, predicate):
    values = [
        row[key] for run in runs for row in run["rows"]
        if predicate(row["input"])
    ]
    return statistics.mean(values) if values else float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", nargs="+", required=True,
                        help="Three Stage3B production output directories")
    parser.add_argument("--silent", help="Optional diagnostic silent-control directory")
    parser.add_argument("--expected-windows", type=int, default=40)
    parser.add_argument("--expected-aux-rows", type=int, default=640)
    parser.add_argument("--evidence", default="results/stage3b_gate_evidence.json")
    args = parser.parse_args()

    runs = [
        read_run(path, args.expected_windows, args.expected_aux_rows)
        for path in args.replicates
    ]
    correlations = [
        pearson(
            [row["input"] for row in run["rows"]],
            [row["response"] for row in run["rows"]],
        )
        for run in runs
    ]
    correlation_passes = [
        math.isfinite(value) and abs(value) > 0.5 for value in correlations
    ]
    low_response = grouped_mean(runs, "response", lambda value: value <= 0.25)
    mid_response = grouped_mean(runs, "response", lambda value: value == 0.5)
    high_response = grouped_mean(runs, "response", lambda value: value >= 0.75)
    mid_occupancy = grouped_mean(runs, "fraction", lambda value: value == 0.5)

    summary_csv_pass = all(
        run["required_columns_pass"]
        and run["summary_csv"]["row_count_pass"]
        and run["summary_csv"]["rectangular_csv_pass"]
        for run in runs
    )
    auxiliary_csv_pass = all(
        validation["row_count_pass"] and validation["rectangular_csv_pass"]
        for run in runs for validation in run["auxiliary_csvs"].values()
    )
    gates = {
        "replicate_count_pass": len(runs) >= 3,
        "summary_csv_validation_pass": summary_csv_pass,
        "auxiliary_csv_validation_pass": auxiliary_csv_pass,
        "every_replicate_absolute_correlation_pass": (
            len(correlation_passes) >= 3 and all(correlation_passes)
        ),
        "mid_occupancy_pass": 0.2 <= mid_occupancy <= 0.8,
        "ordered_response_pass": low_response < mid_response < high_response,
        "response_timescale_pass": RECEIVER_T95_S <= WINDOW_SECONDS,
    }
    evidence = {
        "schema": "BSimReservoirPlanStage3B-gates-v1",
        "receiver": {
            "equation": "dR/dt=(C^2/(K^2+C^2)-R)/tau",
            "k_um": RECEIVER_K_UM,
            "hill_n": RECEIVER_HILL_N,
            "tau_s": RECEIVER_TAU_S,
            "t95_s": RECEIVER_T95_S,
            "t95_semantics": "analytic exact first-order response time, not empirical",
        },
        "gate_rules": {
            "correlation": "|r(AHL_Input, Mean_q)| > 0.5 in every replicate",
            "mid_occupancy": "0.2 <= mean Fraction_q_gt_0_5 at input 0.5 <= 0.8",
            "ordering": "low Mean_q < mid Mean_q < high Mean_q",
            "timescale": "t95 <= 300 s",
            "csv": "40 summary rows and 640 rectangular sample/voxel rows per replicate",
        },
        "replicate_count": len(runs),
        "files": runs,
        "correlations": correlations,
        "correlation_passes": correlation_passes,
        "low_mean_response": low_response,
        "mid_mean_response": mid_response,
        "high_mean_response": high_response,
        "mid_occupancy": mid_occupancy,
        "gates": gates,
    }
    if args.silent:
        silent = read_run(args.silent, args.expected_windows, None)
        silent_responses = [row["response"] for row in silent["rows"]]
        silent_fractions = [row["fraction"] for row in silent["rows"]]
        active_responses = [row["response"] for run in runs for row in run["rows"]]
        active_fractions = [row["fraction"] for run in runs for row in run["rows"]]
        silent_mean_response = statistics.mean(silent_responses)
        silent_mean_fraction = statistics.mean(silent_fractions)
        source_off_pass = (
            silent["required_columns_pass"]
            and silent["summary_csv"]["row_count_pass"]
            and silent["summary_csv"]["rectangular_csv_pass"]
            and silent_mean_response < 0.01
            and silent_mean_fraction < 0.01
            and statistics.mean(active_responses) > silent_mean_response
            and statistics.mean(active_fractions) > silent_mean_fraction
        )
        evidence["source_off_control"] = {
            "file": silent,
            "silent_mean_response": silent_mean_response,
            "silent_mean_fraction": silent_mean_fraction,
            "active_mean_response": statistics.mean(active_responses),
            "active_mean_fraction": statistics.mean(active_fractions),
            "pass": source_off_pass,
        }
        gates["source_off_control_pass"] = source_off_pass
    evidence["overall_gate_pass"] = all(gates.values())

    destination = Path(args.evidence)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(evidence, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(evidence, indent=2, allow_nan=True))
    print(f"evidence_file={destination}")


if __name__ == "__main__":
    main()
