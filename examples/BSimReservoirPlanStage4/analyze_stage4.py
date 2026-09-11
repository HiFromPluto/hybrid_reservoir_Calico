#!/usr/bin/env python3
"""Evaluate Plan Stage 4 death-channel gates and frozen Stage 3B receiver gates."""

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
REQUIRED_COLUMNS = {
    "AHL_Input",
    "Acid_Input",
    "Mean_q",
    "Fraction_q_gt_0_5",
    "Total_Deaths",
    "Clamp_Deaths",
    "Input_Driven_Deaths",
    "OOB_Deaths",
    "pH_Mean",
    "Acid_mM_Mean",
}


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
        "ahl_input": float(row["AHL_Input"]),
        "acid_input": float(row["Acid_Input"]),
        "response": float(row["Mean_q"]),
        "fraction": float(row["Fraction_q_gt_0_5"]),
        "total_deaths": float(row["Total_Deaths"]),
        "clamp_deaths": float(row["Clamp_Deaths"]),
        "input_driven_deaths": float(row["Input_Driven_Deaths"]),
        "oob_deaths": float(row["OOB_Deaths"]),
        "ph_mean": float(row["pH_Mean"]),
        "acid_mm_mean": float(row["Acid_mM_Mean"]),
        "population": float(row["Population"]),
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
        "mean_total_deaths": statistics.mean(
            row["total_deaths"] for row in numeric
        ) if numeric else float("nan"),
        "mean_input_driven_deaths": statistics.mean(
            row["input_driven_deaths"] for row in numeric
        ) if numeric else float("nan"),
        "sum_input_driven_deaths": sum(
            row["input_driven_deaths"] for row in numeric
        ) if numeric else 0.0,
        "mean_ph": statistics.mean(row["ph_mean"] for row in numeric) if numeric else float("nan"),
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
    parser.add_argument("--evidence", default="results/stage4_gate_evidence.json")
    args = parser.parse_args()

    runs = [
        read_run(path, args.expected_windows, args.expected_aux_rows)
        for path in args.replicates
    ]
    silent = read_run(args.silent, args.expected_windows, args.expected_aux_rows)

    death_correlations = [
        pearson(
            [row["acid_input"] for row in run["rows"]],
            [row["input_driven_deaths"] for row in run["rows"]],
        )
        for run in runs
    ]
    death_correlation_passes = [
        math.isfinite(value) and abs(value) > 0.5 for value in death_correlations
    ]
    receiver_correlations = [
        pearson(
            [row["ahl_input"] for row in run["rows"]],
            [row["response"] for row in run["rows"]],
        )
        for run in runs
    ]
    receiver_correlation_passes = [
        math.isfinite(value) and abs(value) > 0.5 for value in receiver_correlations
    ]
    low_response = grouped_mean(runs, "response", "ahl_input", lambda value: value <= 0.25)
    mid_response = grouped_mean(runs, "response", "ahl_input", lambda value: value == 0.5)
    high_response = grouped_mean(runs, "response", "ahl_input", lambda value: value >= 0.75)
    mid_occupancy = grouped_mean(runs, "fraction", "ahl_input", lambda value: value == 0.5)

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
    mean_total_active = statistics.mean(run["mean_total_deaths"] for run in runs)
    death_fold = (
        mean_total_active / silent["mean_total_deaths"]
        if silent["mean_total_deaths"] not in (0, float("nan")) else float("nan")
    )
    source_off_no_input_kill = silent["sum_input_driven_deaths"] == 0

    gates = {
        "replicate_count_pass": len(runs) >= 3,
        "death_correlation_pass": (
            len(death_correlation_passes) >= 3 and all(death_correlation_passes)
        ),
        "death_fold_pass": math.isfinite(death_fold) and death_fold >= 2.0,
        "source_off_no_input_driven_deaths_pass": source_off_no_input_kill,
        "summary_csv_validation_pass": summary_csv_pass,
        "auxiliary_csv_validation_pass": auxiliary_csv_pass,
        "receiver_correlation_pass": (
            len(receiver_correlation_passes) >= 3 and all(receiver_correlation_passes)
        ),
        "receiver_mid_occupancy_pass": 0.2 <= mid_occupancy <= 0.8,
        "receiver_ordered_response_pass": low_response < mid_response < high_response,
        "receiver_timescale_pass": RECEIVER_T95_S <= WINDOW_SECONDS,
    }
    evidence = {
        "schema": "BSimReservoirPlanStage4-gates-v1",
        "receiver": {
            "equation": "dR/dt=(C^2/(K^2+C^2)-R)/tau",
            "k_um": RECEIVER_K_UM,
            "hill_n": RECEIVER_HILL_N,
            "tau_s": RECEIVER_TAU_S,
            "t95_s": RECEIVER_T95_S,
            "status": "frozen Stage 3B; not retuned",
        },
        "gate_rules": {
            "death_correlation": "|r(Acid_Input, Input_Driven_Deaths)| > 0.5 in every replicate",
            "death_fold": "mean Total_Deaths at least 2x source-off baseline",
            "source_off": "source-off Input_Driven_Deaths sum to 0",
            "csv": "40 summary rows and 640 rectangular sample/voxel rows",
            "receiver": "Stage 3B receiver gates still pass on the same protocol",
        },
        "replicate_count": len(runs),
        "files": runs,
        "death_correlations": death_correlations,
        "death_correlation_passes": death_correlation_passes,
        "receiver_correlations": receiver_correlations,
        "receiver_correlation_passes": receiver_correlation_passes,
        "low_mean_response": low_response,
        "mid_mean_response": mid_response,
        "high_mean_response": high_response,
        "mid_occupancy": mid_occupancy,
        "mean_total_deaths_active": mean_total_active,
        "mean_total_deaths_silent": silent["mean_total_deaths"],
        "death_fold_vs_silent": death_fold,
        "source_off_control": {
            "file": silent,
            "sum_input_driven_deaths": silent["sum_input_driven_deaths"],
            "mean_ph": silent["mean_ph"],
        },
        "gates": gates,
        "overall_gate_pass": all(gates.values()),
    }

    destination = Path(args.evidence)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(evidence, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(evidence, indent=2, allow_nan=True))
    print(f"evidence_file={destination}")
    print(f"overall_gate_pass={evidence['overall_gate_pass']}")


if __name__ == "__main__":
    main()
