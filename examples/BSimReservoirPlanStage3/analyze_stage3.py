#!/usr/bin/env python3
"""Aggregate predeclared Plan Stage 3 numeric gates and write evidence JSON."""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def pearson(x, y):
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    mx, my = statistics.mean(x), statistics.mean(y)
    dx, dy = [v - mx for v in x], [v - my for v in y]
    den = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    return sum(a * b for a, b in zip(dx, dy)) / den if den else float("nan")


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


def read_run(path, expected_rows, expected_aux_rows=None):
    path = Path(path)
    csv_path = path if path.is_file() else path / "window_summary.csv"
    with csv_path.open(newline="", encoding="utf-8") as handle:
        raw = list(csv.reader(handle, delimiter=";"))
    widths = [len(row) for row in raw]
    malformed = not raw or any(width != widths[0] for width in widths)
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    numeric = []
    for row in rows:
        numeric.append({
            "window": int(row["Window"]),
            "input": float(row["AHL_Input"]),
            "extracellular_ahl": float(row["Extracellular_AHL_uM_Mean"]),
            "q": float(row["Mean_q"]),
            "fraction": float(row["Fraction_q_gt_0_5"]),
        })
    result = {
        "path": str(csv_path),
        "label": rows[0]["Run_Label"] if rows else "",
        "replicate": rows[0]["Stochastic_Replicate"] if rows else "",
        "rows": numeric,
        "row_count": len(rows),
        "column_count": widths[0] if widths else 0,
        "row_count_pass": len(rows) == expected_rows,
        "rectangular_csv_pass": not malformed,
    }
    if expected_aux_rows is not None:
        result["auxiliary_csvs"] = {
            name: validate_csv(path / name, expected_aux_rows)
            for name in ("results.csv", "voxels.csv")
        }
    return result


def grouped_mean(runs, key, predicate):
    values = [row[key] for run in runs for row in run["rows"] if predicate(row["input"])]
    return statistics.mean(values) if values else float("nan")


def correlations(runs, key, lag):
    values = []
    for run in runs:
        x = [row["input"] for row in run["rows"]]
        y = [row[key] for row in run["rows"]]
        if lag == 0:
            values.append(pearson(x, y))
        else:
            values.append(pearson(x[:-lag], y[lag:]))
    finite = [v for v in values if math.isfinite(v)]
    threshold = 0.5
    per_replicate_pass = [
        math.isfinite(value) and abs(value) > threshold for value in values
    ]
    return {
        "per_replicate": values,
        "mean": statistics.mean(finite) if finite else float("nan"),
        "absolute_threshold": threshold,
        "per_replicate_absolute_pass": per_replicate_pass,
        "all_replicates_absolute_pass": (
            len(values) >= 3 and all(per_replicate_pass)
        ),
    }


def index_correlations(runs, key):
    values = [
        pearson(
            [row["window"] for row in run["rows"]],
            [row[key] for row in run["rows"]],
        )
        for run in runs
    ]
    finite = [value for value in values if math.isfinite(value)]
    return {
        "per_replicate": values,
        "mean": statistics.mean(finite) if finite else float("nan"),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", nargs="+", required=True,
                        help="At least three selected-template output directories/CSVs")
    parser.add_argument("--silent", required=True, help="Silent-control output directory/CSV")
    parser.add_argument("--expected-windows", type=int, default=40)
    parser.add_argument("--expected-aux-rows", type=int, default=640)
    parser.add_argument("--evidence", default="results/stage3_gate_evidence.json")
    args = parser.parse_args()

    runs = [
        read_run(path, args.expected_windows, args.expected_aux_rows)
        for path in args.replicates
    ]
    silent = read_run(args.silent, args.expected_windows)
    low_fraction = grouped_mean(runs, "fraction", lambda u: u <= .25)
    mid_fraction = grouped_mean(runs, "fraction", lambda u: u == .5)
    high_fraction = grouped_mean(runs, "fraction", lambda u: u >= .75)
    active_q = grouped_mean(runs, "q", lambda u: True)
    silent_q = grouped_mean([silent], "q", lambda u: True)
    active_fraction = grouped_mean(runs, "fraction", lambda u: True)
    silent_fraction = grouped_mean([silent], "fraction", lambda u: True)

    zero_q = correlations(runs, "q", 0)
    zero_fraction = correlations(runs, "fraction", 0)
    correlation_gate_pass = (
        zero_q["all_replicates_absolute_pass"]
        or zero_fraction["all_replicates_absolute_pass"]
    )
    evidence = {
        "schema": "BSimReservoirPlanStage3-gates-v1",
        "replicate_semantics": "stochastic_replicate (rng.seed is incomplete)",
        "replicate_count": len(runs),
        "replicate_count_pass": len(runs) >= 3,
        "files": runs + [silent],
        "row_count_gate_pass": all(r["row_count_pass"] and r["rectangular_csv_pass"]
                                   for r in runs + [silent]),
        "auxiliary_csv_gate_pass": all(
            result["row_count_pass"] and result["rectangular_csv_pass"]
            for run in runs for result in run["auxiliary_csvs"].values()
        ),
        "fraction_q_gt_0_5": {
            "low_mean": low_fraction,
            "mid_mean": mid_fraction,
            "high_mean": high_fraction,
            "ordered_separation_pass": low_fraction < mid_fraction < high_fraction,
            "mid_target": [0.2, 0.8],
            "mid_target_pass": 0.2 <= mid_fraction <= 0.8,
        },
        "zero_lag": {
            "mean_q": zero_q,
            "fraction_q_gt_0_5": zero_fraction,
            "gate_rule": (
                "mean_q OR fraction_q_gt_0_5 has |r|>0.5 "
                "in every one of at least three stochastic replicates"
            ),
            "absolute_correlation_gate_pass": correlation_gate_pass,
        },
        "predeclared_lags_not_best_selected": {
            str(lag): {
                "mean_q": correlations(runs, "q", lag),
                "fraction_q_gt_0_5": correlations(runs, "fraction", lag),
            } for lag in (1, 2, 3)
        },
        "diagnostics_not_gate_changing": {
            "input_vs_extracellular_ahl_zero_lag": correlations(
                runs, "extracellular_ahl", 0
            ),
            "window_index_vs_mean_q": index_correlations(runs, "q"),
            "window_index_vs_fraction_q_gt_0_5": index_correlations(
                runs, "fraction"
            ),
        },
        "silent_control_contrast": {
            "active_mean_q": active_q, "silent_mean_q": silent_q,
            "delta_mean_q": active_q - silent_q,
            "active_mean_fraction": active_fraction,
            "silent_mean_fraction": silent_fraction,
            "delta_mean_fraction": active_fraction - silent_fraction,
            "pass": active_q > silent_q and active_fraction > silent_fraction,
        },
    }
    evidence["overall_gate_pass"] = (
        evidence["replicate_count_pass"]
        and evidence["row_count_gate_pass"]
        and evidence["auxiliary_csv_gate_pass"]
        and evidence["fraction_q_gt_0_5"]["ordered_separation_pass"]
        and evidence["fraction_q_gt_0_5"]["mid_target_pass"]
        and evidence["zero_lag"]["absolute_correlation_gate_pass"]
        and evidence["silent_control_contrast"]["pass"]
    )
    destination = Path(args.evidence)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(evidence, indent=2, allow_nan=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(evidence, indent=2, allow_nan=True))
    print(f"evidence_file={destination}")


if __name__ == "__main__":
    main()
