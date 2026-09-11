#!/usr/bin/env python3
"""Summarize Stage 3 source-strength calibration runs."""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def correlation(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    dx, dy = [v - mx for v in x], [v - my for v in y]
    denominator = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denominator if denominator else float("nan")


def summarize(path):
    csv_path = Path(path) / "window_summary.csv"
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    inputs = [float(row["AHL_Input"]) for row in rows]
    extracellular = [float(row["Extracellular_AHL_uM_Mean"]) for row in rows]
    q_values = [float(row["Mean_q"]) for row in rows]
    fractions = [float(row["Fraction_q_gt_0_5"]) for row in rows]
    windows = [float(row["Window"]) for row in rows]

    def grouped(values):
        return {
            str(level): statistics.mean(
                value for stimulus, value in zip(inputs, values) if stimulus == level
            )
            for level in (0.0, 0.25, 0.5, 0.75, 1.0)
        }

    return {
        "path": str(csv_path),
        "rows": len(rows),
        "mean_q_by_input": grouped(q_values),
        "fraction_q_gt_0_5_by_input": grouped(fractions),
        "zero_lag_r_extracellular_ahl": correlation(inputs, extracellular),
        "zero_lag_r_mean_q": correlation(inputs, q_values),
        "zero_lag_r_fraction_q_gt_0_5": correlation(inputs, fractions),
        "r_window_index_mean_q": correlation(windows, q_values),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--evidence")
    args = parser.parse_args()
    evidence = {"sweeps": [summarize(path) for path in args.runs]}
    text = json.dumps(evidence, indent=2) + "\n"
    if args.evidence:
        destination = Path(args.evidence)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
