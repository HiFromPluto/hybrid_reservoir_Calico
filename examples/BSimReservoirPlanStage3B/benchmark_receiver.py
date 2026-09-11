#!/usr/bin/env python3
"""Replay recorded Stage 3 AHL voxel traces through a fast Hill receiver."""

import argparse
import csv
import json
import math
import statistics
from pathlib import Path

MOLECULES_PER_UM3_PER_UM = 602.2
WINDOW_SECONDS = 300.0
HILL_N = 2.0
K_VALUES_UM = (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2)
TAU_VALUES_S = (15.0, 30.0, 60.0, 90.0, 120.0)


def pearson(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    dx, dy = [value - mx for value in x], [value - my for value in y]
    denominator = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denominator if denominator else float("nan")


def hill(concentration_um, k_um):
    concentration_n = concentration_um ** HILL_N
    return concentration_n / (k_um ** HILL_N + concentration_n)


def load_trace(path):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        window_index = header.index("Window")
        time_index = header.index("TimeInWindow_s")
        input_index = header.index("Input_AC1_AHL")
        ahl_indices = [i for i, name in enumerate(header) if name.startswith("AHL_")]
        density_indices = [i for i, name in enumerate(header) if name.startswith("Den_")]
        if len(ahl_indices) != 200 or len(density_indices) != 200:
            raise ValueError(f"{path}: expected 200 AHL and density voxels")
        rows = []
        for raw in reader:
            window = int(raw[window_index])
            offset = float(raw[time_index])
            rows.append({
                "window": window,
                "time": window * WINDOW_SECONDS + offset,
                "input": float(raw[input_index]),
                "ahl_um": [
                    float(raw[index]) / MOLECULES_PER_UM3_PER_UM
                    for index in ahl_indices
                ],
                "density": [float(raw[index]) for index in density_indices],
            })
    if len(rows) != 640:
        raise ValueError(f"{path}: expected 640 samples, found {len(rows)}")
    return {"path": str(path), "rows": rows}


def replay(trace, k_um, tau_s):
    rows = trace["rows"]
    state = [hill(value, k_um) for value in rows[0]["ahl_um"]]
    previous_time = rows[0]["time"]
    windows = {}
    for sample_index, row in enumerate(rows):
        if sample_index:
            dt = row["time"] - previous_time
            decay = math.exp(-dt / tau_s)
            for voxel, concentration in enumerate(row["ahl_um"]):
                target = hill(concentration, k_um)
                state[voxel] = target + (state[voxel] - target) * decay
        previous_time = row["time"]
        total_density = sum(row["density"])
        if total_density > 0:
            mean_response = sum(
                response * density
                for response, density in zip(state, row["density"])
            ) / total_density
            occupancy = sum(
                density
                for response, density in zip(state, row["density"])
                if response > 0.5
            ) / total_density
        else:
            mean_response = statistics.mean(state)
            occupancy = sum(response > 0.5 for response in state) / len(state)
        accumulator = windows.setdefault(row["window"], {
            "input": row["input"], "response": [], "occupancy": []
        })
        accumulator["response"].append(mean_response)
        accumulator["occupancy"].append(occupancy)

    summaries = []
    for window, values in sorted(windows.items()):
        summaries.append({
            "window": window,
            "input": values["input"],
            "mean_response": statistics.mean(values["response"]),
            "occupancy": statistics.mean(values["occupancy"]),
        })
    return summaries


def grouped_mean(rows, key, predicate):
    values = [row[key] for row in rows if predicate(row["input"])]
    return statistics.mean(values) if values else float("nan")


def evaluate(traces, k_um, tau_s):
    replayed = [
        {"path": trace["path"], "windows": replay(trace, k_um, tau_s)}
        for trace in traces
    ]
    correlations = [
        pearson(
            [row["input"] for row in result["windows"]],
            [row["mean_response"] for row in result["windows"]],
        )
        for result in replayed
    ]
    all_rows = [row for result in replayed for row in result["windows"]]
    low = grouped_mean(all_rows, "mean_response", lambda value: value <= 0.25)
    mid = grouped_mean(all_rows, "mean_response", lambda value: value == 0.5)
    high = grouped_mean(all_rows, "mean_response", lambda value: value >= 0.75)
    mid_occupancy = grouped_mean(all_rows, "occupancy", lambda value: value == 0.5)
    t95_s = -math.log(0.05) * tau_s
    gates = {
        "correlation_pass": all(
            math.isfinite(value) and value > 0.5 for value in correlations
        ),
        "mid_occupancy_pass": 0.2 <= mid_occupancy <= 0.8,
        "ordered_response_pass": low < mid < high,
        "timescale_pass": t95_s <= WINDOW_SECONDS,
    }
    return {
        "k_um": k_um,
        "hill_n": HILL_N,
        "tau_s": tau_s,
        "t95_s": t95_s,
        "correlations": correlations,
        "minimum_correlation": min(correlations),
        "low_mean_response": low,
        "mid_mean_response": mid,
        "high_mean_response": high,
        "mid_occupancy": mid_occupancy,
        "gates": gates,
        "overall_pass": all(gates.values()),
        "replayed": replayed,
    }


def candidate_rank(candidate):
    return (
        candidate["overall_pass"],
        sum(candidate["gates"].values()),
        candidate["minimum_correlation"],
        -abs(candidate["mid_occupancy"] - 0.5),
        candidate["high_mean_response"] - candidate["low_mean_response"],
    )


def write_selected_windows(path, selected):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerow([
            "Trace", "Window", "Input", "Mean_Response",
            "Fraction_Response_gt_0_5",
        ])
        for result in selected["replayed"]:
            for row in result["windows"]:
                writer.writerow([
                    result["path"], row["window"], row["input"],
                    row["mean_response"], row["occupancy"],
                ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("traces", nargs="+")
    parser.add_argument("--output-dir", default="results")
    args = parser.parse_args()
    traces = [load_trace(path) for path in args.traces]
    candidates = [
        evaluate(traces, k_um, tau_s)
        for k_um in K_VALUES_UM for tau_s in TAU_VALUES_S
    ]
    selected = max(candidates, key=candidate_rank)
    destination = Path(args.output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    write_selected_windows(destination / "selected_window_summary.csv", selected)
    evidence = {
        "schema": "BSimReservoirPlanStage3B-trace-benchmark-v1",
        "model": "first_order_hill",
        "equation": "dR/dt=(C^n/(K^n+C^n)-R)/tau",
        "hill_n": HILL_N,
        "initialization": "voxel receiver starts at Hill equilibrium of first recorded sample",
        "occupancy_semantics": "density-weighted fraction of recorded state voxels with R>0.5",
        "timescale_gate": "t95=-ln(0.05)*tau <= 300 s",
        "trace_files": [trace["path"] for trace in traces],
        "candidate_count": len(candidates),
        "selected": {key: value for key, value in selected.items() if key != "replayed"},
        "passing_candidate_count": sum(candidate["overall_pass"] for candidate in candidates),
        "candidates": [
            {key: value for key, value in candidate.items() if key != "replayed"}
            for candidate in candidates
        ],
    }
    (destination / "benchmark_evidence.json").write_text(
        json.dumps(evidence, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(evidence["selected"], indent=2, allow_nan=True))
    print(f"passing_candidate_count={evidence['passing_candidate_count']}")


if __name__ == "__main__":
    main()
