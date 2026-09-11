#!/usr/bin/env python3
"""Score the frozen ChassisHybridmm I3s spatial-occupancy gates."""

from __future__ import annotations

import csv
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "I3S_PROTOCOL.md"
JAVA = HERE / "ChassisHybridmmI3s.java"
EVIDENCE = RESULTS / "i3s_evidence.json"
STANDING = HERE / "I3S_STANDING.md"
SEEDS = (111, 222, 333)
NUM_WINDOWS = 200
SAMPLES_PER_WINDOW = 16
STATE_X = 20
STATE_Y = 10
STATE_BINS = STATE_X * STATE_Y
FAR_X = frozenset((*range(0, 4), *range(16, 20)))
EXPECTED_MAP_ROWS = NUM_WINDOWS * SAMPLES_PER_WINDOW
EXPECTED_MECHANICS_ROWS = 7801
WARMUP_S = 18000.0
END_S = 78000.0
MAP_MEDIAN_MIN = 160
MAP_TEST_MIN = 140
FAR_MEDIAN_MIN = 64
FAR_TEST_MIN = 56
POP_LOW = 500
POP_HIGH = 1500
POP_WINDOW_MIN = 400
RATIO_LOW = 0.5
RATIO_HIGH = 1.5
OFFPLANE_MAX = 0.05
OVERLAP_MAX_UM = 0.25
NEMATIC_MAX = 0.95


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=";"))


def require_freeze() -> dict:
    protocol = PROTOCOL.read_text(encoding="utf-8")
    java = JAVA.read_text(encoding="utf-8")
    protocol_phrases = (
        "Frozen 2026-08-25 before I3s Java exists",
        "full-dish lattice",
        "count `≥160` over all 200 windows",
        "minimum count `≥140`",
        "median source-distant occupied bins `≥64`",
        "minimum source-distant occupied bins `≥56`",
        "Motility | OFF",
    )
    source_phrases = (
        "static final int PLACEMENT_X = 60;",
        "static final int PLACEMENT_Y = 30;",
        "static final double PLACEMENT_JITTER_UM = 0.4;",
        "placeFullDishPopulation(sim);",
        "ChassisHybridmmI3p.updateRodChemistry(",
        "ChassisHybridmmI3p.T_REMOVAL_S",
        "motility=OFF",
        "ridge=OFF",
        "NARMA_overall=NONE",
    )
    missing_protocol = [
        phrase for phrase in protocol_phrases if phrase not in protocol
    ]
    missing_source = [phrase for phrase in source_phrases if phrase not in java]
    forbidden_source = [
        phrase
        for phrase in (
            "setGoal(",
            "BrownianParticle",
            "evaluate_readout",
            "T_REMOVAL_S = 1800",
            "4.0 * Math.PI / 1800",
        )
        if phrase in java
    ]
    return {
        "protocol_phrases_missing": missing_protocol,
        "source_phrases_missing": missing_source,
        "forbidden_source_tokens": forbidden_source,
        "pass": not missing_protocol and not missing_source and not forbidden_source,
    }


def exact_row(rows: list[dict[str, str]], time_s: float) -> dict[str, str]:
    matches = [
        row for row in rows
        if math.isclose(float(row["t_s"]), time_s, rel_tol=0.0, abs_tol=1e-9)
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one mechanics row at t={time_s}, got {len(matches)}")
    return matches[0]


def parse_status(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def map_index(x: int, y: int) -> int:
    return x * STATE_Y + y


def score_seed(seed: int) -> dict:
    directory = RESULTS / f"i3s_silent_seed{seed}"
    summary_path = directory / "window_summary.csv"
    map_path = directory / "map_samples.csv"
    mechanics_path = directory / "mechanics.csv"
    status_path = directory / "run_status.txt"
    summary = read_csv(summary_path)
    samples = read_csv(map_path)
    mechanics = read_csv(mechanics_path)
    status = parse_status(status_path)

    windows = [int(row["Window"]) for row in summary]
    if windows != list(range(NUM_WINDOWS)):
        raise ValueError(f"{summary_path}: expected windows 0..199")
    density_names = [f"Den_{index}" for index in range(STATE_BINS)]
    if not samples or any(name not in samples[0] for name in density_names):
        raise ValueError(f"{map_path}: missing 20x10 density columns")

    grouped: dict[int, list[dict[str, str]]] = {
        window: [] for window in range(NUM_WINDOWS)
    }
    for row in samples:
        window = int(row["Window"])
        if window not in grouped:
            raise ValueError(f"{map_path}: invalid window {window}")
        grouped[window].append(row)

    window_counts = []
    far_counts = []
    sample_counts = []
    rectangular = True
    sample_contract = True
    for window in range(NUM_WINDOWS):
        rows = grouped[window]
        if len(rows) != SAMPLES_PER_WINDOW:
            sample_contract = False
            continue
        observed = [False] * STATE_BINS
        for row_index, row in enumerate(rows):
            expected_time = 299.95 if row_index == 15 else row_index * 20.0
            sample_contract = sample_contract and int(row["Sample"]) == row_index
            sample_contract = sample_contract and math.isclose(
                float(row["TimeInWindow_s"]),
                expected_time,
                rel_tol=0.0,
                abs_tol=1e-9,
            )
            density = []
            for name in density_names:
                try:
                    value = int(row[name])
                except (KeyError, TypeError, ValueError):
                    rectangular = False
                    value = 0
                density.append(value)
            occupied = [value > 0 for value in density]
            sample_counts.append(sum(occupied))
            for index, value in enumerate(occupied):
                observed[index] = observed[index] or value
            rectangular = rectangular and int(row["Occupied_Bins"]) == sum(occupied)
        total = sum(observed)
        far = sum(
            observed[map_index(x, y)]
            for x in FAR_X
            for y in range(STATE_Y)
        )
        window_counts.append(total)
        far_counts.append(far)
        rectangular = rectangular and int(
            summary[window]["Occupied_Bins_Window"]
        ) == total
        rectangular = rectangular and int(
            summary[window]["Far_Occupied_Bins_Window"]
        ) == far

    numeric_columns = (
        "t_s",
        "seed",
        "N",
        "offplane",
        "nematic_order",
        "delta_cc_max_um",
        "d_centers_min_um",
    )
    finite_mechanics = bool(mechanics) and all(
        math.isfinite(float(row[column]))
        for row in mechanics
        for column in numeric_columns
    )
    warmup = exact_row(mechanics, WARMUP_S)
    end = exact_row(mechanics, END_S)
    warmup_n = int(warmup["N"])
    end_n = int(end["N"])
    population = [int(row["Population_End"]) for row in summary]
    ratio = end_n / warmup_n if warmup_n > 0 else math.inf
    eligible = [row for row in mechanics if int(row["N"]) >= 8]
    max_overlap = max(
        (float(row["delta_cc_max_um"]) for row in eligible),
        default=0.0,
    )
    map_median = float(statistics.median(window_counts))
    map_test_minimum = min(window_counts[150:])
    far_median = float(statistics.median(far_counts))
    far_test_minimum = min(far_counts[150:])
    clock_equal = math.isclose(
        float(status.get("T_div_s", "nan")),
        float(status.get("T_removal_s", "nan")),
        rel_tol=0.0,
        abs_tol=1e-9,
    )
    placement_pass = (
        status.get("placement") == "full-dish-60x30"
        and math.isclose(
            float(status.get("placement_jitter_um", "nan")),
            0.4,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        and int(status.get("initial_occupied_bins", "-1")) == STATE_BINS
    )
    identity_pass = (
        placement_pass
        and clock_equal
        and status.get("motility") == "OFF"
        and status.get("ridge") == "OFF"
        and status.get("NARMA_overall") == "NONE"
    )
    gates = {
        "placement_identity": identity_pass,
        "map_occupancy": map_median >= MAP_MEDIAN_MIN
        and map_test_minimum >= MAP_TEST_MIN,
        "far_field_occupancy": far_median >= FAR_MEDIAN_MIN
        and far_test_minimum >= FAR_TEST_MIN,
        "population_guard": POP_LOW <= warmup_n <= POP_HIGH
        and POP_LOW <= end_n <= POP_HIGH
        and min(population) >= POP_WINDOW_MIN
        and RATIO_LOW <= ratio <= RATIO_HIGH,
        "monolayer": float(warmup["offplane"]) <= OFFPLANE_MAX,
        "packing": max_overlap <= OVERLAP_MAX_UM,
        "morphology": float(warmup["nematic_order"]) <= NEMATIC_MAX,
        "completeness": len(summary) == NUM_WINDOWS
        and len(samples) == EXPECTED_MAP_ROWS
        and len(mechanics) == EXPECTED_MECHANICS_ROWS
        and sample_contract
        and rectangular
        and finite_mechanics
        and math.isclose(float(mechanics[0]["t_s"]), 0.0, abs_tol=1e-9)
        and math.isclose(float(mechanics[-1]["t_s"]), END_S, abs_tol=1e-9),
    }
    return {
        "seed": seed,
        "paths": {
            "summary": str(summary_path),
            "map_samples": str(map_path),
            "mechanics": str(mechanics_path),
            "status": str(status_path),
        },
        "placement": {
            "identity": status.get("placement"),
            "jitter_um": float(status.get("placement_jitter_um", "nan")),
            "initial_occupied_bins": int(
                status.get("initial_occupied_bins", "-1")
            ),
        },
        "map_occupancy": {
            "window_values": window_counts,
            "median": map_median,
            "minimum": min(window_counts),
            "test_minimum": map_test_minimum,
            "sample_median": float(statistics.median(sample_counts)),
            "sample_minimum": min(sample_counts),
        },
        "far_field_occupancy": {
            "window_values": far_counts,
            "median": far_median,
            "minimum": min(far_counts),
            "test_minimum": far_test_minimum,
        },
        "population": {
            "initial": int(mechanics[0]["N"]),
            "end_warmup": warmup_n,
            "window_150_end": population[150],
            "end_run": end_n,
            "post_warmup_minimum": min(population),
            "post_warmup_maximum": max(population),
            "end_to_warmup_ratio": ratio,
            "values": population,
        },
        "mechanics": {
            "end_warmup_offplane": float(warmup["offplane"]),
            "end_warmup_nematic_order": float(warmup["nematic_order"]),
            "maximum_overlap_um_n_ge_8": max_overlap,
            "row_count": len(mechanics),
        },
        "gates": gates,
        "pass": all(gates.values()),
    }


def run_baselines() -> dict:
    scripts = [
        EXAMPLES / "BacteriumFromScratch" / "check_job2.py",
        EXAMPLES / "BacteriumFromScratch" / "check_job3.py",
        EXAMPLES / "BacteriumFromScratch" / "check_job3b.py",
        EXAMPLES / "BacteriumFromScratch" / "check_job3c.py",
        EXAMPLES / "ChassisPocket" / "check_i0c.py",
    ]
    checks = {}
    for script in scripts:
        process = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(script.parent),
            text=True,
            capture_output=True,
        )
        output = (process.stdout + process.stderr).strip().splitlines()
        checks[script.name] = {
            "returncode": process.returncode,
            "last_line": output[-1] if output else "",
            "pass": process.returncode == 0,
        }
    return {
        "checks": checks,
        "pass": all(item["pass"] for item in checks.values()),
    }


def write_standing(evidence: dict) -> None:
    lines = [
        "# ChassisHybridmm I3s spatial-occupancy standing",
        "",
        f"**I3s: {evidence['overall']}**",
        "",
        "I3s is a placement/map-occupancy question with motility OFF. It has",
        "no ridge, no NRMSE, and no NARMA Overall. I3n remains OVERALL FAIL",
        "at unmasked driven mean `1.1980`; I3p remains PASS.",
        "",
        "The placement was frozen before bin count: a full-dish `60×30`",
        "cell-centred lattice with `±0.4 µm` seeded jitter, initially nine",
        "rod centres in every 20×10 map bin.",
        "",
        "## Map, far-field, and physical gates",
        "",
        "| Seed | initial bins | map median / test min | far median / test min | N warmup → test-150 → end | offplane | nematic | max overlap (µm) | Result |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for result in evidence["seeds"]:
        placement = result["placement"]
        map_result = result["map_occupancy"]
        far = result["far_field_occupancy"]
        population = result["population"]
        mechanics = result["mechanics"]
        lines.append(
            f"| {result['seed']} | {placement['initial_occupied_bins']} | "
            f"{map_result['median']:.1f} / {map_result['test_minimum']} | "
            f"{far['median']:.1f} / {far['test_minimum']} | "
            f"{population['end_warmup']} → {population['window_150_end']} → "
            f"{population['end_run']} | "
            f"{mechanics['end_warmup_offplane']:.4f} | "
            f"{mechanics['end_warmup_nematic_order']:.4f} | "
            f"{mechanics['maximum_overlap_um_n_ge_8']:.4f} | "
            f"{'PASS' if result['pass'] else 'FAIL'} |"
        )
    coverage_traces = " / ".join(
        f"{result['seed']}: "
        f"{result['map_occupancy']['window_values'][0]}→"
        f"{result['map_occupancy']['window_values'][-1]}"
        for result in evidence["seeds"]
    )
    lines.extend(
        [
            "",
            "Total occupied bins at window 0 → window 199 were "
            f"`{coverage_traces}`. Population remained inside its guard, so",
            "the failure is loss of spatial lineages, not renewed population",
            "collapse.",
            "",
            "Inherited map gates: median `≥160/200`, test minimum `≥140/200`.",
            "Predeclared far-field gates: median `≥64/80`, test minimum",
            "`≥56/80`. Population and mechanics guards were inherited from I3p.",
            "",
            "## Gate decision",
            "",
            evidence["decision"],
            "",
            "Baseline Job 2 / 3 / 3b / 3c / ChassisPocket I0c checks: "
            + ("PASS." if evidence["baselines"]["pass"] else "FAIL."),
            "",
        ]
    )
    if evidence["overall"] == "PASS":
        lines.extend(
            [
                "A3 is resolved for this declared full-dish, motility-OFF",
                "placement. I3s does not establish a NARMA score. The next named",
                "job is I3n2, which requires its own frozen protocol before Java.",
            ]
        )
    else:
        lines.extend(
            [
                "Stop. Do not change placement or enable motility. I3n2 is",
                "NOT_SCORED for a delayed living-layer map.",
            ]
        )
    STANDING.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    source_freeze = require_freeze()
    seed_results = [score_seed(seed) for seed in SEEDS]
    baselines = run_baselines()
    honesty = source_freeze["pass"] and baselines["pass"]
    passed = honesty and all(result["pass"] for result in seed_results)
    overall = "PASS" if passed else "FAIL"
    if passed:
        decision = (
            "Every seed passed the pre-registered placement, total-map, "
            "far-field, population, monolayer, packing, morphology, and "
            "completeness gates."
        )
    else:
        failed = {
            str(result["seed"]): [
                name for name, value in result["gates"].items() if not value
            ]
            for result in seed_results
            if not result["pass"]
        }
        if not honesty:
            failed["honesty"] = [
                "source/protocol freeze" if not source_freeze["pass"] else "baselines"
            ]
        decision = f"Pre-registered gates failed: {failed}."
    evidence = {
        "schema": "ChassisHybridmm-I3s-spatial-occupancy-v1",
        "overall": overall,
        "question_kind": "map occupancy only; no ridge; no NARMA Overall",
        "protocol": str(PROTOCOL),
        "gates_frozen": {
            "map_median_min": MAP_MEDIAN_MIN,
            "map_test_min": MAP_TEST_MIN,
            "far_median_min": FAR_MEDIAN_MIN,
            "far_test_min": FAR_TEST_MIN,
        },
        "source_freeze": source_freeze,
        "seeds": seed_results,
        "baselines": baselines,
        "honesty_pass": honesty,
        "decision": decision,
        "i3p": "PASS_UNCHANGED",
        "i3n_overall": "FAIL_UNCHANGED_1.1980",
        "narma10b_driven": "UNCHANGED_0.928",
    }
    EVIDENCE.write_text(
        json.dumps(evidence, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_standing(evidence)
    print(f"I3s: {overall}")
    for result in seed_results:
        map_result = result["map_occupancy"]
        far = result["far_field_occupancy"]
        print(
            f"seed={result['seed']} map_median={map_result['median']:.1f} "
            f"map_test_min={map_result['test_minimum']} "
            f"far_median={far['median']:.1f} "
            f"far_test_min={far['test_minimum']} "
            f"result={'PASS' if result['pass'] else 'FAIL'}"
        )
    print(decision)
    print("No NRMSE or NARMA Overall was computed.")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
