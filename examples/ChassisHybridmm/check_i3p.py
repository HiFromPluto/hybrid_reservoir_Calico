#!/usr/bin/env python3
"""Score the pre-registered ChassisHybridmm I3p population-clock gates."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
RESULTS = HERE / "results"
PROTOCOL = HERE / "I3P_PROTOCOL.md"
JAVA = HERE / "ChassisHybridmmI3p.java"
EVIDENCE = RESULTS / "i3p_evidence.json"
STANDING = HERE / "I3P_STANDING.md"
SEEDS = (111, 222, 333)
WARMUP_S = 18000.0
END_S = 78000.0
EXPECTED_SUMMARY_ROWS = 200
EXPECTED_MECHANICS_ROWS = 7801
POP_WARMUP_LOW = 500
POP_WARMUP_HIGH = 1500
POP_END_LOW = 500
POP_END_HIGH = 1500
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
        "Frozen 2026-08-25 before I3p Java exists",
        "T_removal = T_div",
        "`500 ≤ N(18000 s) ≤ 1500`",
        "`500 ≤ N(78000 s) ≤ 1500`",
        "NARMA Overall. It does not amend I3n",
    )
    source_phrases = (
        "static final double T_DIV_S = ChassisParameters.analyticDivisionTimeS(",
        "static final double T_REMOVAL_S = T_DIV_S;",
        "double pRemoval = (DT / T_REMOVAL_S)",
        "rod.body.elongateCited(ChassisParameters.DT_S,",
        "Double.isFinite(stats.minCentre)",
        "d_centers_no_pair_sentinel_um=%.12g",
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
            "T_REMOVAL_S = 1800",
            "4.0 * Math.PI / 1800",
            "setGoal(",
            "BrownianParticle",
            "evaluate_readout",
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


def score_seed(seed: int) -> dict:
    directory = RESULTS / f"i3p_silent_seed{seed}"
    summary_path = directory / "window_summary.csv"
    mechanics_path = directory / "mechanics.csv"
    status_path = directory / "run_status.txt"
    summary = read_csv(summary_path)
    mechanics = read_csv(mechanics_path)
    windows = [int(row["Window"]) for row in summary]
    if windows != list(range(EXPECTED_SUMMARY_ROWS)):
        raise ValueError(f"{summary_path}: expected windows 0..199")
    if not mechanics:
        raise ValueError(f"{mechanics_path}: empty")

    numeric_columns = (
        "t_s",
        "seed",
        "N",
        "offplane",
        "nematic_order",
        "delta_cc_max_um",
        "d_centers_min_um",
    )
    finite = all(
        math.isfinite(float(row[column]))
        for row in mechanics
        for column in numeric_columns
    )
    warmup = exact_row(mechanics, WARMUP_S)
    end = exact_row(mechanics, END_S)
    warmup_n = int(warmup["N"])
    end_n = int(end["N"])
    window_population = [int(row["Population_End"]) for row in summary]
    ratio = end_n / warmup_n if warmup_n > 0 else math.inf
    eligible = [row for row in mechanics if int(row["N"]) >= 8]
    max_overlap = max(
        (float(row["delta_cc_max_um"]) for row in eligible),
        default=0.0,
    )
    status = status_path.read_text(encoding="utf-8")
    expected_division = 2595.142595
    status_clock_values = {}
    for line in status.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            status_clock_values[key] = value
    clock_equal = math.isclose(
        float(status_clock_values.get("T_div_s", "nan")),
        float(status_clock_values.get("T_removal_s", "nan")),
        rel_tol=0.0,
        abs_tol=1e-9,
    )
    clock_value = float(status_clock_values.get("T_div_s", "nan"))
    no_pair_sentinel = float(
        status_clock_values.get("d_centers_no_pair_sentinel_um", "nan")
    )
    gates = {
        "clock_identity": clock_equal
        and math.isclose(clock_value, expected_division, rel_tol=0.0, abs_tol=0.01)
        and status_clock_values.get("ridge") == "OFF"
        and status_clock_values.get("NARMA_overall") == "NONE",
        "end_warmup_population": POP_WARMUP_LOW <= warmup_n <= POP_WARMUP_HIGH,
        "end_run_population": POP_END_LOW <= end_n <= POP_END_HIGH,
        "post_warmup_hold": min(window_population) >= POP_WINDOW_MIN
        and RATIO_LOW <= ratio <= RATIO_HIGH,
        "monolayer": float(warmup["offplane"]) <= OFFPLANE_MAX,
        "packing": max_overlap <= OVERLAP_MAX_UM,
        "morphology": float(warmup["nematic_order"]) <= NEMATIC_MAX,
        "completeness": len(summary) == EXPECTED_SUMMARY_ROWS
        and len(mechanics) == EXPECTED_MECHANICS_ROWS
        and math.isclose(float(mechanics[0]["t_s"]), 0.0, abs_tol=1e-9)
        and math.isclose(float(mechanics[-1]["t_s"]), END_S, abs_tol=1e-9)
        and math.isclose(no_pair_sentinel, 4.5, abs_tol=1e-12)
        and finite,
    }
    return {
        "seed": seed,
        "paths": {
            "summary": str(summary_path),
            "mechanics": str(mechanics_path),
            "status": str(status_path),
        },
        "population": {
            "initial": int(mechanics[0]["N"]),
            "end_warmup": warmup_n,
            "window_0_end": window_population[0],
            "window_150_end": window_population[150],
            "end_run": end_n,
            "post_warmup_minimum": min(window_population),
            "post_warmup_maximum": max(window_population),
            "end_to_warmup_ratio": ratio,
            "values": window_population,
        },
        "mechanics": {
            "end_warmup_offplane": float(warmup["offplane"]),
            "end_warmup_nematic_order": float(warmup["nematic_order"]),
            "maximum_overlap_um_n_ge_8": max_overlap,
            "row_count": len(mechanics),
        },
        "clock": {
            "T_div_s": clock_value,
            "T_removal_s": float(status_clock_values.get("T_removal_s", "nan")),
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
    overall = evidence["overall"]
    lines = [
        "# ChassisHybridmm I3p population-clock standing",
        "",
        f"**I3p: {overall}**",
        "",
        "I3p is a population-only clock-join question. It has no ridge, no",
        "NRMSE, and no NARMA Overall. I3n remains OVERALL FAIL at unmasked",
        "driven mean `1.1980`; I3n-D remains post-hoc.",
        "",
        "The only model join was declared before population was observed:",
        r"`T_removal = T_div = analyticDivisionTimeS(0.5 mM)` while keeping",
        "the I3n clamp formula, cadence, rods, placement, mechanics, and silent",
        "chemistry scenario.",
        "",
        "## Population and physical gates",
        "",
        "| Seed | N: initial → warmup → window 150 → end | post-warmup min–max | end/warmup | offplane | nematic | max overlap (µm) | Result |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for result in evidence["seeds"]:
        population = result["population"]
        mechanics = result["mechanics"]
        lines.append(
            f"| {result['seed']} | {population['initial']} → "
            f"{population['end_warmup']} → {population['window_150_end']} → "
            f"{population['end_run']} | {population['post_warmup_minimum']}–"
            f"{population['post_warmup_maximum']} | "
            f"{population['end_to_warmup_ratio']:.3f} | "
            f"{mechanics['end_warmup_offplane']:.4f} | "
            f"{mechanics['end_warmup_nematic_order']:.4f} | "
            f"{mechanics['maximum_overlap_um_n_ge_8']:.4f} | "
            f"{'PASS' if result['pass'] else 'FAIL'} |"
        )
    lines.extend(
        [
            "",
            "Pre-registered population bands were `500..1500` at warmup and",
            "end-run, every post-warmup window-end `N≥400`, and end/warmup",
            "ratio `0.5..1.5`. Mechanics retained I3n thresholds.",
            "",
            "## Gate decision",
            "",
            evidence["decision"],
            "",
            "Implementation note: the first execution failed completeness only",
            "because the inherited no-neighbor minimum-distance value was logged",
            "as `Infinity` at `t=0`. The logger now records the finite censored",
            "`4.5 µm` broadphase bound. All three seeds were rerun; population",
            "traces were "
            + (
                "value-for-value identical."
                if evidence["logging_correction"][
                    "population_traces_identical_to_initial_execution"
                ]
                else "not identical (see evidence)."
            ),
            "",
            "Baseline Job 2 / 3 / 3b / 3c / ChassisPocket I0c checks: "
            + ("PASS." if evidence["baselines"]["pass"] else "FAIL."),
            "",
        ]
    )
    if overall == "PASS":
        lines.extend(
            [
                "A2 is resolved for this joined-clock dish. The next named job is",
                "I3s, which requires its own frozen placement protocol before code.",
                "I3p does not establish map occupancy.",
            ]
        )
    else:
        lines.extend(
            [
                "Stop. Do not tune the clock and do not proceed to I3s or I3n2.",
                "I3n2 is NOT_SCORED for a stable living layer.",
            ]
        )
    STANDING.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    source_freeze = require_freeze()
    previous_evidence = (
        json.loads(EVIDENCE.read_text(encoding="utf-8"))
        if EVIDENCE.exists()
        else None
    )
    seed_results = [score_seed(seed) for seed in SEEDS]
    baselines = run_baselines()
    honesty = source_freeze["pass"] and baselines["pass"]
    passed = honesty and all(result["pass"] for result in seed_results)
    overall = "PASS" if passed else "FAIL"
    if passed:
        decision = (
            "Every seed passed every pre-registered population, monolayer, "
            "packing, morphology, completeness, and clock-identity gate."
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

    previous_populations = {}
    if previous_evidence is not None:
        previous_populations = {
            int(result["seed"]): result["population"]["values"]
            for result in previous_evidence.get("seeds", [])
            if "population" in result and "values" in result["population"]
        }
    replay_population_identical = bool(previous_populations) and all(
        previous_populations.get(result["seed"])
        == result["population"]["values"]
        for result in seed_results
    )
    prior_correction = (
        previous_evidence.get("logging_correction", {})
        if previous_evidence is not None
        else {}
    )
    logging_correction = {
        "initial_execution_outcome": prior_correction.get(
            "initial_execution_outcome",
            previous_evidence.get("overall")
            if previous_evidence is not None
            else None,
        ),
        "initial_failure": prior_correction.get(
            "initial_failure",
            (
                "completeness only: t=0 no-pair d_centers_min_um was Infinity"
                if previous_evidence is not None
                else None
            ),
        ),
        "correction": (
            "finite 4.5 um censored broadphase bound; no dynamics or RNG change"
        ),
        "reran_all_seeds_without_editing_result_csvs": True,
        "population_traces_identical_to_initial_execution": (
            prior_correction.get(
                "population_traces_identical_to_initial_execution",
                replay_population_identical,
            )
        ),
    }

    evidence = {
        "schema": "ChassisHybridmm-I3p-population-clock-v1",
        "overall": overall,
        "question_kind": "population-only; no ridge; no NARMA Overall",
        "protocol": str(PROTOCOL),
        "population_bands": {
            "end_warmup": [POP_WARMUP_LOW, POP_WARMUP_HIGH],
            "end_run": [POP_END_LOW, POP_END_HIGH],
            "post_warmup_window_minimum": POP_WINDOW_MIN,
            "end_to_warmup_ratio": [RATIO_LOW, RATIO_HIGH],
        },
        "source_freeze": source_freeze,
        "seeds": seed_results,
        "baselines": baselines,
        "honesty_pass": honesty,
        "decision": decision,
        "logging_correction": logging_correction,
        "i3n_overall": "FAIL_UNCHANGED_1.1980",
        "narma10b_driven": "UNCHANGED_0.928",
    }
    EVIDENCE.write_text(
        json.dumps(evidence, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_standing(evidence)
    print(f"I3p: {overall}")
    for result in seed_results:
        population = result["population"]
        print(
            f"seed={result['seed']} N_warmup={population['end_warmup']} "
            f"N_end={population['end_run']} "
            f"ratio={population['end_to_warmup_ratio']:.3f} "
            f"result={'PASS' if result['pass'] else 'FAIL'}"
        )
    print(decision)
    print("No NRMSE or NARMA Overall was computed.")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
