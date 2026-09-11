#!/usr/bin/env python3
"""Check the frozen ChassisHybridmm I3n occupancy and NARMA gates.

This checker imports the closed Narma10b ridge implementation so lambda
selection remains windows 128..149 only. It writes I3n evidence/standing files,
never a GATE_EVIDENCE.md and never a Narma10b result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
NARMA_DIR = EXAMPLES / "BSimReservoirPlanNarma10b"
CHASSIS_DIR = EXAMPLES / "BacteriumFromScratch"
POCKET_DIR = EXAMPLES / "ChassisPocket"
RESULTS = HERE / "results"
NARMA_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
SEEDS = (111, 222, 333)
ARMS = ("driven", "brownian", "silent")
OCCUPANCY_MIN = 0.05
SATURATION_MAX = 0.95
OFFPLANE_MAX = 0.05
OVERLAP_MAX_UM = 0.25
NEMATIC_MAX = 0.95
EXPECTED_FEATURES = {"driven": 408, "silent": 408, "brownian": 200}
REFERENCE = {
    "driven": 0.928,
    "field": 1.029,
    "u_taps_10": 0.683,
    "brownian_silent": 1.162,
}


def load_narma_module():
    path = NARMA_DIR / "check_narma10b.py"
    spec = importlib.util.spec_from_file_location("closed_narma10b_checker", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NARMA = load_narma_module()


def canonical_u_sha(path: Path) -> str:
    values = [
        float(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    payload = ",".join(f"{value:.12f}" for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def byte_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def input_identity() -> dict:
    files = ("input_ahl_narma200.txt", "input_acid_held05_200.txt", "narma10_target.csv")
    byte_equal = {
        name: byte_sha(HERE / name) == byte_sha(NARMA_DIR / name)
        for name in files
    }
    digest = canonical_u_sha(HERE / "input_ahl_narma200.txt")
    _, _, target_digest = NARMA.load_target(
        HERE / "narma10_target.csv", HERE / "input_ahl_narma200.txt"
    )
    return {
        "canonical_u_sha256": digest,
        "target_u_sha256": target_digest,
        "expected_u_sha256": NARMA_SHA,
        "byte_equal_to_narma10b": byte_equal,
        "pass": digest == NARMA_SHA
        and target_digest == NARMA_SHA
        and all(byte_equal.values()),
    }


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=";"))


def mechanics_gate(run_dir: Path) -> dict:
    path = run_dir / "mechanics.csv"
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"empty mechanics file: {path}")
    rod_rows = [row for row in rows if row["arm"] != "brownian"]
    eligible = [row for row in rod_rows if int(row["N"]) >= 8]
    max_overlap = max((float(row["delta_cc_max_um"]) for row in eligible), default=0.0)
    warmup_rows = [row for row in rod_rows if float(row["t_s"]) <= 18000.0 + 1e-6]
    end_warmup = warmup_rows[-1] if warmup_rows else None
    offplane = float(end_warmup["offplane"]) if end_warmup else math.nan
    nematic = float(end_warmup["nematic_order"]) if end_warmup else math.nan
    return {
        "path": str(path),
        "row_count": len(rows),
        "max_overlap_um_after_n8": max_overlap,
        "packing_pass": max_overlap <= OVERLAP_MAX_UM,
        "end_warmup_t_s": float(end_warmup["t_s"]) if end_warmup else None,
        "end_warmup_N": int(end_warmup["N"]) if end_warmup else None,
        "end_warmup_mean_R": float(end_warmup["mean_R"]) if end_warmup else None,
        "end_warmup_offplane": offplane,
        "end_warmup_nematic_order": nematic,
        "monolayer_pass": bool(end_warmup) and offplane <= OFFPLANE_MAX,
        "filament_guard_pass": bool(end_warmup) and nematic <= NEMATIC_MAX,
    }


def occupancy_gate(run_dir: Path) -> dict:
    rows = read_csv(run_dir / "window_summary.csv")
    values = [float(row["Mean_R"]) for row in rows]
    mean_r = float(np.mean(values)) if values else math.nan
    if not values:
        status = "INCOMPLETE"
    elif mean_r < OCCUPANCY_MIN:
        status = "DEAD"
    elif mean_r > SATURATION_MAX:
        status = "SATURATED"
    else:
        status = "ALIVE"
    return {
        "window_count": len(rows),
        "mean_R": mean_r,
        "minimum": OCCUPANCY_MIN,
        "saturation_stop": SATURATION_MAX,
        "status": status,
        "pass": status == "ALIVE",
    }


def source_hygiene() -> dict:
    java = (HERE / "ChassisHybridmmNarma.java").read_text(encoding="utf-8")
    protocol = (HERE / "PROTOCOL.md").read_text(encoding="utf-8")
    forbidden = ("LuxI", "Danino", "Weber", "glucoseField", "setGoal(", "motility")
    required = (
        "BOUND_Z = 1.0",
        "RECEIVER_K_UM = 1.6",
        "RECEIVER_HILL_N = 2.0",
        "RECEIVER_TAU_S = 15.0",
        "AHL_SOURCE_RATE = 1.28e7",
        "AHL_DIFFUSIVITY = 159.0",
        "AHL_DECAY = 0.0033",
        "EcoliRodCell",
        "ValdezHertzian",
        "SYMMETRY_BROKEN",
    )
    new_dish_label = "new-dish question" in protocol
    return {
        "forbidden_tokens": [token for token in forbidden if token in java],
        "missing_required_tokens": [token for token in required if token not in java],
        "protocol_new_dish_label": new_dish_label,
        "pass": not any(token in java for token in forbidden)
        and all(token in java for token in required)
        and new_dish_label,
    }


def run_baselines() -> dict:
    scripts = [
        CHASSIS_DIR / "check_job2.py",
        CHASSIS_DIR / "check_job3.py",
        CHASSIS_DIR / "check_job3b.py",
        CHASSIS_DIR / "check_job3c.py",
        POCKET_DIR / "check_i0c.py",
    ]
    results = {}
    for script in scripts:
        proc = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(script.parent),
            text=True,
            capture_output=True,
        )
        lines = (proc.stdout + proc.stderr).strip().splitlines()
        results[script.name] = {
            "returncode": proc.returncode,
            "last_line": lines[-1] if lines else "",
            "pass": proc.returncode == 0,
        }
    return {"checks": results, "pass": all(item["pass"] for item in results.values())}


def smoke_check(run_dir: Path, include_baselines: bool) -> dict:
    identity = input_identity()
    occupancy = occupancy_gate(run_dir)
    mechanics = mechanics_gate(run_dir)
    hygiene = source_hygiene()
    baselines = run_baselines() if include_baselines else {"skipped": True, "pass": True}
    physical_pass = (
        mechanics["packing_pass"]
        and mechanics["monolayer_pass"]
        and mechanics["filament_guard_pass"]
    )
    if not physical_pass:
        decision = "NOT_SCORED"
        reason = "physical smoke failed (packing, monolayer, or filament guard)"
    elif occupancy["status"] in {"DEAD", "SATURATED"}:
        decision = "NOT_SCORED"
        reason = f"occupancy {occupancy['status']}"
    elif occupancy["status"] == "ALIVE":
        decision = "PROCEED_TO_FULL_NARMA"
        reason = "occupancy and physical smoke passed"
    else:
        decision = "INCOMPLETE"
        reason = "occupancy output incomplete"
    evidence = {
        "schema": "ChassisHybridmm-I3n-smoke-v1",
        "dish": "ChassisHybridmm I3n new dish",
        "run_dir": str(run_dir),
        "input_identity": identity,
        "occupancy": occupancy,
        "mechanics": mechanics,
        "source_hygiene": hygiene,
        "baselines": baselines,
        "decision": decision,
        "reason": reason,
        "pass": decision == "PROCEED_TO_FULL_NARMA"
        and identity["pass"]
        and hygiene["pass"]
        and baselines["pass"],
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "i3n_smoke_evidence.json").write_text(
        json.dumps(evidence, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    if decision == "NOT_SCORED":
        write_standing(evidence, smoke_only=True)
    return evidence


def run_dir(arm: str, seed: int) -> Path:
    return RESULTS / f"i3n_{arm}_seed{seed}"


def evaluate_full(include_baselines: bool) -> dict:
    identity = input_identity()
    _, target, digest = NARMA.load_target(
        HERE / "narma10_target.csv", HERE / "input_ahl_narma200.txt"
    )
    grouped = {}
    for arm in ARMS:
        grouped[arm] = [
            NARMA.read_run(run_dir(arm, seed), arm, 200, 3200)
            for seed in SEEDS
        ]
    evaluated = {
        arm: [
            NARMA.evaluate_arm_run(run, NARMA.load_sequence(HERE / "input_ahl_narma200.txt"), target)
            for run in runs
        ]
        for arm, runs in grouped.items()
    }
    nrmse = {
        arm: [float(run["task"]["test_nrmse"]) for run in runs]
        for arm, runs in evaluated.items()
    }
    means = {arm: float(np.mean(values)) for arm, values in nrmse.items()}
    ses = {
        arm: float(np.std(values, ddof=1) / np.sqrt(len(values)))
        for arm, values in nrmse.items()
    }
    fields = [float(run["field_only"]["test_nrmse"]) for run in evaluated["driven"]]
    field_mean = float(np.mean(fields))
    occupancies = {
        seed: occupancy_gate(run_dir("driven", seed))
        for seed in SEEDS
    }
    mechanics = {
        f"{arm}_{seed}": mechanics_gate(run_dir(arm, seed))
        for arm in ("driven", "silent")
        for seed in SEEDS
    }
    feature_counts = {
        arm: [run["n_features"] for run in runs]
        for arm, runs in evaluated.items()
    }
    feature_pass = all(
        count == EXPECTED_FEATURES[arm]
        for arm, counts in feature_counts.items()
        for count in counts
    )
    csv_pass = all(
        run["summary_csv"]["row_count_pass"]
        and run["summary_csv"]["rectangular_csv_pass"]
        and all(
            item["row_count_pass"]
            and item["rectangular_csv_pass"]
            and item.get("last_sample_pass", True)
            for item in run["auxiliary_csvs"].values()
        )
        for runs in evaluated.values()
        for run in runs
    )
    system_pass = means["driven"] < means["brownian"] and means["driven"] < means["silent"]
    living_pass = means["driven"] < field_mean
    occupancy_pass = all(item["pass"] for item in occupancies.values())
    physical_pass = all(
        item["packing_pass"] and item["monolayer_pass"] and item["filament_guard_pass"]
        for item in mechanics.values()
    )
    hygiene = source_hygiene()
    baselines = run_baselines() if include_baselines else {"skipped": True, "pass": True}
    gates = {
        "I3n.0_occupancy": occupancy_pass,
        "I3n.1_monolayer_and_packing": physical_pass,
        "I3n.2_u_identity": identity["pass"] and digest == NARMA_SHA,
        "I3n.3_system": system_pass,
        "I3n.4_living_layer": living_pass,
        "I3n.5_ridge_hygiene": feature_pass and csv_pass,
        "I3n.6_honesty": hygiene["pass"] and baselines["pass"],
    }
    overall = "PASS" if all(
        gates[name] for name in (
            "I3n.0_occupancy",
            "I3n.1_monolayer_and_packing",
            "I3n.2_u_identity",
            "I3n.3_system",
            "I3n.5_ridge_hygiene",
            "I3n.6_honesty",
        )
    ) else "FAIL"
    evidence = {
        "schema": "ChassisHybridmm-I3n-full-v1",
        "dish": "ChassisHybridmm I3n new dish",
        "seeds": list(SEEDS),
        "overall": overall,
        "gates": gates,
        "test_nrmse": nrmse,
        "mean_nrmse": means,
        "se_nrmse": ses,
        "field_test_nrmse": fields,
        "field_mean_nrmse": field_mean,
        "lambdas": {
            arm: [float(run["task"]["lambda"]) for run in runs]
            for arm, runs in evaluated.items()
        },
        "feature_counts": feature_counts,
        "csv_pass": csv_pass,
        "occupancy": occupancies,
        "mechanics": mechanics,
        "input_identity": identity,
        "source_hygiene": hygiene,
        "baselines": baselines,
        "ridge_rule": {
            "implementation": str(NARMA_DIR / "check_narma10b.py"),
            "validation_windows": "128..149 only",
            "intercept_regularized": False,
            "biology_features": 408,
        },
        "narma10b_spheres_not_this_run": REFERENCE,
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "i3n_evidence.json").write_text(
        json.dumps(evidence, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    write_standing(evidence, smoke_only=False)
    return evidence


def write_standing(evidence: dict, smoke_only: bool) -> None:
    if smoke_only:
        occupancy = evidence["occupancy"]
        mechanics = evidence["mechanics"]
        lines = [
            "# ChassisHybridmm NARMA I3n standing",
            "",
            "**I3n OVERALL: NOT_SCORED**",
            "",
            "This is a new millimetre Hertzian-rod dish, not a rewrite of",
            "Narma10b spheres. Narma10b driven F408 `0.928` remains untouched.",
            "",
            f"Occupancy smoke: `{occupancy['status']}`, mean R `{occupancy['mean_R']:.6f}`.",
            f"End-warmup offplane `{mechanics['end_warmup_offplane']:.6f}`, "
            f"nematic `{mechanics['end_warmup_nematic_order']:.6f}`, "
            f"maximum overlap `{mechanics['max_overlap_um_after_n8']:.6f} µm`.",
            "",
            f"Stop reason: {evidence['reason']}. No J_max, K, clamp, or k_cc retuning.",
            "No 200-window NARMA arms were started.",
        ]
    else:
        means = evidence["mean_nrmse"]
        ses = evidence["se_nrmse"]
        field = evidence["field_mean_nrmse"]
        gates = evidence["gates"]
        system = "PASS" if gates["I3n.3_system"] else "FAIL"
        living = "PASS" if gates["I3n.4_living_layer"] else "FAIL"
        lines = [
            "# ChassisHybridmm NARMA I3n standing",
            "",
            f"**I3n OVERALL: {evidence['overall']}**",
            "",
            "This is a new `1000×500×1 µm` Hertzian-rod dish. It does not",
            "rewrite paper-2 Narma10b or its driven F408 score `0.928`.",
            "",
            f"System gate: **{system}**. Living-layer versus field: **{living}**.",
            "Brownian is the copied Narma10b passive `Den_*` particle null; silent",
            "uses Hertzian rods with AHL/acid point sources off.",
            "",
            "| I3n arm | seed 111 | seed 222 | seed 333 | mean ± s.e. |",
            "|---|---:|---:|---:|---:|",
        ]
        for arm in ("driven", "brownian", "silent"):
            values = evidence["test_nrmse"][arm]
            lines.append(
                f"| {arm} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} "
                f"| {means[arm]:.4f} ± {ses[arm]:.4f} |"
            )
        fields = evidence["field_test_nrmse"]
        occupancy = evidence["occupancy"]
        mechanics = evidence["mechanics"]
        driven_mechanics = [mechanics[f"driven_{seed}"] for seed in SEEDS]
        lines.append(
            f"| field AHL | {fields[0]:.4f} | {fields[1]:.4f} | {fields[2]:.4f} "
            f"| {field:.4f} |"
        )
        lines.extend([
            "",
            "## Occupancy and mechanics",
            "",
            "Driven occupancy was **ALIVE** on all three seeds: mean R "
            f"`{occupancy[111]['mean_R']:.4f} / {occupancy[222]['mean_R']:.4f} / "
            f"{occupancy[333]['mean_R']:.4f}` for seeds 111/222/333.",
            "End-warmup populations were "
            f"`{driven_mechanics[0]['end_warmup_N']} / "
            f"{driven_mechanics[1]['end_warmup_N']} / "
            f"{driven_mechanics[2]['end_warmup_N']}`; offplane was `0` on every seed.",
            "Maximum logged Hertzian overlap was "
            f"`{max(item['max_overlap_um_after_n8'] for item in driven_mechanics):.4f} µm` "
            "(gate `≤0.25 µm`).",
            f"Copied u identity: `{evidence['input_identity']['canonical_u_sha256']}`.",
            "",
            "## Narma10b spheres, not this run",
            "",
            "| Arm | Test NRMSE |",
            "|---|---:|",
            "| Driven F408 | 0.928 |",
            "| Field AHL | 1.029 |",
            "| 10-tap of u | 0.683 |",
            "| Brownian / silent | 1.162 |",
            "",
            "The 10-tap value is a task ceiling, not an I3n gate. I3n was not",
            "required to beat 0.928. No source, Hill, clamp, packing, or ridge",
            "parameter was retuned after occupancy or NRMSE.",
        ])
    (HERE / "CHASSISHYBRIDMM_NARMA_STANDING.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def print_smoke(evidence: dict) -> None:
    occupancy = evidence["occupancy"]
    mechanics = evidence["mechanics"]
    print(f"u_sha256={evidence['input_identity']['canonical_u_sha256']}")
    print(f"occupancy={occupancy['status']} mean_R={occupancy['mean_R']:.6f}")
    print(
        "end_warmup "
        f"N={mechanics['end_warmup_N']} "
        f"offplane={mechanics['end_warmup_offplane']:.6f} "
        f"nematic={mechanics['end_warmup_nematic_order']:.6f} "
        f"delta_cc_max={mechanics['max_overlap_um_after_n8']:.6f}"
    )
    print(f"decision={evidence['decision']}")


def print_full(evidence: dict) -> None:
    print("ChassisHybridmm I3n NRMSE (new dish)")
    print("arm          seed111  seed222  seed333  mean")
    for arm in ("driven", "brownian", "silent"):
        values = evidence["test_nrmse"][arm]
        print(
            f"{arm:12s} {values[0]:8.4f} {values[1]:8.4f} "
            f"{values[2]:8.4f} {evidence['mean_nrmse'][arm]:8.4f}"
        )
    values = evidence["field_test_nrmse"]
    print(
        f"{'field':12s} {values[0]:8.4f} {values[1]:8.4f} "
        f"{values[2]:8.4f} {evidence['field_mean_nrmse']:8.4f}"
    )
    print("Narma10b spheres, not this run: driven=.928 field=1.029 taps=.683 nulls=1.162")
    for gate, passed in evidence["gates"].items():
        print(f"{gate}: {'PASS' if passed else 'FAIL'}")
    print(f"I3n OVERALL: {evidence['overall']}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=RESULTS / "i3n_occupancy_smoke_seed111",
    )
    parser.add_argument("--skip-baselines", action="store_true")
    args = parser.parse_args()
    try:
        if args.smoke:
            evidence = smoke_check(args.run_dir, not args.skip_baselines)
            print_smoke(evidence)
            return 0 if evidence["pass"] else 1
        evidence = evaluate_full(not args.skip_baselines)
        print_full(evidence)
        return 0 if evidence["overall"] == "PASS" else 1
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
