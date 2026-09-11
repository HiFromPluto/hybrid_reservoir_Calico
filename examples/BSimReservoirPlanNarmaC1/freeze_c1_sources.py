#!/usr/bin/env python3
"""Append on-disk SHA-256 of C1 sources to U_FREEZE.md before any Java job."""

from __future__ import annotations

import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FREEZE = RESULTS / "U_FREEZE.md"

SOURCE_GLOBS = [
    "BSimReservoirPlanNarmaC1.java",
    "VoxelAnalyzer.java",
    "generate_c1_inputs.py",
    "freeze_c1_sources.py",
    "run_c1_jobs.py",
    "check_c1.py",
    "plot_c1.py",
    "PROTOCOL.md",
    "README.md",
    "sim_config_c1.properties",
    "sim_config_c1_driven_traj*_seed*.properties",
    "input_ahl_narma200_traj*.txt",
    "narma10_target_traj*.csv",
    "input_acid_held05_200.txt",
]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if not FREEZE.exists():
        raise SystemExit("ABORT: run generate_c1_inputs.py first")
    text = FREEZE.read_text(encoding="utf-8")
    marker = "## Package source hashes"
    if marker in text:
        text = text.split(marker)[0].rstrip() + "\n\n"
    paths = []
    for pattern in SOURCE_GLOBS:
        paths.extend(sorted(HERE.glob(pattern)))
    seen = []
    lines = [marker, "", "Recorded before the first new Java job.", "", "| File | SHA-256 |", "|---|---|"]
    for path in paths:
        rel = path.relative_to(HERE).as_posix()
        if rel in seen:
            continue
        seen.append(rel)
        lines.append(f"| {rel} | `{sha256_file(path)}` |")
    narma_u = HERE.parent / "BSimReservoirPlanNarma10b" / "input_ahl_narma200.txt"
    lines.append(
        f"| ../BSimReservoirPlanNarma10b/input_ahl_narma200.txt (traj 00) | `{sha256_file(narma_u)}` |"
    )
    lines.append("")
    FREEZE.write_text(text + "\n".join(lines) + "\n", encoding="utf-8")
    print(f"hashed {len(seen) + 1} source files into {FREEZE}")


if __name__ == "__main__":
    main()
