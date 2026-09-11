#!/usr/bin/env python3
"""Generate frozen WaveformS sine / square / triangle sequences.

Track E1. Recognizable 16-sample waveforms on the HybridDish claim
dish. Phase is identically 0. Class vector is frozen; do not reshuffle.
Do not regenerate after seeing reservoir AUC. Do not rewrite templates.
Waveform1 GATE_EVIDENCE.md is not edited. This is not Waveform2d.
"""

from __future__ import annotations

import hashlib
import math
import shutil
from pathlib import Path

import numpy as np

CLASS_SEED = 20260818
NUM_BLOCKS = 25
WINDOWS_PER_BLOCK = 16
NUM_WINDOWS = NUM_BLOCKS * WINDOWS_PER_BLOCK
N_CLASSES = 3
CLASS_NAMES = ("sine", "square", "triangle")

WASHOUT_BLOCKS = 3
TRAIN_BLOCKS = 16
TEST_BLOCKS = 6
WASHOUT_COUNTS = (1, 1, 1)
TRAIN_COUNTS = (6, 5, 5)
TEST_COUNTS = (2, 2, 2)

# Frozen 12-decimal templates. Do not change after seeing AUC.
SINE = (
    0.250000000000,
    0.345670858091,
    0.426776695297,
    0.480969883128,
    0.500000000000,
    0.480969883128,
    0.426776695297,
    0.345670858091,
    0.250000000000,
    0.154329141909,
    0.073223304703,
    0.019030116872,
    0.000000000000,
    0.019030116872,
    0.073223304703,
    0.154329141909,
)
SQUARE = (
    0.500000000000,
    0.500000000000,
    0.500000000000,
    0.500000000000,
    0.500000000000,
    0.500000000000,
    0.500000000000,
    0.500000000000,
    0.000000000000,
    0.000000000000,
    0.000000000000,
    0.000000000000,
    0.000000000000,
    0.000000000000,
    0.000000000000,
    0.000000000000,
)
TRIANGLE = (
    0.000000000000,
    0.062500000000,
    0.125000000000,
    0.187500000000,
    0.250000000000,
    0.312500000000,
    0.375000000000,
    0.437500000000,
    0.500000000000,
    0.437500000000,
    0.375000000000,
    0.312500000000,
    0.250000000000,
    0.187500000000,
    0.125000000000,
    0.062500000000,
)
TEMPLATES = (SINE, SQUARE, TRIANGLE)

CLASS_VECTOR = (
    2, 0, 1, 1, 2, 0, 0, 1, 0, 2, 0, 1, 0, 1, 1, 2, 0, 2, 2, 0, 2, 2, 1, 1, 0
)
CLASS_VECTOR_STR = "2,0,1,1,2,0,0,1,0,2,0,1,0,1,1,2,0,2,2,0,2,2,1,1,0"
CLASS_SHA = "40adec9e1e69314a1250a723173d5d877bc4540e15ad17ae918fac735f081a46"
U_SHA = "12941fff63a4b02641803b3c11fdb89cd812832b54a3e06c6fee5f9a4b9d7c92"

HERE = Path(__file__).resolve().parent
WAVEFORM2C = HERE.parent / "BSimReservoirPlanWaveform2c"


def sha256_u_body(values):
    payload = "\n".join(f"{value:.12f}" for value in values) + "\n"
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def sha256_ascii_ints(values):
    payload = ",".join(str(int(v)) for v in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def block_moments(values):
    arr = np.asarray(values, dtype=float)
    mean = float(arr.mean())
    var = float(arr.var(ddof=0))
    sd = math.sqrt(var)
    power = float(np.mean(arr * arr))
    mn = float(arr.min())
    mx = float(arr.max())
    return {
        "mean": mean,
        "variance": var,
        "pop_sd": sd,
        "power": power,
        "min": mn,
        "max": mx,
        "n": int(arr.size),
    }


def assert_frozen_templates():
    for w, expected in enumerate(SINE):
        computed = 0.25 + 0.25 * math.sin(2.0 * math.pi * w / 16.0)
        if abs(round(computed, 12) - expected) > 1e-12:
            raise RuntimeError(f"sine template writer mismatch at w={w}")
    for w, expected in enumerate(SQUARE):
        computed = 0.5 if w < 8 else 0.0
        if abs(round(computed, 12) - expected) > 1e-12:
            raise RuntimeError(f"square template writer mismatch at w={w}")
    for w, expected in enumerate(TRIANGLE):
        computed = 0.5 * (1.0 - abs(2.0 * (w / 16.0) - 1.0))
        if abs(round(computed, 12) - expected) > 1e-12:
            raise RuntimeError(f"triangle template writer mismatch at w={w}")
    for template in TEMPLATES:
        if any(v < 0.0 - 1e-15 or v > 0.5 + 1e-15 for v in template):
            raise RuntimeError("template value outside [0, 0.5]; clip is unnecessary")
    if ",".join(str(v) for v in CLASS_VECTOR) != CLASS_VECTOR_STR:
        raise RuntimeError("class vector string mismatch")
    if sha256_ascii_ints(CLASS_VECTOR) != CLASS_SHA:
        raise RuntimeError("class vector SHA-256 mismatch; do not reshuffle")


def main():
    directory = HERE
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "results").mkdir(parents=True, exist_ok=True)
    assert_frozen_templates()

    if len(CLASS_VECTOR) != NUM_BLOCKS:
        raise RuntimeError(f"class vector length {len(CLASS_VECTOR)}")
    washout = list(CLASS_VECTOR[:WASHOUT_BLOCKS])
    train = list(CLASS_VECTOR[WASHOUT_BLOCKS:WASHOUT_BLOCKS + TRAIN_BLOCKS])
    test = list(CLASS_VECTOR[WASHOUT_BLOCKS + TRAIN_BLOCKS:])
    for name, split, counts in (
        ("washout", washout, WASHOUT_COUNTS),
        ("train", train, TRAIN_COUNTS),
        ("test", test, TEST_COUNTS),
    ):
        got = (split.count(0), split.count(1), split.count(2))
        if got != counts:
            raise RuntimeError(f"{name} class counts {got} != {counts}")

    phases = [0] * NUM_BLOCKS
    u = []
    y_window = []
    for n in range(NUM_WINDOWS):
        block = n // WINDOWS_PER_BLOCK
        w = n % WINDOWS_PER_BLOCK
        cls = CLASS_VECTOR[block]
        y_window.append(cls)
        u.append(round(TEMPLATES[cls][w], 12))

    body = "\n".join(f"{value:.12f}" for value in u) + "\n"
    digest = sha256_u_body(u)
    if digest != U_SHA:
        raise RuntimeError(
            f"u-body SHA-256 {digest} != frozen {U_SHA}. "
            "Do not fix templates. Debug the writer."
        )
    if [round(v, 12) for v in u[:16]] != list(TRIANGLE):
        raise RuntimeError("first 16 values are not the triangle template")

    ahl_path = directory / "input_ahl_waveforms_400.txt"
    ahl_path.write_text(
        "# Frozen WaveformS AHL. sine/square/triangle, 16 samples/period, phase=0.\n"
        f"# Class seed={CLASS_SEED} (vector frozen; do not reshuffle).\n"
        "# First 16 values are triangle (class vector starts with 2).\n"
        "# Do not regenerate after seeing reservoir AUC.\n"
        + body,
        encoding="utf-8",
    )

    acid_src = WAVEFORM2C / "input_acid_held05_400.txt"
    if not acid_src.exists():
        raise RuntimeError(f"missing Waveform2c acid file {acid_src}")
    shutil.copyfile(acid_src, directory / "input_acid_held05_400.txt")

    rows = ["n;block;class;phase;u"]
    rows.extend(
        f"{n};{n // WINDOWS_PER_BLOCK};{y_window[n]};0;{u[n]:.12f}"
        for n in range(NUM_WINDOWS)
    )
    (directory / "waveforms_labels.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    class_vector = ",".join(str(v) for v in CLASS_VECTOR)
    phase_vector = ",".join(str(v) for v in phases)
    schedule = [
        "# Frozen WaveformS (class, phase) schedule. phase identically 0.",
        f"# class_seed={CLASS_SEED} (vector frozen; do not reshuffle)",
        "# class: 0=sine 1=square 2=triangle",
        f"class={class_vector}",
        f"phase={phase_vector}",
        "",
    ]
    (directory / "class_phase_schedule.txt").write_text("\n".join(schedule), encoding="utf-8")

    class_sha = sha256_ascii_ints(CLASS_VECTOR)
    phase_sha = sha256_ascii_ints(phases)
    moment_rows = [block_moments(t) for t in TEMPLATES]

    print(f"blocks={NUM_BLOCKS} windows={NUM_WINDOWS}")
    print(f"class_seed={CLASS_SEED}")
    print("phase=all_zeros")
    print(f"class_vector={class_vector}")
    print(f"phase_vector={phase_vector}")
    print(f"class_vector_sha256={class_sha}")
    print(f"phase_vector_sha256={phase_sha}")
    print(f"u_sha256={digest}")
    print(
        "washout sine={0} square={1} triangle={2}".format(
            washout.count(0), washout.count(1), washout.count(2)
        )
    )
    print(
        "train sine={0} square={1} triangle={2}".format(
            train.count(0), train.count(1), train.count(2)
        )
    )
    print(
        "test sine={0} square={1} triangle={2}".format(
            test.count(0), test.count(1), test.count(2)
        )
    )
    print("BLOCK_MOMENTS")
    print("class mean variance pop_sd power min max n")
    for cls, row in enumerate(moment_rows):
        print(
            f"{CLASS_NAMES[cls]} {row['mean']:.12f} {row['variance']:.12f} "
            f"{row['pop_sd']:.12f} {row['power']:.12f} {row['min']:.12f} "
            f"{row['max']:.12f} {row['n']}"
        )
    print(f"u_min={min(u):.12f} u_max={max(u):.12f}")
    print("first16_triangle=True")
    print("HASH_OK")


if __name__ == "__main__":
    main()
