#!/usr/bin/env python3
"""Generate the frozen waveform-classification AHL sequence and labels.

This script is the sequence provenance. Re-running it with the same seed
must reproduce the committed input files. Do not retune templates, class
counts, or the shuffle seed after seeing AUC.
"""

from __future__ import annotations

import hashlib
import math
import random
import shutil
from pathlib import Path

LABEL_SEED = 20260822
NUM_BLOCKS = 40
WINDOWS_PER_BLOCK = 5
NUM_WINDOWS = NUM_BLOCKS * WINDOWS_PER_BLOCK
CLASS_SINE, CLASS_SQUARE, CLASS_TRIANGLE = 0, 1, 2
CLASS_NAMES = {0: "sine", 1: "square", 2: "triangle"}
TRIANGLE = (0.05, 0.25, 0.45, 0.25, 0.05)

# Forced class counts inside each split, then shuffled.
WASHOUT_COUNTS = (3, 3, 2)  # sine, square, triangle
TRAIN_COUNTS = (8, 7, 7)
TEST_COUNTS = (3, 3, 4)

HERE = Path(__file__).resolve().parent
STAGE6_ACID = HERE.parent / "BSimReservoirPlanStage6" / "input_acid_held05_200.txt"


def sha256_u(values):
    payload = ",".join(f"{value:.12f}" for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def clip_u(value):
    return min(0.5, max(0.0, value))


def template_u(cls, w):
    if cls == CLASS_SINE:
        return clip_u(0.25 + 0.25 * math.sin(2.0 * math.pi * w / 5.0))
    if cls == CLASS_SQUARE:
        return clip_u(0.45 if w < 3 else 0.05)
    if cls == CLASS_TRIANGLE:
        return clip_u(TRIANGLE[w])
    raise ValueError(f"unknown class {cls}")


def forced_then_shuffled(rng, counts):
    labels = (
        [CLASS_SINE] * counts[0]
        + [CLASS_SQUARE] * counts[1]
        + [CLASS_TRIANGLE] * counts[2]
    )
    rng.shuffle(labels)
    return labels


def main():
    directory = HERE
    rng = random.Random(LABEL_SEED)
    block_y = (
        forced_then_shuffled(rng, WASHOUT_COUNTS)
        + forced_then_shuffled(rng, TRAIN_COUNTS)
        + forced_then_shuffled(rng, TEST_COUNTS)
    )
    if len(block_y) != NUM_BLOCKS:
        raise RuntimeError(f"class vector length {len(block_y)}")
    washout, train, test = block_y[:8], block_y[8:30], block_y[30:]
    for name, split, counts in (
        ("washout", washout, WASHOUT_COUNTS),
        ("train", train, TRAIN_COUNTS),
        ("test", test, TEST_COUNTS),
    ):
        got = (split.count(0), split.count(1), split.count(2))
        if got != counts:
            raise RuntimeError(f"{name} class counts {got} != {counts}")

    class_vector = ",".join(str(v) for v in block_y)
    class_sha = hashlib.sha256(class_vector.encode("ascii")).hexdigest()

    u = []
    y_window = []
    for n in range(NUM_WINDOWS):
        block = n // WINDOWS_PER_BLOCK
        w = n % WINDOWS_PER_BLOCK
        cls = block_y[block]
        y_window.append(cls)
        u.append(template_u(cls, w))

    ahl_path = directory / "input_ahl_waveform200.txt"
    ahl_path.write_text(
        "# Frozen waveform-classification AHL. Deterministic sine/square/triangle\n"
        "# templates on 40 blocks x 5 windows. Label seed=20260822. u clipped to [0, 0.5].\n"
        "# Do not regenerate after seeing AUC.\n"
        + "\n".join(f"{value:.12f}" for value in u)
        + "\n",
        encoding="utf-8",
    )
    shutil.copyfile(STAGE6_ACID, directory / "input_acid_held05_200.txt")

    rows = ["n;block;y;u"]
    rows.extend(
        f"{n};{n // WINDOWS_PER_BLOCK};{y_window[n]};{u[n]:.12f}"
        for n in range(NUM_WINDOWS)
    )
    (directory / "waveform_labels.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    digest = sha256_u(u)
    print(f"blocks={NUM_BLOCKS} windows={NUM_WINDOWS}")
    print(f"label_seed={LABEL_SEED}")
    print(f"class_vector={class_vector}")
    print(f"class_vector_sha256={class_sha}")
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
    print(f"u_sha256={digest}")
    print(f"u_min={min(u):.12f} u_max={max(u):.12f}")
    print("templates: 0=sine 1=square 2=triangle")


if __name__ == "__main__":
    main()
