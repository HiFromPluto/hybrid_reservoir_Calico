#!/usr/bin/env python3
"""Generate frozen Waveform2c aligned equal-histogram sequences.

Track E1. Waveform2 remains TASK_VOID (shared orbit). Waveform2b remains
the disjoint-orbit / random-phase negative (RAW_U_8 0.5298). Waveform2c
uses canonical representatives that all start at MID, with phase p=0 on
every block. Class shuffle uses random.Random(20260819). No phase subseed.
"""

from __future__ import annotations

import hashlib
import math
import random
from pathlib import Path

import numpy as np

CLASS_SEED = 20260819
NUM_BLOCKS = 50
WINDOWS_PER_BLOCK = 8
NUM_WINDOWS = NUM_BLOCKS * WINDOWS_PER_BLOCK
N_CLASSES = 3
CLASS_NAMES = ("ORDER_UNI", "ORDER_DOWN", "ORDER_UP")

WASHOUT_BLOCKS = 8
TRAIN_BLOCKS = 26
TEST_BLOCKS = 16
WASHOUT_COUNTS = (3, 3, 2)
TRAIN_COUNTS = (9, 9, 8)
TEST_COUNTS = (6, 5, 5)

S2 = 0.25 * math.sqrt(2.0) / 2.0
LO = round(0.25 - S2, 12)
MID = round(0.25, 12)
HI = round(0.25 + S2, 12)
MINV = round(0.0, 12)
MAXV = round(0.5, 12)

assert LO == 0.073223304703
assert MID == 0.250000000000
assert HI == 0.426776695297
assert MINV == 0.000000000000
assert MAXV == 0.500000000000

SHARED_MULTISET = (MINV, LO, LO, MID, MID, HI, HI, MAXV)

# Canonical representatives, all starting at MID, phase identically 0.
# ORDER_DOWN / ORDER_UP are Waveform2b DOWN/UP rotated by 3.
TEMPLATES = (
    (MID, HI, MAXV, HI, MID, LO, MINV, LO),
    (MID, MID, LO, LO, MINV, MAXV, HI, HI),
    (MID, MID, HI, HI, MAXV, MINV, LO, LO),
)

HERE = Path(__file__).resolve().parent


def sha256_u(values):
    payload = ",".join(f"{value:.12f}" for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def sha256_ascii_ints(values):
    payload = ",".join(str(int(v)) for v in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def clip_u(value):
    return min(MAXV, max(MINV, value))


def template_u(cls, w):
    return clip_u(TEMPLATES[cls][w])


def sorted_multiset(values):
    return tuple(sorted(round(v, 12) for v in values))


def cyclic_shifts(seq):
    seq = tuple(seq)
    return {tuple(seq[p:] + seq[:p]) for p in range(len(seq))}


def orbit_intersection_matrix(templates):
    orbits = [cyclic_shifts(t) for t in templates]
    matrix = [[len(orbits[i] & orbits[j]) for j in range(len(templates))]
              for i in range(len(templates))]
    return orbits, matrix


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
        "range": mx - mn,
        "n": int(arr.size),
    }


def forced_then_shuffled(rng, counts):
    labels = [0] * counts[0] + [1] * counts[1] + [2] * counts[2]
    rng.shuffle(labels)
    return labels


def prehash_gates():
    for cls, template in enumerate(TEMPLATES):
        if sorted_multiset(template) != sorted_multiset(SHARED_MULTISET):
            raise RuntimeError(f"{CLASS_NAMES[cls]} is not a permutation of the shared multiset")
        if round(template[0], 12) != MID:
            raise RuntimeError(f"{CLASS_NAMES[cls]} does not start at MID")
    if len(set(TEMPLATES)) != N_CLASSES:
        raise RuntimeError("canonical 8-vectors are not pairwise unequal")

    orbits, matrix = orbit_intersection_matrix(TEMPLATES)
    print("ORBIT_INTERSECTION_MATRIX")
    header = "         " + " ".join(f"{name:>12}" for name in CLASS_NAMES)
    print(header)
    for i, name in enumerate(CLASS_NAMES):
        cells = " ".join(f"{matrix[i][j]:12d}" for j in range(N_CLASSES))
        print(f"{name:>8} {cells}")
    for i in range(N_CLASSES):
        for j in range(i + 1, N_CLASSES):
            if orbits[i] & orbits[j]:
                raise RuntimeError(
                    f"ORBIT_GATE_FAIL {CLASS_NAMES[i]} ∩ {CLASS_NAMES[j]} nonempty"
                )
    print("ORBIT_GATE_PASS")

    moment_rows = [block_moments(t) for t in TEMPLATES]
    for row in moment_rows[1:]:
        for key in ("mean", "variance", "pop_sd", "power", "min", "max", "range"):
            if abs(row[key] - moment_rows[0][key]) > 1e-12:
                raise RuntimeError(f"block moment {key} differs across classes")
    print("MOMENTS_MATCH_1e-12")
    return moment_rows


def main():
    directory = HERE
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "results").mkdir(parents=True, exist_ok=True)

    moment_rows = prehash_gates()

    rng = random.Random(CLASS_SEED)
    block_y = (
        forced_then_shuffled(rng, WASHOUT_COUNTS)
        + forced_then_shuffled(rng, TRAIN_COUNTS)
        + forced_then_shuffled(rng, TEST_COUNTS)
    )
    if len(block_y) != NUM_BLOCKS:
        raise RuntimeError(f"class vector length {len(block_y)}")
    washout = block_y[:WASHOUT_BLOCKS]
    train = block_y[WASHOUT_BLOCKS:WASHOUT_BLOCKS + TRAIN_BLOCKS]
    test = block_y[WASHOUT_BLOCKS + TRAIN_BLOCKS:]
    for name, split, counts in (
        ("washout", washout, WASHOUT_COUNTS),
        ("train", train, TRAIN_COUNTS),
        ("test", test, TEST_COUNTS),
    ):
        got = (split.count(0), split.count(1), split.count(2))
        if got != counts:
            raise RuntimeError(f"{name} class counts {got} != {counts}")

    phases = [0] * NUM_BLOCKS
    if any(p != 0 for p in phases):
        raise RuntimeError("phase column is not identically 0")

    u = []
    y_window = []
    for n in range(NUM_WINDOWS):
        block = n // WINDOWS_PER_BLOCK
        w = n % WINDOWS_PER_BLOCK
        cls = block_y[block]
        y_window.append(cls)
        u.append(round(template_u(cls, w), 12))

    ahl_path = directory / "input_ahl_waveform2c_400.txt"
    ahl_path.write_text(
        "# Frozen Waveform2c AHL. Aligned equal-histogram templates. phase=0.\n"
        "# Classes ORDER_UNI/ORDER_DOWN/ORDER_UP. All start at MID. No extra noise.\n"
        f"# Class seed={CLASS_SEED}. No phase subseed. u clipped to [0, 0.5].\n"
        "# Do not regenerate after seeing reservoir AUC.\n"
        + "\n".join(f"{value:.12f}" for value in u)
        + "\n",
        encoding="utf-8",
    )
    (directory / "input_acid_held05_400.txt").write_text(
        "# Frozen Waveform2c acid hold. Every analysis window is 0.5. Off in warmup (Java).\n"
        + "\n".join("0.50" for _ in range(NUM_WINDOWS))
        + "\n",
        encoding="utf-8",
    )

    rows = ["n;block;class;phase;u"]
    rows.extend(
        f"{n};{n // WINDOWS_PER_BLOCK};{y_window[n]};0;{u[n]:.12f}"
        for n in range(NUM_WINDOWS)
    )
    (directory / "waveform2c_labels.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    class_vector = ",".join(str(v) for v in block_y)
    phase_vector = ",".join(str(v) for v in phases)
    schedule = [
        "# Frozen Waveform2c (class, phase) schedule. phase identically 0.",
        f"# class_seed={CLASS_SEED}",
        "# class: 0=ORDER_UNI 1=ORDER_DOWN 2=ORDER_UP",
        f"class={class_vector}",
        f"phase={phase_vector}",
        "",
    ]
    (directory / "class_phase_schedule.txt").write_text("\n".join(schedule), encoding="utf-8")

    class_sha = sha256_ascii_ints(block_y)
    phase_sha = sha256_ascii_ints(phases)
    digest = sha256_u(u)

    print(f"blocks={NUM_BLOCKS} windows={NUM_WINDOWS}")
    print(f"class_seed={CLASS_SEED}")
    print("phase=all_zeros")
    print(f"class_vector={class_vector}")
    print(f"phase_vector={phase_vector}")
    print(f"class_vector_sha256={class_sha}")
    print(f"phase_vector_sha256={phase_sha}")
    print(f"u_sha256={digest}")
    print(
        "washout ORDER_UNI={0} ORDER_DOWN={1} ORDER_UP={2}".format(
            washout.count(0), washout.count(1), washout.count(2)
        )
    )
    print(
        "train ORDER_UNI={0} ORDER_DOWN={1} ORDER_UP={2}".format(
            train.count(0), train.count(1), train.count(2)
        )
    )
    print(
        "test ORDER_UNI={0} ORDER_DOWN={1} ORDER_UP={2}".format(
            test.count(0), test.count(1), test.count(2)
        )
    )
    print("BLOCK_MOMENTS")
    print("class mean variance pop_sd power min max range n")
    for cls, row in enumerate(moment_rows):
        print(
            f"{CLASS_NAMES[cls]} {row['mean']:.12f} {row['variance']:.12f} "
            f"{row['pop_sd']:.12f} {row['power']:.12f} {row['min']:.12f} "
            f"{row['max']:.12f} {row['range']:.12f} {row['n']}"
        )
    print(f"u_min={min(u):.12f} u_max={max(u):.12f}")
    print("all_start_MID=True")


if __name__ == "__main__":
    main()
