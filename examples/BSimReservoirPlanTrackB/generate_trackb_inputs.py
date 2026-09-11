#!/usr/bin/env python3
"""Generate the frozen Track B 5-channel patient classification inputs.

This script is the sequence provenance. Re-running it with the same seeds
must reproduce the committed input files. Do not retune HIGH/LOW, the
product structure, or the seeds after seeing AUC.
"""

from __future__ import annotations

import hashlib
import random
import shutil
from pathlib import Path

LABEL_SEED = 20260820
TRACE_SEED = 20260821
NUM_PATIENTS = 40
WINDOWS_PER_PATIENT = 5
NUM_WINDOWS = NUM_PATIENTS * WINDOWS_PER_PATIENT
HIGH_LOW = (0.30, 0.50)
LOW_RANGE = (0.00, 0.20)
DIST_RANGE = (0.00, 0.50)
WASHOUT_POS, WASHOUT_NEG = 4, 4
TRAIN_POS, TRAIN_NEG = 11, 11
TEST_POS, TEST_NEG = 5, 5

HERE = Path(__file__).resolve().parent
STAGE6_ACID = HERE.parent / "BSimReservoirPlanStage6" / "input_acid_held05_200.txt"


def sha256_csv(values, fmt):
    payload = ",".join(fmt.format(value) for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def sha256_u(u):
    return sha256_csv(u, "{:.12f}")


def balanced_block(rng, n_pos, n_neg):
    labels = [1] * n_pos + [0] * n_neg
    rng.shuffle(labels)
    return labels


def uniform(rng, lo, hi):
    return lo + (hi - lo) * rng.random()


def main():
    directory = HERE
    label_rng = random.Random(LABEL_SEED)
    patient_y = (
        balanced_block(label_rng, WASHOUT_POS, WASHOUT_NEG)
        + balanced_block(label_rng, TRAIN_POS, TRAIN_NEG)
        + balanced_block(label_rng, TEST_POS, TEST_NEG)
    )
    if len(patient_y) != NUM_PATIENTS:
        raise RuntimeError(f"class vector length {len(patient_y)}")
    washout = patient_y[:8]
    train = patient_y[8:30]
    test = patient_y[30:]
    if sum(washout) != WASHOUT_POS or sum(train) != TRAIN_POS or sum(test) != TEST_POS:
        raise RuntimeError("class balance inside a split is wrong")

    class_vector = ",".join(str(v) for v in patient_y)
    class_sha = hashlib.sha256(class_vector.encode("ascii")).hexdigest()

    trace_rng = random.Random(TRACE_SEED)
    u_ahl, u_acid, u_atta, u_attb, u_rep, y_window = [], [], [], [], [], []
    for n in range(NUM_WINDOWS):
        patient = n // WINDOWS_PER_PATIENT
        y = patient_y[patient]
        y_window.append(y)
        # Frozen draw order per window, declared in PROTOCOL.md:
        # DIST AHL, DIST attB, DIST rep, then class-conditional attA/acid.
        u_ahl.append(uniform(trace_rng, *DIST_RANGE))
        u_attb.append(uniform(trace_rng, *DIST_RANGE))
        u_rep.append(uniform(trace_rng, *DIST_RANGE))
        if y == 1:
            u_atta.append(uniform(trace_rng, *HIGH_LOW))
            u_acid.append(uniform(trace_rng, *HIGH_LOW))
        elif trace_rng.random() < 0.5:
            u_atta.append(uniform(trace_rng, *HIGH_LOW))
            u_acid.append(uniform(trace_rng, *LOW_RANGE))
        else:
            u_atta.append(uniform(trace_rng, *LOW_RANGE))
            u_acid.append(uniform(trace_rng, *HIGH_LOW))

    u_singlesite = []
    for n in range(NUM_WINDOWS):
        mean = (u_ahl[n] + u_acid[n] + u_atta[n] + u_attb[n] + u_rep[n]) / 5.0
        u_singlesite.append(min(0.5, max(0.0, mean)))

    files = {
        "input_ahl_trackb200.txt": (
            "# Frozen Track B AHL distractor. DIST ~ Uniform[0.00, 0.50], "
            "trace seed=20260821. Independent of class.\n"
            "# Do not regenerate after seeing AUC.\n",
            u_ahl,
        ),
        "input_acid_trackb200.txt": (
            "# Frozen Track B acid biomarker. Class 1: HIGH; class 0: HIGH xor attA. "
            "Trace seed=20260821. Free input on the driven 5-channel arm.\n"
            "# Do not regenerate after seeing AUC.\n",
            u_acid,
        ),
        "input_atta_trackb200.txt": (
            "# Frozen Track B attA biomarker. Class 1: HIGH; class 0: HIGH xor acid. "
            "Trace seed=20260821.\n"
            "# Do not regenerate after seeing AUC.\n",
            u_atta,
        ),
        "input_attb_trackb200.txt": (
            "# Frozen Track B attB distractor. DIST ~ Uniform[0.00, 0.50], "
            "trace seed=20260821. Independent of class.\n"
            "# Do not regenerate after seeing AUC.\n",
            u_attb,
        ),
        "input_rep_trackb200.txt": (
            "# Frozen Track B repellent distractor. DIST ~ Uniform[0.00, 0.50], "
            "trace seed=20260821. Independent of class.\n"
            "# Do not regenerate after seeing AUC.\n",
            u_rep,
        ),
        "input_ahl_singlesite200.txt": (
            "# Frozen Track B single-site AHL control. Per-window mean of the five "
            "u channels, clipped to [0, 0.5].\n"
            "# H5 control. Acid held at 0.5. attA=attB=rep=0. Do not regenerate after AUC.\n",
            u_singlesite,
        ),
    }
    for name, (header, values) in files.items():
        (directory / name).write_text(
            header + "\n".join(f"{value:.12f}" for value in values) + "\n",
            encoding="utf-8",
        )

    shutil.copyfile(STAGE6_ACID, directory / "input_acid_held05_200.txt")

    rows = ["n;patient;y;u_ahl;u_acid;u_atta;u_attb;u_rep"]
    rows.extend(
        f"{n};{n // WINDOWS_PER_PATIENT};{y_window[n]};"
        f"{u_ahl[n]:.12f};{u_acid[n]:.12f};{u_atta[n]:.12f};"
        f"{u_attb[n]:.12f};{u_rep[n]:.12f}"
        for n in range(NUM_WINDOWS)
    )
    (directory / "trackb_labels.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"patients={NUM_PATIENTS} windows={NUM_WINDOWS}")
    print(f"label_seed={LABEL_SEED}")
    print(f"trace_seed={TRACE_SEED}")
    print(f"class_vector={class_vector}")
    print(f"class_vector_sha256={class_sha}")
    print(f"washout_pos={sum(washout)} washout_neg={len(washout) - sum(washout)}")
    print(f"train_pos={sum(train)} train_neg={len(train) - sum(train)}")
    print(f"test_pos={sum(test)} test_neg={len(test) - sum(test)}")
    print(f"u_ahl_sha256={sha256_u(u_ahl)}")
    print(f"u_acid_sha256={sha256_u(u_acid)}")
    print(f"u_atta_sha256={sha256_u(u_atta)}")
    print(f"u_attb_sha256={sha256_u(u_attb)}")
    print(f"u_rep_sha256={sha256_u(u_rep)}")
    print(f"u_singlesite_sha256={sha256_u(u_singlesite)}")


if __name__ == "__main__":
    main()
