#!/usr/bin/env python3
"""Generate frozen Stage 7 four-source input sequences.

Independent Uniform draws, length 40. AHL in [0, 0.5]; attractant A/B and
acid in [0, 1]. Pairwise |r| must be < 0.3 before any BSim run. Do not
reuse the Stage 6 NARMA file. Do not regenerate after seeing SVD.
"""

from __future__ import annotations

import hashlib
import math
import random
import statistics
from pathlib import Path

NUM_WINDOWS = 40
AHL_HIGH = 0.5
MAX_ABS_R = 0.3
BASE_SEEDS = {
    "att_a": 20260870,
    "ahl": 20260871,
    "acid": 20260872,
    "att_b": 20260873,
}


def pearson(x, y):
    if len(x) != len(y) or len(x) < 2:
        return float("nan")
    mx, my = statistics.mean(x), statistics.mean(y)
    dx = [a - mx for a in x]
    dy = [b - my for b in y]
    denom = math.sqrt(sum(a * a for a in dx) * sum(b * b for b in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denom if denom else float("nan")


def draw(seed, high):
    rng = random.Random(seed)
    return [high * rng.random() for _ in range(NUM_WINDOWS)]


def pairwise(sequences):
    names = list(sequences)
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pairs.append((a, b, pearson(sequences[a], sequences[b])))
    return pairs


def sha256(values):
    return hashlib.sha256(
        (",".join(f"{value:.12f}" for value in values)).encode("ascii")
    ).hexdigest()


def write_sequence(path, header, values):
    path.write_text(
        header + "\n".join(f"{value:.12f}" for value in values) + "\n",
        encoding="utf-8",
    )


def main():
    directory = Path(__file__).resolve().parent
    offset = 0
    while True:
        seeds = {name: seed + offset for name, seed in BASE_SEEDS.items()}
        sequences = {
            "att_a": draw(seeds["att_a"], 1.0),
            "ahl": draw(seeds["ahl"], AHL_HIGH),
            "acid": draw(seeds["acid"], 1.0),
            "att_b": draw(seeds["att_b"], 1.0),
        }
        pairs = pairwise(sequences)
        if all(abs(r) < MAX_ABS_R for _, _, r in pairs):
            break
        offset += 1
        if offset > 10000:
            raise RuntimeError("could not find uncorrelated sequences")

    write_sequence(
        directory / "input_att_a_40.txt",
        "# Frozen Stage 7 Attractant A. Uniform[0,1]. Independent of AHL/acid/AttB.\n",
        sequences["att_a"],
    )
    write_sequence(
        directory / "input_ahl_40.txt",
        "# Frozen Stage 7 AHL. Uniform[0, 0.5]. Not the Stage 6 NARMA file.\n",
        sequences["ahl"],
    )
    write_sequence(
        directory / "input_acid_40.txt",
        "# Frozen Stage 7 acid. Uniform[0,1]. Free input, not held at 0.5.\n",
        sequences["acid"],
    )
    write_sequence(
        directory / "input_att_b_40.txt",
        "# Frozen Stage 7 Attractant B. Uniform[0,1]. Same attractant field, other corner.\n",
        sequences["att_b"],
    )
    report = ["channel;channel;r"]
    report.extend(f"{a};{b};{r:.12f}" for a, b, r in pairs)
    (directory / "input_pairwise_r.csv").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"num_windows={NUM_WINDOWS}")
    print(f"seed_offset={offset}")
    for name, seed in seeds.items():
        print(f"seed_{name}={seed} sha256={sha256(sequences[name])}")
    for a, b, r in pairs:
        print(f"r({a},{b})={r:.6f}")
    print("pairwise_max_abs_r=%.6f" % max(abs(r) for _, _, r in pairs))


if __name__ == "__main__":
    main()
