#!/usr/bin/env python3
"""Validate balance and temporal correlations of Stage 3 input permutations."""

import argparse
import json
import math
import random
import statistics
from collections import Counter
from pathlib import Path


def correlation(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    dx, dy = [value - mx for value in x], [value - my for value in y]
    denominator = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denominator if denominator else float("nan")


def read_values(path):
    return [
        float(line)
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*")
    parser.add_argument("--generate", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260813)
    args = parser.parse_args()
    if args.generate:
        rng = random.Random(args.seed)
        generated = []
        base = [level for level in (0.0, 0.25, 0.5, 0.75, 1.0) for _ in range(8)]
        while len(generated) < args.generate:
            rng.shuffle(base)
            index_r = correlation(list(range(40)), base)
            lag_r = [correlation(base[:-lag], base[lag:]) for lag in (1, 2, 3)]
            if abs(index_r) < 0.05 and max(map(abs, lag_r)) < 0.2:
                generated.append({
                    "values": base.copy(),
                    "r_window_index": index_r,
                    "autocorrelation": lag_r,
                })
        print(json.dumps(generated, indent=2))
        return
    result = {}
    for path in args.files:
        values = read_values(path)
        result[path] = {
            "count": len(values),
            "level_counts": dict(sorted(Counter(values).items())),
            "r_window_index": correlation(list(range(len(values))), values),
            "autocorrelation": {
                str(lag): correlation(values[:-lag], values[lag:])
                for lag in (1, 2, 3)
            },
        }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
