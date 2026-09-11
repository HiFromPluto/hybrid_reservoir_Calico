#!/usr/bin/env python3
"""Generate the frozen Waveform2b distinct-orbit temporal-order sequences.

Track E1 follow-on. Waveform2 remains TASK_VOID (ORDER_A ≡ ORDER_C orbit).
Classes are ORDER_UNI, ORDER_DOWN, ORDER_UP — not sine/square/triangle.
Orbit intersections must be empty before any SHA-256. Class shuffle uses
random.Random(20260818). Phase assignment uses a documented subseed.
"""

from __future__ import annotations

import hashlib
import math
import random
from collections import Counter
from pathlib import Path

import numpy as np

CLASS_SEED = 20260818
NUM_BLOCKS = 50
WINDOWS_PER_BLOCK = 8
NUM_WINDOWS = NUM_BLOCKS * WINDOWS_PER_BLOCK
N_CLASSES = 3
N_PHASES = 8
CLASS_NAMES = ("ORDER_UNI", "ORDER_DOWN", "ORDER_UP")

WASHOUT_BLOCKS = 8
TRAIN_BLOCKS = 26
TEST_BLOCKS = 16
INNER_VAL_BLOCKS = 6
WASHOUT_COUNTS = (3, 3, 2)
TRAIN_COUNTS = (9, 9, 8)
TEST_COUNTS = (6, 5, 5)

RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6)
START_U_MAX_ABS_FROM_CHANCE = 0.15
START_U_LEAK = 0.70
MAX_SUBSEED = 10000

S2 = 0.25 * math.sqrt(2.0) / 2.0
LO = round(0.25 - S2, 12)
MID = round(0.25, 12)
HI = round(0.25 + S2, 12)
MINV = round(0.0, 12)
MAXV = round(0.5, 12)

# Frozen 12-decimal alphabet.
assert LO == 0.073223304703
assert MID == 0.250000000000
assert HI == 0.426776695297
assert MINV == 0.000000000000
assert MAXV == 0.500000000000

SHARED_MULTISET = (MINV, LO, LO, MID, MID, HI, HI, MAXV)

TEMPLATES = (
    (MID, HI, MAXV, HI, MID, LO, MINV, LO),  # ORDER_UNI
    (MAXV, HI, HI, MID, MID, LO, LO, MINV),  # ORDER_DOWN
    (MINV, LO, LO, MID, MID, HI, HI, MAXV),  # ORDER_UP
)

ALPHABET_NAMES = {
    MINV: "MINV",
    LO: "LO",
    MID: "MID",
    HI: "HI",
    MAXV: "MAXV",
}

HERE = Path(__file__).resolve().parent


def sha256_u(values):
    payload = ",".join(f"{value:.12f}" for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def sha256_ascii_ints(values):
    payload = ",".join(str(int(v)) for v in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def clip_u(value):
    return min(MAXV, max(MINV, value))


def template_u(cls, phase, w):
    return clip_u(TEMPLATES[cls][(w + phase) % WINDOWS_PER_BLOCK])


def start_u(cls, phase):
    return template_u(cls, phase, 0)


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


def assert_distinct_orbits(templates):
    """Mandatory pre-hash gate. Pairwise circular-shift intersections empty."""
    orbits, matrix = orbit_intersection_matrix(templates)
    print("ORBIT_INTERSECTION_MATRIX")
    header = "         " + " ".join(f"{name:>12}" for name in CLASS_NAMES)
    print(header)
    for i, name in enumerate(CLASS_NAMES):
        cells = " ".join(f"{matrix[i][j]:12d}" for j in range(N_CLASSES))
        print(f"{name:>8} {cells}")
    uni_rev = tuple(reversed(templates[0]))
    uni_rev_in_uni = uni_rev in orbits[0]
    print(f"ORDER_UNI_REVERSE_IN_UNI_ORBIT={uni_rev_in_uni}")
    if uni_rev == templates[1] or uni_rev == templates[2]:
        raise RuntimeError(
            "ORBIT_GATE_FAIL: reverse of ORDER_UNI is used as another class"
        )
    failed = False
    for i in range(N_CLASSES):
        for j in range(i + 1, N_CLASSES):
            shared = orbits[i] & orbits[j]
            if shared:
                failed = True
                print(f"ORBIT_GATE_FAIL {CLASS_NAMES[i]} ∩ {CLASS_NAMES[j]} nonempty")
            if templates[j] in orbits[i] or templates[i] in orbits[j]:
                failed = True
                print(
                    f"ORBIT_GATE_FAIL {CLASS_NAMES[j]} is a circular shift of "
                    f"{CLASS_NAMES[i]}"
                )
    if failed:
        raise RuntimeError(
            "ORBIT_GATE_FAIL: pairwise cyclic orbits are not disjoint. "
            "STOP before hashes. Do not hand-permute after ridge."
        )
    print("ORBIT_GATE_PASS")
    return orbits, matrix


def forced_then_shuffled(rng, counts):
    labels = (
        [0] * counts[0]
        + [1] * counts[1]
        + [2] * counts[2]
    )
    rng.shuffle(labels)
    return labels


def balanced_phases(n, rng):
    """Spread n assignments across 8 phases as evenly as integers allow."""
    order = list(range(N_PHASES))
    rng.shuffle(order)
    q, r = divmod(n, N_PHASES)
    phases = []
    for i, phase in enumerate(order):
        phases.extend([phase] * (q + (1 if i < r else 0)))
    rng.shuffle(phases)
    if len(phases) != n:
        raise RuntimeError(f"phase assignment length {len(phases)} != {n}")
    return phases


def assign_phases(block_y, subseed):
    rng = random.Random(CLASS_SEED + 1 + subseed)
    phases = [0] * NUM_BLOCKS
    splits = (
        (0, WASHOUT_BLOCKS),
        (WASHOUT_BLOCKS, WASHOUT_BLOCKS + TRAIN_BLOCKS),
        (WASHOUT_BLOCKS + TRAIN_BLOCKS, NUM_BLOCKS),
    )
    for start, end in splits:
        for cls in range(N_CLASSES):
            idxs = [i for i in range(start, end) if block_y[i] == cls]
            assigned = balanced_phases(len(idxs), rng)
            for i, phase in zip(idxs, assigned):
                phases[i] = phase
    return phases


def phase_ok(block_y, phases, start, end, require_multi_phase):
    slice_y = block_y[start:end]
    slice_p = phases[start:end]
    for cls in range(N_CLASSES):
        used = {p for y, p in zip(slice_y, slice_p) if y == cls}
        if not used:
            continue
        if require_multi_phase and len(used) < 2:
            return False
    return True


def start_histogram(block_y, phases, start, end):
    table = {cls: Counter() for cls in range(N_CLASSES)}
    for y, p in zip(block_y[start:end], phases[start:end]):
        table[y][round(start_u(y, p), 12)] += 1
    return table


def start_values_overlap(block_y, phases, start, end):
    """Every class that appears must share at least one start value with another class."""
    hist = start_histogram(block_y, phases, start, end)
    supports = {cls: set(counts) for cls, counts in hist.items() if sum(counts.values())}
    if len(supports) < 2:
        return True
    for cls, support in supports.items():
        others = set()
        for other, other_support in supports.items():
            if other != cls:
                others |= other_support
        if support.isdisjoint(others):
            return False
    return True


def roc_auc(y_true, scores):
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    n_pos = int(np.sum(y_true == 1))
    n_neg = int(np.sum(y_true == 0))
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    sorted_scores = scores[order]
    ranks = np.empty(len(scores), dtype=float)
    i = 0
    n = len(scores)
    while i < n:
        j = i
        while j + 1 < n and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    pos_rank_sum = float(np.sum(ranks[y_true == 1]))
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def macro_ovr_auc(y_true, scores):
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    per_class = [roc_auc((y_true == k).astype(int), scores[:, k]) for k in range(N_CLASSES)]
    return float(np.mean(per_class)), per_class


def standardize(train, *others):
    mean = train.mean(axis=0)
    std = train.std(axis=0)
    std = np.where(std < 1e-12, 1.0, std)
    converted = [(matrix - mean) / std for matrix in (train,) + others]
    return converted, mean, std


def ridge_fit(X, y, lam):
    n_features = X.shape[1]
    Xb = np.column_stack([np.ones(len(X)), X])
    gram = Xb.T @ Xb
    gram[1:, 1:] = gram[1:, 1:] + lam * np.eye(n_features)
    target = Xb.T @ y
    try:
        return np.linalg.solve(gram, target)
    except np.linalg.LinAlgError:
        weights, *_ = np.linalg.lstsq(gram, target, rcond=None)
        return weights


def ridge_predict(X, weights):
    Xb = np.column_stack([np.ones(len(X)), X])
    return Xb @ weights


def ovr_scores(X, y, X_pred, lam):
    columns = []
    for k in range(N_CLASSES):
        yk = (np.asarray(y) == k).astype(float)
        wk = ridge_fit(X, yk, lam)
        columns.append(ridge_predict(X_pred, wk))
    return np.column_stack(columns)


def block_select_lambda(X, y):
    inner_train_end = TRAIN_BLOCKS - INNER_VAL_BLOCKS
    X_inner, y_inner = X[:inner_train_end], y[:inner_train_end]
    X_val, y_val = X[inner_train_end:TRAIN_BLOCKS], y[inner_train_end:TRAIN_BLOCKS]
    (X_inner_z, X_val_z), _, _ = standardize(X_inner, X_val)
    scored = []
    for lam in RIDGE_GRID:
        scores_val = ovr_scores(X_inner_z, y_inner, X_val_z, lam)
        macro, _ = macro_ovr_auc(y_val, scores_val)
        scored.append((-macro, -lam, lam, macro))
    scored.sort()
    return scored[0][2]


def block_test_macro_auc(X, y):
    lam = block_select_lambda(X, y)
    X_train, y_train = X[:TRAIN_BLOCKS], y[:TRAIN_BLOCKS]
    X_test, y_test = X[TRAIN_BLOCKS:], y[TRAIN_BLOCKS:]
    (X_train_z, X_test_z), _, _ = standardize(X_train, X_test)
    test_scores = ovr_scores(X_train_z, y_train, X_test_z, lam)
    macro, per_class = macro_ovr_auc(y_test, test_scores)
    return lam, macro, per_class


def block_features_from_schedule(block_y, phases):
    u_blocks = np.zeros((NUM_BLOCKS, WINDOWS_PER_BLOCK), dtype=float)
    y = np.array(block_y, dtype=int)
    for b in range(NUM_BLOCKS):
        for w in range(WINDOWS_PER_BLOCK):
            u_blocks[b, w] = template_u(block_y[b], phases[b], w)
    start = u_blocks[:, :1]
    return u_blocks, y, start


def start_u_auc(block_y, phases):
    _, y, start = block_features_from_schedule(block_y, phases)
    y_fit = y[WASHOUT_BLOCKS:]
    X_fit = start[WASHOUT_BLOCKS:]
    return block_test_macro_auc(X_fit, y_fit)


def constraints_ok(block_y, phases):
    if not phase_ok(block_y, phases, WASHOUT_BLOCKS, WASHOUT_BLOCKS + TRAIN_BLOCKS, True):
        return False, "train_single_phase"
    if not phase_ok(block_y, phases, WASHOUT_BLOCKS + TRAIN_BLOCKS, NUM_BLOCKS, True):
        return False, "test_single_phase"
    if not start_values_overlap(block_y, phases, WASHOUT_BLOCKS, WASHOUT_BLOCKS + TRAIN_BLOCKS):
        return False, "train_start_separated"
    if not start_values_overlap(block_y, phases, WASHOUT_BLOCKS + TRAIN_BLOCKS, NUM_BLOCKS):
        return False, "test_start_separated"
    _, macro, _ = start_u_auc(block_y, phases)
    if not math.isfinite(macro):
        return False, "start_u_auc_nan"
    if macro >= START_U_LEAK:
        return False, f"start_u_leak_{macro:.4f}"
    if abs(macro - 0.5) > START_U_MAX_ABS_FROM_CHANCE:
        return False, f"start_u_far_{macro:.4f}"
    return True, f"start_u_{macro:.4f}"


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


def format_hist(hist):
    keys = (MINV, LO, MID, HI, MAXV)
    lines = []
    for cls in range(N_CLASSES):
        parts = [f"{ALPHABET_NAMES[k]}={hist[cls].get(k, 0)}" for k in keys]
        lines.append(f"  {CLASS_NAMES[cls]}: " + " ".join(parts))
    return "\n".join(lines)


def main():
    directory = HERE
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "results").mkdir(parents=True, exist_ok=True)

    for cls, template in enumerate(TEMPLATES):
        if sorted_multiset(template) != sorted_multiset(SHARED_MULTISET):
            raise RuntimeError(f"{CLASS_NAMES[cls]} is not a permutation of the shared multiset")

    assert_distinct_orbits(TEMPLATES)

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

    chosen_subseed = None
    phases = None
    constraint_note = None
    for subseed in range(MAX_SUBSEED):
        candidate = assign_phases(block_y, subseed)
        ok, note = constraints_ok(block_y, candidate)
        if ok:
            chosen_subseed = subseed
            phases = candidate
            constraint_note = note
            break
        if subseed < 20 or subseed % 100 == 0:
            print(f"subseed={subseed} rejected ({note})")
    if phases is None:
        raise RuntimeError(
            f"no phase assignment satisfied constraints after {MAX_SUBSEED} subseeds"
        )

    u = []
    y_window = []
    phase_window = []
    for n in range(NUM_WINDOWS):
        block = n // WINDOWS_PER_BLOCK
        w = n % WINDOWS_PER_BLOCK
        cls = block_y[block]
        phase = phases[block]
        y_window.append(cls)
        phase_window.append(phase)
        u.append(round(template_u(cls, phase, w), 12))

    moment_rows = []
    for cls, template in enumerate(TEMPLATES):
        moment_rows.append(block_moments(template))
    for row in moment_rows[1:]:
        for key in ("mean", "variance", "pop_sd", "power", "min", "max", "range"):
            if abs(row[key] - moment_rows[0][key]) > 1e-12:
                raise RuntimeError(f"block moment {key} differs across classes")

    ahl_path = directory / "input_ahl_waveform2b_400.txt"
    ahl_path.write_text(
        "# Frozen Waveform2b AHL. Distinct-orbit equal-histogram templates.\n"
        "# Classes ORDER_UNI/ORDER_DOWN/ORDER_UP. Shared 8-value multiset. No extra noise.\n"
        f"# Class seed={CLASS_SEED}. Phase subseed={chosen_subseed}. u clipped to [0, 0.5].\n"
        "# Do not regenerate after seeing reservoir AUC.\n"
        + "\n".join(f"{value:.12f}" for value in u)
        + "\n",
        encoding="utf-8",
    )
    acid_path = directory / "input_acid_held05_400.txt"
    acid_path.write_text(
        "# Frozen Waveform2b acid hold. Every analysis window is 0.5. Off in warmup (Java).\n"
        + "\n".join("0.50" for _ in range(NUM_WINDOWS))
        + "\n",
        encoding="utf-8",
    )

    rows = ["n;block;class;phase;u"]
    rows.extend(
        f"{n};{n // WINDOWS_PER_BLOCK};{y_window[n]};{phase_window[n]};{u[n]:.12f}"
        for n in range(NUM_WINDOWS)
    )
    (directory / "waveform2b_labels.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    class_vector = ",".join(str(v) for v in block_y)
    phase_vector = ",".join(str(v) for v in phases)
    schedule = [
        "# Frozen Waveform2b (class, phase) schedule.",
        f"# class_seed={CLASS_SEED}",
        f"# phase_subseed={chosen_subseed}",
        "# class: 0=ORDER_UNI 1=ORDER_DOWN 2=ORDER_UP",
        f"# constraint={constraint_note}",
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
    print(f"phase_subseed={chosen_subseed}")
    print(f"constraint={constraint_note}")
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
    print("START_VALUE_HISTOGRAM")
    for name, start, end in (
        ("washout", 0, WASHOUT_BLOCKS),
        ("train", WASHOUT_BLOCKS, WASHOUT_BLOCKS + TRAIN_BLOCKS),
        ("test", WASHOUT_BLOCKS + TRAIN_BLOCKS, NUM_BLOCKS),
    ):
        print(name)
        print(format_hist(start_histogram(block_y, phases, start, end)))
    print("BLOCK_MOMENTS")
    print("class mean variance pop_sd power min max range n")
    for cls, row in enumerate(moment_rows):
        print(
            f"{CLASS_NAMES[cls]} {row['mean']:.12f} {row['variance']:.12f} "
            f"{row['pop_sd']:.12f} {row['power']:.12f} {row['min']:.12f} "
            f"{row['max']:.12f} {row['range']:.12f} {row['n']}"
        )
    print(f"u_min={min(u):.12f} u_max={max(u):.12f}")
    print("ALPHABET LO={:.12f} MID={:.12f} HI={:.12f}".format(LO, MID, HI))


if __name__ == "__main__":
    main()
