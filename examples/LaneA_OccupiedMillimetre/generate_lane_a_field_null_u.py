#!/usr/bin/env python3
"""Write frozen LANE_A_FIELD_NULL_TASK u and labels. Not AUC. Not Java."""

from __future__ import annotations

import csv
import hashlib
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "PROTOCOL_FIELD_NULL.md"
PROTOCOL_JSON = HERE / "configs" / "field_null_protocol.json"
U_FILE = HERE / "input_u_field_null_400.txt"
LABELS = HERE / "field_null_labels.csv"

T = 400
BIT_SEED = 2026083101
U_ONE = 0.5
U_ZERO = 0.0
FORBIDDEN = {
    "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b",
    "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e",
    "1134e3ac2031cf33f699341090b3d186f9a776a81cfe072b18e9d0c9553ba1b5",
    "0979090c051c513dce14632e83ae9bf025c4b666e497dc879ae1db4bba30d465",
}


def sha256_u(u) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def sha256_labels(y) -> str:
    payload = ",".join(str(int(v)) for v in y)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def main() -> None:
    proto = PROTOCOL.read_text(encoding="utf-8")
    js_text = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js_text:
        raise SystemExit("PROTOCOL_FIELD_NULL must exist and be frozen before u files")
    js = json.loads(js_text)
    if js.get("num_windows") != T or js.get("bit_seed") != BIT_SEED:
        raise SystemExit("freeze T/seed drifted before generation")
    if js.get("java_in_this_chat") is not False:
        raise SystemExit("this freeze forbids Java")
    if js.get("two_way") is not False or js.get("extra_ACs") is not False:
        raise SystemExit("field-null freeze drifted into two-way or extra ACs")
    if js.get("lane_a_charc_sweep_started") is not False:
        raise SystemExit("CHARC sweep must stay closed")
    if js.get("raise_J_max") is not False or js.get("second_Hill_K") is not False:
        raise SystemExit("must not raise J_max or add a second K")

    rng = random.Random(BIT_SEED)
    bits = [rng.randint(0, 1) for _ in range(T)]
    u = [U_ONE if b else U_ZERO for b in bits]
    prev = [0] + bits[:-1]
    y = [bits[n] ^ prev[n] for n in range(T)]
    y[0] = 0

    u_sha = sha256_u(u)
    y_sha = sha256_labels(y)
    if u_sha in FORBIDDEN:
        raise SystemExit("XOR stream collided with a forbidden u hash")
    if any(v < 0.0 or v > 0.5 for v in u):
        raise SystemExit("u escaped [0, 0.5]")
    if set(u) - {0.0, 0.5}:
        raise SystemExit("u is not a 0 / 0.5 bit command")

    U_FILE.write_text("".join(f"{value:.12f}\n" for value in u), encoding="utf-8")
    with LABELS.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerow(["n", "bit", "bit_prev", "xor", "u"])
        for n in range(T):
            writer.writerow([n, bits[n], prev[n], y[n], f"{u[n]:.12f}"])

    js["u_sha256"] = u_sha
    js["label_sha256"] = y_sha
    PROTOCOL_JSON.write_text(json.dumps(js, indent=2) + "\n", encoding="utf-8")

    stamp = (
        f"- u SHA-256 (stamped at generation, before Java): `{u_sha}`\n"
        f"- label SHA-256 (stamped at generation, before Java): `{y_sha}`\n"
    )
    marker = (
        "SHA hexes are stamped when the files are written, still **before**\n"
        "Java and **before** AUC.\n"
    )
    if marker not in proto:
        raise SystemExit("PROTOCOL missing SHA stamp marker")
    if "u SHA-256 (stamped at generation" not in proto:
        PROTOCOL.write_text(proto.replace(marker, marker + "\n" + stamp), encoding="utf-8")

    print(f"u {U_FILE.name} sha256={u_sha} n={T}")
    print(f"labels {LABELS.name} sha256={y_sha} n={T}")
    print("u files written. No Java. No AUC.")


if __name__ == "__main__":
    main()
