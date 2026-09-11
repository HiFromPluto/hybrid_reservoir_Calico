#!/usr/bin/env python3
"""Write frozen LANE_A_KR_GR u files. Not a rank. Not Java."""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "PROTOCOL_KR_GR.md"
PROTOCOL_JSON = HERE / "configs" / "kr_gr_protocol.json"
KR_FILE = HERE / "input_u_kr_iid840.txt"
GR_FILE = HERE / "input_u_gr_const840.txt"

T = 840
KR_SEED = 20260831
GR_CONST = 0.25
NARMA_SHA = "2fb692bb50bc8b2a9b8e28d2faf6ee63add7e6c980335c78f4f2c0ef3f07f06b"
NARMA10B_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"


def sha256_u(u) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def write_lines(path: Path, u) -> None:
    path.write_text("".join(f"{value:.12f}\n" for value in u), encoding="utf-8")


def main() -> None:
    proto = PROTOCOL.read_text(encoding="utf-8")
    js_text = PROTOCOL_JSON.read_text(encoding="utf-8")
    if "frozen_before_traces" not in proto or '"frozen_before_traces": true' not in js_text:
        raise SystemExit("PROTOCOL_KR_GR must exist and be frozen before u files")
    js = json.loads(js_text)
    if js.get("T_windows") != T or js.get("kr_u", {}).get("seed") != KR_SEED:
        raise SystemExit("freeze T/seed drifted before generation")
    if js.get("reuse_narma_u_as_kr_stream") is not False:
        raise SystemExit("must not reuse NARMA u")

    rng = random.Random(KR_SEED)
    kr = [rng.uniform(0.0, 0.5) for _ in range(T)]
    gr = [GR_CONST] * T
    kr_sha = sha256_u(kr)
    gr_sha = sha256_u(gr)
    if kr_sha in (NARMA_SHA, NARMA10B_SHA):
        raise SystemExit("KR stream collided with a forbidden NARMA hash")
    if any(v < 0.0 or v > 0.5 for v in kr):
        raise SystemExit("KR u escaped [0, 0.5]")

    write_lines(KR_FILE, kr)
    write_lines(GR_FILE, gr)

    js["kr_u"]["sha256"] = kr_sha
    js["gr_u"]["sha256"] = gr_sha
    PROTOCOL_JSON.write_text(json.dumps(js, indent=2) + "\n", encoding="utf-8")

    stamp = (
        f"- KR SHA-256 (stamped at generation, before Java): `{kr_sha}`\n"
        f"- GR SHA-256 (stamped at generation, before Java): `{gr_sha}`\n"
    )
    marker = "SHA hexes are stamped when the files are written, still **before**\nJava and **before** rank.\n"
    if marker not in proto:
        raise SystemExit("PROTOCOL missing SHA stamp marker")
    if "KR SHA-256 (stamped at generation" not in proto:
        proto = proto.replace(marker, marker + "\n" + stamp)
        PROTOCOL.write_text(proto, encoding="utf-8")

    print(f"KR {KR_FILE.name} sha256={kr_sha} n={T}")
    print(f"GR {GR_FILE.name} sha256={gr_sha} n={T}")
    print("u files written. No Java. No rank.")


if __name__ == "__main__":
    main()
