#!/usr/bin/env python3
"""Generate C1 independent NARMA-10 AHL sequences and freeze hashes.

Traj 00 is the existing Narma10b file and is never rewritten.
Traj 01-10 use predeclared random.Random seeds. Uniform[0,1] is forbidden.
Do not add a trajectory after seeing NRMSE.
"""

from __future__ import annotations

import hashlib
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
NARMA = HERE.parent / "BSimReservoirPlanNarma10b"
RESULTS = HERE / "results"

NUM_WINDOWS = 200
U_LOW = 0.0
U_HIGH = 0.5
TRAJ00_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
TRAJ00_SEED = 20260814
TRAJ00_PATH = NARMA / "input_ahl_narma200.txt"
TRAJ00_TARGET = NARMA / "narma10_target.csv"

# Predeclared. Do not extend after scores.
TRAJ_SEEDS = {
    "00": TRAJ00_SEED,
    "01": 2026081501,
    "02": 2026081502,
    "03": 2026081503,
    "04": 2026081504,
    "05": 2026081505,
    "06": 2026081506,
    "07": 2026081507,
    "08": 2026081508,
    "09": 2026081509,
    "10": 2026081510,
}
INTERACTION_TRAJS = ("01", "02")
PRIMARY_SEED = 111
INTERACTION_SEEDS = (222, 333)


def sha256_u(u) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def narma10(u):
    y = [0.0] * (len(u) + 1)
    for t, u_t in enumerate(u):
        acc = sum(y[t - i] if t - i >= 0 else 0.0 for i in range(10))
        u_lag = u[t - 9] if t >= 9 else 0.0
        y[t + 1] = 0.3 * y[t] + 0.05 * y[t] * acc + 1.5 * u_lag * u_t + 0.1
    return y


def load_sequence(path: Path):
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            values.append(float(line))
    return values


def make_u(seed: int):
    rng = random.Random(seed)
    u = [U_LOW + (U_HIGH - U_LOW) * rng.random() for _ in range(NUM_WINDOWS)]
    if any(value < U_LOW or value > U_HIGH for value in u):
        raise SystemExit("ABORT: u escaped [0, 0.5]; Uniform[0,1] is forbidden")
    return u


def ahl_path(traj_id: str) -> Path:
    if traj_id == "00":
        return TRAJ00_PATH
    return HERE / f"input_ahl_narma200_traj{traj_id}.txt"


def target_path(traj_id: str) -> Path:
    if traj_id == "00":
        return TRAJ00_TARGET
    return HERE / f"narma10_target_traj{traj_id}.csv"


def write_traj(traj_id: str, seed: int, u, y):
    ahl = ahl_path(traj_id)
    target = target_path(traj_id)
    ahl.write_text(
        f"# C1 independent NARMA-10 AHL. u ~ Uniform[0, 0.5], "
        f"random.Random({seed}), traj {traj_id}.\n"
        "# Do not regenerate after seeing NRMSE. Uniform[0,1] is forbidden.\n"
        + "\n".join(f"{value:.12f}" for value in u)
        + "\n",
        encoding="utf-8",
    )
    rows = ["window;u;y_next"]
    rows.extend(
        f"{index};{u[index]:.12f};{y[index + 1]:.12f}" for index in range(NUM_WINDOWS)
    )
    target.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_configs():
    template = HERE / "sim_config_c1.properties"
    template.write_text(
        """# C1 common production template. Frozen HybridDish / Narma10b claim dish.
# Do not retune K, n, tau, ALPHA_LUX, DELTA_LUX, AHL source, acid source,
# K_MAX, GROWTH_RATE, clamp, flow, or layout. FLOW_SPEED is compile-time 0.
# Traj 00 is scored in place from Narma10b. Change only input.ahl.file per run.
dt=0.05
grid.x=50
grid.y=25
grid.z=1
readout.grid.x=20
readout.grid.y=10
readout.grid.z=1
readout.countgrid.x=4
readout.countgrid.y=2
readout.countgrid.z=1
warmup.s=18000
warmup.ahl.input=0.5
window.duration.s=300
pulse.duration.s=75
sampling.duration.s=300
sampling.interval.s=20
num.windows=200
initial.pop=1800
carrying.capacity=2000
field.att.diff=100.0
field.rep.diff=100.0
field.ahl.diff=159.0
field.acid.diff=200.0
field.att.decay=0.0067
field.rep.decay=0.033
field.ahl.decay=0.0033
field.acid.decay=0.0067
field.ahl.source.rate=128000000
field.acid.source.rate=200000000000
field.acid.cell.production.rate=1000000
field.acid.kmax=0.002
luminescence.alpha=0.000666666666667
input.ahl.file=input_ahl_narma200_traj01.txt
input.acid.file=input_acid_held05_200.txt
headless=true
arm=driven
rng.seed=111
output.write.samples=true
output.write.voxels=true
output.run.label=c1_driven
output.stochastic.replicate=driven_traj01_seed111
output.dir=results/c1_driven_traj01_seed111
""",
        encoding="utf-8",
    )
    jobs = [(traj, PRIMARY_SEED) for traj in TRAJ_SEEDS if traj != "00"]
    for traj in INTERACTION_TRAJS:
        for seed in INTERACTION_SEEDS:
            jobs.append((traj, seed))
    for traj, seed in jobs:
        name = f"sim_config_c1_driven_traj{traj}_seed{seed}.properties"
        (HERE / name).write_text(
            f"""config.include=sim_config_c1.properties
arm=driven
rng.seed={seed}
input.ahl.file=input_ahl_narma200_traj{traj}.txt
output.run.label=c1_driven
output.stochastic.replicate=driven_traj{traj}_seed{seed}
output.dir=results/c1_driven_traj{traj}_seed{seed}
""",
            encoding="utf-8",
        )
    return jobs


def write_u_freeze(records):
    RESULTS.mkdir(parents=True, exist_ok=True)
    lines = [
        "# C1 u freeze",
        "",
        "Frozen **before** the first new Java job. Do not add traj 11.",
        "Do not drop a sequence after seeing NRMSE. Traj 00 was not rewritten.",
        "",
        "Payload SHA-256 is the comma-joined 12-decimal helper from",
        "`check_narma10b.py`. File SHA-256 is the on-disk bytes.",
        "",
        "Interaction subset (frozen here, not after scores): traj **01** and **02**.",
        "Those two new sequences, and only those two, also get bacterial seeds",
        "222 and 333. Traj 00 already has 111/222/333 in Narma10b.",
        "",
        "| TrajID | RNG seed | Provenance | payload SHA-256 | file SHA-256 |",
        "|---|---|---|---|---|",
    ]
    for row in records:
        lines.append(
            f"| {row['traj']} | {row['seed']} | {row['provenance']} | "
            f"`{row['payload_sha']}` | `{row['file_sha']}` |"
        )
    lines.extend(
        [
            "",
            f"Collision check: **{records[0]['collision']}**",
            "",
            "Traj 00 path (scored in place, not copied into C1 as a new run):",
            f"`{TRAJ00_PATH.as_posix()}`",
            "",
            "Package source hashes are appended after Java, configs, and the",
            "checker are written, still before the first new Java job.",
            "",
        ]
    )
    (RESULTS / "U_FREEZE.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if not TRAJ00_PATH.exists():
        raise SystemExit(f"ABORT: missing traj 00 file {TRAJ00_PATH}")
    u00 = load_sequence(TRAJ00_PATH)
    if len(u00) != NUM_WINDOWS:
        raise SystemExit("ABORT: traj 00 does not have 200 values")
    digest00 = sha256_u(u00)
    if digest00 != TRAJ00_SHA:
        raise SystemExit(f"ABORT: traj 00 payload hash {digest00} != {TRAJ00_SHA}")
    recomputed00 = make_u(TRAJ00_SEED)
    if sha256_u(recomputed00) != TRAJ00_SHA:
        raise SystemExit("ABORT: random.Random(20260814) does not reproduce traj 00")

    records = [
        {
            "traj": "00",
            "seed": TRAJ00_SEED,
            "provenance": "Existing Narma10b. Not rewritten.",
            "payload_sha": digest00,
            "file_sha": sha256_file(TRAJ00_PATH),
            "collision": "PASS",
        }
    ]
    seen = {digest00: "00"}
    for traj_id, seed in TRAJ_SEEDS.items():
        if traj_id == "00":
            continue
        u = make_u(seed)
        digest = sha256_u(u)
        if digest in seen:
            raise SystemExit(
                f"ABORT: traj {traj_id} payload hash collides with traj {seen[digest]}"
            )
        if digest == TRAJ00_SHA:
            raise SystemExit(f"ABORT: traj {traj_id} collides with traj 00")
        seen[digest] = traj_id
        y = narma10(u)
        write_traj(traj_id, seed, u, y)
        records.append(
            {
                "traj": traj_id,
                "seed": seed,
                "provenance": "new",
                "payload_sha": digest,
                "file_sha": sha256_file(ahl_path(traj_id)),
                "collision": "PASS",
            }
        )
        print(f"traj={traj_id} seed={seed} u_sha256={digest}")

    write_configs()
    write_u_freeze(records)
    print(f"traj=00 seed={TRAJ00_SEED} u_sha256={digest00} (not rewritten)")
    print("collision_check=PASS")
    print(f"wrote {RESULTS / 'U_FREEZE.md'}")
    print("configs=14 driven jobs + template")


if __name__ == "__main__":
    main()
