#!/usr/bin/env python3
"""Generate Track E3 Lorenz L3 development clocks and hash them.

Historical anchors k=1 (BenchA) and k=50 (BenchA2) are verified, not
regenerated. Acid is copied from Stage 6 / BenchA held-0.5. Do not change
sigma, rho, beta, IC, or transient after seeing NRMSE. This script does
not rank baselines.
"""

from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path

import numpy as np

SIGMA = 10.0
RHO = 28.0
BETA = 8.0 / 3.0
DT_L = 0.02
IC = (1.0, 1.0, 1.0)
TRANSIENT_STEPS = 5000
NUM_RECORDED = 201
NUM_WINDOWS = 200
U_HIGH = 0.5

DEVELOPMENT_K = (5, 10, 20, 40)
ANCHOR_K = (1, 50)
ANCHOR_SHA = {
    1: "3bc69bbe32f6f2cf6a8638c4e2775133b4d27026a4223cbe783c057c19232f79",
    50: "35646e7b504940dbc859ec66f1f2c5b49017f25388c28f0d662ffc29a7f9a1e5",
}

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
BENCHA = EXAMPLES / "BSimReservoirPlanBenchA"
BENCHA2 = EXAMPLES / "BSimReservoirPlanBenchA2"
STAGE6 = EXAMPLES / "BSimReservoirPlanStage6"


def lorenz_deriv(state: np.ndarray) -> np.ndarray:
    x, y, z = state
    return np.array([
        SIGMA * (y - x),
        x * (RHO - z) - y,
        x * y - BETA * z,
    ], dtype=float)


def rk4_step(state: np.ndarray, dt: float) -> np.ndarray:
    k1 = lorenz_deriv(state)
    k2 = lorenz_deriv(state + 0.5 * dt * k1)
    k3 = lorenz_deriv(state + 0.5 * dt * k2)
    k4 = lorenz_deriv(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate(skip: int) -> np.ndarray:
    """Sample i is the state after 5000 + i*skip RK4 steps from the IC."""
    state = np.array(IC, dtype=float)
    for _ in range(TRANSIENT_STEPS):
        state = rk4_step(state, DT_L)
    samples = np.empty((NUM_RECORDED, 3), dtype=float)
    samples[0] = state
    for index in range(1, NUM_RECORDED):
        for _ in range(skip):
            state = rk4_step(state, DT_L)
        samples[index] = state
    return samples


def affine_u(x: np.ndarray) -> tuple[np.ndarray, float, float]:
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if xmax == xmin:
        raise SystemExit("generator is wrong: xmax == xmin")
    u = U_HIGH * (x - xmin) / (xmax - xmin)
    if np.any(u < 0.0) or np.any(u > U_HIGH + 1e-15):
        raise SystemExit("u escaped [0, 0.5]")
    u = np.clip(u, 0.0, U_HIGH)
    return u, xmin, xmax


def sha256_u(u: np.ndarray) -> str:
    payload = ",".join(f"{value:.12f}" for value in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def pack_from_samples(samples: np.ndarray) -> dict:
    x_drive = samples[:NUM_WINDOWS, 0]
    x_next = samples[1:NUM_WINDOWS + 1, 0]
    y_next = samples[1:NUM_WINDOWS + 1, 1]
    u, xmin, xmax = affine_u(x_drive)
    return {
        "u": u,
        "x": x_drive,
        "x_next": x_next,
        "y_next": y_next,
        "xmin": xmin,
        "xmax": xmax,
        "digest": sha256_u(u),
    }


def load_stored_u(path: Path) -> np.ndarray:
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        values.append(float(stripped))
    if len(values) != NUM_WINDOWS:
        raise SystemExit(f"{path} has {len(values)} values; expected {NUM_WINDOWS}")
    return np.array(values, dtype=float)


def write_target(path: Path, pack: dict) -> None:
    rows = ["n;u;x;y_next;x_next"]
    for index in range(NUM_WINDOWS):
        rows.append(
            f"{index};{pack['u'][index]:.12f};{pack['x'][index]:.12f};"
            f"{pack['y_next'][index]:.12f};{pack['x_next'][index]:.12f}"
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_ahl(path: Path, pack: dict, skip: int) -> None:
    dt_sample = skip * DT_L
    path.write_text(
        f"# Frozen Track E3 Lorenz L3 AHL input. Affine x[0..199] onto [0, 0.5].\n"
        f"# sigma=10, rho=28, beta=8/3, RK4 dt=0.02, IC=(1,1,1), discard 5000 steps.\n"
        f"# skip={skip} RK4 steps between samples (dt_sample={dt_sample:.2f}).\n"
        f"# Do not regenerate after seeing NRMSE. One AHL channel only.\n"
        + "\n".join(f"{value:.12f}" for value in pack["u"])
        + "\n",
        encoding="utf-8",
    )


def verify_anchor(skip: int, stored_ahl: Path, stored_target: Path, pack: dict) -> None:
    stored_u = load_stored_u(stored_ahl)
    digest = sha256_u(stored_u)
    expected = ANCHOR_SHA[skip]
    if digest != expected:
        raise SystemExit(f"anchor k={skip} stored u hash {digest} != {expected}")
    if pack["digest"] != expected:
        raise SystemExit(
            f"anchor k={skip} recomputed u hash {pack['digest']} != {expected}"
        )
    if not np.allclose(stored_u, pack["u"], rtol=0, atol=1e-9):
        raise SystemExit(f"anchor k={skip} stored u does not match RK4 affine map")
    with stored_target.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))
    if len(rows) != NUM_WINDOWS:
        raise SystemExit(f"anchor k={skip} target has {len(rows)} rows")
    stored_x = np.array([float(row["x"]) for row in rows], dtype=float)
    stored_y = np.array([float(row["y_next"]) for row in rows], dtype=float)
    stored_tu = np.array([float(row["u"]) for row in rows], dtype=float)
    if not np.allclose(stored_x, pack["x"], rtol=0, atol=1e-9):
        raise SystemExit(f"anchor k={skip} stored x does not match RK4")
    if not np.allclose(stored_y, pack["y_next"], rtol=0, atol=1e-9):
        raise SystemExit(f"anchor k={skip} stored y_next does not match RK4")
    if not np.allclose(stored_tu, pack["u"], rtol=0, atol=1e-9):
        raise SystemExit(f"anchor k={skip} stored target u does not match affine map")


def main() -> None:
    stage6_acid = STAGE6 / "input_acid_held05_200.txt"
    if not stage6_acid.exists():
        raise SystemExit(f"missing Stage 6 acid file: {stage6_acid}")
    acid_path = HERE / "input_acid_held05_200.txt"
    shutil.copyfile(stage6_acid, acid_path)

    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    hash_rows = []

    print("LORENZ_L3 generator: hashes before any baseline ranking")
    print(f"sigma={SIGMA} rho={RHO} beta={BETA} dt_L={DT_L} ic={IC}")
    print(f"transient_steps={TRANSIENT_STEPS} recorded={NUM_RECORDED}")

    for skip in ANCHOR_K + DEVELOPMENT_K:
        pack = pack_from_samples(integrate(skip))
        dt_sample = skip * DT_L
        role = "anchor" if skip in ANCHOR_K else "development"
        if skip == 1:
            ahl_path = BENCHA / "input_ahl_lorenz200.txt"
            target_path = BENCHA / "lorenz_target.csv"
            verify_anchor(skip, ahl_path, target_path, pack)
            print(
                f"k={skip} dt={dt_sample:.2f} ANCHOR verified "
                f"sha256={pack['digest']} xmin={pack['xmin']:.12f} xmax={pack['xmax']:.12f}"
            )
        elif skip == 50:
            ahl_path = BENCHA2 / "input_ahl_lorenz_skip50.txt"
            target_path = BENCHA2 / "lorenz_skip50_target.csv"
            verify_anchor(skip, ahl_path, target_path, pack)
            print(
                f"k={skip} dt={dt_sample:.2f} ANCHOR verified "
                f"sha256={pack['digest']} xmin={pack['xmin']:.12f} xmax={pack['xmax']:.12f}"
            )
        else:
            ahl_path = HERE / f"input_ahl_lorenz_k{skip}.txt"
            target_path = HERE / f"lorenz_target_k{skip}.csv"
            write_ahl(ahl_path, pack, skip)
            write_target(target_path, pack)
            print(
                f"k={skip} dt={dt_sample:.2f} WROTE {ahl_path.name} "
                f"sha256={pack['digest']} xmin={pack['xmin']:.12f} xmax={pack['xmax']:.12f}"
            )
        hash_rows.append({
            "k": skip,
            "dt_sample": f"{dt_sample:.2f}",
            "role": role,
            "u_sha256": pack["digest"],
            "xmin": f"{pack['xmin']:.12f}",
            "xmax": f"{pack['xmax']:.12f}",
            "u_min": f"{float(np.min(pack['u'])):.12f}",
            "u_max": f"{float(np.max(pack['u'])):.12f}",
            "ahl_path": str(ahl_path),
            "target_path": str(target_path),
            "regenerated": "FALSE" if skip in ANCHOR_K else "TRUE",
        })

    hash_csv = results / "lorenz_l3_input_hashes.csv"
    with hash_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(hash_rows[0].keys()))
        writer.writeheader()
        writer.writerows(hash_rows)
    print(f"hashes={hash_csv}")
    print(f"acid={acid_path.name} copied from Stage 6; not regenerated")
    print("HASHES_RECORDED_BEFORE_RANKING")


if __name__ == "__main__":
    main()
