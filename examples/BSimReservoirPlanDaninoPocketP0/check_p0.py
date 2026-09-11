#!/usr/bin/env python3
"""Evaluate Gate P0 packed-device mechanics (chemistry OFF, NOT_FIG4B).

Thresholds are frozen in PROTOCOL.md before the production run existed.
Do not retune a threshold, k_cc, Lz, founder, or seed after seeing output.
If a gate fails, the failure is the finding.

Object A remains FAIL. Object B is OFF. This is not Fig. 4b.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
POCKET = ROOT / "examples" / "PocketDish"
CHASSIS = ROOT / "examples" / "BacteriumFromScratch"
PROD = HERE / "results" / "p0_seed101"
SMOKE_A = HERE / "results" / "p0_smoke_a"
SMOKE_B = HERE / "results" / "p0_smoke_b"
CSV_NAME = "colony_timeseries.csv"
CLAIM = POCKET / "CLAIM_FREEZE_DANINO_SI.md"
PROTOCOL = HERE / "PROTOCOL.md"

SEED = 101
T_END_S = 25950.0
T_SMOKE_S = 15570.0
LOG_DT_S = 10.0
N_TARGET = 256
SMOKE_N = 32
OVERLAP_CAP_UM = 0.25
MIN_CENTRE_UM = 0.75
OFFPLANE_MAX = 0.05
Z_EXC_MAX_UM = 0.40
R_RDISC_MIN = 0.7
R_RDISC_MAX = 1.3
D_CENTRES_MAX_UM = 1.2
LX, LY, LZ = 100.0, 100.0, 1.65
FOUNDER_Y = 80.0

NX, NY, NZ = 200, 90, 1
DX, DY, DZ = 2.0, 2.0, 1.65
MASK_X0, MASK_Y0, MASK_Z0 = -150.0, 0.0, 0.0

REQUIRED = [
    "t_s", "seed", "N", "R_um", "R_disc_um", "R_Rdisc", "offplane",
    "z_exc_max_um", "delta_cc_max_um", "d_centers_min_um", "spill_tick",
    "spill_cum", "N_ever", "door_contact", "y_max_um", "wall_leak",
    "n_drop_unexplained",
]

FROZEN_SHA256 = {
    "results/job2_seed101/job2_timeseries.csv":
        "e5ad8c7809d6b1712b47d4f231de59d71698aab888980e644c7604392df746ae",
    "results/job3_seed101/twobody_timeseries.csv":
        "b3a41af567834cc3076056027f9608f46bd541f8d191d88c8119c690b90eb5c9",
    "results/job3_seed101/cluster_timeseries.csv":
        "bceef8bc7a417a069552a82ebeb54efbd3f2ec8df66e9e67760292e512941774",
    "results/job3_seed101/isolated_packing_timeseries.csv":
        "cb21d5ecce63b59813f33ef54af8388ce70fdaf26608356afb37d72a06c6cf9a",
    "results/job3b_seed101/depletion_timeseries.csv":
        "351f2f2ffa5d8dc23b71dd0697eab7dc3056bdd3fcde017925537c6ff311d829",
    "results/job3b_seed101/packing_cluster_timeseries.csv":
        "e6fe563fa6c1b02db5c72c0a97d960c825e2a6658f96348050d48b2ece577a1c",
    "results/job3b_seed101/uniform_timeseries.csv":
        "ebf19ffa47e980f8b04b827ce9bcfd02bb190c500958c2ea717ec55369e74bbc",
    "results/job3c_seed101/colony_timeseries.csv":
        "f1d83f2e87ee36d3217568a78eb3ec9438c94886b712f0fa76e8c66640dc20cc",
}

JOB3B = (
    "results/job3b_seed101/depletion_timeseries.csv",
    "results/job3b_seed101/packing_cluster_timeseries.csv",
    "results/job3b_seed101/uniform_timeseries.csv",
)


def is_fluid(i: int, j: int, k: int) -> bool:
    cx = MASK_X0 + (i + 0.5) * DX
    cy = MASK_Y0 + (j + 0.5) * DY
    cz = MASK_Z0 + (k + 0.5) * DZ
    in_pocket = 0.0 <= cx <= LX and 0.0 <= cy <= LY and 0.0 <= cz <= LZ
    in_bus = MASK_X0 <= cx <= MASK_X0 + 400.0 and LY < cy <= LY + 80.0 and 0.0 <= cz <= LZ
    return in_pocket or in_bus


def mask_sha256() -> str:
    packed = bytearray(NX * NY * NZ)
    p = 0
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                packed[p] = 1 if is_fluid(i, j, k) else 0
                p += 1
    return hashlib.sha256(packed).hexdigest()


def load_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=";")
        if reader.fieldnames is None:
            raise ValueError(f"missing header: {path}")
        missing = [c for c in REQUIRED if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"missing columns: {missing}")
        return list(reader)


def gate1_growth(rows):
    errs = []
    n_ever_final = max(int(r["N_ever"]) for r in rows)
    if n_ever_final < N_TARGET:
        errs.append(f"N_ever never reached {N_TARGET} (max {n_ever_final})")
    prev = None
    for r in rows:
        n_ever = int(r["N_ever"])
        if prev is not None:
            if n_ever < prev:
                errs.append(f"N_ever decreased {prev}->{n_ever} at t={r['t_s']}")
                break
            if n_ever - prev > prev and prev > 0:
                errs.append(f"N_ever grew {prev}->{n_ever} at t={r['t_s']}: not binary")
                break
        prev = n_ever
    last_t = float(rows[-1]["t_s"])
    if abs(last_t - T_END_S) > LOG_DT_S:
        errs.append(f"last t={last_t}, expected {T_END_S}")
    return errs, n_ever_final


def gate2_packing(rows):
    errs = []
    for r in rows:
        if int(r["N"]) < 8:
            continue
        if float(r["delta_cc_max_um"]) > OVERLAP_CAP_UM:
            errs.append(
                f"overlap {r['delta_cc_max_um']} > {OVERLAP_CAP_UM} at t={r['t_s']}"
            )
            break
    for r in rows:
        if int(r["N"]) < 8:
            continue
        if float(r["d_centers_min_um"]) < MIN_CENTRE_UM:
            errs.append(
                f"centres {r['d_centers_min_um']} < {MIN_CENTRE_UM} at t={r['t_s']}"
            )
            break
    return errs


def first_at_target(rows, n=N_TARGET):
    for r in rows:
        if int(r["N"]) >= n:
            return r
    return None


def gate3_monolayer(target):
    if target is None:
        return ["never reached N=256 in pocket, gate not evaluable"]
    errs = []
    off = float(target["offplane"])
    if off > OFFPLANE_MAX:
        errs.append(f"offplane={off} > {OFFPLANE_MAX}")
    zexc = float(target["z_exc_max_um"])
    if zexc > Z_EXC_MAX_UM:
        errs.append(f"z_exc_max={zexc} > {Z_EXC_MAX_UM} um")
    rr = float(target["R_Rdisc"])
    if rr < R_RDISC_MIN or rr > R_RDISC_MAX:
        errs.append(f"R/R_disc={rr:.3f} not in [{R_RDISC_MIN}, {R_RDISC_MAX}]")
        if rr > 2.0:
            errs.append("filament growth — check DivisionMode.SYMMETRY_BROKEN")
    dmin = float(target["d_centers_min_um"])
    if dmin > D_CENTRES_MAX_UM:
        errs.append(
            f"d_centres_min={dmin:.3f} > {D_CENTRES_MAX_UM} um — "
            "lateral Hertzian never engaged"
        )
    return errs


def gate4_retention(rows):
    errs = []
    n_end = int(rows[-1]["N"])
    n_max = max(int(r["N"]) for r in rows)
    if n_end == 0:
        errs.append("N_end=0 EMPTY — T4/T5 class failure")
    if n_end < N_TARGET:
        errs.append(f"N_end={n_end} < {N_TARGET} (drained; N_max={n_max})")
    if any(int(r["wall_leak"]) != 0 for r in rows):
        errs.append("wall_leak > 0 — closed-face Hertzian leaked")
    if not any(int(r["door_contact"]) == 1 for r in rows):
        errs.append(
            "door_contact never fired: P0 did not ask the open-edge question "
            f"(founder_y={FOUNDER_Y})"
        )
    return errs


def gate5_extrusion(rows):
    errs = []
    if any(int(r["n_drop_unexplained"]) != 0 for r in rows):
        errs.append("n_drop_unexplained > 0 — death clamp or hidden deletion")
    first_door = None
    for r in rows:
        if int(r["door_contact"]) == 1:
            first_door = r
            break
    if first_door is None:
        errs.append("no door_contact; cannot test extrusion-only turnover")
        return errs
    t_door = float(first_door["t_s"])
    for r in rows:
        if float(r["t_s"]) < t_door and int(r["spill_cum"]) > 0:
            errs.append(
                f"spill_cum={r['spill_cum']} at t={r['t_s']} before door_contact "
                f"at t={t_door}"
            )
            break
    if int(rows[-1]["spill_cum"]) <= 0:
        errs.append("spill_cum=0 at t_end after door contact — no extrusion")
    return errs


def gate6_repeatability():
    errs = []
    a = SMOKE_A / CSV_NAME
    b = SMOKE_B / CSV_NAME
    if not a.is_file() or not b.is_file():
        errs.append("missing smoke CSVs (need results/p0_smoke_a and p0_smoke_b)")
        return errs
    if a.read_bytes() != b.read_bytes():
        errs.append("smoke CSVs are not byte-identical at seed 101")
        return errs
    try:
        rows = load_csv(a)
    except ValueError as exc:
        errs.append(str(exc))
        return errs
    target = first_at_target(rows, SMOKE_N)
    if target is None:
        errs.append(f"smoke never reached N={SMOKE_N}")
        return errs
    if float(target["offplane"]) > OFFPLANE_MAX:
        errs.append(f"smoke offplane={target['offplane']}")
    if float(target["R_Rdisc"]) >= 2.0:
        errs.append(f"smoke filament R/R_disc={target['R_Rdisc']}")
    if int(target["door_contact"]) != 0:
        errs.append("smoke door_contact fired (unexpected at N=32)")
    if int(target["spill_cum"]) != 0:
        errs.append("smoke spill_cum > 0")
    last_t = float(rows[-1]["t_s"])
    if abs(last_t - T_SMOKE_S) > LOG_DT_S:
        errs.append(f"smoke last t={last_t}, expected {T_SMOKE_S}")
    return errs


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gate_chassis():
    errs = []
    for rel, expected in FROZEN_SHA256.items():
        path = CHASSIS / rel
        if not path.is_file():
            errs.append(f"missing chassis CSV {rel}")
            continue
        got = sha256(path)
        if got != expected:
            kind = "Job 3b byte-identity" if rel in JOB3B else "chassis CSV"
            errs.append(
                f"{kind} moved: {rel}\n    expected {expected}\n    got      {got}"
            )
    return errs


def gate7_honesty(prod_dir: Path):
    errs = []
    if not PROTOCOL.is_file():
        errs.append("PROTOCOL.md missing")
    else:
        proto = PROTOCOL.read_text(encoding="utf-8")
        if "frozen_before_traces:** true" not in proto.replace(" ", ""):
            if "frozen_before_traces:** true" not in proto and "frozen_before_traces: true" not in proto:
                errs.append("PROTOCOL.md must freeze frozen_before_traces")
        if "1.65" not in proto:
            errs.append("PROTOCOL.md must declare Lz=1.65 µm")
        if "BSim.export()" in proto and "Not** `BSim.export()`" not in proto:
            pass
        if "BSimStepScheduler" not in proto:
            errs.append("PROTOCOL.md must name BSimStepScheduler")
    if not CLAIM.is_file():
        errs.append("CLAIM_FREEZE_DANINO_SI.md missing")
    else:
        text = CLAIM.read_text(encoding="utf-8")
        if "NOT_FIG4B" not in text or "DaninoSI_OccupiedDF" not in text:
            errs.append("claim freeze must name Object B NOT_FIG4B")
    freeze_path = prod_dir / "protocol_freeze.json"
    if not freeze_path.is_file():
        errs.append("protocol_freeze.json missing")
        return errs
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("chemistry") != "OFF":
        errs.append("chemistry is not OFF")
    if freeze.get("motility") != "OFF":
        errs.append("motility is not OFF")
    if freeze.get("object_b") != "OFF":
        errs.append("Object B is not OFF")
    if freeze.get("not_fig4b") is not True:
        errs.append("not_fig4b must be true")
    if abs(float(freeze.get("Lz_um", 0.0)) - LZ) > 1e-9:
        errs.append(f"Lz={freeze.get('Lz_um')} is not the SI trap {LZ}")
    if freeze.get("scheduler") != "BSimStepScheduler":
        errs.append("scheduler is not BSimStepScheduler")
    if freeze.get("death_clamp") != "OFF":
        errs.append("death_clamp is not OFF")
    mask_path = prod_dir / "geometry_mask.json"
    if not mask_path.is_file():
        errs.append("geometry_mask.json missing")
    else:
        mask = json.loads(mask_path.read_text(encoding="utf-8"))
        expected = mask_sha256()
        if mask.get("mask_sha256") != expected:
            errs.append(
                f"mask sha256 mismatch: java={mask.get('mask_sha256')} python={expected}"
            )
        if abs(float(mask.get("pocket_um", [0, 0, 0])[2]) - LZ) > 1e-9:
            errs.append("mask height is not 1.65 µm")
        if mask.get("chemistry") != "OFF":
            errs.append("mask chemistry is not OFF")
    return errs


def print_gate(label, errs, failures):
    status = "PASS" if not errs else "FAIL"
    print(f"Gate {label}: {status}")
    for e in errs:
        print(f"    {e}")
    failures.extend(errs)


def main() -> int:
    smoke_only = "--smoke" in sys.argv
    failures = []

    e6 = gate6_repeatability()
    print_gate("P0.6 repeatability", e6, failures)

    if smoke_only:
        print()
        print("  chemistry=OFF  motility=OFF  Object B=OFF  NOT_FIG4B")
        print()
        if failures:
            print("P0 SMOKE: FAIL")
            return 1
        print("P0 SMOKE: PASS")
        return 0

    csv_path = PROD / CSV_NAME
    if not csv_path.is_file():
        print("FAIL: missing", csv_path)
        return 1
    try:
        rows = load_csv(csv_path)
    except ValueError as exc:
        print("FAIL:", exc)
        return 1
    if not rows:
        print("FAIL: empty CSV")
        return 1
    if {int(r["seed"]) for r in rows} != {SEED}:
        print(f"FAIL: seed not uniformly {SEED}")
        return 1

    target = first_at_target(rows)
    e1, n_ever_final = gate1_growth(rows)
    e2 = gate2_packing(rows)
    e3 = gate3_monolayer(target)
    e4 = gate4_retention(rows)
    e5 = gate5_extrusion(rows)
    e7 = gate7_honesty(PROD)
    e_ch = gate_chassis()

    print_gate("P0.1 growth", e1, failures)
    print_gate("P0.2 packing", e2, failures)
    print_gate("P0.3 monolayer", e3, failures)
    print_gate("P0.4 retention", e4, failures)
    print_gate("P0.5 extrusion-only", e5, failures)
    print_gate("P0.7 honesty", e7, failures)
    print_gate("P0 chassis regression (Job 2/3/3b/3c bytes)", e_ch, failures)

    print()
    last = rows[-1]
    print(f"  N_end = {int(last['N'])}   N_ever = {n_ever_final}   "
          f"spill_cum = {int(last['spill_cum'])}   t_end = {last['t_s']} s")
    print(f"  door_contact_any = {max(int(r['door_contact']) for r in rows)}   "
          f"y_max_end = {float(last['y_max_um']):.2f} um   "
          f"wall_leak_max = {max(int(r['wall_leak']) for r in rows)}")
    print(f"  L = {LX:.0f}x{LY:.0f}x{LZ} um   founder_y = {FOUNDER_Y:.0f} um"
          f"   seed = {SEED}   open=+y")
    if target is not None:
        print(f"  at N>=256: t = {target['t_s']} s, R = {float(target['R_um']):.2f} um")
        print(f"    R/R_disc = {float(target['R_Rdisc']):.3f}"
              f"   offplane = {float(target['offplane']):.3f}"
              f"   z_exc_max = {float(target['z_exc_max_um']):.3f} um"
              f"   d_centres_min = {float(target['d_centers_min_um']):.3f} um"
              f"   delta_cc_max = {float(target['delta_cc_max_um']):.4f} um")
    print()
    print("  chemistry=OFF  motility=OFF  Object B=OFF  NOT_FIG4B")
    print("  AHL=OFF  Hill=OFF  NARMA=OFF  AC=OFF  LuxI=OFF  death=OFF")
    print()
    if failures:
        print("P0 OVERALL: FAIL")
        print("C0 may start: NO")
        return 1
    print("P0 OVERALL: PASS")
    print("C0 may start: YES (Object B well-mixed coupling only; NOT_FIG4B; not Fig. 4b)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
