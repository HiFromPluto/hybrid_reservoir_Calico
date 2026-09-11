#!/usr/bin/env python3
"""E4.2 seed-111 scout checker. Closed two-way vs A1/replay.

Does not edit GATE_EVIDENCE. Does not retune K, n, tau, clamp, mortality,
flow, layout, K_AC, Jmax, K_P, k_P, or alpha. Seeds 222/333 are not started.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NARMA = HERE.parent / "BSimReservoirPlanNarma10b"
E4AC = HERE.parent / "BSimReservoirPlanE4AC"
sys.path.insert(0, str(NARMA))
sys.path.insert(0, str(E4AC))
import check_e4ac as E  # noqa: E402
import check_narma10b as N  # noqa: E402

WASHOUT, TRAIN, TEST = 40, 110, 50
NUM_WINDOWS = 200
NARMA_SHA = "d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e"
LENTINI_SHA = "5c4beb444515ffc58ba01d333598ba5ddaa47e75447eaa0e0691ac2549411ab6"
G0, N_AC, K_AC = 0.0, 2.0, 0.25
JMAX_A1 = 74677509.75821304
N_P, K_P, ALPHA = 2.0, 1.6, 0.5
K_P_RATE = 46121.49695387294
WEAK_MARGIN = 0.03
OCC_DEAD_R = 0.05
OCC_DEAD_CORR = 0.5
WARMUP_S = 18_000.0
WINDOW_S = 300.0
PULSE_S = 75.0

A0_DRIVEN = NARMA / "results" / "narma10b_driven_seed111"
A0_BROWN = NARMA / "results" / "narma10b_brownian_seed111"
A0_SILENT = NARMA / "results" / "narma10b_silent_seed111"
A1_DRIVEN = E4AC / "results" / "e4ac_a1_driven_seed111"
A1_BROWN = E4AC / "results" / "e4ac_a1_brownian_seed111"
CLOSED = HERE / "results" / "e42_closed_seed111"
REPLAY = HERE / "results" / "e42_replay_seed111"
SHUFFLE = HERE / "results" / "e42_shuffle_seed111"

FORBIDDEN_JAVA = (
    "glucose", "Glucose", "GLUCOSE",
    "Danino", "danino",
    "CS_in",
    "AiiA",
    "LuxI",
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ac_gate(u):
    u = np.asarray(u, dtype=float)
    un = np.power(np.maximum(u, 0.0), N_AC)
    kn = K_AC ** N_AC
    return G0 + (1.0 - G0) * un / (kn + un)


def ac_h(p, alpha=ALPHA):
    p = np.asarray(p, dtype=float)
    pn = np.power(np.maximum(p, 0.0), N_P)
    kn = K_P ** N_P
    return 1.0 + alpha * pn / (kn + pn)


def load_frozen():
    frozen = json.loads((HERE / "results" / "tw_frozen.json").read_text(encoding="utf-8"))
    if frozen["MODEL_STATUS"] != "HYPOTHETICAL_DESIGN_ENVELOPE":
        raise SystemExit("STOP MODEL_STATUS must remain HYPOTHETICAL_DESIGN_ENVELOPE")
    if frozen.get("digital_twin"):
        raise SystemExit("STOP digital_twin must be false")
    if abs(frozen["Jmax_A1"] - JMAX_A1) > 1e-6:
        raise SystemExit("STOP Jmax_A1 was retuned")
    if abs(frozen["K_AC"] - K_AC) > 1e-15 or abs(frozen["K_P"] - K_P) > 1e-15:
        raise SystemExit("STOP K_AC or K_P was retuned")
    if abs(frozen["k_P"] - K_P_RATE) / K_P_RATE > 1e-9:
        raise SystemExit("STOP k_P was retuned")
    return frozen


def java_ok():
    java = (HERE / "BSimReservoirPlanE42TW.java").read_text(encoding="utf-8")
    if any(token in java for token in FORBIDDEN_JAVA):
        return False
    if "class ArtificialCell" in java:
        return False
    if "FLOW_SPEED = 0.0" not in java:
        return False
    if "RECEIVER_K_UM = 1.6" not in java:
        return False
    if "acFeedback" not in java or "K_P_RATE" not in java:
        return False
    if "package BSimReservoirPlanE42TW;" not in java:
        return False
    if "JMAX_A1_FROZEN" not in java:
        return False
    return True


def load_p_ac(path: Path):
    times, p_ac, h, p_live = [], [], [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        times.append(float(parts[0]))
        p_ac.append(float(parts[1]))
        h.append(float(parts[2]) if len(parts) > 2 else float("nan"))
        p_live.append(float(parts[3]) if len(parts) > 3 else float("nan"))
    return {
        "t": np.array(times),
        "p_ac": np.array(p_ac),
        "h": np.array(h),
        "p_live": np.array(p_live),
    }


def window_series(pfile, u):
    t = pfile["t"]
    p = pfile["p_ac"]
    h = pfile["h"]
    mean_p = np.zeros(NUM_WINDOWS)
    mean_h = np.zeros(NUM_WINDOWS)
    pulse_h = np.zeros(NUM_WINDOWS)
    for window in range(NUM_WINDOWS):
        t0 = WARMUP_S + window * WINDOW_S
        t1 = t0 + WINDOW_S
        tp = t0 + PULSE_S
        mask = (t >= t0) & (t < t1)
        pulse = (t >= t0) & (t < tp)
        if not np.any(mask):
            raise SystemExit(f"STOP P_ac missing window {window}")
        mean_p[window] = float(np.mean(p[mask]))
        mean_h[window] = float(np.mean(h[mask]))
        pulse_h[window] = float(np.mean(h[pulse])) if np.any(pulse) else mean_h[window]
    g = ac_gate(u)
    gh = g * pulse_h
    j_mean = JMAX_A1 * gh * (PULSE_S / WINDOW_S)
    return {
        "mean_P": mean_p,
        "mean_h": mean_h,
        "pulse_h": pulse_h,
        "g": g,
        "gh": gh,
        "j_mean": j_mean,
    }


def write_shuffle(src: Path, dest: Path, seed: int = 111) -> None:
    pfile = load_p_ac(src)
    rng = np.random.default_rng(seed)
    t = pfile["t"]
    warmup = t < WARMUP_S
    blocks = []
    for window in range(NUM_WINDOWS):
        t0 = WARMUP_S + window * WINDOW_S
        t1 = t0 + WINDOW_S
        blocks.append(np.where((t >= t0) & (t < t1))[0])
    order = rng.permutation(NUM_WINDOWS)
    out_t = []
    out_p = []
    out_h = []
    out_live = []
    for idx in np.where(warmup)[0]:
        out_t.append(pfile["t"][idx])
        out_p.append(pfile["p_ac"][idx])
        out_h.append(pfile["h"][idx])
        out_live.append(pfile["p_live"][idx])
    for new_window, old_window in enumerate(order):
        t0 = WARMUP_S + new_window * WINDOW_S
        old = blocks[old_window]
        old_t0 = WARMUP_S + old_window * WINDOW_S
        for idx in old:
            out_t.append(t0 + (pfile["t"][idx] - old_t0))
            out_p.append(pfile["p_ac"][idx])
            out_h.append(pfile["h"][idx])
            out_live.append(pfile["p_live"][idx])
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as handle:
        handle.write("# t_s P_ac_uM h P_live_uM shuffled_from_closed seed=111\n")
        for row in zip(out_t, out_p, out_h, out_live):
            handle.write(f"{row[0]:.2f} {row[1]:.12e} {row[2]:.12e} {row[3]:.12e}\n")
    print(f"wrote shuffle {dest} order={order.tolist()[:8]}...")


def occupancy_from_run(run_dir: Path, u):
    voxels = E.load_voxel_arrays(run_dir / "voxels.csv")
    return E.occupancy(voxels, u), voxels


def gate_evidence_untouched():
    return E.gate_evidence_untouched()


def claim_label(pass_flag: bool) -> str:
    return "PASS" if pass_flag else "FAIL"


def analyze(u, y, frozen):
    alpha = float(frozen["alpha"])
    a0_d, a0_b, a0_s = E.claim_sanity(u, y)
    if not E.csv_ok(A1_DRIVEN) or not E.csv_ok(A1_BROWN):
        raise SystemExit("STOP A1 CSV completeness failed; do not overwrite E4AC")
    if not E.csv_ok(CLOSED):
        raise SystemExit("STOP CLOSED CSV completeness failed")
    a1_driven = N.read_run(A1_DRIVEN, "driven", NUM_WINDOWS, E.EXPECTED_AUX)
    a1_brown = N.read_run(A1_BROWN, "brownian", NUM_WINDOWS, E.EXPECTED_AUX)
    a1_d = N.evaluate_arm_run(a1_driven, u, y)
    a1_b = N.evaluate_arm_run(a1_brown, u, y)
    closed_run = N.read_run(CLOSED, "driven", NUM_WINDOWS, E.EXPECTED_AUX)
    closed_d = N.evaluate_arm_run(closed_run, u, y)
    occ, voxels = occupancy_from_run(CLOSED, u)
    pfile = load_p_ac(CLOSED / "p_ac_timeseries.txt")
    series = window_series(pfile, u)
    sr, sl = E.integrate_surrogate(voxels["absolute_time"], voxels["ahl"])
    mask = (voxels["den"] != 0).astype(float)
    masked = np.column_stack([
        E.window_mean_matrix(sr * mask, voxels["window"]),
        E.window_mean_matrix(sl * mask, voxels["window"]),
    ])
    closed_masked = E.ridge_score(masked, y)
    gh_scalar = E.ridge_score(series["gh"].reshape(-1, 1), y)
    j_mean = E.ridge_score(series["j_mean"].reshape(-1, 1), y)
    gh_taps = E.ridge_score(E.delay_matrix(series["gh"], 10), y)
    g_taps = E.ridge_score(E.delay_matrix(series["g"], 10), y)
    taps_u = E.delay_matrix(u, 10)
    product = np.array([u[n] * u[n - 9] if n >= 9 else 0.0 for n in range(len(u))])
    informed = np.column_stack([taps_u, product])

    replay_exists = (REPLAY / "voxels.csv").exists()
    replay_d = None
    replay_occ = None
    if replay_exists:
        if not E.csv_ok(REPLAY):
            raise SystemExit("STOP REPLAY CSV completeness failed")
        replay_run = N.read_run(REPLAY, "driven", NUM_WINDOWS, E.EXPECTED_AUX)
        replay_d = N.evaluate_arm_run(replay_run, u, y)
        replay_occ, _ = occupancy_from_run(REPLAY, u)

    closed_f408 = closed_d["task"]["test_nrmse"]
    closed_field = closed_d["field_only"]["test_nrmse"]
    a1_f408 = a1_d["task"]["test_nrmse"]
    brown = a1_b["task"]["test_nrmse"]
    silent = a0_s["task"]["test_nrmse"]
    replay_f408 = replay_d["task"]["test_nrmse"] if replay_d else float("nan")
    replay_field = replay_d["field_only"]["test_nrmse"] if replay_d else float("nan")

    system_pass = closed_f408 < brown and closed_f408 < silent
    living_pass = closed_f408 < closed_field
    d_a1 = abs(closed_f408 - a1_f408)
    d_replay = abs(closed_f408 - replay_f408) if math.isfinite(replay_f408) else float("nan")
    ties_replay = math.isfinite(d_replay) and d_replay < WEAK_MARGIN
    ties_a1 = d_a1 < WEAK_MARGIN
    beats_a1 = (not ties_a1) and closed_f408 < a1_f408
    beats_replay = (not ties_replay) and math.isfinite(replay_f408) and closed_f408 < replay_f408
    if ties_a1:
        interaction = "NO_STORY_MOVE"
    elif beats_a1 and ties_replay:
        interaction = "FEEDFORWARD_SCHEDULE"
    elif beats_a1 and beats_replay:
        interaction = "PASS"
    else:
        interaction = "FAIL"
    transducer = "NO_STORY_MOVE" if d_a1 < WEAK_MARGIN else "MOVED"
    shuffle_auth = (not ties_a1) or (math.isfinite(d_replay) and d_replay >= WEAK_MARGIN)
    shuffle_skip = not shuffle_auth
    mean_p = float(np.mean(series["mean_P"]))
    mean_h = float(np.mean(series["mean_h"]))
    r_pu = E.pearson(series["mean_P"], u)
    row = {
        "Seed": 111,
        "MODEL_STATUS": "HYPOTHETICAL_DESIGN_ENVELOPE",
        "digital_twin": False,
        "alpha": alpha,
        "K_P": K_P,
        "k_P": K_P_RATE,
        "occupancy": occ["occupancy"],
        "mean_R": occ["mean_R"],
        "r_meanR_u": occ["r_meanR_u"],
        "mean_P_ac": mean_p,
        "mean_h": mean_h,
        "r_Pac_u": r_pu,
        "mean_AHL": occ["mean_AHL"],
        "A0_F408": a0_d["task"]["test_nrmse"],
        "A0_field": a0_d["field_only"]["test_nrmse"],
        "A0_brownian": a0_b["task"]["test_nrmse"],
        "A0_silent": a0_s["task"]["test_nrmse"],
        "A1_F408": a1_f408,
        "A1_field": a1_d["field_only"]["test_nrmse"],
        "A1_brownian": brown,
        "A1_silent_reused": silent,
        "A1_G_DELAY_10": g_taps["test_nrmse"],
        "CLOSED_F408": closed_f408,
        "CLOSED_field": closed_field,
        "CLOSED_masked": closed_masked["test_nrmse"],
        "REPLAY_F408": replay_f408,
        "REPLAY_field": replay_field,
        "AC_GH_SCALAR": gh_scalar["test_nrmse"],
        "AC_J_WINDOW_MEAN": j_mean["test_nrmse"],
        "AC_GH_DELAY_10": gh_taps["test_nrmse"],
        "LINEAR_U_DELAY_10": E.ridge_score(taps_u, y)["test_nrmse"],
        "NARMA_INFORMED_INPUT": E.ridge_score(informed, y)["test_nrmse"],
        "abs_CLOSED_minus_A1_F408": d_a1,
        "abs_CLOSED_minus_REPLAY_F408": d_replay,
        "system_pass": system_pass,
        "living_layer_pass": living_pass,
        "transducer": transducer,
        "interaction": interaction,
        "replay_tied_closed": ties_replay,
        "shuffle_authorized": shuffle_auth,
        "SKIP_SHUFFLE": shuffle_skip,
        "GATE_EVIDENCE_untouched": gate_evidence_untouched(),
        "java_ok": java_ok(),
        "ahl_finite": occ["ahl_finite"],
        "ahl_nonneg": occ["ahl_nonneg"],
        "p_finite": bool(np.all(np.isfinite(pfile["p_ac"]))),
        "p_nonneg": bool(np.min(pfile["p_ac"]) >= -1e-12),
        "lentini_pdf_sha256": LENTINI_SHA,
        "u_sha256": NARMA_SHA,
        "replay_occupancy": None if replay_occ is None else replay_occ["occupancy"],
    }
    return row, occ, series


def write_report(row):
    system = claim_label(row["system_pass"])
    living = claim_label(row["living_layer_pass"])
    replay = "n/a" if not math.isfinite(row["REPLAY_F408"]) else f"{row['REPLAY_F408']:.4f}"
    lines = [
        "# E4.2 two-way scout (seed 111)",
        "",
        "**MODEL_STATUS = `HYPOTHETICAL_DESIGN_ENVELOPE`.** Not a digital twin.",
        "Module 1 PASSed and Module 1.6 was STABLE_RESPONSIVE before this living",
        "run. Narma10b and E4AC voxels were reused, not overwritten.",
        "`GATE_EVIDENCE.md` was not edited. K, n, tau, clamp, mortality, flow,",
        "layout, `K_AC`, `Jmax_A1`, `K_P`, `k_P`, and `alpha` were not retuned.",
        "Lentini is topology + Turing/replay controls only.",
        "",
        f"u SHA-256 `{NARMA_SHA}`.",
        f"Lentini PDF SHA-256 `{LENTINI_SHA}`.",
        f"A1 `g0=0`, `n_AC=2`, `K_AC=0.25`, `Jmax_A1={JMAX_A1:.12e}`.",
        f"TW `n_P=2`, `K_P={row['K_P']}`, `k_P={row['k_P']:.12e}`, `alpha={row['alpha']}`.",
        "Claim dish CENTER / `FLOW=0`. CSV 200 / 3200 / 3200; last sample",
        "`199;15;299.95`.",
        "",
        f"Occupancy: **{row['occupancy']}** (`mean_R={row['mean_R']:.4f}`, "
        f"`r(mean_R,u)={row['r_meanR_u']:.3f}`).",
        f"`mean_P_ac={row['mean_P_ac']:.4f}`, `mean_h={row['mean_h']:.4f}`, "
        f"`r(P_ac,u)={row['r_Pac_u']:.3f}`.",
        "",
        "## Claims (report, do not retune)",
        "",
        f"- **System** (CLOSED F408 vs Brownian and silent): **{system}**",
        f"  (`{row['CLOSED_F408']:.4f}` vs `{row['A1_brownian']:.4f}` / `{row['A1_silent_reused']:.4f}`).",
        f"- **Living-layer** (CLOSED F408 vs CLOSED field): **{living}**",
        f"  (`{row['CLOSED_F408']:.4f}` vs `{row['CLOSED_field']:.4f}`).",
        f"- **Transducer** |CLOSED F408 − A1 F408| = `{row['abs_CLOSED_minus_A1_F408']:.4f}`.",
        f"  **{row['transducer']}**. Do not promote the two-way dish on a 0.03 tick.",
        f"- **Interaction** CLOSED vs A1 (nonfunctional) AND CLOSED vs REPLAY:",
        f"  **{row['interaction']}**. |CLOSED−A1| = `{row['abs_CLOSED_minus_A1_F408']:.4f}`; "
        f"|CLOSED−REPLAY| = `{row['abs_CLOSED_minus_REPLAY_F408']:.4f}`.",
        f"  Replay tied closed: **{row['replay_tied_closed']}**.",
        "  Interaction PASS only if closed beats both by ≥ 0.03. If closed beats A1",
        "  but replay ties closed: `FEEDFORWARD_SCHEDULE`. If closed ties A1:",
        "  `NO_STORY_MOVE`, two-way untested-as-useful.",
        "- **Ceiling:** LINEAR_U_DELAY_10 remains `0.6828`. Biology is not",
        "  required to beat it. **`0.93` stays a weak predictor.**",
        "",
        "## Reused references",
        "",
        "| Readout | Test NRMSE |",
        "|---|---|",
        f"| A0 F408 | {row['A0_F408']:.4f} |",
        f"| A0 field | {row['A0_field']:.4f} |",
        f"| Brownian (A1 ≡ functional TW) | {row['A1_brownian']:.4f} |",
        f"| silent | {row['A1_silent_reused']:.4f} |",
        f"| A1 F408 (nonfunctional AC) | {row['A1_F408']:.4f} |",
        f"| A1 field | {row['A1_field']:.4f} |",
        f"| A1 10-tap g | {row['A1_G_DELAY_10']:.4f} |",
        "",
        "## Closed loop seed 111",
        "",
        "| Readout | Test NRMSE |",
        "|---|---|",
        f"| CLOSED F408 | {row['CLOSED_F408']:.4f} |",
        f"| CLOSED field | {row['CLOSED_field']:.4f} |",
        f"| occupancy-masked RL surrogate | {row['CLOSED_masked']:.4f} |",
        f"| REPLAY F408 | {replay} |",
        f"| REPLAY field | {row['REPLAY_field']:.4f} |" if math.isfinite(row["REPLAY_field"])
        else "| REPLAY field | n/a |",
        f"| AC-output g[n]*h[n] | {row['AC_GH_SCALAR']:.4f} |",
        f"| AC-output window-mean J | {row['AC_J_WINDOW_MEAN']:.4f} |",
        f"| AC-output 10-tap (g*h) | {row['AC_GH_DELAY_10']:.4f} |",
        "",
        "## Task baselines (unchanged u / target)",
        "",
        "| Baseline | Test NRMSE |",
        "|---|---|",
        f"| LINEAR_U_DELAY_10 | {row['LINEAR_U_DELAY_10']:.4f} |",
        f"| NARMA_INFORMED_INPUT | {row['NARMA_INFORMED_INPUT']:.4f} |",
        "",
        "Direct-input 0.6828 remains the task ceiling. A 10-tap of `g` or `g*h`",
        "is AC-output attribution, not a living-layer win.",
        "",
        "## Shuffle",
        "",
    ]
    if row["SKIP_SHUFFLE"]:
        lines.append(
            f"**SKIP_SHUFFLE.** |CLOSED−A1|={row['abs_CLOSED_minus_A1_F408']:.4f} "
            f"and |CLOSED−REPLAY|={row['abs_CLOSED_minus_REPLAY_F408']:.4f} "
            "are both < 0.03 (`NO_STORY_MOVE` on interaction). Time-shuffle BSim "
            "is not authorized."
        )
    else:
        lines.append(
            "Time-shuffle BSim is authorized "
            f"(|CLOSED−A1|={row['abs_CLOSED_minus_A1_F408']:.4f}, "
            f"|CLOSED−REPLAY|={row['abs_CLOSED_minus_REPLAY_F408']:.4f})."
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "Stop after seed 111. No best seed. Do not retune. Do not promote",
            "the two-way dish as a new claim dish on a 0.03 tick.",
            "Seeds 222/333 are not in this protocol.",
            "",
            f"GATE_EVIDENCE untouched: {row['GATE_EVIDENCE_untouched']}.",
            f"Java constraints: {row['java_ok']}.",
            "",
            "C1 remains **DEFER**. E0.3 remains **NO_STORY_MOVE**. Waveform1 FAIL,",
            "Waveform2c field≥driven, and L3 living-layer FAIL stay.",
            "MODEL_STATUS remains `HYPOTHETICAL_DESIGN_ENVELOPE`.",
            "",
        ]
    )
    (HERE / "results" / "E4_2_SCOUT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    with (HERE / "results" / "e42_scout.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    (HERE / "results" / "e42_scout.json").write_text(
        json.dumps(row, indent=2, default=str) + "\n", encoding="utf-8"
    )


def check_smoke(log_path: Path) -> None:
    text = log_path.read_text(encoding="utf-8")
    needed = (
        "MODEL_STATUS=HYPOTHETICAL_DESIGN_ENVELOPE",
        "digital_twin=false",
        "FLOW_SPEED=0.0",
        "CENTER",
        "h(0)=1.000000000000",
        "h(K_P)=1.250000000000",
        "alpha=0.500000",
        "coupling=closed",
    )
    ok = all(token in text for token in needed)
    print("SMOKE tw_frozen / P=0 identity / h(K_P) / CENTER / FLOW=0:", "PASS" if ok else "FAIL")
    if not ok:
        missing = [token for token in needed if token not in text]
        raise SystemExit(f"STOP smoke print missing {missing}")


def occupancy_gate(u) -> str:
    if not E.csv_ok(CLOSED):
        raise SystemExit("STOP CLOSED CSV incomplete")
    occ, _ = occupancy_from_run(CLOSED, u)
    print(
        "CLOSED occupancy", occ["occupancy"],
        f"mean_R={occ['mean_R']:.4f}",
        f"r(R,u)={occ['r_meanR_u']:.3f}",
        f"AHL_finite={occ['ahl_finite']}",
        f"AHL_nonneg={occ['ahl_nonneg']}",
    )
    if occ["occupancy"] == "DEAD":
        print("DEAD: keep the row; do not raise k_P/alpha/Jmax; do not run replay/shuffle/222/333")
    return occ["occupancy"]


def print_claims(row):
    print("MODEL_STATUS=HYPOTHETICAL_DESIGN_ENVELOPE digital_twin=false")
    print(f"occupancy={row['occupancy']} mean_R={row['mean_R']:.4f} mean_h={row['mean_h']:.4f}")
    print(f"system={'PASS' if row['system_pass'] else 'FAIL'} "
          f"living={'PASS' if row['living_layer_pass'] else 'FAIL'} "
          f"transducer={row['transducer']} interaction={row['interaction']}")
    print(f"replay_tied_closed={row['replay_tied_closed']} "
          f"SKIP_SHUFFLE={row['SKIP_SHUFFLE']}")
    print("REFUSE_RETUNE: K, n, tau, clamp, K_AC, Jmax, K_P, k_P, alpha are frozen")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke-log", default="")
    parser.add_argument("--occupancy-gate", action="store_true")
    parser.add_argument("--write-shuffle", action="store_true")
    args = parser.parse_args()
    frozen = load_frozen()
    if args.smoke_log:
        check_smoke(Path(args.smoke_log))
        return
    u, y, digest = N.load_target(HERE / "narma10_target.csv", HERE / "input_ahl_narma200.txt")
    if digest != NARMA_SHA:
        raise SystemExit(f"STOP u hash {digest}")
    pdf = HERE.parent / "HybridDish" / "external_analysis" / "oc6b00330.pdf"
    if sha256_file(pdf) != LENTINI_SHA:
        raise SystemExit("STOP Lentini PDF hash mismatch")
    if not java_ok():
        raise SystemExit("STOP Java constraints failed")
    if args.occupancy_gate:
        occupancy_gate(u)
        return
    if args.write_shuffle:
        write_shuffle(CLOSED / "p_ac_timeseries.txt", CLOSED / "p_ac_shuffled.txt")
        return
    row, occ, _series = analyze(u, y, frozen)
    write_report(row)
    print_claims(row)
    print(
        "CLOSED", occ["occupancy"],
        f"F408={row['CLOSED_F408']:.4f}",
        f"field={row['CLOSED_field']:.4f}",
        f"A1={row['A1_F408']:.4f}",
        f"REPLAY={row['REPLAY_F408']:.4f}",
    )


if __name__ == "__main__":
    main()
