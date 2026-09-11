#!/usr/bin/env python3
"""SweepMGE03: Mackey-Glass three-geometry scout. Claim sanity, then occupancy/ridge."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCHA = HERE.parent / "BSimReservoirPlanBenchA"
sys.path.insert(0, str(BENCHA))
import check_bencha as A  # noqa: E402

WASHOUT, TRAIN, TEST = 40, 110, 50
NUM_WINDOWS = 200
MG_SHA = "e06810ea70f61bbe4baf61939f93b1e1cae9d92d9eee8d757865c4935dbac780"
RX, RY = 20, 10
K_HILL = 1.6
D0, F0, G0 = 0.3306, 0.0229, 0.3077
SANITY_TOL = 1e-3

# Piecewise-constant AHL on the 20 s voxel samples: between sample t_{k-1} and
# t_k the field is held at AHL[t_{k-1}]. Stage 3B / Stage 5 ODEs:
#   dR/dt = (H(C) - R) / tau_R,  tau_R = 15 s, H = C^2 / (K^2 + C^2), K = 1.6
#   dL/dt = (R - L) / tau_L,     tau_L = 1500 s
# Exact update on a constant-H interval of length delta.


def hill(c):
    c = np.asarray(c, dtype=float)
    return c * c / (K_HILL * K_HILL + c * c)


def pearson(a, b):
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return 0.0 if den < 1e-15 else float(np.sum(a * b) / den)


def constant_nrmse(y, constant):
    y_test = y[WASHOUT + TRAIN :]
    pred = np.full_like(y_test, constant)
    return A.nrmse(y_test, pred)


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_voxel_arrays(path):
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    ahl_i = [header.index(f"AHL_uM_{i}") for i in range(200)]
    den_i = [header.index(f"Den_{i}") for i in range(200)]
    r_i = [header.index(f"Receiver_R_{i}") for i in range(200)]
    l_i = [header.index(f"Lum_Mean_{i}") for i in range(200)]
    win = header.index("Window")
    samp = header.index("Sample")
    tcol = header.index("TimeInWindow_s")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return {
        "window": data[:, win].astype(int),
        "sample": data[:, samp].astype(int),
        "t_in": data[:, tcol],
        "ahl": data[:, ahl_i],
        "den": data[:, den_i],
        "receiver": data[:, r_i],
        "lum": data[:, l_i],
        "absolute_time": data[:, win] * 300.0 + data[:, tcol],
    }


def window_mean_matrix(values, windows):
    out = []
    for window in range(NUM_WINDOWS):
        mask = windows == window
        out.append(np.mean(values[mask], axis=0))
    return np.vstack(out)


def occupancy(voxels, u):
    r_win = window_mean_matrix(voxels["receiver"], voxels["window"])
    a_win = window_mean_matrix(voxels["ahl"], voxels["window"])
    mean_r_w = np.mean(r_win, axis=1)
    mean_a_w = np.mean(a_win, axis=1)
    mean_r = float(np.mean(mean_r_w))
    r_u = pearson(mean_r_w, u)
    if abs(r_u) < 0.5 or mean_r < 0.05:
        label = "DEAD"
    elif mean_r > 0.8 or float(np.mean(np.mean(r_win > 0.5, axis=1))) > 0.8:
        label = "SATURATED"
    else:
        label = "ALIVE"
    return {
        "occupancy": label,
        "r_meanR_u": r_u,
        "mean_R": mean_r,
        "mean_AHL": float(np.mean(mean_a_w)),
        "mean_L": float(np.mean(voxels["lum"])),
    }


def integrate_surrogate(absolute_time, ahl):
    """Exact R/L update assuming AHL is piecewise-constant on 20 s samples."""
    R = np.zeros(ahl.shape[1], dtype=float)
    L = np.zeros(ahl.shape[1], dtype=float)
    out_r = np.empty_like(ahl)
    out_l = np.empty_like(ahl)
    out_r[0], out_l[0] = R, L
    for index in range(1, len(absolute_time)):
        delta = absolute_time[index] - absolute_time[index - 1]
        C = ahl[index - 1]
        h = hill(C)
        er = math.exp(-delta / 15.0)
        el = math.exp(-delta / 1500.0)
        old_r, old_l = R, L
        R = h + (old_r - h) * er
        L = h + (old_l - h) * el + (old_r - h) * 15.0 / (15.0 - 1500.0) * (er - el)
        out_r[index], out_l[index] = R, L
    return out_r, out_l


def ridge_score(X, y):
    metrics = A.evaluate_readout(np.asarray(X), np.asarray(y), list(range(NUM_WINDOWS)))
    return {
        "lambda": metrics["lambda"],
        "test_nrmse": metrics["test_nrmse"],
        "test_r2": metrics["test_r2"],
        "n_features": np.asarray(X).shape[1],
    }


def summary_stats(run_dir):
    path = Path(run_dir) / "window_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))

    def col(name):
        return np.array([float(row[name]) for row in rows], dtype=float)

    return {
        "mean_pop": float(np.mean(col("Population"))),
        "births": float(np.sum(col("Births"))),
        "clamp_deaths": float(np.sum(col("Clamp_Deaths"))),
        "acid_deaths": float(np.sum(col("Input_Driven_Deaths"))),
        "oob_deaths": float(np.sum(col("OOB_Deaths"))),
    }


def csv_ok(run_dir, expected_windows=200, expected_aux=3200):
    d = Path(run_dir)
    s = A.validate_csv(d / "window_summary.csv", expected_windows)
    r = A.validate_csv(d / "results.csv", expected_aux, check_last_sample=True)
    v = A.validate_csv(d / "voxels.csv", expected_aux, check_last_sample=True)
    return (
        s["row_count_pass"]
        and r["row_count_pass"]
        and v["row_count_pass"]
        and r.get("last_sample_pass")
        and v.get("last_sample_pass")
    )


def story_rules(occupancy_label, d, f, brownian, d0=D0, f0=F0, g0=G0):
    g = d - f
    fired = []
    if occupancy_label != "ALIVE":
        return fired, g
    if abs(f - f0) >= 0.020:
        fired.append("Carrier moved")
    if abs(d - d0) >= 0.050:
        fired.append("Living moved")
    if abs(g - g0) >= 0.050:
        fired.append("Gap moved")
    if d < f:
        fired.append("Living-layer flip")
    if d >= brownian:
        fired.append("System lost")
    return fired, g


def baselines(y):
    intercept = constant_nrmse(y, float(np.mean(y[WASHOUT:WASHOUT + TRAIN])))
    persistence = A.nrmse(
        y[WASHOUT + TRAIN :],
        y[WASHOUT + TRAIN - 1 : WASHOUT + TRAIN + TEST - 1],
    )
    return intercept, persistence


def sanity():
    u, y = A.load_target(HERE / "mg_target.csv", "mg")
    digest = A.sha256_u(u)
    if digest != MG_SHA:
        raise SystemExit(f"ABORT u hash {digest}")
    ahl_lines = [
        line.strip()
        for line in (HERE / "input_ahl_mg200.txt").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    ahl = np.array([float(line) for line in ahl_lines], dtype=float)
    if not np.allclose(ahl[:NUM_WINDOWS], u, rtol=0, atol=1e-12):
        raise SystemExit("ABORT copied AHL file does not match mg_target u")
    driven_run = A.read_run(BENCHA / "results" / "mg_driven_seed101", "driven", NUM_WINDOWS, 3200)
    brown_run = A.read_run(BENCHA / "results" / "mg_brownian_seed101", "brownian", NUM_WINDOWS, 3200)
    silent_run = A.read_run(BENCHA / "results" / "mg_silent_seed101", "silent", NUM_WINDOWS, 3200)
    d = A.evaluate_readout(driven_run["matrix"]["X"], y, driven_run["matrix"]["windows"])
    f = A.evaluate_readout(driven_run["field"]["X"], y, driven_run["field"]["windows"])
    b = A.evaluate_readout(brown_run["matrix"]["X"], y, brown_run["matrix"]["windows"])
    s = A.evaluate_readout(silent_run["matrix"]["X"], y, silent_run["matrix"]["windows"])
    intercept, persistence = baselines(y)
    checks = {
        "driven_F408": (d["test_nrmse"], 0.3306),
        "field": (f["test_nrmse"], 0.0229),
        "brownian": (b["test_nrmse"], 1.0415),
        "silent": (s["test_nrmse"], 1.0416),
    }
    print("CLAIM SANITY seed 101 (reused BenchA MG voxels)")
    failed = False
    for name, (got, expected) in checks.items():
        ok = abs(got - expected) <= SANITY_TOL
        print(f"  {name:12} {got:.4f} expected {expected:.4f} {'PASS' if ok else 'FAIL'}")
        failed = failed or not ok
    print(f"  persistence  {persistence:.4f}  (BenchA pack ~0.1799)")
    print(f"  intercept    {intercept:.4f}  (train-mean predictor ~1.0416)")
    print("Do not report BenchA driven MC 15.331 as CHARC memory capacity.")
    if failed:
        raise SystemExit("STOP claim sanity failed; do not interpret new geometries")
    print("SANITY PASS")
    return u, y, intercept, persistence


GEOMETRIES = [
    ("PRI_CENTER_F0p0", 500.0, 0.0, True),
    ("PRI_CENTER_F0p25", 500.0, 0.25, False),
    ("PRI_UPSTREAM_CENTER_F0p0", 200.0, 0.0, False),
]
CONFIRM_SEEDS = (202, 303)


def driven_dir(cid, seed):
    if cid == "PRI_CENTER_F0p0":
        return BENCHA / "results" / f"mg_driven_seed{seed}"
    if cid == "PRI_CENTER_F0p25":
        return HERE / "results" / f"mge03_center_f0p25_driven_seed{seed}"
    if cid == "PRI_UPSTREAM_CENTER_F0p0":
        return HERE / "results" / f"mge03_upstream_f0p0_driven_seed{seed}"
    raise ValueError(cid)


def null_dir(flow, arm, seed):
    if flow == 0.0:
        return BENCHA / "results" / f"mg_{arm}_seed{seed}"
    return HERE / "results" / f"mge03_center_f{str(flow).replace('.', 'p')}_{arm}_seed{seed}"


def score_condition(cid, x, flow, reused, seed, u, y, intercept, persistence, ref=None):
    ddir = Path(driven_dir(cid, seed))
    sdir = Path(null_dir(flow, "silent", seed))
    bdir = Path(null_dir(flow, "brownian", seed))
    for path in (ddir, sdir, bdir):
        if not (path / "voxels.csv").exists():
            return None, str(path)
        if not csv_ok(path):
            raise SystemExit(f"STOP CSV completeness failed for {path}")
    voxels = load_voxel_arrays(ddir / "voxels.csv")
    occ = occupancy(voxels, u)
    stats = summary_stats(ddir)
    driven_run = A.read_run(ddir, "driven", NUM_WINDOWS, 3200)
    silent_run = A.read_run(sdir, "silent", NUM_WINDOWS, 3200)
    brown_run = A.read_run(bdir, "brownian", NUM_WINDOWS, 3200)
    driven = A.evaluate_readout(driven_run["matrix"]["X"], y, driven_run["matrix"]["windows"])
    field = A.evaluate_readout(driven_run["field"]["X"], y, driven_run["field"]["windows"])
    silent = A.evaluate_readout(silent_run["matrix"]["X"], y, silent_run["matrix"]["windows"])
    brown = A.evaluate_readout(brown_run["matrix"]["X"], y, brown_run["matrix"]["windows"])
    sr, sl = integrate_surrogate(voxels["absolute_time"], voxels["ahl"])
    mask = (voxels["den"] != 0).astype(float)
    unmasked = np.column_stack([
        window_mean_matrix(sr, voxels["window"]),
        window_mean_matrix(sl, voxels["window"]),
    ])
    masked = np.column_stack([
        window_mean_matrix(sr * mask, voxels["window"]),
        window_mean_matrix(sl * mask, voxels["window"]),
    ])
    d = driven["test_nrmse"]
    f = field["test_nrmse"]
    brown_nrmse = brown["test_nrmse"]
    if ref is None:
        d0, f0, g0 = D0, F0, G0
        ref_name = "claim_seed101"
    else:
        d0, f0, g0 = ref
        ref_name = "own_CENTER_F0p0"
    fired, g = story_rules(occ["occupancy"], d, f, brown_nrmse, d0, f0, g0)
    transit = float("inf") if flow == 0 else 1000.0 / flow
    row = {
        "ConditionID": cid,
        "Seed": seed,
        "ReusedClaimRun": reused and cid == "PRI_CENTER_F0p0",
        "AHL_x": x,
        "AHL_y": 250.0,
        "AHL_z": 5.0,
        "Flow_um_s": flow,
        "Boundary": "NO_FLUX" if flow == 0 else "OUTFLOW",
        "Transit_s": transit,
        "Courant": flow * 0.05 / 20.0,
        **occ,
        **stats,
        "driven_F408_nrmse": d,
        "driven_F408_lambda": driven["lambda"],
        "field_nrmse": f,
        "field_lambda": field["lambda"],
        "gap_G": g,
        "unmasked_RL_nrmse": ridge_score(unmasked, y)["test_nrmse"],
        "masked_RL_nrmse": ridge_score(masked, y)["test_nrmse"],
        "silent_nrmse": silent["test_nrmse"],
        "brownian_nrmse": brown_nrmse,
        "intercept_nrmse": intercept,
        "persistence_nrmse": persistence,
        "story_ref": ref_name,
        "story_move_rules": "|".join(fired) if fired else "NO_STORY_MOVE",
    }
    print(
        f"seed {seed}",
        cid,
        occ["occupancy"],
        f"R={occ['mean_R']:.4f}",
        f"r(R,u)={occ['r_meanR_u']:.4f}",
        f"D={d:.4f}",
        f"F={f:.4f}",
        f"G={g:.4f}",
        row["story_move_rules"],
    )
    return row, None


def score_seed(seed, u, y, intercept, persistence, geometries=None):
    rows = []
    missing = []
    own_ref = None
    for cid, x, flow, reused in (geometries or GEOMETRIES):
        use_ref = None if seed == 101 else own_ref
        row, miss = score_condition(
            cid, x, flow, reused, seed, u, y, intercept, persistence, ref=use_ref
        )
        if miss:
            missing.append(miss)
            print(f"missing {miss}")
            continue
        if cid == "PRI_CENTER_F0p0":
            own_ref = (row["driven_F408_nrmse"], row["field_nrmse"], row["gap_G"])
            if seed != 101:
                row["story_move_rules"] = "REFERENCE"
                row["story_ref"] = "own_CENTER_F0p0"
        rows.append(row)
    return rows, missing


def analyze(u, y, intercept, persistence):
    rows_101, missing_101 = score_seed(101, u, y, intercept, persistence)
    if missing_101 or len(rows_101) < 3:
        raise SystemExit("STOP missing production voxels; do not interpret incomplete geometries")
    any_story = any(row["story_move_rules"] != "NO_STORY_MOVE" for row in rows_101)
    verdict = "STORY_MOVE" if any_story else "NO_STORY_MOVE"
    all_rows = list(rows_101)
    missing = list(missing_101)
    if any_story:
        fired_ids = [
            row["ConditionID"] for row in rows_101 if row["story_move_rules"] != "NO_STORY_MOVE"
        ]
        confirm_geoms = [g for g in GEOMETRIES if g[0] in fired_ids or g[0] == "PRI_CENTER_F0p0"]
        for seed in CONFIRM_SEEDS:
            extra, miss = score_seed(seed, u, y, intercept, persistence, confirm_geoms)
            all_rows.extend(extra)
            missing.extend(miss)
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    csv_path = out / "mg_geom_scout.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    print("wrote", csv_path)
    write_scout(all_rows, intercept, persistence, verdict, missing)
    if any_story:
        print("STORY_MOVE on seed 101. Confirmation seeds 202/303 are authorized")
        print("only for the conditions that fired, plus claim CENTER/0 already on disk.")
        print("Do not freeze a new claim dish in this job. No best-seed.")
        if any("seed202" in path or "seed303" in path for path in missing):
            print("Confirmation voxels still missing; run run_mge03_confirm_jobs.py")
    else:
        print("NO_STORY_MOVE. Do not start seeds 202/303. Do not promote a new claim dish.")
        print("0.26 stays a weak one-step predictor vs field 0.023 and persistence ~0.18.")
    return all_rows, verdict


def fmt(value, digits=4):
    if isinstance(value, float) and not math.isfinite(value):
        return "inf"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def table_row(row):
    return (
        "| {ConditionID} | {occupancy} | {mean_R} | {r_meanR_u} | {mean_AHL} | {mean_pop} | "
        "{driven_F408_nrmse} | {field_nrmse} | {gap_G} | {brownian_nrmse} | "
        "{silent_nrmse} | {unmasked_RL_nrmse} | {masked_RL_nrmse} | {story_move_rules} |".format(
            ConditionID=row["ConditionID"],
            occupancy=row["occupancy"],
            mean_R=fmt(row["mean_R"]),
            r_meanR_u=fmt(row["r_meanR_u"]),
            mean_AHL=fmt(row["mean_AHL"]),
            mean_pop=fmt(row["mean_pop"], 1),
            driven_F408_nrmse=fmt(row["driven_F408_nrmse"]),
            field_nrmse=fmt(row["field_nrmse"]),
            gap_G=fmt(row["gap_G"]),
            brownian_nrmse=fmt(row["brownian_nrmse"]),
            silent_nrmse=fmt(row["silent_nrmse"]),
            unmasked_RL_nrmse=fmt(row["unmasked_RL_nrmse"]),
            masked_RL_nrmse=fmt(row["masked_RL_nrmse"]),
            story_move_rules=row["story_move_rules"],
        )
    )


def seed101_interpretation(rows, persistence):
    lines = ["## Seed 101 reading", ""]
    claim = next(row for row in rows if row["ConditionID"] == "PRI_CENTER_F0p0")
    lines.append(
        f"Claim CENTER/0: D={fmt(claim['driven_F408_nrmse'])}, "
        f"F={fmt(claim['field_nrmse'])}, G={fmt(claim['gap_G'])}, "
        f"mean_R={fmt(claim['mean_R'])}."
    )
    lines.append("")
    field_moved = False
    living_moved = False
    for row in rows:
        if row["ConditionID"] == "PRI_CENTER_F0p0":
            continue
        d_f = abs(row["field_nrmse"] - F0)
        d_d = abs(row["driven_F408_nrmse"] - D0)
        field_moved = field_moved or d_f >= 0.020
        living_moved = living_moved or d_d >= 0.050
        vs_p = (
            "still worse than persistence"
            if row["driven_F408_nrmse"] > persistence
            else "beats persistence"
        )
        lines.append(
            f"- {row['ConditionID']}: occupancy {row['occupancy']}, "
            f"mean_R={fmt(row['mean_R'])}, D={fmt(row['driven_F408_nrmse'])} "
            f"({vs_p} {fmt(persistence)}), F={fmt(row['field_nrmse'])} "
            f"(|F-F0|={d_f:.4f}), G={fmt(row['gap_G'])}. "
            f"Unmasked surrogate {fmt(row['unmasked_RL_nrmse'])}; "
            f"occupancy-masked {fmt(row['masked_RL_nrmse'])}."
        )
    lines.append("")
    if living_moved and not field_moved:
        lines.extend([
            "Field did not meet the carrier-moved threshold. Driven NRMSE moved.",
            "The gap shrank because the living layer changed, not because mild",
            "flow smeared the spatial delay line. That is a result. It is not",
            "a living-layer flip (D remains >> F). Unmasked kinetic surrogate",
            "stays with the field; occupancy-masked stays with driven. Driven",
            "still does not beat persistence or field. Do not promote.",
            "",
        ])
    elif field_moved and living_moved:
        lines.extend([
            "Both field and driven moved. Layout is affecting transport and the",
            "living readout together. Recorded; not a new claim dish.",
            "",
        ])
    elif field_moved and not living_moved:
        lines.extend([
            "Field moved while driven held. Consistent with smearing the spatial",
            "delay line more than L. Recorded; not a new claim dish.",
            "",
        ])
    return lines


def write_scout(rows, intercept, persistence, verdict, missing):
    lines = [
        "# Mackey–Glass three-geometry scout (SweepMGE03)",
        "",
        "Dated 2026-08-18. Development only. Not a new claim dish. Does not",
        "rewrite BenchA Mackey–Glass Overall PASS. Not E0.3 NARMA and not the",
        "nine-condition factorial. C1 stays DEFER. Waveform1 FAIL and WaveformS",
        "stay as written. Two-way is future work.",
        "",
        "Protocol (frozen before scores):",
        "`examples/BSimReservoirPlanSweepMGE03/PROTOCOL.md`",
        "",
        "CSV:",
        "`examples/BSimReservoirPlanSweepMGE03/results/mg_geom_scout.csv`",
        "",
        "## Scientific question",
        "",
        "At fixed commanded MG payload, do source position and mild +x flow",
        "change occupancy and the driven–field gap relative to the BenchA",
        "claim dish?",
        "",
        "Hypothesis (predeclared): MG is carrier-rich (field NRMSE 0.0229 vs",
        "driven 0.3306 on seed 101). Mild flow may smear the spatial delay",
        "line more than it smears L. If field gets worse while driven holds,",
        "the gap shrinks. If both arms move together, layout is transport,",
        "not living computation. Either outcome is a result.",
        "",
        "Field-only is a diagnostic, not a kill switch. Persistence of x",
        "(~0.180) is a mandatory trivial baseline. Do not report BenchA",
        "driven MC 15.331 as CHARC memory capacity (smooth correlated drive).",
        "",
        "## Frozen geometries",
        "",
        "| ConditionID | AHL position | Flow | Boundary |",
        "|---|---|---|---|",
        "| PRI_CENTER_F0p0 | (500, 250, 5) µm | 0 µm/s | NO_FLUX |",
        "| PRI_CENTER_F0p25 | (500, 250, 5) µm | 0.25 µm/s | OUTFLOW |",
        "| PRI_UPSTREAM_CENTER_F0p0 | (200, 250, 5) µm | 0 µm/s | NO_FLUX |",
        "",
        "Acid stays at (300, 375, 5). Attractant/repellent stay silent.",
        "Kinetics were not retuned. Rate was not scaled with flow or position.",
        "",
        "## Target-only baselines (identical for every geometry)",
        "",
        f"- persistence x[n] → x[n+1] on the closed split: {persistence:.4f}",
        f"- intercept / train-mean predictor: {intercept:.4f}",
        "",
        "Teacher-forced one-step. Not Jaeger free-run. NRMSE = RMSE / pop-std",
        "of the test target. Independent lambda per condition×arm.",
        "",
        "## Surrogate approximation",
        "",
        "Unmasked and occupancy-masked kinetic surrogates integrate the frozen",
        "Stage 3B / Stage 5 ODEs (tau_R=15, K=1.6, n=2, tau_L=1500) on driven",
        "voxel AHL. Between 20 s samples the AHL concentration is held piecewise",
        "constant at the previous sample. Occupancy mask zeros R and L in empty",
        "density voxels. No extra BSim.",
        "",
        "## Seed 101 rows",
        "",
        "| ConditionID | occ | mean_R | r(R,u) | mean_AHL | mean_pop | D | F | G | Brownian | silent | unmasked | masked | rules |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    seed101 = [row for row in rows if int(row["Seed"]) == 101]
    confirm = [row for row in rows if int(row["Seed"]) != 101]
    for row in seed101:
        lines.append(table_row(row))
    lines.extend(["", f"## Verdict: {verdict}", ""])
    if verdict == "NO_STORY_MOVE":
        lines.extend([
            "No ALIVE condition met a predeclared story-move rule on seed 101.",
            "Seeds 202/303 were not started. No new claim dish. Ticks smaller",
            "than the thresholds are the same story, different 0.26.",
            "",
        ])
    else:
        fired = [
            row["ConditionID"] for row in seed101 if row["story_move_rules"] != "NO_STORY_MOVE"
        ]
        lines.extend([
            "Story-move fired on seed 101: " + ", ".join(fired) + ".",
            "Seeds 202/303 are confirmation only for those conditions, plus the",
            "claim CENTER/0 ridge already on disk. No best-seed. Do not freeze",
            "a new claim dish in this job.",
            "",
        ])
        lines.extend(seed101_interpretation(seed101, persistence))
    if confirm:
        lines.extend([
            "## Confirmation seeds 202 / 303 (not a claim freeze)",
            "",
            "Reported because seed 101 fired. Confirmation story-move rules are",
            "versus that seed's own CENTER/0, not versus seed-101 D0=0.3306.",
            "Not a best-seed. Not a new dish.",
            "",
            "| Seed | ConditionID | occ | D | F | G | Brownian | silent | unmasked | masked | ref | rules |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|",
        ])
        for row in confirm:
            lines.append(
                "| {Seed} | {ConditionID} | {occupancy} | {D} | {F} | {G} | {B} | {S} | {U} | {M} | {ref} | {R} |".format(
                    Seed=row["Seed"],
                    ConditionID=row["ConditionID"],
                    occupancy=row["occupancy"],
                    D=fmt(row["driven_F408_nrmse"]),
                    F=fmt(row["field_nrmse"]),
                    G=fmt(row["gap_G"]),
                    B=fmt(row["brownian_nrmse"]),
                    S=fmt(row["silent_nrmse"]),
                    U=fmt(row["unmasked_RL_nrmse"]),
                    M=fmt(row["masked_RL_nrmse"]),
                    ref=row.get("story_ref", ""),
                    R=row["story_move_rules"],
                )
            )
        lines.append("")
        claim_rows = [row for row in confirm if row["ConditionID"] == "PRI_CENTER_F0p0"]
        new_rows = [row for row in confirm if row["ConditionID"] != "PRI_CENTER_F0p0"]
        if claim_rows:
            ds = ", ".join(f"{row['Seed']}={fmt(row['driven_F408_nrmse'])}" for row in claim_rows)
            lines.extend([
                "Claim CENTER/0 driven NRMSE on confirmation seeds: " + ds + ".",
                "Seed 101 claim D=0.3306 is the high seed of the BenchA trio.",
                "New geometries on 202/303 are compared to that seed's own CENTER/0.",
                "",
            ])
        if new_rows and all(row["story_move_rules"] == "NO_STORY_MOVE" for row in new_rows):
            lines.extend([
                "Confirmation: no new geometry meets a story-move threshold versus",
                "its own CENTER/0. The seed-101 living/gap move does not survive",
                "seed-matched claim comparison. Do not promote.",
                "",
            ])
    if missing:
        lines.extend(["Missing production paths (not interpreted):", ""])
        lines.extend(f"- `{path}`" for path in missing)
        lines.append("")
    lines.extend([
        "Do not select on NRMSE alone. Do not edit BenchA Overall PASS.",
        "Do not promote CENTER/0.25 or UPSTREAM to the claim dish.",
        "Do not quote MG MC 15.331 as a capacity result.",
        "",
    ])
    path = HERE / "results" / "MG_GEOM_SCOUT.md"
    path.parent.mkdir(exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    print("wrote", path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sanity", action="store_true")
    args = parser.parse_args()
    u, y, intercept, persistence = sanity()
    if args.sanity:
        return
    analyze(u, y, intercept, persistence)


if __name__ == "__main__":
    main()
