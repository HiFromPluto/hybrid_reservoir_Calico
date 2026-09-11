#!/usr/bin/env python3
"""E5.2 Module 0: occupancy-first encoding freeze. No F408. No Java."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MAT = REPO / "examples" / "HybridDish" / "external_analysis" / "synthetic_biomarker_BSim_inputs.mat"
PREFLIGHT = REPO / "examples" / "HybridDish" / "design_space_preflight"
E51 = REPO / "examples" / "BSimReservoirPlanE51CRP"
INPUT = HERE / "input"
RESULTS = HERE / "results"
MAPS_DIR = INPUT / "maps"

DEV = [1, 4, 6, 7, 9, 13, 16, 18, 20, 22, 24, 25]
CONF = [8, 12, 14, 15, 17, 19, 21, 23]
DEV_SHA = "83acc7102e532112c961bb45a1327504290ffddea1c06e659bfb5e68be032915"
CONF_SHA = "c90cca8f5627c737b38024b4863415e46a384d8715e6d06a5823e8cb93ba3787"
LAMBDA_GRID = np.array([1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6])
T = 48
K_HILL = 1.6
OCC_DEAD_R = 0.05
OCC_DEAD_CORR = 0.5
SAT_MEDIAN_R = 0.40
SAT_FRAC_R = 0.50
JMAX = 1.28e8

MAPS = (
    {"id": "M_E51", "u_floor": 0.00, "g": 0.50, "notes": "E5.1 replay (mostly DEAD). Reference."},
    {"id": "M_F15", "u_floor": 0.15, "g": 0.35, "notes": "Lift floor, keep contrast"},
    {"id": "M_F20", "u_floor": 0.20, "g": 0.30, "notes": "Mid"},
    {"id": "M_F25", "u_floor": 0.25, "g": 0.25, "notes": "Highest floor, least CRP range"},
)


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_u12(u) -> str:
    payload = ",".join(f"{float(v):.12f}" for v in u)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def write_lines(path: Path, values, fmt="{:.12f}"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(fmt.format(float(v)) + "\n" for v in values), encoding="utf-8")


def clip01_half(u):
    return np.clip(np.asarray(u, dtype=float), 0.0, 0.5)


def encode(x0, u_floor, g):
    return clip01_half(u_floor + g * np.asarray(x0, dtype=float))


def pearson(a, b):
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return 0.0 if den < 1e-15 else float(np.sum(a * b) / den)


def nrmse_pop(ytrue, ypred):
    ytrue = np.asarray(ytrue, dtype=float).ravel()
    ypred = np.asarray(ypred, dtype=float).ravel()
    rmse = np.sqrt(np.mean((ytrue - ypred) ** 2))
    denom = ytrue.std(ddof=0)
    return float(rmse / denom) if denom > 0 else float("inf")


def patient_pooled(ytrue, ypred, pids):
    per = []
    for p in np.unique(pids):
        m = pids == p
        per.append(nrmse_pop(ytrue[m], ypred[m]))
    return float(np.mean(per)), [float(v) for v in per]


def standardize_fit(X):
    mu = X.mean(axis=0)
    sg = X.std(axis=0, ddof=0)
    sg = np.where(sg < 1e-12, 1.0, sg)
    return (X - mu) / sg, mu, sg


def ridge_fit(X, y, lam):
    n = X.shape[1]
    Xb = np.column_stack([np.ones(len(X)), X])
    G = Xb.T @ Xb
    G[1:, 1:] = G[1:, 1:] + lam * np.eye(n)
    return np.linalg.solve(G, Xb.T @ y)


def ridge_predict(X, w):
    return np.column_stack([np.ones(len(X)), X]) @ w


def delay_rows(Xnorm, patients, t_lo, t_hi, channels):
    feats, targets, pids = [], [], []
    n_lag = 10
    for p in patients:
        for t in range(t_lo, t_hi + 1):
            row = []
            for b in channels:
                for k in range(n_lag):
                    row.append(Xnorm[p, t - k, b])
            feats.append(row)
            targets.append(Xnorm[p, t + 1, 0])
            pids.append(p)
    return np.asarray(feats, float), np.asarray(targets, float), np.asarray(pids, int)


def persist_rows(Xnorm, patients, t_lo, t_hi):
    x, y, pids = [], [], []
    for p in patients:
        for t in range(t_lo, t_hi + 1):
            x.append(Xnorm[p, t, 0])
            y.append(Xnorm[p, t + 1, 0])
            pids.append(p)
    return np.asarray(x, float), np.asarray(y, float), np.asarray(pids, int)


def loo_lambda(X, y, pids, patients):
    best = np.inf
    best_lam = LAMBDA_GRID[-1]
    grid = []
    for lam in LAMBDA_GRID:
        scores = []
        for held in patients:
            tr = pids != held
            va = pids == held
            Xtr_z, mu, sg = standardize_fit(X[tr])
            w = ridge_fit(Xtr_z, y[tr], lam)
            pred = ridge_predict((X[va] - mu) / sg, w)
            scores.append(nrmse_pop(y[va], pred))
        pooled = float(np.mean(scores))
        grid.append({"lambda": float(lam), "loo_pooled": pooled})
        if pooled < best - 1e-15 or (abs(pooled - best) <= 1e-15 and lam > best_lam):
            best = pooled
            best_lam = lam
    return float(best_lam), float(best), grid


def fit_score(Xtr, ytr, ptr, Xte, yte, pte, lam):
    Xtr_z, mu, sg = standardize_fit(Xtr)
    w = ridge_fit(Xtr_z, ytr, lam)
    pred_tr = ridge_predict(Xtr_z, w)
    pred_te = ridge_predict((Xte - mu) / sg, w)
    tr_p, _ = patient_pooled(ytr, pred_tr, ptr)
    te_p, te_per = patient_pooled(yte, pred_te, pte)
    return {
        "lambda": float(lam),
        "development_pooled": tr_p,
        "confirmation_pooled": te_p,
        "confirmation_concat": nrmse_pop(yte, pred_te),
        "confirmation_per_patient": te_per,
    }


def persist_score(Xnorm, patients, t_lo, t_hi):
    x, y, pids = persist_rows(Xnorm, patients, t_lo, t_hi)
    pooled, per = patient_pooled(y, x, pids)
    return {
        "pooled": pooled,
        "concat": nrmse_pop(y, x),
        "per_patient": per,
        "n_pairs": int(len(y)),
    }


def algebraic_hill(c):
    c = np.asarray(c, dtype=float)
    return c * c / (K_HILL * K_HILL + c * c)


def load_e51_window_fields(patient: int):
    path = E51 / "results" / f"e51_driven_p{patient:04d}_seed111" / "voxels.csv"
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    den_i = [header.index(f"Den_{i}") for i in range(200)]
    ahl_i = [header.index(f"AHL_uM_{i}") for i in range(200)]
    win_i = header.index("Window")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    windows = data[:, win_i].astype(int)
    den = np.zeros((T, 200))
    ahl = np.zeros(T)
    for window in range(T):
        mask = windows == window
        den[window] = data[mask][:, den_i].mean(axis=0)
        ahl[window] = float(np.mean(data[mask][:, ahl_i]))
    return {"den": den, "ahl_mean": ahl, "path": str(path)}


def occupancy_worker(payload: dict) -> dict:
    """One map × patient occupancy screen. Imported under a process pool."""
    sys.path.insert(0, str(PREFLIGHT))
    import run_transport_screen as ts

    u = np.asarray(payload["u"], dtype=float)
    den = np.asarray(payload["den"], dtype=float)
    ahl_e51 = np.asarray(payload["ahl_e51"], dtype=float)
    u_e51 = np.asarray(payload["u_e51"], dtype=float)
    map_id = payload["map_id"]
    patient = payload["patient"]
    path_used = "E0.2_TRANSPORT_STAGE3B"

    try:
        condition = next(c for c in ts.make_conditions() if c.condition_id == "PRI_CENTER_F0p0")
        state = ts.initial_state()
        u0 = float(u[0])
        for _ in range(int(ts.WARMUP_S / ts.WINDOW_S)):
            state, _ = ts.simulate_window(state, condition, u0, collect=False)
        _, samples = ts.simulate_drive(state, condition, u, clear_budget=True, collect=True)
        assert samples is not None
        n_per = len(ts.SAMPLE_OFFSETS)
        r_win = np.zeros(T)
        r_uniform = np.zeros(T)
        for window in range(T):
            chunk = samples[window * n_per : (window + 1) * n_per]
            r_readout = ts.readout_state(chunk, ts.N).mean(axis=0)
            denom = float(np.sum(den[window]))
            r_win[window] = float(np.sum(r_readout * den[window]) / denom) if denom > 0 else float(np.mean(r_readout))
            r_uniform[window] = float(np.mean(r_readout))
    except Exception as exc:  # noqa: BLE001 — documented fallback
        path_used = f"FALLBACK_AHL_RESCALE:{type(exc).__name__}:{exc}"
        ahl_new = ahl_e51 * u / np.maximum(u_e51, 1e-6)
        r_win = algebraic_hill(ahl_new)
        r_uniform = r_win.copy()

    mean_r = float(np.mean(r_win))
    r_u = pearson(r_win, u)
    alive = mean_r >= OCC_DEAD_R and abs(r_u) >= OCC_DEAD_CORR
    frac_high = float(np.mean(r_win > 0.5))
    return {
        "map_id": map_id,
        "patient": patient,
        "path": path_used,
        "mean_R": mean_r,
        "median_R_win": float(np.median(r_win)),
        "r_meanR_u": r_u,
        "frac_R_gt_0_5": frac_high,
        "mean_R_uniform": float(np.mean(r_uniform)),
        "occupancy": "ALIVE" if alive else "DEAD",
        "u0": float(u[0]),
        "u_min": float(np.min(u)),
        "u_max": float(np.max(u)),
        "R_win": [float(v) for v in r_win],
    }


def run_occupancy_screen(u_by_map: dict, fields: dict, u_e51: dict) -> list[dict]:
    jobs = []
    for spec in MAPS:
        for p in DEV:
            jobs.append(
                {
                    "map_id": spec["id"],
                    "patient": p,
                    "u": u_by_map[spec["id"]][p].tolist(),
                    "den": fields[p]["den"].tolist(),
                    "ahl_e51": fields[p]["ahl_mean"].tolist(),
                    "u_e51": u_e51[p].tolist(),
                }
            )
    rows = []
    workers = min(8, max(1, len(jobs)))
    print(f"occupancy screen jobs={len(jobs)} workers={workers}", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(occupancy_worker, job): job for job in jobs}
        for fut in as_completed(futures):
            row = fut.result()
            rows.append(row)
            print(
                f"  {row['map_id']} p{row['patient']:04d} {row['occupancy']} "
                f"mean_R={row['mean_R']:.4f} r={row['r_meanR_u']:.3f} path={row['path'][:40]}",
                flush=True,
            )
    order = {(spec["id"], p): i for i, spec in enumerate(MAPS) for p in DEV}
    rows.sort(key=lambda r: order[(r["map_id"], r["patient"])])
    return rows


def summarize_maps(rows: list[dict], u_by_map: dict) -> list[dict]:
    summaries = []
    for spec in MAPS:
        sub = sorted((r for r in rows if r["map_id"] == spec["id"]), key=lambda r: r["patient"])
        mean_rs = [r["mean_R"] for r in sub]
        all_win = np.concatenate([np.asarray(r["R_win"], float) for r in sub])
        u_cat = np.concatenate([u_by_map[spec["id"]][p] for p in DEV])
        n_alive = sum(1 for r in sub if r["occupancy"] == "ALIVE")
        saturated = (float(np.median(mean_rs)) > SAT_MEDIAN_R) or (float(np.mean(all_win > 0.5)) > SAT_FRAC_R)
        summaries.append(
            {
                **spec,
                "n_ALIVE": n_alive,
                "n_DEAD": 12 - n_alive,
                "median_mean_R": float(np.median(mean_rs)),
                "frac_windows_R_gt_0_5": float(np.mean(all_win > 0.5)),
                "saturated": bool(saturated),
                "u_contrast": float(np.max(u_cat) - np.min(u_cat)),
                "u_min": float(np.min(u_cat)),
                "u_max": float(np.max(u_cat)),
                "path": sub[0]["path"] if sub else "",
                "patients": sub,
            }
        )
    return summaries


def select_winner(summaries: list[dict]) -> dict | None:
    eligible = [s for s in summaries if not s["saturated"]]
    if not eligible:
        return None
    best_n = max(s["n_ALIVE"] for s in eligible)
    tied = [s for s in eligible if s["n_ALIVE"] == best_n]
    if best_n < 10:
        return None
    tied.sort(key=lambda s: (-s["u_contrast"], s["u_floor"], s["id"]))
    return tied[0]


def write_encoding_freeze(summaries, winner, hashes, screen_meta):
    lines = [
        "# E5.2 encoding freeze",
        "",
        "Occupancy screen only. No F408. No living ridge. `COHORT_STATUS: SYNTHETIC_COHORT`.",
        "Not clinical. `Jmax` stays `1.28e8`. No fifth map. E5.1 `NOT_SCORED` is unchanged.",
        "",
        f"Screen path: `{screen_meta['path_summary']}`",
        f"Density mask: `{screen_meta['density_mask']}`",
        f"Transport status: `{screen_meta['transport_status']}`",
        "",
        "Formula:",
        "",
        "```",
        "u[n] = clip(u_floor + g * Xnorm[p, n, 0], 0, 0.5)",
        "```",
        "",
        "ALIVE if `mean_R ≥ 0.05` and `|r(R_win, u)| ≥ 0.5`.",
        "SATURATED if median `mean_R > 0.40` or fraction of windows with `R > 0.5` exceeds 0.50.",
        "",
        "## Map table (development 12)",
        "",
        "| MapID | u_floor | g | n_ALIVE | median mean_R | frac R>0.5 | SATURATED | u contrast | Notes |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for s in summaries:
        lines.append(
            f"| {s['id']} | {s['u_floor']:.2f} | {s['g']:.2f} | {s['n_ALIVE']}/12 | "
            f"{s['median_mean_R']:.4f} | {s['frac_windows_R_gt_0_5']:.3f} | "
            f"{'yes' if s['saturated'] else 'no'} | {s['u_contrast']:.4f} | {s['notes']} |"
        )
    lines.extend(["", "## Per-patient occupancy", ""])
    lines.append("| MapID | Patient | Occupancy | mean_R | r(R,u) | u[0] | uniform mean_R |")
    lines.append("|---|---|---|---|---|---|---|")
    for s in summaries:
        for r in s["patients"]:
            lines.append(
                f"| {r['map_id']} | {r['patient']} | {r['occupancy']} | {r['mean_R']:.4f} | "
                f"{r['r_meanR_u']:.3f} | {r['u0']:.4f} | {r['mean_R_uniform']:.4f} |"
            )
    lines.extend(["", "## Selection", ""])
    if winner is None:
        lines.extend(
            [
                "**ENCODING_FAIL.** No non-saturated map reached `n_ALIVE ≥ 10`.",
                "Do not raise `Jmax`. Do not invent Map 5. Do not run confirmation.",
                "E5.1 stays `NOT_SCORED`. Living BSim for a winner is not started.",
                "",
            ]
        )
        status = "ENCODING_FAIL"
        winner_id = None
    else:
        lines.extend(
            [
                f"**Winner: `{winner['id']}`** (`u_floor={winner['u_floor']:.2f}`, `g={winner['g']:.2f}`).",
                f"`n_ALIVE={winner['n_ALIVE']}/12`, contrast `{winner['u_contrast']:.4f}`.",
                "This map is frozen. Do not add a map after this file.",
                "",
            ]
        )
        status = "ENCODING_FROZEN"
        winner_id = winner["id"]
    lines.extend(
        [
            "## Development u SHA-256 (before F408)",
            "",
            "| MapID | File | SHA-256 | u12 | u[0] |",
            "|---|---|---|---|---|",
        ]
    )
    for spec in MAPS:
        for p in DEV:
            rec = hashes["maps"][spec["id"]][f"input_ahl_{spec['id']}_p{p:04d}.txt"]
            lines.append(
                f"| {spec['id']} | `maps/{spec['id']}/input_ahl_p{p:04d}.txt` | "
                f"`{rec['file_sha256']}` | `{rec['u12_sha256']}` | {rec['u0']:.12f} |"
            )
    lines.extend(
        [
            "",
            f"Scout development SHA-256: `{DEV_SHA}`",
            f"Scout confirmation SHA-256: `{CONF_SHA}`",
            f"Vendored mat SHA-256: `{hashes['mat_sha256']}`",
            "",
            "Fallback formula (not used unless the preferred path failed):",
            "",
            "```",
            "AHL_new[n] = AHL_E51[n] * u_new[n] / max(u_E51[n], 1e-6)",
            "R_proxy[n] = AHL_new[n]^2 / (K^2 + AHL_new[n]^2)   # K=1.6",
            "```",
            "",
        ]
    )
    (RESULTS / "ENCODING_FREEZE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    payload = {
        "cohort_status": "SYNTHETIC_COHORT",
        "clinical": False,
        "status": status,
        "winner": winner_id,
        "jmax": JMAX,
        "density_mask": screen_meta["density_mask"],
        "screen_path": screen_meta["path_summary"],
        "transport_status": screen_meta["transport_status"],
        "maps": [
            {k: v for k, v in s.items() if k != "patients"} | {
                "patients": [{kk: vv for kk, vv in r.items() if kk != "R_win"} for r in s["patients"]]
            }
            for s in summaries
        ],
        "hashes": hashes,
    }
    (RESULTS / "e52_encoding_freeze.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    with (RESULTS / "e52_occupancy.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "map_id", "patient", "split", "occupancy", "mean_R", "r_meanR_u",
                "median_R_win", "frac_R_gt_0_5", "mean_R_uniform", "u0", "u_min",
                "u_max", "saturated_map", "n_ALIVE_map", "path",
            ],
        )
        writer.writeheader()
        for s in summaries:
            for r in s["patients"]:
                writer.writerow(
                    {
                        "map_id": r["map_id"],
                        "patient": r["patient"],
                        "split": "development",
                        "occupancy": r["occupancy"],
                        "mean_R": f"{r['mean_R']:.10g}",
                        "r_meanR_u": f"{r['r_meanR_u']:.10g}",
                        "median_R_win": f"{r['median_R_win']:.10g}",
                        "frac_R_gt_0_5": f"{r['frac_R_gt_0_5']:.10g}",
                        "mean_R_uniform": f"{r['mean_R_uniform']:.10g}",
                        "u0": f"{r['u0']:.12f}",
                        "u_min": f"{r['u_min']:.12f}",
                        "u_max": f"{r['u_max']:.12f}",
                        "saturated_map": s["saturated"],
                        "n_ALIVE_map": s["n_ALIVE"],
                        "path": r["path"],
                    }
                )
    return status, winner_id


def write_ceilings(Xnorm):
    persist_all_dev = persist_score(Xnorm, DEV, 0, T - 2)
    persist_all_conf = persist_score(Xnorm, CONF, 0, T - 2)
    persist_al_dev = persist_score(Xnorm, DEV, 9, T - 2)
    persist_al_conf = persist_score(Xnorm, CONF, 9, T - 2)
    Xtr_c, ytr_c, ptr_c = delay_rows(Xnorm, DEV, 9, T - 2, (0,))
    Xte_c, yte_c, pte_c = delay_rows(Xnorm, CONF, 9, T - 2, (0,))
    lam_c, loo_c, grid_c = loo_lambda(Xtr_c, ytr_c, ptr_c, DEV)
    crp10 = fit_score(Xtr_c, ytr_c, ptr_c, Xte_c, yte_c, pte_c, lam_c)
    Xtr_5, ytr_5, ptr_5 = delay_rows(Xnorm, DEV, 9, T - 2, (0, 1, 2, 3, 4))
    Xte_5, yte_5, pte_5 = delay_rows(Xnorm, CONF, 9, T - 2, (0, 1, 2, 3, 4))
    lam_5, loo_5, grid_5 = loo_lambda(Xtr_5, ytr_5, ptr_5, DEV)
    ch5 = fit_score(Xtr_5, ytr_5, ptr_5, Xte_5, yte_5, pte_5, lam_5)
    ceilings = {
        "cohort_status": "SYNTHETIC_COHORT",
        "clinical": False,
        "note": "Target is Xnorm[:,:,0], not u. 5-channel delay-10 is a ceiling, not a living-layer bar.",
        "development": DEV,
        "confirmation": CONF,
        "persistence_t0_46": {
            "development_pooled": persist_all_dev["pooled"],
            "confirmation_pooled": persist_all_conf["pooled"],
            "confirmation_concat_diagnostic": persist_all_conf["concat"],
        },
        "persistence_aligned_t9_46": {
            "development_pooled": persist_al_dev["pooled"],
            "confirmation_pooled": persist_al_conf["pooled"],
            "confirmation_concat_diagnostic": persist_al_conf["concat"],
        },
        "crp_delay10": {**crp10, "loo_pooled_development": loo_c, "loo_grid": grid_c},
        "five_channel_delay10": {**ch5, "loo_pooled_development": loo_5, "loo_grid": grid_5},
    }
    (RESULTS / "e52_input_ceilings.json").write_text(json.dumps(ceilings, indent=2) + "\n", encoding="utf-8")
    return ceilings


def write_configs(winner: dict, u_winner: dict, hashes: dict):
    base = """# E5.2 common template. Claim dish copied from E5.1 / Narma10b. Do not retune.
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
num.windows=48
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
headless=true
arm=driven
rng.seed=111
output.write.samples=true
output.write.voxels=true
"""
    (HERE / "sim_config_e52.properties").write_text(base, encoding="utf-8")
    map_id = winner["id"]

    def driven_cfg(p):
        u0 = float(u_winner[p][0])
        tag = f"p{p:04d}"
        return (
            f"config.include=sim_config_e52.properties\n"
            f"arm=driven\n"
            f"patient.id={p}\n"
            f"encoding.map.id={map_id}\n"
            f"rng.seed=111\n"
            f"warmup.ahl.input={u0:.12f}\n"
            f"input.ahl.file=input/input_ahl_{tag}.txt\n"
            f"input.acid.file=input/input_acid_held05_48.txt\n"
            f"output.run.label=e52_driven_{tag}\n"
            f"output.stochastic.replicate=driven_{tag}_seed111\n"
            f"output.dir=results/e52_driven_{tag}_seed111\n"
        )

    for p in DEV + CONF:
        (HERE / f"sim_config_e52_driven_p{p:04d}_seed111.properties").write_text(
            driven_cfg(p), encoding="utf-8"
        )
    (HERE / "sim_config_e52_silent_seed111.properties").write_text(
        "config.include=sim_config_e52.properties\n"
        "arm=silent\n"
        "patient.id=-1\n"
        "encoding.map.id=SILENT\n"
        "rng.seed=111\n"
        "warmup.ahl.input=0\n"
        "input.ahl.file=input/input_ahl_zero48.txt\n"
        "input.acid.file=input/input_acid_zero48.txt\n"
        "output.run.label=e52_silent\n"
        "output.stochastic.replicate=silent_seed111\n"
        "output.dir=results/e52_silent_seed111\n",
        encoding="utf-8",
    )
    (HERE / "sim_config_e52_brownian_seed111.properties").write_text(
        "config.include=sim_config_e52.properties\n"
        "arm=brownian\n"
        "patient.id=-1\n"
        "encoding.map.id=BROWNIAN\n"
        "rng.seed=111\n"
        "warmup.ahl.input=0\n"
        "input.ahl.file=input/input_ahl_zero48.txt\n"
        "input.acid.file=input/input_acid_held05_48.txt\n"
        "output.run.label=e52_brownian\n"
        "output.stochastic.replicate=brownian_seed111\n"
        "output.dir=results/e52_brownian_seed111\n",
        encoding="utf-8",
    )
    hashes["winner_configs"] = {
        "sim_config_e52.properties": sha256_file(HERE / "sim_config_e52.properties"),
    }
    for p in DEV + CONF:
        name = f"sim_config_e52_driven_p{p:04d}_seed111.properties"
        hashes["winner_configs"][name] = sha256_file(HERE / name)
    print("winner configs written", map_id)


def main():
    INPUT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    MAPS_DIR.mkdir(parents=True, exist_ok=True)
    assert sha256_text(",".join(str(i) for i in DEV)) == DEV_SHA, "development scout hash mismatch"
    assert sha256_text(",".join(str(i) for i in CONF)) == CONF_SHA, "confirmation scout hash mismatch"

    mat = loadmat(str(MAT))
    Xnorm = np.asarray(mat["Xnorm"], dtype=float)
    assert Xnorm.shape[0] >= 26 and Xnorm.shape[1] == T and Xnorm.shape[2] == 5

    transport_status = "UNKNOWN"
    status_line = (PREFLIGHT / "TRANSPORT_MODEL_VALIDATION.md").read_text(encoding="utf-8")
    if "VALIDATED_FOR_SCREENING" in status_line:
        transport_status = "VALIDATED_FOR_SCREENING"

    hashes = {
        "cohort_status": "SYNTHETIC_COHORT",
        "mat_sha256": sha256_file(MAT),
        "development_index_sha256": DEV_SHA,
        "confirmation_index_sha256": CONF_SHA,
        "development": DEV,
        "confirmation": CONF,
        "maps": {},
        "files": {},
    }

    acid = np.full(T, 0.5)
    zeros = np.zeros(T)
    acid_path = INPUT / "input_acid_held05_48.txt"
    zero_ahl = INPUT / "input_ahl_zero48.txt"
    zero_acid = INPUT / "input_acid_zero48.txt"
    write_lines(acid_path, acid)
    write_lines(zero_ahl, zeros)
    write_lines(zero_acid, zeros)
    for path in (acid_path, zero_ahl, zero_acid):
        hashes["files"][path.name] = sha256_file(path)

    for p in DEV + CONF:
        target = Xnorm[p, 1:, 0]
        tgt_path = INPUT / f"target_crp_p{p:04d}.txt"
        write_lines(tgt_path, target)
        hashes["files"][tgt_path.name] = {
            "file_sha256": sha256_file(tgt_path),
            "n": int(len(target)),
            "patient": p,
        }

    u_by_map = {spec["id"]: {} for spec in MAPS}
    u_e51 = {}
    for spec in MAPS:
        hashes["maps"][spec["id"]] = {}
        mdir = MAPS_DIR / spec["id"]
        mdir.mkdir(parents=True, exist_ok=True)
        for p in DEV + CONF:
            u = encode(Xnorm[p, :, 0], spec["u_floor"], spec["g"])
            u_by_map[spec["id"]][p] = u
            if spec["id"] == "M_E51" and p in DEV + CONF:
                u_e51[p] = u
            path = mdir / f"input_ahl_p{p:04d}.txt"
            write_lines(path, u)
            if p in DEV:
                hashes["maps"][spec["id"]][f"input_ahl_{spec['id']}_p{p:04d}.txt"] = {
                    "file_sha256": sha256_file(path),
                    "u12_sha256": sha256_u12(u),
                    "u0": float(u[0]),
                    "patient": p,
                }

    fields = {}
    for p in DEV:
        loaded = load_e51_window_fields(p)
        if loaded is None:
            raise SystemExit(f"E5.1 driven density missing for development patient {p} (read-only reuse)")
        fields[p] = loaded

    rows = run_occupancy_screen(u_by_map, fields, u_e51)
    summaries = summarize_maps(rows, u_by_map)
    winner = select_winner(summaries)
    paths = {r["path"] for r in rows}
    path_summary = "E0.2_TRANSPORT_STAGE3B" if all(p == "E0.2_TRANSPORT_STAGE3B" for p in paths) else ";".join(sorted(paths))
    screen_meta = {
        "density_mask": "FROZEN_E51_DRIVEN_DEN_WINDOWMEAN",
        "path_summary": path_summary,
        "transport_status": transport_status,
    }
    status, winner_id = write_encoding_freeze(summaries, winner, hashes, screen_meta)
    ceilings = write_ceilings(Xnorm)
    (RESULTS / "e52_input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")

    if winner is not None:
        u_winner = u_by_map[winner["id"]]
        for p in DEV + CONF:
            src = MAPS_DIR / winner["id"] / f"input_ahl_p{p:04d}.txt"
            dst = INPUT / f"input_ahl_p{p:04d}.txt"
            dst.write_bytes(src.read_bytes())
            hashes["files"][dst.name] = {
                "file_sha256": sha256_file(dst),
                "u12_sha256": sha256_u12(u_winner[p]),
                "u0": float(u_winner[p][0]),
                "patient": p,
                "map_id": winner["id"],
            }
        write_configs(winner, u_winner, hashes)
        (RESULTS / "e52_input_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
        freeze = json.loads((RESULTS / "e52_encoding_freeze.json").read_text(encoding="utf-8"))
        freeze["hashes"] = hashes
        (RESULTS / "e52_encoding_freeze.json").write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")

    print("DEV_SHA_OK", DEV_SHA)
    print("CONF_SHA_OK", CONF_SHA)
    print("STATUS", status)
    print("WINNER", winner_id)
    print("persist_conf", ceilings["persistence_aligned_t9_46"]["confirmation_pooled"])
    print("crp10_conf", ceilings["crp_delay10"]["confirmation_pooled"])
    print("ch5_conf", ceilings["five_channel_delay10"]["confirmation_pooled"])
    if status == "ENCODING_FAIL":
        raise SystemExit(0)


if __name__ == "__main__":
    main()
