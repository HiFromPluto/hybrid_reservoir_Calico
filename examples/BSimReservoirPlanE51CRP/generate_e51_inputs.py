#!/usr/bin/env python3
"""E5.1 Module 0: hashed CRP inputs and CRP-only ceilings. No Java."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MAT = REPO / "examples" / "HybridDish" / "external_analysis" / "synthetic_biomarker_BSim_inputs.mat"
INPUT = HERE / "input"
RESULTS = HERE / "results"

DEV = [1, 4, 6, 7, 9, 13, 16, 18, 20, 22, 24, 25]
CONF = [8, 12, 14, 15, 17, 19, 21, 23]
DEV_SHA = "83acc7102e532112c961bb45a1327504290ffddea1c06e659bfb5e68be032915"
CONF_SHA = "c90cca8f5627c737b38024b4863415e46a384d8715e6d06a5823e8cb93ba3787"
LAMBDA_GRID = np.array([1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6])
T = 48


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


def concat_nrmse(ytrue, ypred):
    return nrmse_pop(ytrue, ypred)


def standardize_fit(X):
    mu = X.mean(axis=0)
    sg = X.std(axis=0, ddof=0)
    sg = np.where(sg < 1e-12, 1.0, sg)
    return (X - mu) / sg, mu, sg


def standardize_apply(X, mu, sg):
    return (X - mu) / sg


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
    n_ch = len(channels)
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
            pred = ridge_predict(standardize_apply(X[va], mu, sg), w)
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
    pred_te = ridge_predict(standardize_apply(Xte, mu, sg), w)
    tr_p, _ = patient_pooled(ytr, pred_tr, ptr)
    te_p, te_per = patient_pooled(yte, pred_te, pte)
    return {
        "lambda": float(lam),
        "development_pooled": tr_p,
        "confirmation_pooled": te_p,
        "confirmation_concat": concat_nrmse(yte, pred_te),
        "confirmation_per_patient": te_per,
    }


def persist_score(Xnorm, patients, t_lo, t_hi):
    x, y, pids = persist_rows(Xnorm, patients, t_lo, t_hi)
    pooled, per = patient_pooled(y, x, pids)
    return {
        "pooled": pooled,
        "concat": concat_nrmse(y, x),
        "per_patient": per,
        "n_pairs": int(len(y)),
    }


def main():
    INPUT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    assert ",".join(str(i) for i in DEV)
    assert sha256_text(",".join(str(i) for i in DEV)) == DEV_SHA, "development scout hash mismatch"
    assert sha256_text(",".join(str(i) for i in CONF)) == CONF_SHA, "confirmation scout hash mismatch"

    mat = loadmat(str(MAT))
    Xnorm = np.asarray(mat["Xnorm"], dtype=float)
    assert Xnorm.shape[0] >= 26 and Xnorm.shape[1] == T and Xnorm.shape[2] == 5

    hashes = {
        "cohort_status": "SYNTHETIC_COHORT",
        "mat_sha256": sha256_file(MAT),
        "development_index_sha256": DEV_SHA,
        "confirmation_index_sha256": CONF_SHA,
        "development": DEV,
        "confirmation": CONF,
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
        u = 0.5 * Xnorm[p, :, 0]
        target = Xnorm[p, 1:, 0]
        ahl_path = INPUT / f"input_ahl_p{p:04d}.txt"
        tgt_path = INPUT / f"target_crp_p{p:04d}.txt"
        write_lines(ahl_path, u)
        write_lines(tgt_path, target)
        hashes["files"][ahl_path.name] = {
            "file_sha256": sha256_file(ahl_path),
            "u12_sha256": sha256_u12(u),
            "u0": float(u[0]),
            "patient": p,
        }
        hashes["files"][tgt_path.name] = {
            "file_sha256": sha256_file(tgt_path),
            "n": int(len(target)),
            "patient": p,
        }

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
        "note": "5-channel delay-10 is an analysis ceiling, not a living-layer bar.",
    }
    (RESULTS / "e51_input_hashes.json").write_text(
        json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
    )
    (RESULTS / "e51_input_ceilings.json").write_text(
        json.dumps(ceilings, indent=2) + "\n", encoding="utf-8"
    )

    md = []
    md.append("# E5.1 Module 0 — hashed inputs and CRP-only ceilings")
    md.append("")
    md.append("No Java. No BSim. `COHORT_STATUS: SYNTHETIC_COHORT`.")
    md.append("Not a clinical dataset. Task A / `Y` were not used.")
    md.append("Scout hashes verified against E5.0. Do not redraw.")
    md.append("")
    md.append("| Scout | Indices | SHA-256 |")
    md.append("|---|---|---|")
    md.append(f"| development | {','.join(str(i) for i in DEV)} | `{DEV_SHA}` |")
    md.append(f"| confirmation | {','.join(str(i) for i in CONF)} | `{CONF_SHA}` |")
    md.append("")
    md.append(f"Vendored mat SHA-256: `{hashes['mat_sha256']}`")
    md.append("")
    md.append("Map: `u[n] = 0.5 * Xnorm[p, n, 0]`. Target: `Xnorm[p, t+1, 0]`.")
    md.append("Warmup command is `u[0]` for that patient.")
    md.append("")
    md.append("## File hashes")
    md.append("")
    md.append("| File | SHA-256 |")
    md.append("|---|---|")
    for name, rec in hashes["files"].items():
        digest = rec if isinstance(rec, str) else rec["file_sha256"]
        md.append(f"| `{name}` | `{digest}` |")
    md.append("")
    md.append("## Patient-pooled ceilings on the 12+8 (before Java)")
    md.append("")
    md.append("Lambda by leave-one-patient-out on the 12, refit 12, score 8 once.")
    md.append("Concat NRMSE is a diagnostic only. No new 0.35 void.")
    md.append("")
    md.append("| Baseline | Dev pooled | Conf pooled | Conf concat (diagnostic) | λ |")
    md.append("|---|---|---|---|---|")
    md.append(
        f"| persistence `t=0..46` | {persist_all_dev['pooled']:.4f} | "
        f"{persist_all_conf['pooled']:.4f} | {persist_all_conf['concat']:.4f} | — |"
    )
    md.append(
        f"| persistence aligned `t=9..46` | {persist_al_dev['pooled']:.4f} | "
        f"{persist_al_conf['pooled']:.4f} | {persist_al_conf['concat']:.4f} | — |"
    )
    md.append(
        f"| CRP-only delay-10 (10-D) | {crp10['development_pooled']:.4f} | "
        f"{crp10['confirmation_pooled']:.4f} | {crp10['confirmation_concat']:.4f} | {lam_c:g} |"
    )
    md.append(
        f"| 5-channel delay-10 (50-D, ceiling) | {ch5['development_pooled']:.4f} | "
        f"{ch5['confirmation_pooled']:.4f} | {ch5['confirmation_concat']:.4f} | {lam_5:g} |"
    )
    md.append("")
    md.append("5-channel delay-10 sees IL6/ferritin/lactate/glucose. The dish")
    md.append("does not. Fair input ceilings are persistence and CRP-only")
    md.append("delay-10. E5.0 full-cohort 5-ch 0.805 is a different patient")
    md.append("set; this table is the 12+8 scout.")
    md.append("")
    (RESULTS / "e51_input_ceilings.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("DEV_SHA_OK", DEV_SHA)
    print("CONF_SHA_OK", CONF_SHA)
    print("persist_conf", persist_all_conf["pooled"])
    print("persist_aligned_conf", persist_al_conf["pooled"])
    print("crp10_conf", crp10["confirmation_pooled"], "lam", lam_c)
    print("ch5_conf", ch5["confirmation_pooled"], "lam", lam_5)
    print("files", len(hashes["files"]))
    write_configs(hashes)


def write_configs(hashes):
    base = """# E5.1 common template. Claim dish copied from Narma10b. Do not retune.
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
    (HERE / "sim_config_e51.properties").write_text(base, encoding="utf-8")

    def driven_cfg(p):
        rec = hashes["files"][f"input_ahl_p{p:04d}.txt"]
        u0 = rec["u0"]
        tag = f"p{p:04d}"
        return (
            f"config.include=sim_config_e51.properties\n"
            f"arm=driven\n"
            f"patient.id={p}\n"
            f"rng.seed=111\n"
            f"warmup.ahl.input={u0:.12f}\n"
            f"input.ahl.file=input/input_ahl_{tag}.txt\n"
            f"input.acid.file=input/input_acid_held05_48.txt\n"
            f"output.run.label=e51_driven_{tag}\n"
            f"output.stochastic.replicate=driven_{tag}_seed111\n"
            f"output.dir=results/e51_driven_{tag}_seed111\n"
        )

    for p in DEV + CONF:
        (HERE / f"sim_config_e51_driven_p{p:04d}_seed111.properties").write_text(
            driven_cfg(p), encoding="utf-8"
        )
    (HERE / "sim_config_e51_silent_seed111.properties").write_text(
        "config.include=sim_config_e51.properties\n"
        "arm=silent\n"
        "patient.id=-1\n"
        "rng.seed=111\n"
        "warmup.ahl.input=0\n"
        "input.ahl.file=input/input_ahl_zero48.txt\n"
        "input.acid.file=input/input_acid_zero48.txt\n"
        "output.run.label=e51_silent\n"
        "output.stochastic.replicate=silent_seed111\n"
        "output.dir=results/e51_silent_seed111\n",
        encoding="utf-8",
    )
    (HERE / "sim_config_e51_brownian_seed111.properties").write_text(
        "config.include=sim_config_e51.properties\n"
        "arm=brownian\n"
        "patient.id=-1\n"
        "rng.seed=111\n"
        "warmup.ahl.input=0\n"
        "input.ahl.file=input/input_ahl_zero48.txt\n"
        "input.acid.file=input/input_acid_held05_48.txt\n"
        "output.run.label=e51_brownian\n"
        "output.stochastic.replicate=brownian_seed111\n"
        "output.dir=results/e51_brownian_seed111\n",
        encoding="utf-8",
    )
    print("configs written")


if __name__ == "__main__":
    main()
