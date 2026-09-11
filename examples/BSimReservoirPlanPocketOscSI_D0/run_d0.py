#!/usr/bin/env python3
"""D0 evidence suite: SI transcription, traces, period vs mu, checks.

Frozen before traces. Does not retune. Does not see NARMA.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from danino_si_dde import (  # noqa: E402
    PRODUCTION_POCKETOSCSI_SPLIT,
    PRODUCTION_SI_HILL,
    DaninoSIDDE,
    load_json,
    production,
)
from period_check import extract_period, refuse_narma  # noqa: E402

RESULTS = HERE / "results"
FIGURES = RESULTS / "figures"
FIXTURES = RESULTS / "fixtures"


def _arm_state(arm: dict) -> tuple[np.ndarray, float]:
    y0 = np.array([arm["A"], arm["I"], arm["H_i"], arm["H_e"]], dtype=float)
    return y0, float(arm["H_i_history"])


def _means(t: np.ndarray, y: np.ndarray, t_discard: float) -> dict:
    z = y[t >= t_discard - 1e-12]
    return {
        "mean_A": float(np.mean(z[:, 0])),
        "mean_I": float(np.mean(z[:, 1])),
        "mean_Hi": float(np.mean(z[:, 2])),
        "mean_He": float(np.mean(z[:, 3])),
        "min_state": float(np.min(y)),
    }


def equation_diff(model: DaninoSIDDE) -> dict:
    hs = [0.0, 0.01, 0.1, 1.0, 3.0, 10.0]
    rows = []
    for h in hs:
        si = production(h, model.delta, model.alpha, model.k1, PRODUCTION_SI_HILL)
        split = production(h, model.delta, model.alpha, model.k1, PRODUCTION_POCKETOSCSI_SPLIT)
        rows.append(
            {
                "H_tau": h,
                "P_si_hill": si,
                "P_pocketoscsi_split": split,
                "abs_diff": abs(si - split),
            }
        )
    return {
        "si_quote": (
            "P(alpha, tau) = (delta + alpha H_tau^2) / (1 + k1 H_tau^2) "
            "as typeset on SI p.5 with a single fraction bar."
        ),
        "pocketoscsi": "P = delta + alpha H_tau^2 / (1 + k1 H_tau^2)",
        "rows": rows,
    }


def print_transcription_table(diff: dict) -> None:
    print("=== EQUATION TRANSCRIPTION ===", flush=True)
    print("SI P:        (delta + alpha H_tau^2) / (1 + k1 H_tau^2)", flush=True)
    print("D0 oracle:   production_form=si_hill", flush=True)
    print("PocketOsc-SI: delta + alpha H_tau^2 / (1 + k1 H_tau^2)", flush=True)
    print(f"{'H_tau':>8} {'P_SI':>14} {'P_prior':>14} {'|diff|':>14}", flush=True)
    for row in diff["rows"]:
        print(
            f"{row['H_tau']:8.3f} {row['P_si_hill']:14.6g} "
            f"{row['P_pocketoscsi_split']:14.6g} {row['abs_diff']:14.6g}",
            flush=True,
        )


def save_fixture(name: str, traj: dict, peaks: dict) -> None:
    FIXTURES.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        FIXTURES / f"{name}.npz",
        t=traj["t"],
        Y=traj["Y"],
        mu=traj["mu"],
        d=traj["d"],
        y0=traj["y0"],
        H_i_history=traj["H_i_history"],
        peak_times=np.array(peaks.get("peak_times", []), dtype=float),
    )
    with (FIXTURES / f"{name}.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "mu": traj["mu"],
                "d": traj["d"],
                "production_form": traj["production_form"],
                "flag": peaks["flag"],
                "period": peaks["period"],
                "n_peaks": peaks["n_peaks"],
                "y0": traj["y0"].tolist(),
                "H_i_history": traj["H_i_history"],
            },
            handle,
            indent=2,
        )


def plot_trace(path: Path, traj: dict, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    t = traj["t"]
    y = traj["Y"]
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.0), sharex=True)
    axes[0].plot(t, y[:, 1], label="I (LuxI)", color="tab:cyan")
    axes[0].plot(t, y[:, 0], label="A (AiiA)", color="tab:blue", alpha=0.8)
    axes[0].set_ylabel("protein")
    axes[0].legend(loc="upper right")
    axes[1].plot(t, y[:, 2], label="Hi", color="tab:green")
    axes[1].plot(t, y[:, 3], label="He", color="tab:red")
    axes[1].set_ylabel("AHL")
    axes[1].set_xlabel("SI-scaled time")
    axes[1].legend(loc="upper right")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def integrate_arm(
    model: DaninoSIDDE,
    ens: dict,
    mu: float,
    protocol: dict,
) -> tuple[dict, dict, dict]:
    y0, hist = _arm_state(ens)
    traj = model.integrate(
        d=0.5,
        mu=mu,
        y0=y0,
        hi_history=hist,
        t_end=float(protocol["t_end"]),
        sample_dt=float(protocol["sample_dt"]),
    )
    peaks = extract_period(traj["t"], traj["Y"][:, 1], protocol)
    if traj["negative_state"]:
        peaks = dict(peaks)
        peaks["flag"] = "NEGATIVE_STATE"
        peaks["period"] = float("nan")
    stats = _means(traj["t"], traj["Y"], float(protocol["t_discard"]))
    return traj, peaks, stats


def spearman(xs: list[float], ys: list[float]) -> float:
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def main() -> None:
    refuse_narma()
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true", help="primary mu=1.5 only")
    args = parser.parse_args()

    params = load_json("params.json")
    protocol = load_json("protocol.json")
    grid = load_json("mu_grid.json")
    ensemble = load_json("ensemble.json")
    ens_by_id = {row["id"]: row for row in ensemble["arms"]}
    primary = ens_by_id[ensemble["primary_id"]]
    prior = ens_by_id["FAILED_PRIOR_IVP"]

    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    FIXTURES.mkdir(parents=True, exist_ok=True)

    model = DaninoSIDDE(params=params, production_form=PRODUCTION_SI_HILL)
    diff = equation_diff(model)
    print_transcription_table(diff)
    (RESULTS / "equation_diff.json").write_text(json.dumps(diff, indent=2), encoding="utf-8")

    print("=== ENSEMBLE (frozen) ===", flush=True)
    print(f"label={ensemble['label']} primary={ensemble['primary_id']}", flush=True)
    for row in ensemble["arms"]:
        print(
            f"  {row['id']:18} A={row['A']} I={row['I']} Hi={row['H_i']} "
            f"He={row['H_e']} hist={row['H_i_history']} role={row['role']}",
            flush=True,
        )

    mu_arms = grid["arms"]
    if args.smoke:
        mu_arms = [row for row in mu_arms if row["id"] == "MU_1p5"]

    period_rows = []
    osc_fixture_saved = False
    for row in mu_arms:
        print(f"PRIMARY {row['id']} mu={row['mu']} ...", flush=True)
        traj, peaks, stats = integrate_arm(model, primary, float(row["mu"]), protocol)
        rec = {
            "arm": row["id"],
            "role": row["role"],
            "ensemble": primary["id"],
            "mu": row["mu"],
            "d": 0.5,
            "production_form": PRODUCTION_SI_HILL,
            "flag": peaks["flag"],
            "period": peaks["period"],
            "n_peaks": peaks["n_peaks"],
            **stats,
            "rel_amplitude": peaks["rel_amplitude"],
            "amp_persist": peaks["amp_persist"],
        }
        period_rows.append(rec)
        print(
            f"  flag={peaks['flag']} T={peaks['period']} n_peaks={peaks['n_peaks']} "
            f"mean_I={stats['mean_I']:.4g} mean_Hi={stats['mean_Hi']:.4g}",
            flush=True,
        )
        if row["id"] == "MU_1p5":
            save_fixture("primary_mu_1p5", traj, peaks)
            plot_trace(
                FIGURES / "primary_mu_1p5.png",
                traj,
                f"PRIMARY {primary['id']} mu=1.5 si_hill",
            )
            osc_fixture_saved = True
        elif peaks["flag"] == "OSC" and not osc_fixture_saved:
            save_fixture("primary_osc", traj, peaks)
            plot_trace(
                FIGURES / "primary_osc.png",
                traj,
                f"PRIMARY {primary['id']} mu={row['mu']} si_hill",
            )
            osc_fixture_saved = True

    print("=== FAILED_PRIOR_IVP replay ===", flush=True)
    prior_mu = 1.5
    prior_records = []
    for form in (PRODUCTION_SI_HILL, PRODUCTION_POCKETOSCSI_SPLIT):
        replay = DaninoSIDDE(params=params, production_form=form)
        traj, peaks, stats = integrate_arm(replay, prior, prior_mu, protocol)
        rec = {
            "arm": "FAILED_PRIOR_IVP",
            "ensemble": "FAILED_PRIOR_IVP",
            "mu": prior_mu,
            "d": 0.5,
            "production_form": form,
            "flag": peaks["flag"],
            "period": peaks["period"],
            "n_peaks": peaks["n_peaks"],
            **stats,
        }
        prior_records.append(rec)
        print(
            f"  form={form} flag={peaks['flag']} T={peaks['period']} "
            f"n_peaks={peaks['n_peaks']} mean_I={stats['mean_I']:.4g} "
            f"mean_Hi={stats['mean_Hi']:.4g}",
            flush=True,
        )
        if form == PRODUCTION_SI_HILL:
            save_fixture("failed_prior_mu_1p5", traj, peaks)
            plot_trace(
                FIGURES / "failed_prior_mu_1p5.png",
                traj,
                "FAILED_PRIOR_IVP mu=1.5 si_hill",
            )

    print("=== HISTORY SENSITIVITY at mu=1.5 ===", flush=True)
    hist_rows = []
    if not args.smoke:
        for ens in ensemble["arms"]:
            traj, peaks, stats = integrate_arm(model, ens, 1.5, protocol)
            rec = {
                "ensemble": ens["id"],
                "role": ens["role"],
                "mu": 1.5,
                "flag": peaks["flag"],
                "period": peaks["period"],
                "n_peaks": peaks["n_peaks"],
                **stats,
            }
            hist_rows.append(rec)
            print(
                f"  {ens['id']:18} flag={peaks['flag']} T={peaks['period']} "
                f"n_peaks={peaks['n_peaks']}",
                flush=True,
            )
            if ens["id"] != "SI_BASAL_PERTURB":
                plot_trace(
                    FIGURES / f"ensemble_{ens['id'].lower()}_mu_1p5.png",
                    traj,
                    f"{ens['id']} mu=1.5 si_hill",
                )

    print("=== TOLERANCE CHECK primary mu=1.5 ===", flush=True)
    tight = DaninoSIDDE(
        params=params,
        production_form=PRODUCTION_SI_HILL,
        rtol=float(params["integrator"]["tight_rtol"]),
        atol=float(params["integrator"]["tight_atol"]),
    )
    traj_t, peaks_t, _ = integrate_arm(tight, primary, 1.5, protocol)
    base = next(r for r in period_rows if r["arm"] == "MU_1p5")
    t_base = float(base["period"])
    t_tight = float(peaks_t["period"])
    if math_nan(t_base) or math_nan(t_tight):
        rel = float("nan")
        period_ok = base["flag"] == peaks_t["flag"]
    else:
        rel = abs(t_tight - t_base) / t_base
        period_ok = rel <= float(params["integrator"]["period_rel_tol"])
    flag_ok = base["flag"] == peaks_t["flag"]
    tol = {
        "nominal_flag": base["flag"],
        "tight_flag": peaks_t["flag"],
        "nominal_period": t_base,
        "tight_period": t_tight,
        "rel_period_change": rel,
        "flag_match": flag_ok,
        "period_ok": period_ok,
        "pass": bool(flag_ok and period_ok),
    }
    print(json.dumps(tol, indent=2), flush=True)

    identity = [r for r in period_rows if r["role"] == "IDENTITY"]
    osc_id = [r for r in identity if r["flag"] == "OSC"]
    trend = {
        "applicable": bool(identity),
        "n_identity": len(identity),
        "n_identity_osc": len(osc_id),
        "all_identity_osc": False,
        "consecutive_increase": False,
        "spearman": None,
        "T_mu_1": None,
        "T_mu_2": None,
        "T2_gt_T1": False,
        "period_class_ok": False,
        "pass": False,
    }
    if osc_id:
        mus = [r["mu"] for r in osc_id]
        ts = [r["period"] for r in osc_id]
        consecutive = all(ts[i + 1] > ts[i] for i in range(len(ts) - 1)) if len(ts) >= 2 else False
        rho = spearman(mus, ts) if len(ts) >= 2 else float("nan")
        t1 = next((r["period"] for r in osc_id if r["mu"] == 1.0), float("nan"))
        t2 = next((r["period"] for r in osc_id if r["mu"] == 2.0), float("nan"))
        class_ok = all(
            protocol["period_class_min"] <= r["period"] <= protocol["period_class_max"]
            for r in osc_id
        )
        all_identity_osc = len(osc_id) == len(identity) and len(identity) > 0
        trend = {
            "applicable": True,
            "n_identity": len(identity),
            "n_identity_osc": len(osc_id),
            "all_identity_osc": all_identity_osc,
            "consecutive_increase": consecutive,
            "spearman": rho,
            "T_mu_1": t1,
            "T_mu_2": t2,
            "T2_gt_T1": (t2 > t1) if not (math_nan(t1) or math_nan(t2)) else False,
            "period_class_ok": class_ok,
            "pass": bool(
                all_identity_osc
                and consecutive
                and (rho >= protocol["spearman_min"] if len(ts) >= 2 else False)
                and (t2 > t1 if not (math_nan(t1) or math_nan(t2)) else False)
                and class_ok
            ),
        }
    print("=== TREND ===", flush=True)
    print(json.dumps(trend, indent=2), flush=True)

    extras = [r for r in period_rows if r["role"] == "FINITE_INTERVAL_EXTRA"]
    finite = {
        "n_extra": len(extras),
        "extra_flags": {r["arm"]: r["flag"] for r in extras},
        "note": (
            "Fig.4c finite band is primarily in cell density d. "
            "Extras are reported, not a hard PASS requirement."
        ),
    }

    write_csv(RESULTS / "period_vs_mu.csv", period_rows)
    write_csv(RESULTS / "failed_prior_replay.csv", prior_records)
    write_csv(RESULTS / "history_sensitivity.csv", hist_rows)

    summary = {
        "primary": ensemble["primary_id"],
        "ics_label": ensemble["label"],
        "smoke": bool(args.smoke),
        "trend": trend,
        "tolerance": tol,
        "finite_interval": finite,
        "period_rows": period_rows,
        "failed_prior": prior_records,
        "history_sensitivity": hist_rows,
    }
    if not args.smoke:
        d0_pass = bool(trend.get("pass") and tol.get("pass"))
        summary["D0"] = "PASS" if d0_pass else "FAIL"
    else:
        summary["D0"] = "SMOKE"
    (RESULTS / "d0_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2),
        encoding="utf-8",
    )
    print(f"D0={summary['D0']}", flush=True)


def math_nan(x: float) -> bool:
    return x != x


def json_safe(obj):
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, float) and obj != obj:
        return None
    if isinstance(obj, np.ndarray):
        return json_safe(obj.tolist())
    if isinstance(obj, (np.floating,)):
        value = float(obj)
        return None if value != value else value
    return obj


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key, value in out.items():
                if isinstance(value, float) and value != value:
                    out[key] = ""
            writer.writerow(out)


if __name__ == "__main__":
    main()
