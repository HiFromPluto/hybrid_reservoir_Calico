#!/usr/bin/env python3
"""D0b evidence suite: AHL-kick traces, period vs mu, identity gates.

PROTOCOL.md must exist with frozen_before_traces=true before this
identity-grid run. Does not retune. Does not hunt ICs. Does not see NARMA.
Does not import D0 diagnosis scripts.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
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


def require_frozen_protocol() -> None:
    proto_md = HERE / "PROTOCOL.md"
    if not proto_md.exists():
        raise SystemExit("PROTOCOL.md must exist before traces.")
    text = proto_md.read_text(encoding="utf-8")
    if "frozen_before_traces:** true" not in text and "frozen_before_traces: true" not in text:
        raise SystemExit("PROTOCOL.md must declare frozen_before_traces true before traces.")
    protocol = load_json("protocol.json")
    if protocol.get("frozen_before_traces") is not True:
        raise SystemExit("configs/protocol.json must have frozen_before_traces=true before traces.")
    ensemble = load_json("ensemble.json")
    if ensemble.get("frozen_before_traces") is not True:
        raise SystemExit("configs/ensemble.json must have frozen_before_traces=true before traces.")
    grid = load_json("mu_grid.json")
    if grid.get("frozen_before_traces") is not True:
        raise SystemExit("configs/mu_grid.json must have frozen_before_traces=true before traces.")


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
        "I_end": float(y[-1, 1]),
        "Hi_end": float(y[-1, 2]),
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
    print("D0b oracle:  production_form=si_hill (copied from D0)", flush=True)
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
    period = peaks["period"]
    with (FIXTURES / f"{name}.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "mu": traj["mu"],
                "d": traj["d"],
                "production_form": traj["production_form"],
                "flag": peaks["flag"],
                "period": None if math_nan(float(period)) else float(period),
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


def plot_period_vs_mu(path: Path, rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    ident = [r for r in rows if r["role"] == "IDENTITY"]
    extra = [r for r in rows if r["role"] == "FINITE_INTERVAL_EXTRA"]
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    for group, marker, label in (
        (ident, "o", "identity"),
        (extra, "x", "extra"),
    ):
        osc = [r for r in group if r["flag"] == "OSC"]
        nop = [r for r in group if r["flag"] != "OSC"]
        if osc:
            ax.plot(
                [r["mu"] for r in osc],
                [r["period"] for r in osc],
                marker=marker,
                linestyle="-" if group is ident else "None",
                label=f"{label} OSC",
            )
        if nop:
            ax.scatter(
                [r["mu"] for r in nop],
                [0.0 for _ in nop],
                marker=marker,
                facecolors="none",
                edgecolors="0.5",
                label=f"{label} NO_PERIOD (plotted at 0)",
            )
    ax.set_xlabel(r"$\mu$ (SI-scaled)")
    ax.set_ylabel("period (SI-scaled)")
    ax.set_title("D0b primary AHL_KICK_005 period vs mu")
    ax.legend(loc="best")
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


def record_row(
    ens: dict,
    grid_row: dict | None,
    mu: float,
    role: str,
    peaks: dict,
    stats: dict,
    production_form: str = PRODUCTION_SI_HILL,
) -> dict:
    return {
        "arm": None if grid_row is None else grid_row["id"],
        "role": role,
        "ensemble": ens["id"],
        "mu": float(mu),
        "d": 0.5,
        "production_form": production_form,
        "flag": peaks["flag"],
        "period": peaks["period"],
        "n_peaks": peaks["n_peaks"],
        **stats,
        "rel_amplitude": peaks["rel_amplitude"],
        "amp_persist": peaks["amp_persist"],
    }


def spearman(xs: list[float], ys: list[float]) -> float:
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def math_nan(x: float) -> bool:
    return x != x


def fmt_period(value: float) -> str:
    if value is None or math_nan(float(value)):
        return "—"
    return f"{float(value):.3f}"


def write_period_table(path: Path, rows: list[dict], primary_id: str) -> None:
    lines = [
        f"# D0b period vs μ (primary arm `{primary_id}`, d=0.5, D1=0, si_hill)",
        "",
        "Frozen grid. Do not add points after occupancy.",
        "",
        "| Arm | role | μ | flag | period | n_peaks | mean I | mean Hi |",
        "|---|---|---:|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        mu = r["mu"]
        lines.append(
            f"| {r['arm']} | {r['role']} | {mu:.2f} | {r['flag']} | "
            f"{fmt_period(r['period'])} | {r['n_peaks']} | "
            f"{r['mean_I']:.4g} | {r['mean_Hi']:.4g} |"
        )
    lines.append("")
    lines.append("Source CSV: `period_vs_mu.csv`.")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def score_d0b(
    protocol: dict,
    period_rows: list[dict],
    control_rows: list[dict],
    tol: dict,
) -> dict:
    identity = [r for r in period_rows if r["role"] == "IDENTITY"]
    osc_id = [r for r in identity if r["flag"] == "OSC"]
    extras = [r for r in period_rows if r["role"] == "FINITE_INTERVAL_EXTRA"]
    extras_ge = [r for r in extras if r["mu"] >= float(protocol["extras_off_mu_min"]) - 1e-12]
    extras_off_ok = all(r["flag"] in ("NO_PERIOD", "NEGATIVE_STATE") for r in extras_ge) and bool(extras_ge)

    occupancy_ok = len(osc_id) >= int(protocol["identity_min_osc"])
    class_ok = bool(osc_id) and all(
        protocol["period_class_min"] <= r["period"] <= protocol["period_class_max"]
        for r in osc_id
    )

    consecutive = False
    rho = float("nan")
    trend_ok = False
    mus: list[float] = []
    ts: list[float] = []
    if len(osc_id) >= 2:
        osc_sorted = sorted(osc_id, key=lambda r: r["mu"])
        mus = [r["mu"] for r in osc_sorted]
        ts = [r["period"] for r in osc_sorted]
        consecutive = all(ts[i + 1] > ts[i] for i in range(len(ts) - 1))
        rho = spearman(mus, ts)
        trend_ok = bool(consecutive and rho >= float(protocol["spearman_min"]))

    wrong = next(
        (
            r
            for r in control_rows
            if r["ensemble"] == protocol["wrong_kick_id"]
            and abs(r["mu"] - float(protocol["wrong_kick_mu"])) < 1e-12
        ),
        None,
    )
    wrong_ok = bool(wrong) and wrong["flag"] == protocol["wrong_kick_required_flag"]
    tol_ok = bool(tol.get("pass"))

    gates = {
        "occupancy_ok": occupancy_ok,
        "n_identity": len(identity),
        "n_identity_osc": len(osc_id),
        "identity_min_osc": int(protocol["identity_min_osc"]),
        "period_class_ok": class_ok,
        "consecutive_increase": consecutive,
        "spearman": None if math_nan(rho) else rho,
        "trend_ok": trend_ok,
        "osc_identity_mu": mus,
        "osc_identity_period": ts,
        "extras_mu_ge_1p2_off_ok": extras_off_ok,
        "extras_mu_ge_1p2_flags": {r["arm"]: r["flag"] for r in extras_ge},
        "wrong_kick_ok": wrong_ok,
        "wrong_kick_flag": None if wrong is None else wrong["flag"],
        "tolerance_ok": tol_ok,
    }

    if not occupancy_ok:
        status = "FAIL"
        reason = "primary arm did not occupy the identity line (<6/8 OSC)"
    elif not trend_ok:
        status = "FAIL_NO_IDENTITY"
        reason = "primary occupied the identity line but period did not increase with mu"
    elif not (class_ok and extras_off_ok and wrong_ok and tol_ok):
        status = "FAIL"
        reason = "occupancy and trend held but another predeclared gate failed"
    else:
        status = "PASS"
        reason = "all six identity gates held on the primary arm"

    gates["D0b"] = status
    gates["reason"] = reason
    gates["d1_may_start"] = status == "PASS"
    return gates


def evaluate_tolerance(
    params: dict,
    protocol: dict,
    ensemble: dict,
    ens_by_id: dict,
    period_rows: list[dict],
    control_rows: list[dict],
) -> dict:
    tight_model = DaninoSIDDE(
        params=params,
        production_form=PRODUCTION_SI_HILL,
        rtol=float(params["integrator"]["tight_rtol"]),
        atol=float(params["integrator"]["tight_atol"]),
    )
    rel_tol = float(params["integrator"]["period_rel_tol"])
    items = []
    all_ok = True
    for spec in protocol["tolerance_subset"]:
        ens = ens_by_id[spec["ensemble"]]
        mu = float(spec["mu"])
        print(f"  tight {ens['id']} mu={mu} ...", flush=True)
        _, peaks_t, _ = integrate_arm(tight_model, ens, mu, protocol)
        if ens["id"] == ensemble["primary_id"]:
            base = next(r for r in period_rows if abs(r["mu"] - mu) < 1e-12)
        else:
            base = next(
                r
                for r in control_rows
                if r["ensemble"] == ens["id"] and abs(r["mu"] - mu) < 1e-12
            )
        t_base = float(base["period"])
        t_tight = float(peaks_t["period"])
        flag_ok = base["flag"] == peaks_t["flag"]
        if math_nan(t_base) or math_nan(t_tight):
            rel = float("nan")
            period_ok = flag_ok
        else:
            rel = abs(t_tight - t_base) / t_base
            period_ok = rel <= rel_tol
        item_pass = bool(flag_ok and period_ok)
        all_ok = all_ok and item_pass
        items.append(
            {
                "ensemble": ens["id"],
                "mu": mu,
                "nominal_flag": base["flag"],
                "tight_flag": peaks_t["flag"],
                "nominal_period": None if math_nan(t_base) else t_base,
                "tight_period": None if math_nan(t_tight) else t_tight,
                "rel_period_change": None if math_nan(rel) else rel,
                "flag_match": flag_ok,
                "period_ok": period_ok,
                "pass": item_pass,
            }
        )
    return {"items": items, "pass": all_ok}


def main() -> None:
    refuse_narma()
    require_frozen_protocol()
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true", help="primary mu=0.40 only")
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
        mu_arms = [row for row in mu_arms if abs(row["mu"] - 0.4) < 1e-12]

    period_rows = []
    osc_fixture_saved = False
    primary_040_traj = None
    primary_040_peaks = None
    for row in mu_arms:
        print(f"PRIMARY {row['id']} mu={row['mu']} ...", flush=True)
        traj, peaks, stats = integrate_arm(model, primary, float(row["mu"]), protocol)
        rec = record_row(primary, row, float(row["mu"]), row["role"], peaks, stats)
        period_rows.append(rec)
        print(
            f"  flag={peaks['flag']} T={peaks['period']} n_peaks={peaks['n_peaks']} "
            f"mean_I={stats['mean_I']:.4g} mean_Hi={stats['mean_Hi']:.4g}",
            flush=True,
        )
        if abs(row["mu"] - 0.4) < 1e-12:
            primary_040_traj, primary_040_peaks = traj, peaks
            plot_trace(
                FIGURES / "primary_mu_0p40.png",
                traj,
                f"PRIMARY {primary['id']} mu=0.40 si_hill",
            )
        if peaks["flag"] == "OSC" and not osc_fixture_saved:
            save_fixture("primary_osc", traj, peaks)
            plot_trace(
                FIGURES / "primary_osc.png",
                traj,
                f"PRIMARY {primary['id']} mu={row['mu']} si_hill",
            )
            osc_fixture_saved = True

    if primary_040_traj is not None:
        save_fixture("primary_mu_0p40", primary_040_traj, primary_040_peaks)

    print("=== CONTROLS at mu=0.40 ===", flush=True)
    control_rows = []
    if not args.smoke:
        for ens in ensemble["arms"]:
            if ens["id"] == primary["id"]:
                continue
            traj, peaks, stats = integrate_arm(model, ens, 0.4, protocol)
            rec = record_row(ens, None, 0.4, ens["role"], peaks, stats)
            rec["arm"] = f"{ens['id']}_MU_0p40"
            control_rows.append(rec)
            print(
                f"  {ens['id']:18} flag={peaks['flag']} T={peaks['period']} "
                f"n_peaks={peaks['n_peaks']} mean_I={stats['mean_I']:.4g}",
                flush=True,
            )
            plot_trace(
                FIGURES / f"control_{ens['id'].lower()}_mu_0p40.png",
                traj,
                f"{ens['id']} mu=0.40 si_hill",
            )
            if ens["id"] == "SI_BASAL_PERTURB":
                save_fixture("si_basal_perturb_mu_0p40", traj, peaks)

    print("=== FAILED_PRIOR_IVP replay mu=1.5 ===", flush=True)
    prior_records = []
    if not args.smoke:
        prior_mu = float(protocol["failed_prior_replay_mu"])
        for form in (PRODUCTION_SI_HILL, PRODUCTION_POCKETOSCSI_SPLIT):
            replay = DaninoSIDDE(params=params, production_form=form)
            traj, peaks, stats = integrate_arm(replay, prior, prior_mu, protocol)
            rec = record_row(prior, None, prior_mu, "FAILED_PRIOR_IVP", peaks, stats, form)
            rec["arm"] = f"FAILED_PRIOR_IVP_MU_{prior_mu}"
            prior_records.append(rec)
            print(
                f"  form={form} flag={peaks['flag']} T={peaks['period']} "
                f"n_peaks={peaks['n_peaks']} mean_I={stats['mean_I']:.4g}",
                flush=True,
            )
            if form == PRODUCTION_SI_HILL:
                save_fixture("failed_prior_mu_1p5", traj, peaks)
                plot_trace(
                    FIGURES / "failed_prior_mu_1p5.png",
                    traj,
                    "FAILED_PRIOR_IVP mu=1.5 si_hill",
                )

    print("=== TOLERANCE CHECK predeclared subset ===", flush=True)
    if args.smoke:
        tol = {"items": [], "pass": False, "smoke": True}
    else:
        tol = evaluate_tolerance(
            params, protocol, ensemble, ens_by_id, period_rows, control_rows
        )
    print(json.dumps(json_safe(tol), indent=2), flush=True)

    if args.smoke:
        trend = {"applicable": False, "smoke": True}
        extras = {}
        status = "SMOKE"
        reason = "Smoke only. Full grid not scored."
        gates = {
            "D0b": status,
            "reason": reason,
            "d1_may_start": False,
        }
    else:
        gates = score_d0b(protocol, period_rows, control_rows, tol)
        trend = gates
        extras = {
            "n_extra": sum(1 for r in period_rows if r["role"] == "FINITE_INTERVAL_EXTRA"),
            "extra_flags": {
                r["arm"]: r["flag"]
                for r in period_rows
                if r["role"] == "FINITE_INTERVAL_EXTRA"
            },
        }
        status = gates["D0b"]
        reason = gates["reason"]
        plot_period_vs_mu(FIGURES / "period_vs_mu.png", period_rows)
        write_period_table(RESULTS / "PERIOD_TABLE.md", period_rows, primary["id"])

    print("=== TREND / GATES ===", flush=True)
    print(json.dumps(json_safe(trend), indent=2), flush=True)

    write_csv(RESULTS / "period_vs_mu.csv", period_rows)
    write_csv(RESULTS / "controls.csv", control_rows)
    write_csv(RESULTS / "failed_prior_replay.csv", prior_records)

    summary = {
        "primary": ensemble["primary_id"],
        "ics_label": ensemble["label"],
        "smoke": bool(args.smoke),
        "gates": json_safe(gates),
        "trend": json_safe(trend),
        "tolerance": json_safe(tol),
        "finite_interval": extras,
        "period_rows": json_safe(period_rows),
        "controls": json_safe(control_rows),
        "failed_prior": json_safe(prior_records),
        "D0b": status,
        "reason": reason,
        "d1_may_start": bool(gates.get("d1_may_start")),
    }
    (RESULTS / "d0b_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2),
        encoding="utf-8",
    )
    print(f"D0b={summary['D0b']}", flush=True)
    print(summary["reason"], flush=True)


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
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
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
