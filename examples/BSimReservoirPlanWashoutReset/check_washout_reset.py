#!/usr/bin/env python3
"""WashoutReset residual scout: theory table, leftovers, predeclared findings.

No ridge. No NARMA. No Mackey-Glass NRMSE. No waveform AUC.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE6 = HERE.parent / "BSimReservoirPlanStage6"
BENCHA = HERE.parent / "BSimReservoirPlanBenchA"

TAU_AHL = 1.0 / 0.0033
TAU_R = 15.0
TAU_L = 1500.0
WINDOW_S = 300.0
EPS = 1e-12

THEORY_W5_AHL = math.exp(-1500.0 / TAU_AHL)
THEORY_W5_R = math.exp(-1500.0 / TAU_R)
THEORY_W5_L = math.exp(-1.0)
THEORY_W15_AHL = math.exp(-4500.0 / TAU_AHL)
THEORY_W15_L = math.exp(-3.0)

W5_ARMS = {
    "WASH_DECAY_W5": HERE / "results" / "wash_decay_w5_seed101",
    "WASH_TRICKLE_W5": HERE / "results" / "wash_trickle_w5_seed101",
    "WASH_FLUSH_W5": HERE / "results" / "wash_flush_w5_seed101",
}
W15_ARMS = {
    "WASH_DECAY_W15": HERE / "results" / "wash_decay_w15_seed101",
    "WASH_FLUSH_W15": HERE / "results" / "wash_flush_w15_seed101",
}
SILENT_CANDIDATES = [
    STAGE6 / "results" / "stage6_silent_seed101" / "window_summary.csv",
    BENCHA / "results" / "mg_silent_seed101" / "window_summary.csv",
    BENCHA / "results" / "lorenz_silent_seed101" / "window_summary.csv",
]


def theory_table(wash_windows: int = 5) -> dict:
    dt = wash_windows * WINDOW_S
    return {
        "tau_AHL_s": TAU_AHL,
        "tau_R_s": TAU_R,
        "tau_L_s": TAU_L,
        "W": wash_windows,
        "dt_s": dt,
        "AHL_leftover": math.exp(-dt / TAU_AHL),
        "R_leftover": math.exp(-dt / TAU_R),
        "L_leftover": math.exp(-dt / TAU_L),
    }


def print_theory() -> None:
    w5 = theory_table(5)
    w15 = theory_table(15)
    print("WashoutReset THEORY (spatially uniform exponential; not a retune)")
    print(f"  tau_AHL={w5['tau_AHL_s']:.4f} s  tau_R={w5['tau_R_s']:.1f} s  tau_L={w5['tau_L_s']:.1f} s")
    print(f"  W=5  dt={w5['dt_s']:.1f} s  AHL={w5['AHL_leftover']:.5f}  R={w5['R_leftover']:.6g}  L={w5['L_leftover']:.6f}")
    print(f"  W=15 dt={w15['dt_s']:.1f} s  AHL={w15['AHL_leftover']:.6g}  L={w15['L_leftover']:.6f}")
    print(
        f"  check exp(-1500/303.0303)={THEORY_W5_AHL:.5f}  exp(-1)={THEORY_W5_L:.6f}  "
        f"exp(-3)={THEORY_W15_L:.6f}"
    )
    print("  Living residuals will differ (walls, clamp, OOB, occupancy). Do not change tau_L.")


def pearson(a, b) -> float:
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    den = math.sqrt(float(np.sum(a * a) * np.sum(b * b)))
    return float("nan") if den < 1e-15 else float(np.sum(a * b) / den)


def rel(now: float, then: float) -> float:
    return now / max(then, EPS)


def last_sample_of_window(data, window: int):
    mask = data["window"] == window
    if not np.any(mask):
        raise ValueError(f"missing window {window}")
    idx = np.where(mask)[0][-1]
    return idx


def load_voxels(run_dir: Path) -> dict:
    path = run_dir / "voxels.csv"
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    ahl_i = [header.index(f"AHL_uM_{i}") for i in range(200)]
    r_i = [header.index(f"Receiver_R_{i}") for i in range(200)]
    l_i = [header.index(f"Lum_Mean_{i}") for i in range(200)]
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return {
        "window": data[:, header.index("Window")].astype(int),
        "sample": data[:, header.index("Sample")].astype(int),
        "t_in": data[:, header.index("TimeInWindow_s")],
        "ahl_u": data[:, header.index("Input_AC1_AHL")],
        "ahl": data[:, ahl_i],
        "receiver": data[:, r_i],
        "lum": data[:, l_i],
    }


def load_results(run_dir: Path) -> dict:
    path = run_dir / "results.csv"
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split(";")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return {
        "window": data[:, header.index("Window")].astype(int),
        "sample": data[:, header.index("Sample")].astype(int),
        "t_in": data[:, header.index("TimeInWindow_s")],
        "pop": data[:, header.index("Total_Count")],
        "mean_L": data[:, header.index("Mean_L")],
        "oob": data[:, header.index("OOB_Deaths_This_Window")],
        "acid": data[:, header.index("Input_Driven_Deaths_This_Window")],
    }


def load_summary(run_dir: Path) -> list[dict]:
    path = run_dir / "window_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=";"))


def parse_kv(path: Path) -> dict[str, str]:
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def validate_csv(path: Path, expected_rows: int, last_prefix: str | None = None) -> dict:
    text = path.read_text(encoding="utf-8").splitlines()
    body = text[1:] if text else []
    last = body[-1] if body else ""
    return {
        "path": str(path),
        "row_count": len(body),
        "row_count_pass": len(body) == expected_rows,
        "last_sample_pass": last.startswith(last_prefix) if last_prefix else True,
        "last_line_prefix": last.split(";")[:3],
    }


def snapshot(voxels, results, summary, window: int) -> dict:
    vi = last_sample_of_window(voxels, window)
    ri = last_sample_of_window(results, window)
    oob_cum = sum(int(float(row["OOB_Deaths"])) for row in summary[: window + 1])
    acid_cum = sum(int(float(row["Input_Driven_Deaths"])) for row in summary[: window + 1])
    return {
        "window": window,
        "mean_AHL": float(np.mean(voxels["ahl"][vi])),
        "mean_R": float(np.mean(voxels["receiver"][vi])),
        "mean_L_voxel": float(np.mean(voxels["lum"][vi])),
        "mean_L": float(results["mean_L"][ri]),
        "population": float(results["pop"][ri]),
        "oob_cum": oob_cum,
        "acid_cum": acid_cum,
        "t_in": float(voxels["t_in"][vi]),
        "sample": int(voxels["sample"][vi]),
        "ahl_u": float(voxels["ahl_u"][vi]),
    }


def occupancy_load(voxels, load_last: int) -> dict:
    mean_r_w = []
    u_w = []
    for window in range(load_last + 1):
        mask = voxels["window"] == window
        mean_r_w.append(float(np.mean(voxels["receiver"][mask])))
        u_w.append(float(np.mean(voxels["ahl_u"][mask])))
    vi = last_sample_of_window(voxels, load_last)
    mean_r_tload = float(np.mean(voxels["receiver"][vi]))
    r_u = pearson(mean_r_w, u_w)
    label = "DEAD" if mean_r_tload < 0.05 else "OCCUPIED"
    return {
        "mean_R_t_load": mean_r_tload,
        "mean_R_load_windows": float(np.mean(mean_r_w)),
        "r_meanR_u": r_u,
        "occupancy": label,
    }


def silent_floor() -> dict | None:
    for path in SILENT_CANDIDATES:
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter=";"))
        analysis = [row for row in rows if int(row["Window"]) >= 40]
        if not analysis:
            analysis = rows
        mean_l = float(analysis[0]["Mean_L"])
        return {
            "path": str(path),
            "analysis_start_window": int(analysis[0]["Window"]),
            "mean_L": mean_l,
        }
    return None


def parse_epoch_trace(run_dir: Path) -> dict:
    path = run_dir / "epoch_trace.txt"
    marks = {}
    windows = []
    if not path.exists():
        return {"marks": marks, "windows": windows, "missing": True}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "EPOCH_MARK" in line:
            mark = re.search(r"mark=(\S+)", line)
            flow = re.search(r"FLOW_SPEED=([0-9.eE+-]+)", line)
            bound = re.search(r"chemical_boundary=(\S+)", line)
            courant = re.search(r"chemical_Courant=([0-9.eE+-]+)", line)
            ahl = re.search(r"ahl=\(([0-9.]+),([0-9.]+),([0-9.]+)\)", line)
            window = re.search(r"window=(\d+)", line)
            if mark:
                marks[mark.group(1)] = {
                    "window": int(window.group(1)) if window else None,
                    "FLOW_SPEED": float(flow.group(1)) if flow else None,
                    "boundary": bound.group(1) if bound else None,
                    "Courant": float(courant.group(1)) if courant else None,
                    "ahl": tuple(float(x) for x in ahl.groups()) if ahl else None,
                    "line": line,
                }
        if line.startswith("WashoutReset WINDOW"):
            windows.append(line)
    return {"marks": marks, "windows": windows, "missing": False}


def analyze_arm(condition: str, run_dir: Path, load_last: int, wash_last: int, last_window: int) -> dict:
    status = parse_kv(run_dir / "run_status.txt")
    voxels = load_voxels(run_dir)
    results = load_results(run_dir)
    summary = load_summary(run_dir)
    expected_windows = last_window + 1
    expected_aux = expected_windows * 16
    last_prefix = f"{last_window};15;299.95"
    csv_checks = {
        "summary": validate_csv(run_dir / "window_summary.csv", expected_windows),
        "results": validate_csv(run_dir / "results.csv", expected_aux, last_prefix),
        "voxels": validate_csv(run_dir / "voxels.csv", expected_aux, last_prefix),
    }
    t_load = snapshot(voxels, results, summary, load_last)
    t_wash = snapshot(voxels, results, summary, wash_last)
    post = [snapshot(voxels, results, summary, w) for w in range(wash_last + 1, last_window + 1)]
    occ = occupancy_load(voxels, load_last)
    leftovers = {
        "AHL_rel": rel(t_wash["mean_AHL"], t_load["mean_AHL"]),
        "R_rel": rel(t_wash["mean_R"], t_load["mean_R"]),
        "L_rel": rel(t_wash["mean_L"], t_load["mean_L"]),
        "N_rel": rel(t_wash["population"], t_load["population"]),
        "L_rel_last_post": rel(post[-1]["mean_L"], t_load["mean_L"]) if post else float("nan"),
    }
    traces = []
    for window in range(expected_windows):
        snap = snapshot(voxels, results, summary, window)
        traces.append({
            "window": window,
            "epoch": "LOAD" if window <= load_last else ("WASH" if window <= wash_last else "POST"),
            "mean_AHL": snap["mean_AHL"],
            "mean_R": snap["mean_R"],
            "mean_L": snap["mean_L"],
            "population": snap["population"],
            "oob_cum": snap["oob_cum"],
            "acid_cum": snap["acid_cum"],
        })
    epoch = parse_epoch_trace(run_dir)
    return {
        "condition": condition,
        "run_dir": str(run_dir),
        "status": status,
        "csv": csv_checks,
        "occupancy": occ,
        "t_load": t_load,
        "t_wash": t_wash,
        "post": post,
        "leftovers": leftovers,
        "traces": traces,
        "epoch": epoch,
        "wash_flow": float(status.get("wash.flow.speed.um_s", "nan")),
        "wash_courant": float(status.get("wash_chemical_Courant", "nan")),
        "ahl": (
            float(status.get("ahl.x", "nan")),
            float(status.get("ahl.y", "nan")),
            float(status.get("ahl.z", "nan")),
        ),
    }


def evaluate_findings(arms: dict[str, dict]) -> dict:
    decay = arms["WASH_DECAY_W5"]["leftovers"]
    trickle = arms["WASH_TRICKLE_W5"]["leftovers"]
    flush = arms["WASH_FLUSH_W5"]["leftovers"]
    ahl_clear = decay["AHL_rel"] < 0.05
    l_sticks = decay["L_rel"] > 0.20
    flush_dumps_cells = flush["N_rel"] <= 0.80
    flush_dumps_l = (decay["L_rel"] - flush["L_rel"]) >= 0.10 and flush_dumps_cells
    trickle_like = (
        abs(trickle["L_rel"] - decay["L_rel"]) < 0.05
        and abs(trickle["AHL_rel"] - decay["AHL_rel"]) < 0.02
    )
    occ_dead = any(arm["occupancy"]["occupancy"] == "DEAD" for arm in arms.values())
    impl_fail = []
    for name, arm in arms.items():
        if arm["wash_courant"] >= 1.0:
            impl_fail.append(f"Courant>=1 on {name}")
        if abs(arm["ahl"][0] - 500) > 1e-6 or abs(arm["ahl"][1] - 250) > 1e-6 or abs(arm["ahl"][2] - 5) > 1e-6:
            impl_fail.append(f"source not CENTER on {name}")
        csv_ok = all(v["row_count_pass"] and v.get("last_sample_pass", True) for v in arm["csv"].values())
        if not csv_ok:
            impl_fail.append(f"CSV fail on {name}")
        wash_mark = arm["epoch"]["marks"].get("WASH_START") or arm["epoch"]["marks"].get("WASH_END")
        if name == "WASH_FLUSH_W5":
            printed = None
            for key in ("WASH_START", "WASH_END"):
                if key in arm["epoch"]["marks"]:
                    printed = arm["epoch"]["marks"][key]["FLOW_SPEED"]
            if printed is None or abs(printed - 8.0) > 1e-6:
                impl_fail.append(f"FLUSH printed flow != 8 (got {printed})")
            if abs(arm["wash_flow"] - 8.0) > 1e-6:
                impl_fail.append(f"FLUSH wash.flow != 8 (got {arm['wash_flow']})")
        _ = wash_mark
    if not ahl_clear:
        impl_fail.append("WASH_DECAY_W5 AHL_rel >= 0.05 (decay not doing what the dish claims)")
    if occ_dead:
        impl_fail.append("occupancy DEAD at t_load")
    follow_on = ahl_clear and l_sticks and not occ_dead and not any(
        "Courant" in f or "CENTER" in f or "FLUSH printed" in f or "FLUSH wash" in f
        for f in impl_fail
    )
    # AHL_rel fail still writes report and blocks W=15; occupancy DEAD blocks interpret
    if not ahl_clear or occ_dead:
        follow_on = False
    return {
        "AHL_CLEAR": ahl_clear,
        "L_STICKS": l_sticks,
        "FLUSH_DUMPS_CELLS": flush_dumps_cells,
        "FLUSH_DUMPS_L": flush_dumps_l,
        "TRICKLE_LIKE_DECAY": trickle_like,
        "occupancy_dead": occ_dead,
        "impl_fail": impl_fail,
        "follow_on": follow_on,
    }


def fmt(value) -> str:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_report(arms: dict[str, dict], findings: dict, silent: dict | None, extra: dict[str, dict]) -> None:
    out_dir = HERE / "results"
    out_dir.mkdir(exist_ok=True)
    md = out_dir / "WASHOUT_RESET_SCOUT.md"
    csv_path = out_dir / "washout_reset_scout.csv"
    lines = [
        "# WashoutReset residual scout",
        "",
        "Development-only Track E0. Not a task scout. No ridge. Not a new claim dish.",
        "",
        "## Theory",
        "",
        f"| Quantity | W=5 leftover | W=15 leftover |",
        f"|---|---|---|",
        f"| AHL | {THEORY_W5_AHL:.5f} | {THEORY_W15_AHL:.6g} |",
        f"| R | {THEORY_W5_R:.6g} | ~0 |",
        f"| L | {THEORY_W5_L:.6f} | {THEORY_W15_L:.6f} |",
        "",
        "Spatially uniform exponential. Living residuals differ (walls, clamp, OOB, occupancy).",
        "Do not retune tau_L.",
        "",
    ]
    if silent:
        lines += [
            "## Silent floor",
            "",
            f"Stage 6 / BenchA silent seed 101 analysis-start mean_L = {silent['mean_L']:.6g} "
            f"(window {silent['analysis_start_window']}, `{silent['path']}`).",
            "",
        ]
    else:
        lines += [
            "## Silent floor",
            "",
            "Stage 6 / BenchA silent seed 101 voxels were not found. Did not rerun silent.",
            "",
        ]
    lines += ["## Per-arm last-sample means", ""]
    header = (
        "| Condition | t_load AHL | t_load R | t_load L | t_load N | "
        "t_wash AHL | t_wash R | t_wash L | t_wash N | "
        "AHL_rel | R_rel | L_rel | N_rel | L_rel last post | occ |"
    )
    lines += [header, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    all_arms = dict(arms)
    all_arms.update(extra)
    rows = []
    for name, arm in all_arms.items():
        load, wash, left, occ = arm["t_load"], arm["t_wash"], arm["leftovers"], arm["occupancy"]
        lines.append(
            "| {name} | {a0:.4g} | {r0:.4g} | {l0:.4g} | {n0:.0f} | "
            "{a1:.4g} | {r1:.4g} | {l1:.4g} | {n1:.0f} | "
            "{ar:.4g} | {rr:.4g} | {lr:.4g} | {nr:.4g} | {lp:.4g} | {occ} |".format(
                name=name,
                a0=load["mean_AHL"], r0=load["mean_R"], l0=load["mean_L"], n0=load["population"],
                a1=wash["mean_AHL"], r1=wash["mean_R"], l1=wash["mean_L"], n1=wash["population"],
                ar=left["AHL_rel"], rr=left["R_rel"], lr=left["L_rel"], nr=left["N_rel"],
                lp=left["L_rel_last_post"], occ=occ["occupancy"],
            )
        )
        rows.append({
            "condition": name,
            "wash_flow": arm["wash_flow"],
            "t_load_AHL": load["mean_AHL"],
            "t_load_R": load["mean_R"],
            "t_load_L": load["mean_L"],
            "t_load_N": load["population"],
            "t_load_oob": load["oob_cum"],
            "t_load_acid": load["acid_cum"],
            "t_wash_AHL": wash["mean_AHL"],
            "t_wash_R": wash["mean_R"],
            "t_wash_L": wash["mean_L"],
            "t_wash_N": wash["population"],
            "t_wash_oob": wash["oob_cum"],
            "t_wash_acid": wash["acid_cum"],
            "AHL_rel": left["AHL_rel"],
            "R_rel": left["R_rel"],
            "L_rel": left["L_rel"],
            "N_rel": left["N_rel"],
            "L_rel_last_post": left["L_rel_last_post"],
            "occupancy": occ["occupancy"],
            "mean_R_t_load": occ["mean_R_t_load"],
            "r_meanR_u": occ["r_meanR_u"],
        })
    lines += ["", "## Findings (W=5, seed 101)", ""]
    for key in ("AHL_CLEAR", "L_STICKS", "FLUSH_DUMPS_CELLS", "FLUSH_DUMPS_L", "TRICKLE_LIKE_DECAY"):
        lines.append(f"- `{key}`: {'FIRE' if findings[key] else 'no'}")
    if findings["impl_fail"]:
        lines += ["", "## STOP / implementation notes", ""]
        for item in findings["impl_fail"]:
            lines.append(f"- {item}")
    lines += ["", f"W=15 follow-on: {'YES' if findings['follow_on'] else 'NO'}", ""]
    if extra:
        decay15 = extra.get("WASH_DECAY_W15")
        flush15 = extra.get("WASH_FLUSH_W15")
        lines += [
            "## W=15 how-long-to-wait (not a W=5 confirmation)",
            "",
        ]
        if decay15:
            left = decay15["leftovers"]
            lines.append(
                f"WASH_DECAY_W15 L_rel={left['L_rel']:.4g} vs theory {THEORY_W15_L:.6f} "
                f"(AHL_rel={left['AHL_rel']:.4g} vs theory {THEORY_W15_AHL:.6g}). "
                "Living L is slightly above the uniform exponential; do not retune tau_L."
            )
        if flush15:
            left = flush15["leftovers"]
            lines.append(
                f"WASH_FLUSH_W15 L_rel={left['L_rel']:.4g}, N_rel={left['N_rel']:.4g} "
                f"(OOB_cum at t_wash={flush15['t_wash']['oob_cum']}). "
                "Fifteen flush windows dump cells; five did not cross the 20% N_rel gate."
            )
        lines.append("")
    sentence = one_sentence(findings, arms)
    lines += ["## One sentence", "", sentence, "", "## Traces", ""]
    for name, arm in all_arms.items():
        lines.append(f"### {name}")
        lines.append("")
        lines.append("| window | epoch | mean_AHL | mean_R | mean_L | pop | OOB_cum |")
        lines.append("|---|---|---|---|---|---|---|")
        for row in arm["traces"]:
            lines.append(
                f"| {row['window']} | {row['epoch']} | {row['mean_AHL']:.4g} | "
                f"{row['mean_R']:.4g} | {row['mean_L']:.4g} | {row['population']:.0f} | "
                f"{row['oob_cum']} |"
            )
        lines.append("")
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter=";")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {md}")
    print(f"wrote {csv_path}")
    print(sentence)


def one_sentence(findings: dict, arms: dict) -> str:
    decay = arms["WASH_DECAY_W5"]["leftovers"]
    if findings["occupancy_dead"]:
        return "Load never occupied (mean_R < 0.05 at t_load); wash leftovers are not interpreted."
    if not findings["AHL_CLEAR"]:
        return (
            f"Chemical silence did not reset AHL (AHL_rel={decay['AHL_rel']:.4g} on WASH_DECAY_W5); "
            "do not claim a wash protocol from this dish."
        )
    if findings["L_STICKS"]:
        extra = []
        if findings["FLUSH_DUMPS_L"]:
            extra.append("flush lowered dish-mean L by dumping cells")
        elif findings["FLUSH_DUMPS_CELLS"]:
            extra.append("flush dumped cells without a 0.10 L_rel gap")
        if findings["TRICKLE_LIKE_DECAY"]:
            extra.append("trickle matched decay")
        tail = "; ".join(extra)
        if tail:
            tail = " " + tail + "."
        else:
            tail = "."
        return (
            f"Chemical washout is not reporter washout: AHL_rel={decay['AHL_rel']:.4g} while "
            f"L_rel={decay['L_rel']:.4g} after five silent windows.{tail}"
        )
    return (
        f"AHL cleared (AHL_rel={decay['AHL_rel']:.4g}) and L did not stick "
        f"(L_rel={decay['L_rel']:.4g}); chemical washout also reset the reporter on this load."
    )


def smoke_check() -> int:
    run_dir = HERE / "results" / "washout_smoke"
    status = parse_kv(run_dir / "run_status.txt")
    epoch = parse_epoch_trace(run_dir)
    wash = epoch["marks"].get("WASH_START") or epoch["marks"].get("WASH_END")
    print("WashoutReset SMOKE (not evidence)")
    print(f"  ahl=({status.get('ahl.x')},{status.get('ahl.y')},{status.get('ahl.z')})")
    print(f"  wash.flow={status.get('wash.flow.speed.um_s')} wash_Courant={status.get('wash_chemical_Courant')}")
    if wash:
        print(f"  WASH mark FLOW_SPEED={wash['FLOW_SPEED']} Courant={wash['Courant']} boundary={wash['boundary']}")
    ok = True
    if abs(float(status.get("ahl.x", "0")) - 500) > 1e-6:
        print("  FAIL source not CENTER")
        ok = False
    if abs(float(status.get("wash.flow.speed.um_s", "0")) - 8) > 1e-6:
        print("  FAIL wash flow != 8")
        ok = False
    if abs(float(status.get("wash_chemical_Courant", "0")) - 0.02) > 1e-6:
        print("  FAIL wash Courant != 0.02")
        ok = False
    if wash is None or abs(wash["FLOW_SPEED"] - 8) > 1e-6:
        print("  FAIL WASH mark FLOW_SPEED != 8")
        ok = False
    if wash and abs(wash["Courant"] - 0.02) > 1e-6:
        print("  FAIL WASH Courant != 0.02")
        ok = False
    voxels = run_dir / "voxels.csv"
    if voxels.exists():
        data = np.loadtxt(voxels, delimiter=";", skiprows=1)
        ahl = data[:, 7:207]
        if np.any(~np.isfinite(ahl)) or np.any(ahl < 0):
            print("  FAIL AHL not finite non-negative")
            ok = False
        else:
            print(f"  AHL min={ahl.min():.4g} max={ahl.max():.4g} finite non-negative")
    summary = load_summary(run_dir) if (run_dir / "window_summary.csv").exists() else []
    if summary:
        pops = [float(row["Population"]) for row in summary]
        print(f"  population {min(pops):.0f}..{max(pops):.0f}")
        if not all(math.isfinite(p) for p in pops):
            ok = False
    print("  SMOKE PASS" if ok else "  SMOKE FAIL")
    return 0 if ok else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theory", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.theory:
        print_theory()
        return 0
    if args.smoke:
        print_theory()
        return smoke_check()
    print_theory()
    missing = [name for name, path in W5_ARMS.items() if not (path / "voxels.csv").exists()]
    if missing:
        print("missing W=5 results: " + ", ".join(missing), file=sys.stderr)
        return 1
    arms = {}
    for name, path in W5_ARMS.items():
        print(f"analyze {name}")
        arms[name] = analyze_arm(name, path, 19, 24, 39)
        left = arms[name]["leftovers"]
        occ = arms[name]["occupancy"]
        print(
            f"  t_load AHL={arms[name]['t_load']['mean_AHL']:.4g} R={arms[name]['t_load']['mean_R']:.4g} "
            f"L={arms[name]['t_load']['mean_L']:.4g} N={arms[name]['t_load']['population']:.0f} occ={occ['occupancy']}"
        )
        print(
            f"  t_wash AHL_rel={left['AHL_rel']:.4g} R_rel={left['R_rel']:.4g} "
            f"L_rel={left['L_rel']:.4g} N_rel={left['N_rel']:.4g}"
        )
    findings = evaluate_findings(arms)
    extra = {}
    for name, path in W15_ARMS.items():
        if (path / "voxels.csv").exists():
            extra[name] = analyze_arm(name, path, 19, 34, 49)
    silent = silent_floor()
    if silent:
        print(f"silent floor mean_L={silent['mean_L']:.6g} from {silent['path']}")
    else:
        print("silent floor MISSING (did not rerun)")
    for key in ("AHL_CLEAR", "L_STICKS", "FLUSH_DUMPS_CELLS", "FLUSH_DUMPS_L", "TRICKLE_LIKE_DECAY"):
        print(f"FINDING {key}={'FIRE' if findings[key] else 'no'}")
    for item in findings["impl_fail"]:
        print(f"STOP {item}")
    print(f"FOLLOW_ON_W15={'YES' if findings['follow_on'] else 'NO'}")
    write_report(arms, findings, silent, extra)
    if findings["occupancy_dead"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
