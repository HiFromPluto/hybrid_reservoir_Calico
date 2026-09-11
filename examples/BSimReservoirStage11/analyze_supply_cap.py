#!/usr/bin/env python3
"""Supply-cap confirmation + input-correlation matrix from Stage11 sweep CSVs."""
import csv, json, math, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
SX, SY = 20, 10
CX, CY = 4, 2
RX, RY = SX // CX, SY // CY  # 5 x 5

SWEEPS = {
    "1x": os.path.join(ROOT, "results", "sweep_g1x"),
    "4x": os.path.join(ROOT, "results", "sweep_g4x"),
    "8x": os.path.join(ROOT, "results", "sweep_g8x"),
}

def pearson(x, y):
    n = min(len(x), len(y))
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    num = dx = dy = 0.0
    for i in range(n):
        num += (x[i] - mx) * (y[i] - my)
        dx += (x[i] - mx) ** 2
        dy += (y[i] - my) ** 2
    if dx <= 0 or dy <= 0:
        return 0.0
    return num / math.sqrt(dx * dy)

def t_from_r(r, n):
    if n < 3 or abs(r) >= 1:
        return float("nan")
    return r * math.sqrt((n - 2) / (1 - r * r))

def last_samples(path):
    """Last sample row per window."""
    by_win = {}
    with open(path, newline="") as f:
        r = csv.reader(f, delimiter=";")
        header = next(r)
        for row in r:
            if not row:
                continue
            w, s = int(row[0]), int(row[1])
            by_win[w] = (s, row, header)
    return [by_win[w] for w in sorted(by_win)]

def col_idx(header, prefix):
    return [i for i, h in enumerate(header) if h.startswith(prefix)]

def aggregate_den(den200):
    bins = [0] * (CX * CY)
    for sid, d in enumerate(den200):
        xi, yi = sid // SY, sid % SY
        cid = (xi // RX) * CY + (yi // RY)
        bins[cid] += int(float(d))
    return bins

def load_inputs(n):
    out = []
    for i in range(3):
        p = os.path.join(ROOT, f"input_AC{i}.txt")
        vals = []
        with open(p) as f:
            for line in f:
                line = line.strip()
                if line:
                    vals.append(float(line))
        out.append(vals[:n])
    return out

def analyze_results(dirpath):
    rows = last_samples(os.path.join(dirpath, "stage11_results.csv"))
    n = len(rows)
    pops = [int(r[1][3]) for r in rows]
    deaths = [int(r[1][4]) for r in rows]
    births = [int(r[1][5]) for r in rows]
    u = load_inputs(n)
    metrics = {}
    for name, series in [("N", pops), ("births", births), ("deaths", deaths)]:
        for ch, uu in enumerate(u):
            r = pearson(series, uu)
            metrics[f"r({name},AC{ch})"] = round(r, 3)
            metrics[f"t({name},AC{ch})"] = round(t_from_r(r, n), 2)
        metrics[f"mean_{name}"] = round(sum(series) / n, 1)
    metrics["n_windows"] = n
    metrics["pops"] = pops
    metrics["births"] = births
    metrics["deaths"] = deaths
    metrics["u0"] = u[0]
    return metrics

def analyze_voxels(dirpath):
    rows = last_samples(os.path.join(dirpath, "stage11_voxels.csv"))
    header = rows[0][2]
    den_cols = col_idx(header, "Den_")
    birth_cols = col_idx(header, "Birth_")
    pairs = []  # (den_agg, births) per bin per window
    for _, row, _ in rows:
        den200 = [row[i] for i in den_cols]
        den_agg = aggregate_den(den200)
        births = [int(float(row[i])) for i in birth_cols]
        for d, b in zip(den_agg, births):
            pairs.append((d, b))
    dens = [p[0] for p in pairs]
    brs = [p[1] for p in pairs]
    r = pearson(dens, brs)
    # equal-count bins on density
    order = sorted(range(len(dens)), key=lambda i: dens[i])
    nb = 8
    binned = []
    for k in range(nb):
        lo = k * len(order) // nb
        hi = (k + 1) * len(order) // nb
        idx = order[lo:hi]
        md = sum(dens[i] for i in idx) / max(1, len(idx))
        mb = sum(brs[i] for i in idx) / max(1, len(idx))
        binned.append({
            "bin": k,
            "n_points": len(idx),
            "mean_den": round(md, 2),
            "mean_births": round(mb, 2),
            "den_lo": dens[idx[0]],
            "den_hi": dens[idx[-1]],
        })
    # subsample scatter (every 4th) for canvas
    scatter = [{"den": dens[i], "births": brs[i]} for i in range(0, len(pairs), 4)]
    return {
        "n_pairs": len(pairs),
        "r(births,den_agg)": round(r, 3),
        "t(births,den_agg)": round(t_from_r(r, len(pairs)), 2),
        "mean_den": round(sum(dens) / len(dens), 2),
        "mean_births": round(sum(brs) / len(brs), 2),
        "binned": binned,
        "scatter": scatter,
        "den_min": min(dens),
        "den_max": max(dens),
        "births_min": min(brs),
        "births_max": max(brs),
    }

def main():
    report = {}
    for tag, d in SWEEPS.items():
        res_csv = os.path.join(d, "stage11_results.csv")
        vox_csv = os.path.join(d, "stage11_voxels.csv")
        if not os.path.exists(res_csv):
            continue
        block = {"results": analyze_results(d)}
        if os.path.exists(vox_csv):
            block["voxels"] = analyze_voxels(d)
        report[tag] = block

    out_json = os.path.join(ROOT, "results", "supply_cap_analysis.json")
    # drop bulky series from file used by humans; keep scatter in json for canvas
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2)

    txt = os.path.join(ROOT, "results", "supply_cap_analysis.txt")
    lines = ["Supply-cap + input-coupling matrix", "=" * 40, ""]
    lines.append("Per-bin Birth vs aggregated Den (20x10 -> 4x2), last sample/window")
    for tag, block in report.items():
        v = block.get("voxels", {})
        lines.append(
            f"  {tag}: r(births,den)={v.get('r(births,den_agg)')}  "
            f"t={v.get('t(births,den_agg)')}  n={v.get('n_pairs')}  "
            f"den[{v.get('den_min')},{v.get('den_max')}]  "
            f"births_mean={v.get('mean_births')}"
        )
        for b in v.get("binned", []):
            lines.append(
                f"      bin{b['bin']}: den={b['mean_den']:.1f} "
                f"[{b['den_lo']}-{b['den_hi']}]  mean_births={b['mean_births']:.2f}"
            )
    lines.append("")
    lines.append("Input correlations (last sample / window, n=40)")
    hdr = f"{'':4} {'r(N,AC0)':>10} {'r(N,AC1)':>10} {'r(N,AC2)':>10} " \
          f"{'r(B,AC0)':>10} {'r(D,AC0)':>10} {'r(B,AC1)':>10} {'r(D,AC1)':>10} " \
          f"{'r(B,AC2)':>10} {'r(D,AC2)':>10}"
    lines.append(hdr)
    for tag, block in report.items():
        m = block["results"]
        lines.append(
            f"{tag:4} {m['r(N,AC0)']:10.3f} {m['r(N,AC1)']:10.3f} {m['r(N,AC2)']:10.3f} "
            f"{m['r(births,AC0)']:10.3f} {m['r(deaths,AC0)']:10.3f} "
            f"{m['r(births,AC1)']:10.3f} {m['r(deaths,AC1)']:10.3f} "
            f"{m['r(births,AC2)']:10.3f} {m['r(deaths,AC2)']:10.3f}"
        )
    text = "\n".join(lines) + "\n"
    with open(txt, "w") as f:
        f.write(text)
    print(text)
    print("wrote", out_json)

if __name__ == "__main__":
    main()
