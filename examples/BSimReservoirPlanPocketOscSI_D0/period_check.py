#!/usr/bin/env python3
"""NARMA-blind period extractor for the D0 bulk oracle.

Peak-to-peak of intracellular LuxI I, frozen in configs/protocol.json.
"""

from __future__ import annotations

import math
import sys

import numpy as np
from scipy.signal import find_peaks


def refuse_narma(argv: list[str] | None = None) -> None:
    blob = " ".join(sys.argv if argv is None else argv).lower()
    if "narma" in blob or "mackey" in blob:
        raise SystemExit("D0 period checker cannot see a NARMA target.")


def _window(t: np.ndarray, y: np.ndarray, t_discard: float) -> tuple[np.ndarray, np.ndarray]:
    mask = t >= t_discard - 1e-12
    return t[mask], y[mask]


def extract_period(
    t: np.ndarray,
    lux_i: np.ndarray,
    protocol: dict,
) -> dict:
    t_discard = float(protocol["t_discard"])
    sample_dt = float(protocol["sample_dt"])
    spacing = float(protocol["peak_min_spacing"])
    rel_prom = float(protocol["peak_rel_prominence"])
    min_peaks = int(protocol["min_peaks"])
    persist_frac = float(protocol["persist_last_fraction"])
    amp_frac = float(protocol["amp_last_fraction"])
    amp_ratio = float(protocol["amp_persist_ratio"])
    rel_amp_min = float(protocol["rel_amplitude_min"])

    t_w, y_w = _window(t, lux_i, t_discard)
    empty = {
        "flag": "NO_PERIOD",
        "n_peaks": 0,
        "period": float("nan"),
        "peak_times": [],
        "rel_amplitude": float("nan"),
        "amp_persist": float("nan"),
        "mean_I": float(np.mean(y_w)) if len(y_w) else float("nan"),
    }
    if len(y_w) < 8:
        return empty

    p5, p95 = np.percentile(y_w, [5.0, 95.0])
    amp = float(p95 - p5)
    median = float(np.median(y_w))
    rel_amp = amp / median if median > 0.0 else 0.0
    if not math.isfinite(amp) or amp <= 0.0:
        return empty

    prominence = rel_prom * amp
    distance = max(1, int(math.ceil(spacing / sample_dt)))
    idx, _ = find_peaks(y_w, prominence=prominence, distance=distance)
    times = [float(t_w[j]) for j in idx]
    n = len(times)

    t_lo, t_hi = float(t_w[0]), float(t_w[-1])
    persist_cut = t_lo + (1.0 - persist_frac) * (t_hi - t_lo)
    n_persist = sum(1 for tau in times if tau >= persist_cut)

    n_amp = max(8, int(amp_frac * len(y_w)))
    amp_first = float(np.percentile(y_w[:n_amp], 95) - np.percentile(y_w[:n_amp], 5))
    amp_last = float(np.percentile(y_w[-n_amp:], 95) - np.percentile(y_w[-n_amp:], 5))
    persist = amp_last / amp_first if amp_first > 0.0 else 0.0

    osc = (
        n >= min_peaks
        and n_persist >= 1
        and persist >= amp_ratio
        and rel_amp >= rel_amp_min
    )
    if osc:
        period = float(np.mean(np.diff(np.array(times))))
        flag = "OSC"
    else:
        period = float("nan")
        flag = "NO_PERIOD"

    return {
        "flag": flag,
        "n_peaks": n,
        "n_persist_peaks": n_persist,
        "period": period,
        "peak_times": times,
        "rel_amplitude": rel_amp,
        "amp_persist": persist,
        "mean_I": float(np.mean(y_w)),
        "p5_I": float(p5),
        "p95_I": float(p95),
    }
