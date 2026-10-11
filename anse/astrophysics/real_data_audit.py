#!/usr/bin/env python3
"""
Independent real-data checks for TIC 4206066 (no synthetic inputs, no constants
copied from the claims being tested).

1. verify_consumed_zenodo_files   - SHA-256 of every Zenodo file we read vs the
                                    deposit's own SHA256SUMS.
2. leave_one_sector_out_timing    - genuine hold-out: fit a weighted linear
                                    ephemeris on two sectors, predict the third.
3. photometric_depth_at_ephemeris - fold the real MAST/QLP photometry on the
                                    published ephemeris and measure the depth.
4. weighted_odd_even              - inverse-variance odd/even comparison using
                                    per-event errors (depth / S/N).
"""

from __future__ import annotations

import csv
import hashlib
import math
import os
from typing import Any

import numpy as np

from anse.astrophysics.ephemeris import Ephemeris, load_ephemeris_a

ZENODO_PART1 = "data/astrophysics/tic4206066/zenodo/part1"
CONSUMED_FILES = (
    "02_transits_A/per_transit_A.csv",
    "02_transits_A/ephemeris_A.json",
)


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _find_part1_dir(preferred: str = ZENODO_PART1) -> str:
    for candidate in [
        preferred,
        "data/zenodo/part1",
        os.path.join(os.path.dirname(__file__), "..", "..", "data", "zenodo", "part1"),
        os.path.join(os.path.dirname(__file__), "..", "..", "data", "astrophysics", "tic4206066", "zenodo", "part1"),
    ]:
        if os.path.exists(os.path.join(candidate, "SHA256SUMS")):
            return candidate
    return preferred


def verify_consumed_zenodo_files(part1_dir: str = ZENODO_PART1) -> dict[str, Any]:
    """Compare the SHA-256 of each file this pipeline reads to the deposit's SHA256SUMS."""
    part1_dir = _find_part1_dir(part1_dir)
    sums: dict[str, str] = {}
    with open(os.path.join(part1_dir, "SHA256SUMS")) as f:
        for line in f:
            digest, _, rel = line.strip().partition("  ")
            sums[os.path.normpath(rel)] = digest
    report: dict[str, Any] = {}
    for rel in CONSUMED_FILES:
        expected = sums.get(os.path.normpath(rel))
        actual = sha256_file(os.path.join(part1_dir, rel))
        report[rel] = {
            "expected": expected,
            "actual": actual,
            "match": expected is not None and expected == actual,
            "listed_in_sha256sums": expected is not None,
        }
    report["all_match"] = all(v["match"] for v in report.values() if isinstance(v, dict))
    return report


def load_events(part1_dir: str = ZENODO_PART1) -> list[dict[str, Any]]:
    """Per-event table including the injection-calibrated timing error (minutes)."""
    path = os.path.join(part1_dir, "02_transits_A/per_transit_A.csv")
    events = []
    with open(path) as f:
        for row in csv.DictReader(f):
            depth = float(row["depth_ppm_verifier_event"])
            snr = float(row["snr_approx_sqrt_dchi2"])
            events.append(
                {
                    "sector": int(row["sector"]),
                    "epoch": int(row["epoch"]),
                    "t_mid_bjd": float(row["t_mid_BJD_TDB"]),
                    "err_min": float(row["err_min_injection_sd"]),
                    "depth_ppm": depth,
                    "depth_err_ppm": depth / snr if snr > 0 else float("inf"),
                    "used_in_timing": row["used_in_19_timing_fit"].strip().lower() == "true",
                }
            )
    return events


def _weighted_line(epochs: np.ndarray, t: np.ndarray, sig_d: np.ndarray):
    w = 1.0 / sig_d**2
    A = np.vstack([np.ones_like(epochs, dtype=float), epochs.astype(float)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    beta = cov @ (A.T @ (w * t))
    return beta, cov


def leave_one_sector_out_timing(events: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """
    For each sector S: fit T0, P on the timing events of the *other* sectors,
    predict the measured mid-times in S and report residuals. The prediction
    uncertainty combines the fit covariance and the per-event measurement error.
    """
    events = load_events() if events is None else events
    used = [e for e in events if e["used_in_timing"]]
    out: dict[str, Any] = {}
    for held in sorted({e["sector"] for e in used}):
        train = [e for e in used if e["sector"] != held]
        test = [e for e in used if e["sector"] == held]
        ep = np.array([e["epoch"] for e in train])
        t = np.array([e["t_mid_bjd"] for e in train])
        sig = np.array([e["err_min"] for e in train]) / 1440.0
        beta, cov = _weighted_line(ep, t, sig)
        res_min, z = [], []
        for e in test:
            pred = beta[0] + beta[1] * e["epoch"]
            var_fit = cov[0, 0] + 2 * e["epoch"] * cov[0, 1] + e["epoch"] ** 2 * cov[1, 1]
            sigma = math.sqrt(var_fit + (e["err_min"] / 1440.0) ** 2)
            res_min.append((e["t_mid_bjd"] - pred) * 1440.0)
            z.append((e["t_mid_bjd"] - pred) / sigma)
        out[f"S{held}"] = {
            "n_train": len(train),
            "n_test": len(test),
            "residual_min": [round(r, 2) for r in res_min],
            "rms_residual_min": round(float(np.sqrt(np.mean(np.square(res_min)))), 2),
            "chi2": round(float(np.sum(np.square(z))), 2),
            "chi2_per_event": round(float(np.mean(np.square(z))), 2),
        }
    out["max_chi2_per_event"] = max(v["chi2_per_event"] for v in out.values())
    out["passed"] = out["max_chi2_per_event"] < 3.0
    return out


def photometric_depth_at_ephemeris(
    mast_lcs: dict[int, dict[str, np.ndarray]],
    eph: Ephemeris | None = None,
    duration_h: float = 2.1,
) -> dict[int, dict[str, float]]:
    """
    Mean in-transit depth of the *real* light curve folded on the published
    ephemeris, with an error from the out-of-transit scatter. Verifies the
    signal is actually present at the stated epoch (and, via the phase control,
    absent half an orbit away).
    """
    eph = eph or load_ephemeris_a()
    P, t0 = eph.period_d, eph.t0_btjd
    half = duration_h / 24.0 / 2.0 * 0.8  # central 80% of the transit (flat bottom)
    res: dict[int, dict[str, float]] = {}
    for sec, d in mast_lcs.items():
        t, f = d["times"], d["fluxes"]

        def box(center_phase_days: float):
            ph = ((t - t0 - center_phase_days + P / 2) % P) - P / 2
            inn = np.abs(ph) < half
            out = (np.abs(ph) > duration_h / 24.0) & (np.abs(ph) < 0.3 * P)
            return inn, out

        inn, out = box(0.0)
        inn_c, out_c = box(P / 2.0)  # phase control: secondary-eclipse position
        if inn.sum() < 5 or out.sum() < 20:
            continue
        base = np.median(f[out])
        depth = (base - np.mean(f[inn])) * 1e6
        err = np.std(f[out]) / math.sqrt(inn.sum()) * 1e6
        base_c = np.median(f[out_c]) if out_c.sum() >= 20 else base
        depth_c = (base_c - np.mean(f[inn_c])) * 1e6 if inn_c.sum() >= 5 else float("nan")
        res[sec] = {
            "n_in_transit": int(inn.sum()),
            "depth_ppm": round(float(depth), 1),
            "depth_err_ppm": round(float(err), 1),
            "snr": round(float(depth / err), 2),
            "control_phase_0.5_depth_ppm": round(float(depth_c), 1),
        }
    return res


def weighted_odd_even(events: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """Inverse-variance mean depth of odd vs even epochs using per-event errors."""
    events = load_events() if events is None else events

    def wmean(sel):
        d = np.array([e["depth_ppm"] for e in sel])
        s = np.array([e["depth_err_ppm"] for e in sel])
        w = 1.0 / s**2
        return float(np.sum(w * d) / np.sum(w)), float(1.0 / math.sqrt(np.sum(w)))

    odd = [e for e in events if e["epoch"] % 2 != 0]
    even = [e for e in events if e["epoch"] % 2 == 0]
    mo, so = wmean(odd)
    me, se = wmean(even)
    sig = abs(mo - me) / math.hypot(so, se)
    return {
        "n_odd": len(odd),
        "n_even": len(even),
        "odd_ppm": round(mo, 1),
        "odd_err_ppm": round(so, 1),
        "even_ppm": round(me, 1),
        "even_err_ppm": round(se, 1),
        "disparity_sigma": round(sig, 2),
        "eclipsing_binary_ruled_out": sig < 3.0,
    }
