#!/usr/bin/env python3
"""
Single source of truth for the TIC 4206066 ephemerides.

Values are read from the Zenodo deposit (doi:10.5281/zenodo.22967099), never
hardcoded, so that predictions, vetting windows and the web interface cannot
drift apart (an earlier revision carried a reference epoch 0.70 d away from
the published fit, which put every "in-transit" window out of phase).

Conventions: times are BJD_TDB; BTJD = BJD_TDB - 2457000.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass

BTJD_OFFSET = 2457000.0
ZENODO_DIR = "data/astrophysics/tic4206066/zenodo"
EPHEMERIS_A_PATH = "part1/02_transits_A/ephemeris_A.json"
EPHEMERIS_B_PATH = "part2/07_signal_B/ephemeris_B.json"


@dataclass(frozen=True)
class Ephemeris:
    """Linear ephemeris T(n) = t0_bjd + n * period_d with 1-sigma errors."""

    name: str
    period_d: float
    period_err_d: float
    t0_bjd: float
    t0_err_d: float
    n_events: int
    rms_oc_min: float | None
    source: str

    @property
    def t0_btjd(self) -> float:
        return self.t0_bjd - BTJD_OFFSET

    def transit_time(self, epoch: int) -> float:
        return self.t0_bjd + epoch * self.period_d

    def sigma_minutes(self, epoch: int) -> float:
        """Propagated 1-sigma timing error at `epoch` (T0 and P assumed independent)."""
        return math.hypot(self.t0_err_d, epoch * self.period_err_d) * 1440.0

    def epochs_in_window(self, start_bjd: float, end_bjd: float) -> list[int]:
        n_lo = math.ceil((start_bjd - self.t0_bjd) / self.period_d)
        n_hi = math.floor((end_bjd - self.t0_bjd) / self.period_d)
        return list(range(n_lo, n_hi + 1))


def _find_zenodo_dir(preferred: str = ZENODO_DIR) -> str:
    for candidate in [
        preferred,
        "data/zenodo",
        os.path.join(os.path.dirname(__file__), "..", "..", "data", "zenodo"),
        os.path.join(os.path.dirname(__file__), "..", "..", "data", "astrophysics", "tic4206066", "zenodo"),
    ]:
        if os.path.exists(os.path.join(candidate, EPHEMERIS_A_PATH)):
            return candidate
    return preferred


def _read(rel: str, zenodo_dir: str) -> dict:
    zenodo_dir = _find_zenodo_dir(zenodo_dir)
    with open(os.path.join(zenodo_dir, rel)) as f:
        return json.load(f)


def load_ephemeris_a(zenodo_dir: str = ZENODO_DIR) -> Ephemeris:
    zenodo_dir = _find_zenodo_dir(zenodo_dir)
    fit = _read(EPHEMERIS_A_PATH, zenodo_dir)["fit"]
    return Ephemeris(
        name="A",
        period_d=float(fit["P_d"]),
        period_err_d=float(fit["eP_d"]),
        t0_bjd=float(fit["T0_BJD_TDB"]),
        t0_err_d=float(fit["eT0_d"]),
        n_events=int(fit["n"]),
        rms_oc_min=float(fit["rms_oc_min"]),
        source=f"zenodo:{EPHEMERIS_A_PATH}",
    )


def load_ephemeris_b(zenodo_dir: str = ZENODO_DIR) -> Ephemeris:
    """Signal B: tentative (~3 sigma post hoc). Fit of the 4 archival events."""
    zenodo_dir = _find_zenodo_dir(zenodo_dir)
    fit = _read(EPHEMERIS_B_PATH, zenodo_dir)["fits_4_events_manuscript_like"]
    return Ephemeris(
        name="B",
        period_d=float(fit["P"]),
        period_err_d=float(fit["eP"]),
        t0_bjd=float(fit["T0"]) + BTJD_OFFSET,
        t0_err_d=float(fit["eT0"]),
        n_events=int(fit["n"]),
        rms_oc_min=None,
        source=f"zenodo:{EPHEMERIS_B_PATH}",
    )

