"""Sector 110 predictions must be derived from the ephemeris, never hand-typed."""
import numpy as np

from anse.astrophysics.ephemeris import load_ephemeris_a, load_ephemeris_b
from subprojects.tic4206066_sector110_runner import (
    SECTOR_110_END_BJD,
    SECTOR_110_START_BJD,
    _predict,
)


def _check(eph):
    pred = _predict(eph, 500.0)
    t = np.array(pred["predicted_mid_transit_bjd"])
    assert pred["expected_transits_count"] == len(t) >= 2
    assert t.min() >= SECTOR_110_START_BJD - 1e-3 and t.max() <= SECTOR_110_END_BJD + 1e-3
    assert np.allclose(np.diff(t), eph.period_d, atol=1e-3)  # rounding to 1e-4 d
    for n, ti in zip(pred["epochs"], t):
        assert abs(ti - eph.transit_time(n)) < 1e-3
    # timing uncertainty grows with epoch (propagated, not constant)
    assert pred["sigma_min_per_transit"][-1] >= pred["sigma_min_per_transit"][0]


def test_signal_a_predictions_reproduce():
    _check(load_ephemeris_a())


def test_signal_b_predictions_reproduce():
    _check(load_ephemeris_b())
