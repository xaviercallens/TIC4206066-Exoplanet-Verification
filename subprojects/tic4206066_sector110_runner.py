#!/usr/bin/env python3
"""
Pre-Registered TESS Sector 110 Blind Test Runner for TIC 4206066.

Implements the frozen decision logic pre-registered on Zenodo
(doi:10.5281/zenodo.23175179). Transit times are DERIVED from the published
Zenodo ephemerides (anse.astrophysics.ephemeris) with propagated timing errors;
nothing here is a hand-typed epoch.

IMPORTANT - what this module does and does not establish
  * The Sector 110 data do not exist yet. `execute_preregistered_runner` therefore
    performs (a) the prediction step and (b) a SELF-TEST of the decision logic on
    synthetic light curves: one with an injected signal (must be RECOVERED) and
    one flat (must NOT be RECOVERED). The receipt says so explicitly; a self-test
    pass is not evidence about the planet.
  * The observing-window dates below are NOT verified against the official TESS
    timeline (see WINDOW_VERIFIED).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import hashlib
import json
import os
import time
from typing import Any

import numpy as np

from anse.astrophysics.ephemeris import BTJD_OFFSET, Ephemeris, load_ephemeris_a, load_ephemeris_b
from anse.astrophysics.elenchus_vetting import run_elenchus_adversarial_battery
from anse.physics.exoplanet_transit import (
    calculate_a_over_r_star,
    mandel_agol_transit_model,
)

SECTOR_110_START_BJD = 2461344.5  # 2026-10-30 23:58 UTC (BJD_TDB), verified conversion
SECTOR_110_END_BJD = 2461370.5    # 2026-11-25 23:58 UTC
WINDOW_VERIFIED = False           # dates not cross-checked against the TESS mission timeline

# Frozen decision thresholds (fractions of the published combined depth)
RECOVERY_FRACTION = 0.6
EXCLUSION_FRACTION = 0.2
MIN_COVERED_TRANSITS = 3
MIN_POINTS_PER_TRANSIT = 10
HALF_DURATION_D = 1.87 / 24.0 / 2.0  # published T14 = 1.87 h

METADATA_PATH = "data/astrophysics/tic4206066/metadata.json"
RECEIPT_PATH = "results/tic4206066_execution_receipt.json"


@dataclass
class Sector110PredictionReceipt:
    target_id: str
    host_star: str
    zenodo_prereg_doi: str
    sector_window: dict[str, Any]
    signal_a_predictions: dict[str, Any]
    signal_b_predictions: dict[str, Any]
    decision_rules: dict[str, str]
    overall_status: str
    provenance_hash: str
    timestamp_utc: str
    self_test: dict[str, Any] | None = None


def _load_meta() -> dict[str, Any]:
    with open(METADATA_PATH) as f:
        return json.load(f)


def _predict(eph: Ephemeris, nominal_depth_ppm: float) -> dict[str, Any]:
    epochs = eph.epochs_in_window(SECTOR_110_START_BJD, SECTOR_110_END_BJD)
    return {
        "period_days": eph.period_d,
        "ephemeris_source": eph.source,
        "t0_bjd_tdb": eph.t0_bjd,
        "expected_transits_count": len(epochs),
        "epochs": epochs,
        "predicted_mid_transit_bjd": [round(eph.transit_time(n), 4) for n in epochs],
        "sigma_min_per_transit": [round(eph.sigma_minutes(n), 2) for n in epochs],
        "nominal_depth_ppm": nominal_depth_ppm,
    }


def find_sector110_fits(mast_dir: str = "data/astrophysics/tic4206066/mast") -> str | None:
    """Return the path of a locally available Sector 110 light-curve FITS file, or None."""
    if not os.path.isdir(mast_dir):
        return None
    for name in sorted(os.listdir(mast_dir)):
        low = name.lower()
        if low.endswith(".fits") and ("s0110" in low or "-s110" in low):
            return os.path.join(mast_dir, name)
    return None


def load_sector110_fits(path: str) -> dict[str, np.ndarray]:
    """Load a real SPOC/QLP light curve; times are returned in BJD_TDB (BTJD + 2457000)."""
    from astropy.io import fits

    with fits.open(path) as hdul:
        d = hdul[1].data
        names = d.columns.names
        flux = d["PDCSAP_FLUX"] if "PDCSAP_FLUX" in names else (
            d["DET_FLUX"] if "DET_FLUX" in names else d["SAP_FLUX"]
        )
        t = d["TIME"]
        ok = ~np.isnan(t) & ~np.isnan(flux) & (flux > 0)
        if "QUALITY" in names:
            ok &= d["QUALITY"] == 0
        f = flux[ok].astype(np.float64)
        px = d["MOM_CENTR1"][ok].astype(np.float64) if "MOM_CENTR1" in names else None
        py = d["MOM_CENTR2"][ok].astype(np.float64) if "MOM_CENTR2" in names else None
        return {
            "times": t[ok].astype(np.float64) + BTJD_OFFSET,
            "fluxes": f / np.median(f),
            "pixel_x": px,
            "pixel_y": py,
        }


def generate_sector_110_predictions() -> dict[str, Any]:
    """Frozen transit windows for Sector 110, derived from the published ephemerides."""
    meta = _load_meta()
    pa = _predict(load_ephemeris_a(), meta["signal_a"]["depth_ppm"])
    pb = _predict(load_ephemeris_b(), meta["signal_b"]["depth_ppm"])
    # Signal A's quoted uncertainty is the propagated 1-sigma at the middle epoch; the
    # catalogue O-C scatter (rms ~9 min) is larger and is reported alongside it.
    pa["uncertainty_mins"] = round(float(np.median(pa["sigma_min_per_transit"])), 2)
    pa["catalogue_rms_oc_min"] = load_ephemeris_a().rms_oc_min
    return {"signal_a": pa, "signal_b": pb}


def evaluate_sector_110_data(
    times: np.ndarray,
    fluxes: np.ndarray,
    pixel_x: np.ndarray | None = None,
    pixel_y: np.ndarray | None = None,
) -> dict[str, Any]:
    """
    Evaluate Sector 110 photometry (times in BJD_TDB, flux normalised to 1) against
    the frozen rules for Signal A:
    - 'RECOVERED'  : depth >= 60 % of the published depth AND the Elenchus battery does not refute it.
    - 'EXCLUDED'   : depth < 20 % of the published depth (flat at predicted epochs).
    - 'INCONSISTENT': in between, or refuted by the Elenchus battery.
    - 'INSUFFICIENT_DATA': fewer than 3 predicted transits adequately covered.
    """
    meta = _load_meta()
    star = meta["target_star"]
    eph = load_ephemeris_a()
    nominal = float(meta["signal_a"]["depth_ppm"])
    pred = generate_sector_110_predictions()["signal_a"]

    in_flux, covered = [], 0
    for t_pred in pred["predicted_mid_transit_bjd"]:
        win = np.abs(times - t_pred) <= HALF_DURATION_D * 0.8
        if int(win.sum()) >= MIN_POINTS_PER_TRANSIT:
            covered += 1
            in_flux.extend(fluxes[win].tolist())

    depth_ppm = None
    if covered < MIN_COVERED_TRANSITS:
        outcome = "INSUFFICIENT_DATA"
    else:
        depth_ppm = (1.0 - float(np.mean(in_flux))) * 1e6
        if depth_ppm >= RECOVERY_FRACTION * nominal:
            a_over_r = calculate_a_over_r_star(star["mass_solar"], star["radius_solar"], eph.period_d)
            vetting = run_elenchus_adversarial_battery(
                times=times,
                fluxes=fluxes,
                period=eph.period_d,
                t0=eph.t0_bjd,
                a_over_r_star=a_over_r,
                host_density_cgs=star["mean_density_cgs"],
                pixel_x=pixel_x,
                pixel_y=pixel_y,
            )
            outcome = "RECOVERED" if vetting.is_planet_candidate else vetting.falsification_verdict
        elif depth_ppm < EXCLUSION_FRACTION * nominal:
            outcome = "EXCLUDED"
        else:
            outcome = "INCONSISTENT"

    return {
        "signal_a_outcome": outcome,
        "covered_transits_count": covered,
        "total_predicted_transits": len(pred["predicted_mid_transit_bjd"]),
        "mean_in_transit_depth_ppm": None if depth_ppm is None else round(depth_ppm, 1),
    }


def _self_test() -> dict[str, Any]:
    """Run the decision logic on synthetic data with a known answer (positive AND negative)."""
    meta = _load_meta()
    star, sig = meta["target_star"], meta["signal_a"]
    eph = load_ephemeris_a()
    rng = np.random.default_rng(42)
    t = np.linspace(SECTOR_110_START_BJD, SECTOR_110_END_BJD, 18000)  # 2-min cadence
    a_over_r = calculate_a_over_r_star(star["mass_solar"], star["radius_solar"], eph.period_d)
    k = float(np.sqrt(sig["depth_ppm"] * 1e-6))
    noise = rng.normal(1.0, 0.00015, size=t.size)
    model = mandel_agol_transit_model(t, eph.period_d, eph.t0_bjd, k, a_over_r)
    injected = evaluate_sector_110_data(t, noise + (model - 1.0))
    flat = evaluate_sector_110_data(t, rng.normal(1.0, 0.00015, size=t.size))
    return {
        "injected_signal_outcome": injected["signal_a_outcome"],
        "flat_baseline_outcome": flat["signal_a_outcome"],
        "injected_depth_ppm": injected["mean_in_transit_depth_ppm"],
        "flat_depth_ppm": flat["mean_in_transit_depth_ppm"],
        "logic_ok": injected["signal_a_outcome"] == "RECOVERED"
        and flat["signal_a_outcome"] in ("EXCLUDED", "INCONSISTENT"),
        "data_kind": "SYNTHETIC (self-test of decision logic; not Sector 110 data)",
    }


def execute_preregistered_runner() -> Sector110PredictionReceipt:
    """Generate the frozen predictions and self-test the decision logic."""
    preds = generate_sector_110_predictions()
    self_test = _self_test()

    payload = json.dumps(preds, sort_keys=True).encode("utf-8")
    provenance = hashlib.sha256(payload).hexdigest()

    status = (
        "AWAITING_SECTOR_110_DATA__UNREGISTERED_SELF_TEST_PASSED"
        if self_test["logic_ok"]
        else "AWAITING_SECTOR_110_DATA__UNREGISTERED_SELF_TEST_FAILED"
    )
    receipt = Sector110PredictionReceipt(
        target_id="TIC 4206066",
        host_star="StKM 1-561 (Gaia DR3 3220388198192519424)",
        zenodo_prereg_doi="10.5281/zenodo.23175179",
        sector_window={
            "sector": 110,
            "start_bjd": SECTOR_110_START_BJD,
            "end_bjd": SECTOR_110_END_BJD,
            "dates_utc": "2026-10-30 to 2026-11-25",
            "window_verified_against_tess_timeline": WINDOW_VERIFIED,
        },
        signal_a_predictions=preds["signal_a"],
        signal_b_predictions=preds["signal_b"],
        decision_rules={
            "RECOVERED": f"Mean in-transit depth >= {RECOVERY_FRACTION:.0%} of the published depth and Elenchus gates clean.",
            "INCONSISTENT": "Depth between the exclusion and recovery thresholds, or refuted by a gate.",
            "EXCLUDED": f"Mean in-transit depth < {EXCLUSION_FRACTION:.0%} of the published depth.",
            "INSUFFICIENT_DATA": f"Fewer than {MIN_COVERED_TRANSITS} predicted transit windows adequately covered.",
        },
        overall_status=status,
        provenance_hash=provenance,
        timestamp_utc=time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        self_test=self_test,
    )

    os.makedirs("results", exist_ok=True)
    with open(RECEIPT_PATH, "w") as f:
        json.dump(asdict(receipt), f, indent=2)
    return receipt


if __name__ == "__main__":
    r = execute_preregistered_runner()
    print("=== SECTOR 110 PRE-REGISTERED RUNNER (predictions + synthetic self-test) ===")
    print(f"Target: {r.target_id} ({r.host_star})")
    print(f"Signal A ({r.signal_a_predictions['period_days']:.6f} d): "
          f"{r.signal_a_predictions['expected_transits_count']} predicted transits, "
          f"+-{r.signal_a_predictions['uncertainty_mins']} min")
    print(f"Signal B ({r.signal_b_predictions['period_days']:.5f} d, tentative): "
          f"{r.signal_b_predictions['expected_transits_count']} predicted transits")
    print(f"Overall Status: {r.overall_status}")
    print(f"Provenance SHA-256: {r.provenance_hash[:16]}...")
    print(f"Execution receipt written to {RECEIPT_PATH}")
