#!/usr/bin/env python3
"""
Socratic Elenchus Falsification Battery for Exoplanet Transit Vetting.
Applies rigorous cross-examination to candidate transit signals:
- Mandatory Positive Control (synthetic planet injection & recovery)
- Negative Control 1: Odd-Even Transit Depth Parity (Eclipsing Binary check)
- Negative Control 2: Secondary Eclipse Search at Phase 0.5 (Occultation check)
- Negative Control 3: Centroid Pixel Motion Decorrelation (Background EB check)
- Negative Control 4: Stellar Density Consistency (Keplerian physical prior)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from anse.physics.exoplanet_transit import (
    derive_stellar_density,
    mandel_agol_transit_model,
    phase_fold,
)


@dataclass
class VettingResult:
    is_planet_candidate: bool
    energy: float
    falsification_verdict: str  # "NOT_REFUTED", "REFUTED_ODD_EVEN", "REFUTED_SECONDARY", etc.
    odd_even_sigma: float
    secondary_eclipse_depth_ppm: float
    centroid_offset_sigma: float
    stellar_density_ratio: float
    positive_control_passed: bool
    negative_controls_rejected: bool
    diagnostic_telemetry: dict[str, Any] = field(default_factory=dict)


def evaluate_odd_even_parity(
    times: np.ndarray,
    fluxes: np.ndarray,
    period: float,
    t0: float,
    duration_days: float = 2.0 / 24.0,
) -> tuple[float, float, float]:
    """
    Evaluates depth difference between odd and even transits.
    Returns: (depth_odd_ppm, depth_even_ppm, delta_sigma)
    """
    epoch_numbers = np.round((times - t0) / period).astype(int)
    phase = phase_fold(times, period, t0)
    in_transit = np.abs(phase * period) <= (duration_days / 2.0)
    out_of_transit = np.abs(phase * period) > (duration_days * 1.5)

    base_flux = np.median(fluxes[out_of_transit]) if np.any(out_of_transit) else 1.0

    even_mask = in_transit & (epoch_numbers % 2 == 0)
    odd_mask = in_transit & (epoch_numbers % 2 == 1)

    if np.sum(even_mask) < 3 or np.sum(odd_mask) < 3:
        return 0.0, 0.0, 0.0

    flux_even = fluxes[even_mask]
    flux_odd = fluxes[odd_mask]

    depth_even = (base_flux - np.mean(flux_even)) * 1e6
    depth_odd = (base_flux - np.mean(flux_odd)) * 1e6

    sem_even = (np.std(flux_even) / math.sqrt(len(flux_even))) * 1e6 + 1e-6
    sem_odd = (np.std(flux_odd) / math.sqrt(len(flux_odd))) * 1e6 + 1e-6

    delta = abs(depth_even - depth_odd)
    combined_sigma = math.sqrt(sem_even ** 2 + sem_odd ** 2)
    delta_sigma = delta / combined_sigma

    return float(depth_odd), float(depth_even), float(delta_sigma)


def search_secondary_eclipse(
    times: np.ndarray,
    fluxes: np.ndarray,
    period: float,
    t0: float,
    duration_days: float = 2.0 / 24.0,
) -> tuple[float, float]:
    """
    Searches for an inverted flux dip at orbital phase phi = 0.50 (secondary occultation).
    Returns: (secondary_depth_ppm, secondary_sigma)
    """
    phase = phase_fold(times, period, t0)
    # Phase around 0.5 (or -0.5)
    sec_mask = np.abs(np.abs(phase) - 0.5) <= (duration_days / (2.0 * period))
    out_mask = (np.abs(phase * period) > (duration_days * 1.5)) & (
        np.abs(np.abs(phase) - 0.5) > (duration_days / period)
    )

    if np.sum(sec_mask) < 3 or np.sum(out_mask) < 10:
        return 0.0, 0.0

    base_flux = np.median(fluxes[out_mask])
    sec_flux = fluxes[sec_mask]

    sec_depth = (base_flux - np.mean(sec_flux)) * 1e6
    sec_sem = (np.std(sec_flux) / math.sqrt(len(sec_flux))) * 1e6 + 1e-6
    sec_sigma = max(0.0, sec_depth / sec_sem)

    return float(sec_depth), float(sec_sigma)


def evaluate_centroid_decorrelation(
    pixel_x: np.ndarray,
    pixel_y: np.ndarray,
    in_transit_mask: np.ndarray,
) -> float:
    """
    Measures center-of-light spatial shift (x, y) during in-transit vs out-of-transit frames.
    Returns: centroid_offset_sigma (standard deviations of displacement).
    """
    if np.sum(in_transit_mask) < 3 or np.sum(~in_transit_mask) < 10:
        return 0.0

    x_in = pixel_x[in_transit_mask]
    y_in = pixel_y[in_transit_mask]
    x_out = pixel_x[~in_transit_mask]
    y_out = pixel_y[~in_transit_mask]

    mean_x_in, mean_y_in = np.mean(x_in), np.mean(y_in)
    mean_x_out, mean_y_out = np.mean(x_out), np.mean(y_out)

    dx = mean_x_in - mean_x_out
    dy = mean_y_in - mean_y_out
    dist = math.sqrt(dx ** 2 + dy ** 2)

    var_x = np.var(x_out) / len(x_out) + np.var(x_in) / len(x_in)
    var_y = np.var(y_out) / len(y_out) + np.var(y_in) / len(y_in)
    sigma_dist = math.sqrt(var_x + var_y) + 1e-9

    return float(dist / sigma_dist)


def run_elenchus_adversarial_battery(
    times: np.ndarray,
    fluxes: np.ndarray,
    period: float,
    t0: float,
    a_over_r_star: float,
    host_density_cgs: float,
    pixel_x: np.ndarray | None = None,
    pixel_y: np.ndarray | None = None,
    duration_days: float = 2.0 / 24.0,
) -> VettingResult:
    """
    Runs the complete Elenchus falsification battery against a candidate transit signal.
    Calculates objective physical Energy E:
    E = 0.0 - 5.0 (good candidate)
    E >= 1000.0 (refuted / false positive / physical invariant breach)
    """
    # 1. Odd-Even Parity Check
    depth_odd, depth_even, odd_even_sigma = evaluate_odd_even_parity(
        times, fluxes, period, t0, duration_days
    )

    # 2. Secondary Eclipse Search
    sec_depth, sec_sigma = search_secondary_eclipse(
        times, fluxes, period, t0, duration_days
    )

    # 3. Centroid Decorrelation Check
    phase = phase_fold(times, period, t0)
    in_transit = np.abs(phase * period) <= (duration_days / 2.0)

    centroid_checked = pixel_x is not None and pixel_y is not None
    if centroid_checked:
        centroid_sigma = evaluate_centroid_decorrelation(pixel_x, pixel_y, in_transit)
    else:
        centroid_sigma = 0.0  # gate NOT evaluated; flagged in diagnostic_telemetry

    # 4. Stellar Density Consistency
    derived_rho = derive_stellar_density(period, a_over_r_star)
    density_ratio = derived_rho / host_density_cgs if host_density_cgs > 0 else 1.0
    log_ratio = abs(math.log10(max(1e-4, density_ratio)))

    # Compute base Chi-Square on out-of-transit baseline
    out_of_transit = ~in_transit
    baseline_noise = np.std(fluxes[out_of_transit]) if np.any(out_of_transit) else 0.001
    chi2_red = float(np.mean(((fluxes - 1.0) / baseline_noise) ** 2))

    energy = chi2_red
    falsification_verdict = "NOT_REFUTED"
    is_candidate = True

    # Elenchus Hard Rejection Gates
    if odd_even_sigma >= 3.0:
        energy += 1000.0
        falsification_verdict = "REFUTED_ODD_EVEN_PARITY"
        is_candidate = False

    elif sec_sigma >= 3.5 and sec_depth >= 150.0:
        energy += 1000.0
        falsification_verdict = "REFUTED_SECONDARY_ECLIPSE"
        is_candidate = False

    elif centroid_sigma >= 3.0:
        energy += 1000.0
        falsification_verdict = "REFUTED_CENTROID_SHIFT"
        is_candidate = False

    elif log_ratio >= 0.6:  # Factor of 4 discrepancy in host density
        energy += 1000.0
        falsification_verdict = "REFUTED_STELLAR_DENSITY_INCONSISTENCY"
        is_candidate = False

    # 5. Mandatory Controls (local seeded RNG: no global-state side effects)
    rng = np.random.default_rng(20260101)
    inj_p, inj_t0, inj_rp, inj_ar = 3.18, 0.5, 0.022, 11.5
    syn_times = np.linspace(0.0, 10.0, 1000)
    syn_noise = rng.normal(1.0, 0.0001, size=len(syn_times))
    syn_model = mandel_agol_transit_model(syn_times, inj_p, inj_t0, inj_rp, inj_ar)
    syn_flux = syn_noise + (syn_model - 1.0)

    # Positive control: the injected planet must be RECOVERED (depth within 30 % of
    # truth at >5 sigma) and must not be flagged as an odd/even false positive.
    syn_phase = phase_fold(syn_times, inj_p, inj_t0)
    syn_in = np.abs(syn_phase * inj_p) <= (duration_days / 2.0) * 0.8
    syn_out = np.abs(syn_phase * inj_p) > (duration_days * 1.5)
    true_depth = float(np.mean(1.0 - syn_model[syn_in]))
    rec_depth = float(np.median(syn_flux[syn_out]) - np.mean(syn_flux[syn_in]))
    rec_err = float(np.std(syn_flux[syn_out]) / math.sqrt(max(1, int(syn_in.sum()))))
    _, _, syn_oe_sig = evaluate_odd_even_parity(syn_times, syn_flux, inj_p, inj_t0)
    depth_recovered = abs(rec_depth - true_depth) / true_depth < 0.30
    significant = rec_err > 0 and rec_depth / rec_err > 5.0
    positive_control_passed = bool(depth_recovered and significant and syn_oe_sig < 2.5)

    # Negative control: an injected eclipsing binary (alternate transits twice as deep)
    # must be REFUTED by the odd/even gate.
    eb_flux = np.copy(syn_noise)
    for i, t in enumerate(syn_times):
        epoch = round((t - inj_t0) / inj_p)
        d_factor = 2.0 if (epoch % 2 == 1) else 1.0
        eb_flux[i] += (syn_model[i] - 1.0) * d_factor
    _, _, eb_oe_sig = evaluate_odd_even_parity(syn_times, eb_flux, inj_p, inj_t0)
    negative_control_rejected = bool(eb_oe_sig >= 3.0)

    return VettingResult(
        is_planet_candidate=is_candidate,
        energy=round(energy, 4),
        falsification_verdict=falsification_verdict,
        odd_even_sigma=round(odd_even_sigma, 3),
        secondary_eclipse_depth_ppm=round(sec_depth, 1),
        centroid_offset_sigma=round(centroid_sigma, 3),
        stellar_density_ratio=round(density_ratio, 3),
        positive_control_passed=positive_control_passed,
        negative_controls_rejected=negative_control_rejected,
        diagnostic_telemetry={
            "depth_odd_ppm": round(depth_odd, 1),
            "depth_even_ppm": round(depth_even, 1),
            "derived_rho_cgs": round(derived_rho, 3),
            "host_density_cgs": round(host_density_cgs, 3),
            "chi2_red": round(chi2_red, 4),
            "centroid_checked": bool(centroid_checked),
            "positive_control_recovered_depth_ppm": round(rec_depth * 1e6, 1),
            "positive_control_true_depth_ppm": round(true_depth * 1e6, 1),
            "negative_control_odd_even_sigma": round(float(eb_oe_sig), 2),
        },
    )
