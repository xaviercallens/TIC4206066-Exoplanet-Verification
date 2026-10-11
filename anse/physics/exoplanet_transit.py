#!/usr/bin/env python3
"""
High-Performance Exoplanet Transit & Keplerian Dynamics Engine.
Pure NumPy / SciPy implementation of Keplerian orbital mechanics,
analytical quadratic limb-darkening transit profiles, and fast BLS/TLS period search.
"""

from __future__ import annotations

import math
from typing import Any
import numpy as np


# Physical constants (CGS & SI)
G_SI = 6.67430e-11          # m^3 kg^-1 s^-2
SOLAR_MASS_KG = 1.98847e30  # kg
SOLAR_RADIUS_M = 6.957e8    # m
EARTH_RADIUS_M = 6.371e6    # m
DAY_SEC = 86400.0           # seconds


def kepler_semi_major_axis(m_star_solar: float, period_days: float) -> tuple[float, float]:
    """
    Computes semi-major axis 'a' in AU and meters from Kepler's 3rd law.
    Returns: (a_au, a_meters)
    """
    m_star_kg = m_star_solar * SOLAR_MASS_KG
    p_sec = period_days * DAY_SEC
    # a^3 = G * M * P^2 / (4 * pi^2)
    a_m = (G_SI * m_star_kg * (p_sec ** 2) / (4.0 * math.pi ** 2)) ** (1.0 / 3.0)
    a_au = a_m / 1.495978707e11
    return a_au, a_m


def calculate_a_over_r_star(m_star_solar: float, r_star_solar: float, period_days: float) -> float:
    """Computes normalized orbital distance a / R_*."""
    _, a_m = kepler_semi_major_axis(m_star_solar, period_days)
    r_star_m = r_star_solar * SOLAR_RADIUS_M
    return a_m / r_star_m


def phase_fold(times: np.ndarray, period: float, t0: float) -> np.ndarray:
    """
    Folds time array into orbital phase in range [-0.5, 0.5).
    Phase 0.0 corresponds to mid-transit (t0).
    """
    phase = ((times - t0) % period) / period
    phase = np.where(phase >= 0.5, phase - 1.0, phase)
    return phase


def mandel_agol_transit_model(
    times: np.ndarray,
    period: float,
    t0: float,
    r_p_over_r_star: float,
    a_over_r_star: float,
    inclination_deg: float = 90.0,
    u1: float = 0.35,
    u2: float = 0.25,
) -> np.ndarray:
    """
    Quadratic limb-darkened transit light curve. Fully-in-transit points use the
    limb-darkening intensity at the planet-centre separation (error ~1 ppm for a
    500 ppm transit); ingress/egress points integrate the intensity over the planet
    disk numerically. This is NOT the closed-form Mandel & Agol (2002) solution; its
    residual versus an independent brute-force integration is bounded in
    tests/astrophysics/test_audit_regressions.py. Out-of-transit flux is 1.0.
    """
    p = r_p_over_r_star
    inc_rad = math.radians(inclination_deg)
    cos_i = math.cos(inc_rad)

    phase = phase_fold(times, period, t0)
    theta = 2.0 * math.pi * phase

    # Projected separation on sky plane in units of R_*
    # z^2 = (a/R_*)^2 * (sin^2(theta) + cos^2(i)*cos^2(theta))
    z = a_over_r_star * np.sqrt(np.sin(theta) ** 2 + (cos_i * np.cos(theta)) ** 2)

    flux = np.ones_like(times, dtype=np.float64)

    # In-transit mask: only when projected center distance < 1 + p and planet in front
    in_transit_geometry = (z < (1.0 + p)) & (np.cos(theta) > 0.0)

    if not np.any(in_transit_geometry):
        return flux

    z_sub = z[in_transit_geometry]
    depth_base = p ** 2

    # Quadratic limb darkening normalization factor
    omega = 1.0 - u1 / 3.0 - u2 / 6.0

    # Approximate limb-darkening profile: I(mu) = 1 - u1*(1 - mu) - u2*(1 - mu)^2
    # mu = sqrt(1 - z_sub^2)
    mu_sq = np.maximum(0.0, 1.0 - np.minimum(1.0, z_sub ** 2))
    mu = np.sqrt(mu_sq)
    limb_factor = (1.0 - u1 * (1.0 - mu) - u2 * ((1.0 - mu) ** 2)) / omega

    # Exact fraction of the planet disk lying inside the stellar disk
    # (circle-circle lens area / pi p^2). Full transit: z <= 1 - p.
    overlap = np.ones_like(z_sub)
    partial_mask = z_sub > (1.0 - p)
    if np.any(partial_mask):
        zp = z_sub[partial_mask]
        k0 = np.arccos(np.clip((p ** 2 + zp ** 2 - 1.0) / (2.0 * p * zp), -1.0, 1.0))
        k1 = np.arccos(np.clip((1.0 - p ** 2 + zp ** 2) / (2.0 * zp), -1.0, 1.0))
        root = 0.5 * np.sqrt(np.maximum(4.0 * zp ** 2 - (1.0 + zp ** 2 - p ** 2) ** 2, 0.0))
        area = p ** 2 * k0 + k1 - root
        overlap[partial_mask] = np.clip(area / (math.pi * p ** 2), 0.0, 1.0)

    delta = depth_base * limb_factor * overlap

    # Ingress/egress: the intensity varies across the planet disk, so integrate it
    # numerically (midpoint rule in polar coordinates centred on the planet) instead of
    # evaluating it at the planet centre. Only these few cadences pay the cost.
    if np.any(partial_mask):
        zp = z_sub[partial_mask]
        n_r, n_t = 48, 128
        rr = (np.arange(n_r) + 0.5) / n_r * p
        th = (np.arange(n_t) + 0.5) / n_t * 2.0 * math.pi
        x = zp[:, None, None] + rr[None, :, None] * np.cos(th)[None, None, :]
        y = rr[None, :, None] * np.sin(th)[None, None, :]
        rho2 = x ** 2 + y ** 2
        m = np.sqrt(np.clip(1.0 - rho2, 0.0, 1.0))
        inten = (1.0 - u1 * (1.0 - m) - u2 * (1.0 - m) ** 2) * (rho2 < 1.0)
        d_area = (rr * (p / n_r))[None, :, None] * (2.0 * math.pi / n_t)
        blocked = np.sum(inten * d_area, axis=(1, 2))
        delta[partial_mask] = blocked / (math.pi * omega)

    flux[in_transit_geometry] = 1.0 - delta

    return flux


def transit_duration_t14(
    period_days: float, a_over_r_star: float, k: float, impact_b: float
) -> float:
    """
    Total transit duration T14 in hours (Winn 2010, eq. 14, circular orbit):
    T14 = P/pi * asin( sqrt((1+k)^2 - b^2) / (a/R* sin i) ),  cos i = b / (a/R*).
    """
    sin_i = math.sqrt(max(0.0, 1.0 - (impact_b / a_over_r_star) ** 2))
    arg = math.sqrt(max(0.0, (1.0 + k) ** 2 - impact_b ** 2)) / (a_over_r_star * sin_i)
    return period_days / math.pi * math.asin(min(1.0, arg)) * 24.0


def derive_stellar_density(period_days: float, a_over_r_star: float) -> float:
    """
    Derives mean stellar density rho_* in g/cm^3 from transit parameters (Kepler's 3rd law).
    rho_* = (3 * pi / (G * P^2)) * (a / R_*)^3
    """
    p_sec = period_days * DAY_SEC
    rho_si = (3.0 * math.pi / (G_SI * (p_sec ** 2))) * (a_over_r_star ** 3)
    # Convert kg/m^3 to g/cm^3 (1 kg/m^3 = 1e-3 g/cm^3)
    rho_cgs = rho_si * 1e-3
    return rho_cgs


def fast_box_least_squares(
    times: np.ndarray,
    fluxes: np.ndarray,
    period_grid: np.ndarray,
    duration_days: float = 2.0 / 24.0,
) -> dict[str, Any]:
    """
    Vectorized Box Least Squares (BLS) periodogram search.
    Returns best period, best power, best t0, and full power spectrum.
    """
    n_points = len(times)
    mean_flux = np.mean(fluxes)
    var_flux = np.var(fluxes) + 1e-12

    powers = np.zeros_like(period_grid)
    best_t0s = np.zeros_like(period_grid)
    best_depths = np.zeros_like(period_grid)

    # Test phase bins
    n_bins = 50
    phase_bins = np.linspace(-0.5, 0.5, n_bins + 1)
    bin_centers = 0.5 * (phase_bins[:-1] + phase_bins[1:])

    for i, p in enumerate(period_grid):
        phases = ((times - times[0]) % p) / p
        phases = np.where(phases >= 0.5, phases - 1.0, phases)

        # Ingress/egress box width in phase
        q = duration_days / p
        half_q = q / 2.0

        best_p_power = -1.0
        best_p_t0 = times[0]
        best_p_depth = 0.0

        for phi_center in bin_centers:
            in_box = np.abs(phases - phi_center) <= half_q
            n_in = np.sum(in_box)
            n_out = n_points - n_in

            if n_in < 3 or n_out < 10:
                continue

            flux_in = np.mean(fluxes[in_box])
            flux_out = np.mean(fluxes[~in_box])
            depth = flux_out - flux_in

            if depth <= 0:
                continue

            # BLS power statistic
            power = (depth ** 2) * (n_in * n_out) / (n_points * var_flux)
            if power > best_p_power:
                best_p_power = power
                best_p_t0 = times[0] + phi_center * p
                best_p_depth = depth

        powers[i] = max(0.0, best_p_power)
        best_t0s[i] = best_p_t0
        best_depths[i] = best_p_depth

    best_idx = int(np.argmax(powers))
    return {
        "best_period": float(period_grid[best_idx]),
        "best_power": float(powers[best_idx]),
        "best_t0": float(best_t0s[best_idx]),
        "best_depth_ppm": float(best_depths[best_idx] * 1e6),
        "period_grid": period_grid.tolist(),
        "powers": powers.tolist(),
    }


def leave_one_epoch_out_predict(
    train_sectors_transits: list[float],
    withheld_window: tuple[float, float],
    period: float,
    t0: float,
) -> list[float]:
    """
    Calculates predicted transit times falling inside the withheld window based on fitted ephemeris.
    """
    t_start, t_end = withheld_window
    predicted_times: list[float] = []

    # Find integer transit number range
    n_min = math.floor((t_start - t0) / period) - 1
    n_max = math.ceil((t_end - t0) / period) + 1

    for n in range(n_min, n_max + 1):
        t_pred = t0 + n * period
        if t_start <= t_pred <= t_end:
            predicted_times.append(t_pred)

    return predicted_times
