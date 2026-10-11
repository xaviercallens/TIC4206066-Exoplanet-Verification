#!/usr/bin/env python3
"""
Frontier Domain: The Autonomous Astrophysicist.
Evaluates empirical exoplanetary transit candidates against physical invariants,
Keplerian orbital conservation, and Socratic Elenchus falsification gates.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from anse.astrophysics.elenchus_vetting import VettingResult, run_elenchus_adversarial_battery


@dataclass
class ExoplanetDomainResult:
    target_id: str
    period_days: float
    energy: float
    is_valid_candidate: bool
    verdict: str
    duration_ms: float
    controls_passed: bool
    telemetry: dict[str, Any]


class AutonomousAstrophysicist:
    """
    Evaluates exoplanet candidates against physical laws and falsification batteries.
    Energy E = 0.0 - 5.0 on robust candidates; E >= 1000.0 on falsification or unphysical states.
    """

    def evaluate_transit_candidate(
        self,
        target_id: str,
        times: np.ndarray,
        fluxes: np.ndarray,
        period: float,
        t0: float,
        a_over_r_star: float,
        host_density_cgs: float,
        pixel_x: np.ndarray | None = None,
        pixel_y: np.ndarray | None = None,
    ) -> ExoplanetDomainResult:
        start_t = time.perf_counter()

        # Execute Elenchus adversarial falsification battery
        vetting: VettingResult = run_elenchus_adversarial_battery(
            times=times,
            fluxes=fluxes,
            period=period,
            t0=t0,
            a_over_r_star=a_over_r_star,
            host_density_cgs=host_density_cgs,
            pixel_x=pixel_x,
            pixel_y=pixel_y,
        )

        dur_ms = (time.perf_counter() - start_t) * 1000.0
        controls_ok = vetting.positive_control_passed and vetting.negative_controls_rejected

        # If controls failed, the test suite itself is an unverified wind-egg
        if not controls_ok:
            energy = max(vetting.energy, 1000.0)
            is_valid = False
            verdict = "REFUTED_INVALID_CONTROLS"
        else:
            energy = vetting.energy
            is_valid = vetting.is_planet_candidate
            verdict = vetting.falsification_verdict

        return ExoplanetDomainResult(
            target_id=target_id,
            period_days=period,
            energy=energy,
            is_valid_candidate=is_valid,
            verdict=verdict,
            duration_ms=round(dur_ms, 2),
            controls_passed=controls_ok,
            telemetry={
                "odd_even_sigma": vetting.odd_even_sigma,
                "secondary_eclipse_ppm": vetting.secondary_eclipse_depth_ppm,
                "centroid_offset_sigma": vetting.centroid_offset_sigma,
                "stellar_density_ratio": vetting.stellar_density_ratio,
                **vetting.diagnostic_telemetry,
            },
        )
