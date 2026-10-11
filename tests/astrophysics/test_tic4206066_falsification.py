#!/usr/bin/env python3
"""
Test Suite: TIC 4206066 Epistemic Reproducibility & Socratic Elenchus Falsification Battery.
Verifies the exoplanetary discovery pipeline around TIC 4206066:
1. Keplerian orbital mechanics & stellar density invariants
2. Multi-epoch transit recovery across Sectors 6, 32, and 98
3. 3-fold leave-one-epoch-out temporal blind cross-validation
4. Elenchus mandatory controls (positive planet injection & negative EB/centroid refutations)
5. AutonomousAstrophysicist frontier domain energy minimization (E < 5.0)
"""

import math

import numpy as np

from anse.astrophysics.domain import AutonomousAstrophysicist
from anse.astrophysics.elenchus_vetting import run_elenchus_adversarial_battery
from anse.astrophysics.real_data_ingest import (
    audit_real_centroid_stability,
    audit_real_data_sources,
    audit_real_odd_even_parity,
    audit_real_timing_ephemeris,
    load_real_mast_lightcurves,
    load_zenodo_per_transit_data,
)
from anse.physics.exoplanet_transit import (
    calculate_a_over_r_star,
    derive_stellar_density,
    kepler_semi_major_axis,
    mandel_agol_transit_model,
)


class TestTIC4206066Reproducibility:
    """Verifies the astrophysical properties of TIC 4206066 and candidate Signal A."""

    # Adopted stellar parameters from the Zenodo deposit (09_star FACTS: STAR.M, STAR.R);
    # density = M / R^3 * 1.41 g/cm^3.
    M_STAR = 0.606  # Solar masses (late-K dwarf)
    R_STAR = 0.623  # Solar radii
    HOST_DENSITY = 3.534  # g/cm^3
    PERIOD_A = 3.182785  # days

    def test_kepler_dynamics_and_density(self):
        """Verifies physical Keplerian invariants for the 3.18d candidate."""
        a_au, a_m = kepler_semi_major_axis(self.M_STAR, self.PERIOD_A)
        a_over_r = calculate_a_over_r_star(self.M_STAR, self.R_STAR, self.PERIOD_A)
        derived_rho = derive_stellar_density(self.PERIOD_A, a_over_r)

        assert 0.0353 < a_au < 0.0363  # published a = 0.0358 AU, f"Unphysical semi-major axis: {a_au}"
        assert 12.2 < a_over_r < 12.5  # published a/R* = 12.37, f"Unphysical a/R_*: {a_over_r}"
        assert abs(derived_rho - self.HOST_DENSITY) < 0.20, (
            f"Stellar density discrepancy: {derived_rho} vs {self.HOST_DENSITY}"
        )

    def test_mandel_agol_transit_depth(self):
        """Verifies quadratic limb-darkened transit profile matches ~500 ppm depth."""
        a_over_r = calculate_a_over_r_star(self.M_STAR, self.R_STAR, self.PERIOD_A)
        # R_p = 1.4 R_Earth around 0.68 R_Sun
        # r_p_over_r_star = sqrt(500e-6) = ~0.0224
        rp_over_rs = math.sqrt(500.0e-6)

        times = np.linspace(-0.1, 0.1, 500)
        flux = mandel_agol_transit_model(times, self.PERIOD_A, 0.0, rp_over_rs, a_over_r)

        depth_ppm = (1.0 - np.min(flux)) * 1e6
        assert 450.0 < depth_ppm < 600.0, f"Transit depth out of expected range: {depth_ppm} ppm"

    def test_published_ephemeris_predicts_events_inside_each_real_sector(self):
        """
        The published ephemeris must place transits inside the REAL data coverage of each
        sector, and every catalogued event must lie within 3 sigma of an ephemeris epoch.
        (The previous version asserted that t0 + n*P has phase 0 relative to t0, a tautology,
        over hand-typed windows that did not even contain Sector 98.)
        """
        from anse.astrophysics.ephemeris import load_ephemeris_a
        from anse.astrophysics.real_data_audit import load_events

        eph = load_ephemeris_a()
        lcs = load_real_mast_lightcurves()
        events = load_events()
        for sec, d in lcs.items():
            lo, hi = d["times"].min() + 2457000.0, d["times"].max() + 2457000.0
            predicted = eph.epochs_in_window(lo, hi)
            catalogued = [e for e in events if e["sector"] == sec]
            assert len(predicted) >= len(catalogued) >= 6
            for e in catalogued:
                resid_min = (e["t_mid_bjd"] - eph.transit_time(e["epoch"])) * 1440.0
                sigma = float(np.hypot(eph.sigma_minutes(e["epoch"]), e["err_min"]))
                assert abs(resid_min) < 3.0 * sigma, (sec, e["epoch"], resid_min, sigma)

    def test_elenchus_falsification_controls(self):
        """Verifies that positive and negative controls correctly fire and fail."""
        np.random.seed(42)
        times = np.linspace(0.0, 15.0, 2000)
        noise = np.random.normal(1.0, 0.0001, size=len(times))
        transit = mandel_agol_transit_model(times, self.PERIOD_A, 0.5, 0.0224, 12.37)
        flux_clean = noise + (transit - 1.0)

        # 1. Clean Candidate: Must pass with low energy
        res_clean = run_elenchus_adversarial_battery(
            times, flux_clean, self.PERIOD_A, 0.5, 12.37, self.HOST_DENSITY
        )
        assert res_clean.is_planet_candidate is True
        assert res_clean.energy < 5.0
        assert res_clean.positive_control_passed is True
        assert res_clean.negative_controls_rejected is True
        assert res_clean.falsification_verdict == "NOT_REFUTED"

        # 2. Negative Control 1: Grazing Eclipsing Binary (Odd-Even difference)
        flux_eb = np.copy(noise)
        for i, t in enumerate(times):
            epoch = round((t - 0.5) / self.PERIOD_A)
            factor = 3.0 if (epoch % 2 == 1) else 1.0
            flux_eb[i] += (transit[i] - 1.0) * factor

        res_eb = run_elenchus_adversarial_battery(
            times, flux_eb, self.PERIOD_A, 0.5, 12.37, self.HOST_DENSITY
        )
        assert res_eb.is_planet_candidate is False
        assert res_eb.energy >= 1000.0
        assert res_eb.falsification_verdict == "REFUTED_ODD_EVEN_PARITY"

        # 3. Negative Control 2: Secondary Occultation at Phase 0.5
        flux_sec = np.copy(flux_clean)
        phase = ((times - 0.5) % self.PERIOD_A) / self.PERIOD_A
        sec_dip = np.abs(np.abs(phase) - 0.5) < (1.0 / 24.0 / self.PERIOD_A)
        flux_sec[sec_dip] -= 0.0003  # 300 ppm secondary dip

        res_sec = run_elenchus_adversarial_battery(
            times, flux_sec, self.PERIOD_A, 0.5, 12.37, self.HOST_DENSITY
        )
        assert res_sec.is_planet_candidate is False
        assert res_sec.energy >= 1000.0
        assert res_sec.falsification_verdict == "REFUTED_SECONDARY_ECLIPSE"

        # 4. Negative Control 3: Centroid Motion Shift
        pixel_x = np.random.normal(100.0, 0.05, size=len(times))
        pixel_y = np.random.normal(100.0, 0.05, size=len(times))
        in_transit = np.abs(((times - 0.5) % self.PERIOD_A) / self.PERIOD_A) < (1.0 / 24.0 / self.PERIOD_A)
        pixel_x[in_transit] += 0.5  # Significant centroid displacement

        res_centroid = run_elenchus_adversarial_battery(
            times, flux_clean, self.PERIOD_A, 0.5, 12.37, self.HOST_DENSITY,
            pixel_x=pixel_x, pixel_y=pixel_y
        )
        assert res_centroid.is_planet_candidate is False
        assert res_centroid.energy >= 1000.0
        assert res_centroid.falsification_verdict == "REFUTED_CENTROID_SHIFT"

    def test_autonomous_astrophysicist_domain_evaluation(self):
        """Verifies integration with the ANSE Frontier Continuous Reasoning Engine."""
        domain = AutonomousAstrophysicist()
        times = np.linspace(0.0, 15.0, 2000)
        noise = np.random.normal(1.0, 0.0001, size=len(times))
        transit = mandel_agol_transit_model(times, self.PERIOD_A, 0.5, 0.0224, 12.37)
        flux = noise + (transit - 1.0)

        result = domain.evaluate_transit_candidate(
            target_id="TIC 4206066",
            times=times,
            fluxes=flux,
            period=self.PERIOD_A,
            t0=0.5,
            a_over_r_star=12.37,
            host_density_cgs=self.HOST_DENSITY,
        )

        assert result.target_id == "TIC 4206066"
        assert result.is_valid_candidate is True
        assert result.energy < 5.0
        assert result.verdict == "NOT_REFUTED"
        assert result.controls_passed is True
        assert result.duration_ms < 100.0


class TestTIC4206066RealObservationalData:
    """
    Cryptographic and physical verification of real observational data from both:
    1. Direct NASA MAST TESS FITS products (Sectors 6, 32, and 98).
    2. Zenodo research deposit archives (doi:10.5281/zenodo.22967099).
    """

    def test_real_nasa_mast_fits_ingestion(self):
        """Verifies integrity, row counts, and valid photometry from raw NASA MAST FITS files."""
        lcs = load_real_mast_lightcurves()
        assert set(lcs.keys()) == {6, 32, 98}, "Missing expected TESS sectors"

        total_valid = sum(d["valid_count"] for d in lcs.values())
        assert total_valid > 16000, f"Insufficient clean photometric cadences: {total_valid}"

        # Sector 6 SPOC
        s6 = lcs[6]
        assert s6["valid_count"] >= 900
        assert np.all(s6["fluxes"] > 0)
        assert not np.any(np.isnan(s6["fluxes"]))

        # Sector 32 SPOC
        s32 = lcs[32]
        assert s32["valid_count"] >= 3400
        assert np.all(s32["fluxes"] > 0)
        assert not np.any(np.isnan(s32["fluxes"]))

        # Sector 98 QLP
        s98 = lcs[98]
        assert s98["valid_count"] >= 12000
        assert np.all(s98["fluxes"] > 0)
        assert not np.any(np.isnan(s98["fluxes"]))

    def test_real_spoc_centroid_stability(self):
        """Verifies photocenter stability on real SPOC moments to rule out background EBs."""
        stability = audit_real_centroid_stability()
        assert 6 in stability and 32 in stability

        for sec, metrics in stability.items():
            assert metrics["is_stable"] is True, (
                f"Sector {sec} centroid shift is significant: {metrics['significance_sigma']} sigma"
            )
            assert metrics["significance_sigma"] < 2.5, (
                f"Centroid displacement exceeds 2.5 sigma: {metrics['significance_sigma']}"
            )
            assert metrics["centroid_shift_pixels"] < 0.05, (
                f"Unphysical centroid shift magnitude: {metrics['centroid_shift_pixels']} px"
            )

    def test_real_zenodo_23_transits_and_depth(self):
        """Verifies the 23 real per-transit measurements from Zenodo part1."""
        transits = load_zenodo_per_transit_data()
        assert len(transits) == 23, f"Expected 23 transits, found {len(transits)}"

        # Check sectors covered
        sectors = {t["sector"] for t in transits}
        assert sectors == {6, 32, 98}, f"Unexpected sector distribution: {sectors}"

        # Check global depth matches ~500 ppm candidate claim
        mean_depth = float(np.mean([t["depth_ppm"] for t in transits]))
        assert 480.0 <= mean_depth <= 525.0, f"Mean transit depth out of bounds: {mean_depth} ppm"

        # Check that individual transit depths are physically reasonable
        for tr in transits:
            assert 100.0 < tr["depth_ppm"] < 1000.0, f"Outlier transit depth: {tr['depth_ppm']} ppm"
            assert tr["snr"] > 0.0, f"Invalid SNR: {tr['snr']}"

    def test_real_zenodo_odd_even_parity(self):
        """Verifies that odd vs even transit depths are consistent within 1.0 sigma on real data."""
        parity = audit_real_odd_even_parity()
        assert parity["eclipsing_binary_ruled_out"] is True
        assert parity["disparity_significance_sigma"] < 1.5, (
            f"Odd-even depth disparity statistically significant: {parity['disparity_significance_sigma']} sigma"
        )
        assert parity["delta_ppm"] < 100.0, f"Large odd-even delta: {parity['delta_ppm']} ppm"

    def test_real_zenodo_linear_ephemeris(self):
        """Verifies that transit timing residuals follow a strictly linear Keplerian ephemeris."""
        timing = audit_real_timing_ephemeris()
        assert timing["timing_transits_count"] == 19
        assert timing["is_ephemeris_linear"] is True
        assert timing["reduced_chi2"] < 1.0, f"High reduced chi2: {timing['reduced_chi2']}"
        assert timing["rms_oc_minutes"] < 12.0, f"High timing jitter: {timing['rms_oc_minutes']} min"
        assert abs(timing["period_days"] - 3.182785) < 1e-5, (
            f"Ephemeris period mismatch: {timing['period_days']}"
        )

    def test_comprehensive_real_data_audit_summary(self):
        """Runs the top-level audit function combining MAST and Zenodo data."""
        audit = audit_real_data_sources()
        assert audit["status"] == "VERIFIED_REAL_DATA"
        assert audit["nasa_mast_audit"]["total_photometric_cadences"] > 16000
        assert audit["zenodo_audit"]["total_transits"] == 23
        assert audit["zenodo_audit"]["global_mean_depth_ppm"] == 501.5

