#!/usr/bin/env python3
"""
Automated End-to-End Workflow: TIC 4206066 Epistemic Reproducibility & Falsification.
Orchestrates all 7 steps of the ExoElenchus-TIC4206066 pipeline:
1. Ingestion & Stellar Parameter Verification
2. High-Performance Keplerian Transit Modeling
3. Socratic Elenchus Falsification Battery (with Positive & Negative Controls)
4. 3-Fold Temporal Holdout Cross-Validation
5. Pre-Registered NASA TESS Sector 110 Blind Runner & Receipt Signing
6. SMCCC Audience Reward & Quality Gate Audit
7. Pytest Test Suite & Attestation Execution
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from anse.astrophysics.elenchus_vetting import run_elenchus_adversarial_battery
from anse.astrophysics.ephemeris import load_ephemeris_a
from anse.astrophysics.real_data_audit import (
    leave_one_sector_out_timing,
    load_events,
    photometric_depth_at_ephemeris,
    verify_consumed_zenodo_files,
    weighted_odd_even,
)
from anse.astrophysics.real_data_ingest import (
    audit_real_centroid_stability,
    audit_real_data_sources,
    audit_real_timing_ephemeris,
    load_real_mast_lightcurves,
)
from anse.community.audience_reward import compute_audience_reward
from anse.physics.exoplanet_transit import (
    calculate_a_over_r_star,
    derive_stellar_density,
    kepler_semi_major_axis,
    mandel_agol_transit_model,
    transit_duration_t14,
)
from subprojects.tic4206066_sector110_runner import execute_preregistered_runner

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("ExoElenchusWorkflow")


class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for NumPy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


@dataclass
class StepResult:
    step_id: int
    name: str
    passed: bool
    duration_ms: float
    details: dict[str, Any]


class TIC4206066Workflow:
    """Executes and verifies the full TIC 4206066 scientific reproduction pipeline."""

    def __init__(self, metadata_path: str | None = None):
        if metadata_path is None:
            for p in ["data/metadata.json", "data/astrophysics/tic4206066/metadata.json"]:
                if os.path.exists(p):
                    metadata_path = p
                    break
            else:
                metadata_path = "data/metadata.json"
        self.metadata_path = metadata_path
        with open(metadata_path) as f:
            self.meta = json.load(f)
        self.results: list[StepResult] = []

    def step1_ingestion(self) -> StepResult:
        """Step 1: Ingestion & real-data audit (MAST photometry + Zenodo deposit, checksummed)."""
        start = time.perf_counter()
        star = self.meta["target_star"]
        sig_a = self.meta["signal_a"]
        eph = load_ephemeris_a()

        real_audit = audit_real_data_sources()
        integrity = verify_consumed_zenodo_files()
        events = load_events()
        ev_depth = np.array([e["depth_ppm"] for e in events])
        ev_err = np.array([e["depth_err_ppm"] for e in events])
        w = 1.0 / ev_err**2
        comb_depth = float(np.sum(w * ev_depth) / np.sum(w))
        comb_err = float(1.0 / np.sqrt(np.sum(w)))

        # Cross-check the metadata against the deposit (independent sources, not self-comparison)
        checks = {
            "zenodo_files_sha256_match": integrity["all_match"],
            "metadata_period_within_3sigma_of_ephemeris": abs(sig_a["period_days"] - eph.period_d)
            < 3.0 * eph.period_err_d,
            "metadata_t0_matches_ephemeris": abs(sig_a["t0_bjd"] - eph.t0_bjd) < eph.t0_err_d,
            "metadata_depth_within_3sigma_of_event_mean": abs(sig_a["depth_ppm"] - comb_depth)
            < 3.0 * max(comb_err, 1.0),
            "event_weighted_mean_depth_ppm": round(comb_depth, 1),
            "event_weighted_mean_depth_err_ppm": round(comb_err, 1),
            "gaia_id_present": "3220388198192519424" in star["gaia_dr3_id"],
            "transits_count_in_catalogue": len(events),
            "transits_count_matches_metadata": len(events) == sig_a["total_transits_observed"],
            "real_mast_valid_cadences": real_audit["nasa_mast_audit"]["total_photometric_cadences"],
            "real_mast_cadences_clean": real_audit["nasa_mast_audit"]["total_photometric_cadences"] > 16000,
            "real_data_sources_verified": real_audit["status"] == "VERIFIED_REAL_DATA",
        }
        passed = all(v is True for v in checks.values() if isinstance(v, bool))
        dur = (time.perf_counter() - start) * 1000.0
        return StepResult(1, "Ingestion & Real Observational Data Audit", passed, round(dur, 2), checks)

    def step2_transit_kernel(self) -> StepResult:
        """Step 2: Keplerian geometry & transit model vs INDEPENDENT published values."""
        start = time.perf_counter()
        star = self.meta["target_star"]
        P = load_ephemeris_a().period_d

        a_au, _ = kepler_semi_major_axis(star["mass_solar"], P)
        a_over_r = calculate_a_over_r_star(star["mass_solar"], star["radius_solar"], P)
        derived_rho = derive_stellar_density(P, a_over_r)

        # Published (Zenodo FACTS): k = 0.0204, b = 0.35 (Mann-prior fit), T14 = 1.87 +- 0.05 h,
        # a = 0.0358 AU, a/R* = 12.37, combined depth 498 +- 29 ppm, rho_free 97.7% UL 3.4-3.7 rho_sun.
        k_pub, b_pub = 0.0204, 0.35
        t14_b0 = transit_duration_t14(P, a_over_r, k_pub, 0.0)
        t14_b062 = transit_duration_t14(P, a_over_r, k_pub, 0.62)
        t14_pub = transit_duration_t14(P, a_over_r, k_pub, b_pub)
        central = mandel_agol_transit_model(np.array([0.0]), P, 0.0, k_pub, a_over_r)
        model_depth = float((1.0 - central[0]) * 1e6)
        rho_sun_cgs = 1.41
        checks = {
            "semi_major_axis_au": round(a_au, 4),
            "a_over_r_star": round(a_over_r, 3),
            "derived_rho_cgs": round(derived_rho, 3),
            "derived_rho_rho_sun": round(derived_rho / rho_sun_cgs, 3),
            "t14_hours_b0": round(t14_b0, 3),
            "t14_hours_b0.35": round(t14_pub, 3),
            "t14_hours_b0.62": round(t14_b062, 3),
            "model_central_depth_ppm": round(model_depth, 1),
            "a_au_matches_published_0.0358": abs(a_au - 0.0358) < 0.0005,
            "a_over_r_matches_published_12.37": abs(a_over_r - 12.37) < 0.15,
            "published_t14_inside_b_range": t14_b062 <= 1.87 <= t14_b0,
            "t14_matches_published_1.87h": abs(t14_pub - 1.87) < 0.10,
            "model_depth_matches_published_498ppm": abs(model_depth - 498.0) < 0.10 * 498.0,
            "density_below_free_fit_upper_limit": derived_rho / rho_sun_cgs < 3.4,
        }
        passed = all(v is True for v in checks.values() if isinstance(v, bool))
        dur = (time.perf_counter() - start) * 1000.0
        return StepResult(2, "Keplerian Dynamics & Transit Modeling", passed, round(dur, 2), checks)

    def step3_elenchus_falsification(self) -> StepResult:
        """Step 3: Elenchus battery (synthetic controls) + falsification tests on REAL data."""
        start = time.perf_counter()
        star = self.meta["target_star"]
        eph = load_ephemeris_a()
        P = eph.period_d
        a_over_r = calculate_a_over_r_star(star["mass_solar"], star["radius_solar"], P)

        # Synthetic injection: a self-test of the battery, NOT evidence about the star.
        rng = np.random.default_rng(42)
        times = np.linspace(0.0, 15.0, 2000)
        noise = rng.normal(1.0, 0.0001, size=len(times))
        transit = mandel_agol_transit_model(times, P, 0.5, 0.0204, a_over_r)
        flux = noise + (transit - 1.0)
        vetting = run_elenchus_adversarial_battery(
            times, flux, P, 0.5, a_over_r, star["mean_density_cgs"]
        )

        # Real-data falsification tests
        mast = load_real_mast_lightcurves()
        real_centroid = audit_real_centroid_stability(mast)
        real_parity = weighted_odd_even()
        real_timing = audit_real_timing_ephemeris()
        real_depth = photometric_depth_at_ephemeris(mast, eph)

        w = np.array([1.0 / v["depth_err_ppm"] ** 2 for v in real_depth.values()])
        d = np.array([v["depth_ppm"] for v in real_depth.values()])
        comb_depth = float(np.sum(w * d) / np.sum(w))
        comb_snr = comb_depth * float(np.sqrt(np.sum(w)))

        checks = {
            "synthetic_selftest_candidate_valid": vetting.is_planet_candidate,
            "synthetic_selftest_energy": vetting.energy,
            "positive_control_recovered": vetting.positive_control_passed,
            "negative_controls_rejected": vetting.negative_controls_rejected,
            "real_photometry_depth_by_sector": real_depth,
            "real_photometry_combined_depth_ppm": round(comb_depth, 1),
            "real_photometry_combined_snr": round(comb_snr, 2),
            "real_photometry_note": "unsmoothed boxes, no detrending: phase-0.5 controls show red noise "
            "at ~2-3 sigma in S32/S98, so only the combined S/N is gated",
            "real_centroid_by_sector": real_centroid,
            "real_centroid_stable": all(m["is_stable"] for m in real_centroid.values()),
            "real_odd_even": real_parity,
            "real_eb_ruled_out": real_parity["eclipsing_binary_ruled_out"],
            "real_ephemeris_reduced_chi2": real_timing["reduced_chi2"],
            "real_ephemeris_linear": real_timing["is_ephemeris_linear"],
        }
        passed = bool(
            vetting.is_planet_candidate
            and vetting.positive_control_passed
            and vetting.negative_controls_rejected
            and comb_snr >= 5.0
            and all(v["depth_ppm"] > 0 for v in real_depth.values())
            and checks["real_centroid_stable"]
            and checks["real_eb_ruled_out"]
            and checks["real_ephemeris_linear"]
        )
        dur = (time.perf_counter() - start) * 1000.0
        return StepResult(3, "Socratic Elenchus Falsification Battery & Real Vetting", passed, round(dur, 2), checks)

    def step4_temporal_holdout(self) -> StepResult:
        """Step 4: genuine leave-one-sector-out timing hold-out on the 19 timing events."""
        start = time.perf_counter()
        res = leave_one_sector_out_timing()
        passed = bool(res["passed"])
        dur = (time.perf_counter() - start) * 1000.0
        return StepResult(4, "Leave-One-Sector-Out Timing Hold-Out", passed, round(dur, 2), res)

    def step5_sector110_preregistration(self) -> StepResult:
        """Step 5: frozen Sector 110 predictions + synthetic self-test of the decision logic."""
        start = time.perf_counter()
        receipt = execute_preregistered_runner()
        st = receipt.self_test or {}
        checks = {
            "target": receipt.target_id,
            "sector": receipt.sector_window["sector"],
            "dates_utc": receipt.sector_window["dates_utc"],
            "window_verified_against_tess_timeline": receipt.sector_window[
                "window_verified_against_tess_timeline"
            ],
            "signal_a_expected_transits": receipt.signal_a_predictions["expected_transits_count"],
            "signal_a_sigma_min": receipt.signal_a_predictions["uncertainty_mins"],
            "signal_b_expected_transits": receipt.signal_b_predictions["expected_transits_count"],
            "status": receipt.overall_status,
            "provenance_hash": receipt.provenance_hash,
            "selftest_injected_outcome": st.get("injected_signal_outcome"),
            "selftest_flat_outcome": st.get("flat_baseline_outcome"),
            "selftest_logic_ok": st.get("logic_ok"),
            "receipt_file_exists": os.path.exists("results/tic4206066_execution_receipt.json"),
            "note": "predictions + SYNTHETIC self-test only; no Sector 110 data evaluated",
        }
        passed = bool(checks["selftest_logic_ok"] and checks["receipt_file_exists"])
        dur = (time.perf_counter() - start) * 1000.0
        return StepResult(5, "Pre-Registered Sector 110 Predictions & Decision-Logic Self-Test", passed, round(dur, 2), checks)

    def step6_smccc_audit(self) -> StepResult:
        """Step 6: SMCCC Dissemination & Audience Reward Audit."""
        start = time.perf_counter()
        eval_path = "results/reddit_evaluation_claudeai_exoplanet_audit.json"

        if os.path.exists(eval_path):
            with open(eval_path) as f:
                saved = json.load(f)
            res = saved["audit_metrics"]
        else:
            title = "How to audit agent-written analysis when a result can look completely convincing and still be wrong?"
            sub = "r/ClaudeAI"
            stmt = "Context: Response to r/ClaudeAI exoplanet discovery thread detailing falsification harnesses, negative controls, and zero-trust verification."
            body_snippet = (
                "In our research lab Socrate AI Lab, we formalize this into an open-source verification "
                "harness called Elenchus. Never treat a green test suite as a scientific license. "
                "Every automated vetting script must ship with a positive control demonstrated to fire "
                "and a negative control demonstrated to fail."
            )
            res = compute_audience_reward(
                subreddit=sub,
                title=title,
                body=body_snippet,
                submission_statement=stmt,
            )

        checks = {
            "total_reward": res["total_reward"],
            "grounding_score": res["grounding_score"],
            "modesty_score": res["modesty_score"],
            "moderation_risk": res["moderation_risk"],
            "is_recommended": res["is_recommended"],
        }
        passed = bool(res["is_recommended"] and res["moderation_risk"] == 0.0 and res["modesty_score"] >= 0.8)
        dur = (time.perf_counter() - start) * 1000.0
        return StepResult(6, "SMCCC Dissemination & Quality Audit", passed, round(dur, 2), checks)

    def step7_pytest_suite(self) -> StepResult:
        """Step 7: Automated Pytest Suite Execution."""
        start = time.perf_counter()
        test_paths = [
            "tests/test_tic4206066_falsification.py",
            "tests/astrophysics/test_tic4206066_falsification.py",
        ]
        test_target = next((p for p in test_paths if os.path.exists(p)), "tests/test_tic4206066_falsification.py")
        cmd = [sys.executable, "-m", "pytest", test_target, "-q"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        dur = (time.perf_counter() - start) * 1000.0

        passed = proc.returncode == 0
        checks = {
            "exit_code": proc.returncode,
            "stdout": proc.stdout.strip(),
            "stderr": proc.stderr.strip()[:300],
            "test_target": test_target,
        }
        return StepResult(7, "Automated Pytest Suite Execution", passed, round(dur, 2), checks)

    def step8_lean4_formal_verification(self) -> StepResult:
        """Step 8: Lean 4 Formal Verification (Keplerian Invariants & Socratic Elenchus Null Hypothesis)."""
        start = time.perf_counter()
        formal_dir = os.path.join(os.getcwd(), "formal")
        lean_file = os.path.join(formal_dir, "ANSE", "ExoplanetTransit.lean")
        olean_file = os.path.join(formal_dir, ".lake", "build", "lib", "lean", "ANSE", "ExoplanetTransit.olean")

        if not os.path.exists(lean_file):
            return StepResult(8, "Lean 4 Formal Verification", False, 0.0, {"error": "Lean file missing"})

        # Zero sorry check
        with open(lean_file, encoding="utf-8") as f:
            lines = f.readlines()
        has_sorry = False
        for line in lines:
            code_part = line.split("--")[0].strip()
            if re.search(r"\bsorry\b", code_part):
                has_sorry = True
                break

        theorems = sorted(
            set(re.findall(r"^\s*(?:theorem|lemma)\s+([A-Za-z0-9_'.]+)", "".join(lines), re.M))
        )

        import shutil
        lake_bin = shutil.which("lake")
        clean_elaboration = False
        exit_code = -1  # fail closed: -1 means Lean never ran
        stderr = ""

        # Optional: reuse a prebuilt Lake/Mathlib cache (set LAKE_CACHE_DIR); no machine-specific default.
        cache_source = os.environ.get("LAKE_CACHE_DIR", "")
        lake_dir = os.path.join(formal_dir, ".lake")
        if cache_source and not os.path.exists(lake_dir) and os.path.exists(cache_source):
            try:
                os.symlink(cache_source, lake_dir)
            except OSError:
                pass

        has_mathlib = os.path.exists(os.path.join(formal_dir, ".lake", "packages", "mathlib")) or os.path.exists(os.path.join(formal_dir, ".lake", "build"))
        if lake_bin and has_mathlib and os.path.exists(os.path.join(formal_dir, "lakefile.toml")):
            try:
                cmd = [lake_bin, "env", "lean", "ANSE/ExoplanetTransit.lean"]
                proc = subprocess.run(cmd, cwd=formal_dir, capture_output=True, text=True, timeout=900)
                exit_code = proc.returncode
                stderr = proc.stderr.strip()
                clean_elaboration = (proc.returncode == 0 and len(stderr) == 0)
            except Exception as exc:
                logger.warning("Lake elaboration check: %s", exc)

        checks = {
            "lean_file_exists": True,
            "olean_compiled": os.path.exists(olean_file),
            "zero_sorry_soundness": not has_sorry,
            "lean_elaboration_exit_code": exit_code,
            "lean_elaboration_clean": clean_elaboration,
            "formalized_theorems": theorems,
            "theorem_count": len(theorems),
        }
        passed = (
            checks["lean_file_exists"]
            and len(checks["formalized_theorems"]) >= 6
            and checks["zero_sorry_soundness"]
            and checks["lean_elaboration_clean"]
        )
        dur = (time.perf_counter() - start) * 1000.0
        return StepResult(8, "Lean 4 Formal Verification (Keplerian Invariants & Socratic Elenchus)", passed, round(dur, 2), checks)

    def run_all(self) -> dict[str, Any]:
        """Executes all 8 workflow steps sequentially."""
        logger.info("Starting ExoElenchus-TIC4206066 Workflow...")
        self.results = [
            self.step1_ingestion(),
            self.step2_transit_kernel(),
            self.step3_elenchus_falsification(),
            self.step4_temporal_holdout(),
            self.step5_sector110_preregistration(),
            self.step6_smccc_audit(),
            self.step7_pytest_suite(),
            self.step8_lean4_formal_verification(),
        ]

        all_passed = all(r.passed for r in self.results)
        import random
        seed = 42
        random.seed(seed)
        np.random.seed(seed)
        git_hash = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        summary = {
            "workflow_name": "ExoElenchus-TIC4206066",
            "target": "TIC 4206066 (StKM 1-561)",
            "all_passed": all_passed,
            "total_duration_ms": round(sum(r.duration_ms for r in self.results), 2),
            "seed": seed,
            "git_hash": git_hash,
            "steps": [asdict(r) for r in self.results],
            "timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        }

        os.makedirs("results", exist_ok=True)
        with open("results/tic4206066_workflow_summary.json", "w") as f:
            json.dump(summary, f, indent=2, cls=NumpyEncoder)

        os.makedirs("data", exist_ok=True)
        with open("data/tic4206066_workflow_summary.json", "w") as f:
            json.dump(summary, f, indent=2, cls=NumpyEncoder)

        logger.info(f"Workflow Complete. Status: {'ALL PASSED' if all_passed else 'FAILED'}")
        return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TIC 4206066 ExoElenchus Workflow Runner")
    parser.add_argument("--step", choices=["all", "1", "2", "3", "4", "5", "6", "7", "8"], default="all")
    args = parser.parse_args()

    wf = TIC4206066Workflow()

    if args.step == "all":
        summary = wf.run_all()
        print("\n" + "=" * 60)
        print("         TIC 4206066 WORKFLOW EXECUTION SUMMARY")
        print("=" * 60)
        for s in summary["steps"]:
            status_str = "✅ PASS" if s["passed"] else "❌ FAIL"
            print(f"[{status_str}] Step {s['step_id']}: {s['name']} ({s['duration_ms']} ms)")
        print("-" * 60)
        print(f"OVERALL STATUS: {'✅ ALL 8 STEPS PASSED' if summary['all_passed'] else '❌ WORKFLOW FAILED'}")
        print(f"Total Duration: {summary['total_duration_ms']} ms")
        print("Summary logged to results/tic4206066_workflow_summary.json")
        print("=" * 60)
        sys.exit(0 if summary["all_passed"] else 1)
    else:
        step_map = {
            "1": wf.step1_ingestion,
            "2": wf.step2_transit_kernel,
            "3": wf.step3_elenchus_falsification,
            "4": wf.step4_temporal_holdout,
            "5": wf.step5_sector110_preregistration,
            "6": wf.step6_smccc_audit,
            "7": wf.step7_pytest_suite,
            "8": wf.step8_lean4_formal_verification,
        }
        res = step_map[args.step]()
        status_str = "✅ PASS" if res.passed else "❌ FAIL"
        print(f"[{status_str}] Step {res.step_id}: {res.name} ({res.duration_ms} ms)")
        print(json.dumps(res.details, indent=2, cls=NumpyEncoder))
        sys.exit(0 if res.passed else 1)
