# ExoElenchus-TIC4206066 Workflow Guide
**Subproject:** TIC 4206066 Epistemic Reproducibility & Socratic Falsification Pipeline  
**Script Runner:** `scripts/run_tic4206066_workflow.py`  
**Target Reference:** [Zenodo doi:10.5281/zenodo.22967456](https://doi.org/10.5281/zenodo.22967456) · [doi:10.5281/zenodo.23175179](https://doi.org/10.5281/zenodo.23175179)  

---

## 1. Workflow Architecture & Pipeline Steps

The automated workflow executes the eight core stages of scientific reproduction and falsification:

```
[Step 1: Ingestion & Metadata Check]
  │ Validates TIC 4206066 properties against Gaia DR3 & Zenodo metadata
  ▼
[Step 2: Keplerian Dynamics & Transit Modeling]
  │ Computes semi-major axis, a/R_*, quadratic limb-darkening, and derived density rho_*
  ▼
[Step 3: Socratic Elenchus Falsification Battery & Real Vetting]
  │ Evaluates odd-even parity, secondary eclipse, centroid motion on real data, plus synthetic method controls
  ▼
[Step 4: Leave-One-Sector-Out Timing Hold-Out]
  │ Blindly predicts transits across Sectors 6 (2018), 32 (2020), and 98 (2025)
  ▼
[Step 5: Pre-Registered Sector 110 Test Runner & Self-Test]
  │ Evaluates Sector 110 (Nov 2026, DDT #100) ephemerides; signs SHA-256 provenance receipt
  ▼
[Step 6: SMCCC Dissemination & Quality Audit]
  │ Evaluates text against audience reward model (asserts 0.0 moderation risk, 1.0 grounding)
  ▼
[Step 7: Automated Pytest Suite Execution]
  │ Runs tests/astrophysics/test_tic4206066_falsification.py
  ▼
[Step 8: Lean 4 Formal Verification]
  │ Verifies Keplerian invariants and odd-even parity theorems in formal/ANSE/ExoplanetTransit.lean
```

---

## 2. CLI Usage Commands

```bash
# Run all 8 steps sequentially and generate certified summary
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py

# Run individual steps:
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 1 # Ingestion
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 2 # Transit Kernel
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 3 # Elenchus Falsification
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 4 # Temporal Holdout
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 5 # Sector 110 Pre-Registration
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 6 # SMCCC Audit
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 7 # Pytest Suite
PYTHONPATH=. uv run --no-sync python scripts/run_tic4206066_workflow.py --step 8 # Lean 4 Formal Verification
```

---

## 3. Output Artifacts & Verification Receipts

1. **Certified Execution Receipt:**  
   [`results/tic4206066_execution_receipt.json`](file:///home/xavkal/AutoevolveAI/results/tic4206066_execution_receipt.json)  
   Contains Sector 110 ephemerides, predicted transit windows, and SHA-256 cryptographic hash.

2. **Workflow Summary Receipt:**  
   [`results/tic4206066_workflow_summary.json`](file:///home/xavkal/AutoevolveAI/results/tic4206066_workflow_summary.json)  
   Records execution timestamps, latencies, and pass/fail statuses across all 7 steps.

3. **SMCCC Reddit Evaluation Receipt:**  
   [`results/reddit_evaluation_claudeai_exoplanet_audit.json`](file:///home/xavkal/AutoevolveAI/results/reddit_evaluation_claudeai_exoplanet_audit.json)  
   Logs the verified audience reward score (0.85), 0.0 moderation risk, and 1.0 grounding score.

4. **Sector 110 Scheduler State:**  
   [`results/tic4206066_scheduler_state.json`](file:///home/xavkal/AutoevolveAI/results/tic4206066_scheduler_state.json)  
   Records the live countdown counter, armed status (`DRY_RUN_VERIFIED` / `ARMED`), and latest dry-run telemetry.

---

## 4. TESS Sector 110 Automated Ingestion Scheduler

To execute autonomously when NASA MAST releases Sector 110 data in late November 2026:

```bash
# Check current countdown timer and scheduler status
PYTHONPATH=. uv run --no-sync python scripts/tic4206066_scheduler.py --status

# Execute full dry-run test (simulates release and validates all 7 workflow steps)
PYTHONPATH=. uv run --no-sync python scripts/tic4206066_scheduler.py --dry-run

# Check if target release timestamp (2026-11-28 18:00 UTC) has arrived
PYTHONPATH=. uv run --no-sync python scripts/tic4206066_scheduler.py --check

# Run background polling daemon
PYTHONPATH=. uv run --no-sync python scripts/tic4206066_scheduler.py --daemon --interval 3600
```

