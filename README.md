# TIC 4206066 Exoplanet Verification & Socratic Elenchus Harness

[![CI Verification Battery](https://github.com/xaviercallens/TIC4206066-Exoplanet-Verification/actions/workflows/ci.yml/badge.svg)](https://github.com/xaviercallens/TIC4206066-Exoplanet-Verification/actions/workflows/ci.yml)
[![Lean 4 Soundness](https://img.shields.io/badge/Lean_4-Formal_Proof_0_Sorry-blue.svg)](formal/ANSE/ExoplanetTransit.lean)
[![GCP Cloud Run](https://img.shields.io/badge/GCP_Cloud_Run-Ready-success.svg)](docs/GCP_DEPLOYMENT.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zenodo Preprint](https://img.shields.io/badge/Zenodo-Preprint%2010.5281%2Fzenodo.22967456-blue)](https://doi.org/10.5281/zenodo.22967456)
[![Zenodo Data & Code](https://img.shields.io/badge/Zenodo-Data%20%26%20Code%2010.5281%2Fzenodo.22967099-green)](https://doi.org/10.5281/zenodo.22967099)
[![Zenodo Pre-Registration](https://img.shields.io/badge/Zenodo-Pre--Registration%2010.5281%2Fzenodo.23175179-orange)](https://doi.org/10.5281/zenodo.23175179)

Independent, zero-trust scientific reproduction, adversarial falsification battery, Lean 4 formal verification, and Google Cloud Run discovery portal for candidate exoplanet **TIC 4206066 b** discovered around late-K dwarf **TIC 4206066** (StKM 1-561) by independent researcher **Pavel Rabtsevich**.

---

## 1. Executive Summary & Astrophysical Grounding

Following the core principle of the **Socratic Elenchus** verification methodology:
> *"The useful way to use an agent for science isn't to ask: 'Is this a planet?'  
> It's to keep giving it ways to prove that it isn't."*

This repository packages the complete verification pipeline, from raw NASA MAST TESS FITS photometry (Sectors 6, 32, 98) and Zenodo research deposits to Lean 4 machine-checked mathematical proofs and a Google Cloud Run production dashboard.

```
Target Host Star:
- ID:             TIC 4206066 (StKM 1-561, Gaia DR3 3220388198192519424)
- Spectral Type:  Late-K dwarf (T_eff ≈ 4020 K)
- Stellar Mass:   0.606 ± 0.015 M_Sun (Mann+2019)
- Stellar Radius: 0.623 ± 0.019 R_Sun (Mann+2015)
- Host Density:   3.534 g/cm^3
- Distance:       ~116 light-years

Primary Candidate (Signal A / TIC 4206066.01):
- Orbital Period:   3.182785 ± 0.000006 days
- Transit Depth:    498.0 ppm (~0.05%)
- Transit Duration: 1.87 hours
- Planetary Radius: 1.4 ± 0.1 R_Earth (Super-Earth regime)
- Semi-Major Axis:  0.0358 AU (a/R_* = 12.37)
- Observations:     23 transits across TESS Sectors 6 (2018), 32 (2020), 98 (2025)

Secondary Candidate (Signal B / TIC 4206066.02):
- Orbital Period:   11.13274 days
- Transit Depth:    ~1000 ppm (~0.10%)
- Planetary Radius: ~2.2 R_Earth (Sub-Neptune regime, tentative ~3-sigma)
```

---

## 2. Architecture & Components

```
TIC4206066-Exoplanet-Verification/
├── app.py                      # FastAPI ASGI backend (GCP Cloud Run production service)
├── stitch_template.html        # Stitch-designed astrophysics dashboard (GCP front end)
├── Dockerfile                  # Production container definition for Google Cloud Run
├── cloudbuild.yaml             # Google Cloud Build deployment pipeline
├── deploy_cloud_run.sh         # One-command GCP Cloud Run deployment script
├── requirements.txt            # Python dependencies
├── pyproject.toml              # Project metadata & pytest configuration
├── anse/                       # Core astrophysics & verification engine
│   ├── physics/
│   │   └── exoplanet_transit.py   # Analytical Mandel-Agol limb-darkened transit kernel
│   ├── astrophysics/
│   │   ├── elenchus_vetting.py    # Socratic Elenchus adversarial falsification battery
│   │   ├── ephemeris.py           # Propagated Keplerian ephemeris modeling
│   │   ├── real_data_ingest.py    # NASA MAST FITS ingestion (S06, S32, S98)
│   │   ├── real_data_audit.py     # Checksum verification & leave-one-out cross-validation
│   │   ├── sector110_scheduler.py # Live countdown ticker & autonomous ingestion scheduler
│   │   └── domain.py              # Autonomous astrophysicist reasoning & energy gates
│   └── community/                 # Audience reward model & dissemination audit
├── subprojects/
│   └── tic4206066_sector110_runner.py # Pre-registered Sector 110 blind test runner
├── scripts/
│   ├── run_tic4206066_workflow.py # 8-step end-to-end automated verification runner
│   └── tic4206066_scheduler.py    # Real-time countdown & scheduler CLI
├── formal/                     # Formal verification in Lean 4
│   ├── lean-toolchain          # Lean 4 version (v4.34.0-rc2)
│   ├── lakefile.toml           # Lake build manifest with Mathlib dependency
│   ├── ANSE/
│   │   └── ExoplanetTransit.lean # 6 Formal proofs with 0 sorry axioms
│   └── README.md               # Mathematical proof specifications
├── data/                       # Ground-truth observational data & FITS products
│   ├── metadata.json           # Target stellar & planetary parameters
│   ├── mast/                   # Authentic NASA MAST TESS FITS (Sectors 6, 32, 98)
│   └── zenodo/                 # Zenodo research deposits (part1, part2, part3)
├── results/                    # Certified execution receipts & provenance hashes
├── tests/                      # Automated test suite (28 tests)
└── docs/                       # Technical & scientific documentation
```

---

## 3. Socratic Elenchus Falsification Battery

To prevent confirmation bias, candidate transit signals are subjected to strict mechanical refutation gates with mandatory controls:

| Gate / Invariant | Physical Threat Ruled Out | Test Metric & Threshold | Epistemic Verdict |
| :--- | :--- | :--- | :--- |
| **Odd-Even Parity** | Grazing Eclipsing Binary (EB) with alternating temperatures | $\|\delta_{\text{odd}} - \delta_{\text{even}}\| / \sigma < 2.5$ | **PASSED** ($\Delta = 60.7\text{ ppm}, 0.87\sigma$) |
| **Secondary Eclipse** | Occulting secondary companion at phase $\phi = 0.5$ | $\delta_{\text{sec}} < 150\text{ ppm}$ or $\sigma < 3.5$ | **PASSED** (no secondary eclipse detected) |
| **Centroid Decorrelation** | Blended background eclipsing binary (BGEB) within TESS PSF | In-transit photocenter shift $\Delta r < 2.5\sigma$ | **PASSED** (S6: $0.03\sigma$, S32: $0.52\sigma$) |
| **Keplerian Density Prior** | Unphysical orbital fitting violating host star structure | $\rho_{*,\text{transit}} / \rho_{*,\text{host}} \in [0.25, 4.0]$ | **PASSED** ($\rho_* = 3.534\text{ g/cm}^3$, ratio $= 1.006$) |

### Mandatory Controls:
- **Positive Control:** Injected limb-darkened synthetic planet ($P=3.18\text{d}, \delta=500\text{ ppm}$) into Gaussian noise $\rightarrow$ **100% Recovered** (`NOT_REFUTED`).
- **Negative Control 1 (Eclipsing Binary):** Injected 3x depth asymmetry into alternate transits $\rightarrow$ **REJECTED** (`REFUTED_ODD_EVEN_PARITY`).
- **Negative Control 2 (Secondary Occultation):** Injected $300\text{ ppm}$ dip at phase $0.5$ $\rightarrow$ **REJECTED** (`REFUTED_SECONDARY_ECLIPSE`).
- **Negative Control 3 (Centroid Motion):** Injected $0.5\text{ px}$ centroid shift $\rightarrow$ **REJECTED** (`REFUTED_CENTROID_SHIFT`).

---

## 4. Formal Verification in Lean 4

The orbital dynamics and Elenchus vetting criteria are formally formalized and proved in [`formal/ANSE/ExoplanetTransit.lean`](formal/ANSE/ExoplanetTransit.lean) using **Lean 4** and **Mathlib**:

- **Theorem 1 (`stellar_density_pos`):** Mean spherical stellar density $\rho_* = \frac{M_*}{\frac{4}{3}\pi R_*^3} > 0$ for all physical stars.
- **Theorem 2 (`kepler_transit_stellar_density_invariant`):** Under Kepler's Third Law ($a^3 = \frac{G M_* P^2}{4\pi^2}$), the transit-derived density $\rho_{\text{transit}} = \frac{3\pi}{G P^2}\left(\frac{a}{R_*}\right)^3$ identically equals the physical stellar density $\rho_*$.
- **Theorems 3 & 4 (`transit_depth_nonneg`, `transit_depth_pos`):** Planetary transit depth $\delta = (R_p / R_*)^2 \ge 0$, and strictly positive for non-zero planet radius.
- **Theorem 5 (`odd_even_transit_depth_symmetric`):** **Socratic Elenchus Null Hypothesis**: For genuine planetary transits, depth disparity $\delta(2k+1) - \delta(2k) \equiv 0$.
- **Theorem 6 (`eclipsing_binary_odd_even_disparity`):** Eclipsing binaries with alternating depths $\delta_{\text{prim}} \ne \delta_{\text{sec}}$ produce non-zero odd-even disparity ($\delta(2k) \ne \delta(2k+1)$), enabling mathematical falsification.

All proofs compile with **0 axioms** and **0 `sorry` placeholders**.

---

## 5. Pre-Registered NASA TESS Sector 110 Blind Test

The upcoming NASA TESS observation:
- **Sector 110 Window:** October 30 to November 25, 2026 (Director's Discretionary Time DDT #100).
- **Nominal Data Release:** Late November 2026 (~2026-11-28 18:00 UTC).
- **Signal A Predictions:** 8 predicted transits (BJD 2461345.66 to 2461367.94), mid-transit uncertainty $\pm 4.25\text{ min}$.
- **Signal B Predictions:** 3 predicted transits (BJD 2461345.98, 2461357.11, 2461368.25).
- **Execution Receipt:** Certified with SHA-256 provenance hash: `5a23f5dcf0b1535820bfa1a0478756f74ee893f80fd2e944f0d1178bbf2ebccd`.
- **Status:** `AWAITING_SECTOR_110_DATA__UNREGISTERED_SELF_TEST_PASSED` (Fails closed until authentic Sector 110 FITS are released).

---

## 6. Quick Start & Execution

### 6.1 Install Dependencies
```bash
pip install -r requirements.txt
```

### 6.2 Run Complete 8-Step Verification Workflow
```bash
PYTHONPATH=. python scripts/run_tic4206066_workflow.py
```
Expected output:
```
============================================================
         TIC 4206066 WORKFLOW EXECUTION SUMMARY
============================================================
[✅ PASS] Step 1: Ingestion & Real Observational Data Audit
[✅ PASS] Step 2: Keplerian Dynamics & Transit Modeling
[✅ PASS] Step 3: Socratic Elenchus Falsification Battery & Real Vetting
[✅ PASS] Step 4: Leave-One-Sector-Out Timing Hold-Out
[✅ PASS] Step 5: Pre-Registered Sector 110 Predictions & Decision-Logic Self-Test
[✅ PASS] Step 6: SMCCC Dissemination & Quality Audit
[✅ PASS] Step 7: Automated Pytest Suite Execution
[✅ PASS] Step 8: Lean 4 Formal Verification (Keplerian Invariants & Socratic Elenchus)
------------------------------------------------------------
OVERALL STATUS: ✅ ALL 8 STEPS PASSED
```

### 6.3 Run Test Suite
```bash
PYTHONPATH=. pytest -v
```
Runs all unit, physical invariant, real-data ingestion, and backend API tests.

### 6.4 Sector 110 Scheduler CLI
```bash
# Check real-time countdown to Sector 110 release
PYTHONPATH=. python scripts/tic4206066_scheduler.py --status

# Execute full dry-run simulation
PYTHONPATH=. python scripts/tic4206066_scheduler.py --dry-run
```

### 6.5 Launch Local Web Portal
```bash
uvicorn app:app --host 0.0.0.0 --port 8080
```
Open `http://localhost:8080` to interact with the Stitch-designed dashboard, live countdown ticker, and verification console.

---

## 7. Google Cloud Run Deployment

Deploy the dashboard to Google Cloud Run in one command:

```bash
# Using automated script
./deploy_cloud_run.sh

# Or directly using gcloud CLI
gcloud run deploy tic4206066-portal \
  --source . \
  --region us-central1 \
  --platform managed \
  --allow-unauthenticated \
  --memory 1Gi \
  --cpu 1 \
  --port 8080
```

See [`docs/GCP_DEPLOYMENT.md`](docs/GCP_DEPLOYMENT.md) for detailed deployment options and Cloud Build CI/CD triggers.

---

## 8. Epistemic Provenance & Citations

- **Preprint:** Pavel Rabtsevich (2026), *Candidate Exoplanet TIC 4206066 b: A 1.4 Earth-Radius Planet in a 3.18-Day Orbit Around Late-K Dwarf StKM 1-561*, [doi:10.5281/zenodo.22967456](https://doi.org/10.5281/zenodo.22967456).
- **Data & Code Deposit:** [doi:10.5281/zenodo.22967099](https://doi.org/10.5281/zenodo.22967099).
- **Pre-Registration:** [doi:10.5281/zenodo.23175179](https://doi.org/10.5281/zenodo.23175179).
- **Observational Data:** NASA TESS Mission / Mikulski Archive for Space Telescopes (MAST STScI).

---

## 9. License

This project is licensed under the [MIT License](LICENSE).
Observational datasets retain their original open-access scientific licenses (CC-BY 4.0 / Public Domain NASA MAST).

## 10. Sector 110 blind-test procedure
```bash
python scripts/fetch_sector110.py --list      # exit 2 = AWAITING_DATA (nothing published yet)
python scripts/fetch_sector110.py --fetch     # downloads ONE file by frozen pipeline priority, records SHA-256
python scripts/fetch_sector110.py --evaluate  # frozen rules on the hashed file -> results/sector110_outcome.json
```
Thresholds (`subprojects/tic4206066_sector110_runner.py`) and pipeline priority (`PIPELINE_PRIORITY`) are frozen before data exists; do not change them after viewing the data. Open items: the observing-window and release dates are unverified (`WINDOW_VERIFIED = False`), and the pipeline priority must be confirmed against the Zenodo pre-registration.
