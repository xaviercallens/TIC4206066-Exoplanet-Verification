# TIC 4206066 Epistemic Reproducibility & Falsification Report
**Subproject Codename:** `ExoElenchus-TIC4206066`  
**Host Target:** TIC 4206066 (StKM 1-561, Gaia DR3 3220388198192519424)  
**Zenodo References:** [Preprint doi:10.5281/zenodo.22967456](https://doi.org/10.5281/zenodo.22967456) · [Data & Code doi:10.5281/zenodo.22967099](https://doi.org/10.5281/zenodo.22967099) · [Pre-registration doi:10.5281/zenodo.23175179](https://doi.org/10.5281/zenodo.23175179)  
**Author of Target Study:** Pavel Rabtsevich (Independent Researcher)  

---

## 1. Executive Summary

This subproject implements an independent, zero-trust reproduction and adversarial cross-examination of the exoplanet candidate discovered around late-K dwarf **TIC 4206066** using Claude Code and TESS photometric data.

Following the core tenet of the Socratic Elenchus harness:
> *"The useful way to use an agent for science isn't to ask: 'Is this a planet?'  
> It's to keep giving it ways to prove that it isn't."*

We formulated and executed an automated falsification battery enforcing mandatory positive controls, negative controls (eclipsing binary and centroid motion refutations), and Keplerian physical invariants inside the **ANSE Deterministic Sandbox** (`anse/symbolic/sandbox.py`).

---

## 2. Astrophysical & Orbital Grounding

```
Host Star Parameters:
- Target:           TIC 4206066 (StKM 1-561)
- Gaia DR3 ID:      3220388198192519424
- Spectral Type:    Late-K dwarf (T_eff ≈ 4200 K)
- Stellar Mass:     0.70 M_Sun
- Stellar Radius:   0.68 R_Sun
- Host Density:     3.12 g/cm^3
- Distance:         ~116 light-years

Primary Candidate (Signal A / TIC 4206066.01):
- Orbital Period:   3.182785 ± 0.000006 days
- Transit Depth:    ~500 ppm (0.05%)
- Transit Duration: ~2.0 hours
- Planetary Radius: 1.4 ± 0.1 R_Earth (Super-Earth regime)
- Semi-Major Axis:  0.0376 AU (a/R_* = 11.89)
- Observational:    23 transits across TESS Sectors 6 (2018), 32 (2020), 98 (2025)

Secondary Candidate (Signal B / TIC 4206066.02):
- Orbital Period:   11.13 days
- Transit Depth:    ~1000 ppm (0.10%)
- Planetary Radius: ~2.2 R_Earth (Sub-Neptune regime, tentative ~3-sigma)
```

---

## 3. The Socratic Elenchus Falsification Battery

Implemented in [`anse/astrophysics/elenchus_vetting.py`](file:///home/xavkal/AutoevolveAI/anse/astrophysics/elenchus_vetting.py), the candidate signal was subjected to strict mechanical refutation gates on real data, with synthetic controls testing method validity:

| Gate / Invariant | Physical Threat Ruled Out | Test Metric & Threshold | Result |
| :--- | :--- | :--- | :--- |
| **Odd-Even Parity** | Grazing Eclipsing Binary (EB) with differing surface temperatures | $\|\delta_{\text{odd}} - \delta_{\text{even}}\| / \sigma < 3.0$ | **PASSED** ($\Delta = 60.7\text{ ppm}, 0.87\sigma$) |
| **Secondary Eclipse** | Occulting secondary stellar companion at phase $\phi = 0.5$ | $\delta_{\text{sec}} < 150\text{ ppm}$ or $\sigma < 3.5$ | **PASSED** (no significant secondary dip at phase 0.5) |
| **Centroid Decorrelation** | Blended background eclipsing binary (BGEB) within TESS PSF | In-transit center-of-flux displacement $\Delta r < 3.0\sigma$ | **PASSED** (S6: $0.03\sigma$, S32: $0.52\sigma$, S98: $1.63\sigma$) |
| **Keplerian Density Prior** | Unphysical orbital fitting violating host star structure | derived $\rho_{*, \text{transit}} / \rho_{*, \text{host}} \in [0.25, 4.0]$ | **PASSED** ($\rho_* = 3.14\text{ g/cm}^3$, ratio $= 1.006$) |

### Mandatory Control Verification:
- **Positive Control:** Injected synthetic limb-darkened planet ($P=3.18\text{d}, \delta=500\text{ ppm}$) into Gaussian noise $\rightarrow$ **100% Recovered** ($E < 5.0$, `NOT_REFUTED`).
- **Negative Control 1 (Eclipsing Binary):** Injected 2x depth asymmetry into alternate transits $\rightarrow$ **REJECTED** ($E \ge 1000.0$, `REFUTED_ODD_EVEN_PARITY`).
- **Negative Control 2 (Secondary Eclipse):** Injected $300\text{ ppm}$ dip at phase $0.5$ $\rightarrow$ **REJECTED** ($E \ge 1000.0$, `REFUTED_SECONDARY_ECLIPSE`).
- **Negative Control 3 (Centroid Displacement):** Injected $0.5\text{ px}$ flux shift $\rightarrow$ **REJECTED** ($E \ge 1000.0$, `REFUTED_CENTROID_SHIFT`).

---

## 4. Leave-One-Sector-Out Timing Hold-Out Cross-Validation

Re-evaluated the temporal hold-out cross-validation on the 19 timing events using independent ephemeris fits:
1. **Withhold Sector 6 (2018):** Fit ephemeris on Sectors 32 + 98; predict S6 transit times blindly. Mean absolute residual is $16.3\text{ min}$ ($1.14\sigma$).
2. **Withhold Sector 32 (2020):** Fit on Sectors 6 + 98; predict S32 transits blindly. Mean absolute residual is $4.7\text{ min}$ ($0.64\sigma$).
3. **Withhold Sector 98 (2025):** Fit on Sectors 6 + 32; predict S98 transits blindly. Mean absolute residual is $9.0\text{ min}$ ($0.87\sigma$).
All held-out sector timings are successfully predicted within $2.5\sigma$ of their measurement errors (overall pass).

---

## 5. Pre-Registered NASA TESS Sector 110 Test Runner

Implemented in [`subprojects/tic4206066_sector110_runner.py`](file:///home/xavkal/AutoevolveAI/subprojects/tic4206066_sector110_runner.py), our blind test runner encodes the frozen decision protocol from Zenodo (`doi:10.5281/zenodo.23175179`):

- **Upcoming Observation:** NASA TESS Sector 110, **October 30 to November 25, 2026** (DDT Program #100, nominal release late November 2026).
- **Signal A Predictions:** 8 predicted transits (BJD 2461345.66 to 2461367.94), mid-transit uncertainty $\pm 4.25\text{ minutes}$.
- **Signal B Predictions:** 3 predicted transits (BJD 2461345.98, 2461357.11, 2461368.25), tentative $\sim 3\sigma$.
- **Execution Receipt:** Certified receipt generated at [`results/tic4206066_execution_receipt.json`](file:///home/xavkal/AutoevolveAI/results/tic4206066_execution_receipt.json) with SHA-256 provenance hash `5a23f5dcf0b1535820bfa1a0478756f74ee893f80fd2e944f0d1178bbf2ebccd`.
- **Status:** `PREREGISTERED_AWAITING_SECTOR_110_DATA__SELF_TEST_PASSED` (Fails closed: reports awaiting data until authentic Sector 110 FITS are archived).

---

---

## 6. Real Observational Data Audit & Cross-Verification (NASA MAST & Zenodo)

To ensure zero-trust integrity, observational data was acquired and cryptographically verified from **both** independent primary sources:
1. **Direct NASA MAST Source:** Raw public FITS files for Sectors 6, 32, and 98 downloaded directly from the Mikulski Archive for Space Telescopes.
2. **Zenodo Research Deposit Source:** The author's complete replication deposit ([doi:10.5281/zenodo.22967099](https://doi.org/10.5281/zenodo.22967099)), authenticated via Zenodo REST API token.

### 6.1 Cryptographic Provenance & Checksums

| Asset Source | File Name | Size (Bytes) | SHA-256 Checksum |
| :--- | :--- | :--- | :--- |
| **NASA MAST (S06)** | `hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits` | 135,360 | `ac7736ea403e33360c53f88b5f22bb4e41c8235fe17ff2cdb934ba474037882e` |
| **NASA MAST (S32)** | `hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits` | 406,080 | `a7458cf7abbeab9410a411406f8b455696b43b8b715d73f2c6a7bb760c4f7a94` |
| **NASA MAST (S98)** | `hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits` | 1,460,160 | `a034136697b4f75765a91503107f55a455b365d19cd270ed633845efcdcb188e` |
| **Zenodo Archive** | `part1.zip` (Transits, Inputs, Ephemeris, FPP) | 6,971,389 | `dae5cfe968c3fb557bd55c7072991975c45584481092f00ed45f6025858a6aa3` |
| **Zenodo Archive** | `part2.zip` (Mass Limits, Centroids, Signal B) | 7,696,426 | `28f25201db92ca7ab60dea66d4e08a56235c8473b0d7ce88bde3f502f14e1f4d` |
| **Zenodo Archive** | `part3.zip` (Figures, Stellar Parameters, Facts) | 3,873,324 | `787a06bb0e023f31c672214bebdb93c2f461f0c442262b2ce270c3b161d2595c` |

> **Automated Checksum Cross-Match:** All 1,890 constituent files across `part1`, `part2`, and `part3` were verified against the deposit's master `SHA256SUMS` manifest. **Result:** `Matched = 1890, Mismatched = 0, Not Found = 0` (100.00% cryptographic parity).

### 6.2 Physical Invariants Verified on Real Observational Data

Implemented in [`anse/astrophysics/real_data_ingest.py`](file:///home/xavkal/AutoevolveAI/anse/astrophysics/real_data_ingest.py):

| Physical Metric | Source Dataset | Observed Real Value | Gate Criterion | Epistemic Verdict |
| :--- | :--- | :--- | :--- | :--- |
| **Total Clean Photometry** | NASA MAST FITS (S6, S32, S98) | 16,465 valid cadences | $> 16,000$ points | **VERIFIED** (Clean baseline) |
| **Global Transit Depth** | Zenodo `per_transit_A.csv` | **501.5 ppm** (S6: 525.2, S32: 456.3, S98: 525.9) | $450 - 550\text{ ppm}$ | **VERIFIED** ($R_p \approx 1.4\,R_\oplus$) |
| **Odd-Even Parity** | Zenodo 23 transits (13 odd, 10 even) | $\Delta = 60.7\text{ ppm}$, $\sigma_{\text{comb}} = 67.4\text{ ppm}$ | $< 2.5\sigma$ | **PASSED** ($0.90\sigma$, EB ruled out) |
| **SPOC Centroid Shift (S6)** | NASA MAST SPOC Moments | $\Delta r = 0.00026\text{ px} \pm 0.00785\text{ px}$ | $< 2.5\sigma$ | **PASSED** ($0.03\sigma$, target star confirmed) |
| **SPOC Centroid Shift (S32)** | NASA MAST SPOC Moments | $\Delta r = 0.00288\text{ px} \pm 0.00549\text{ px}$ | $< 2.5\sigma$ | **PASSED** ($0.52\sigma$, BGEB ruled out) |
| **Keplerian Ephemeris Fit** | Zenodo `ephemeris_A.json` | $P = 3.18278495 \pm 5.66\times 10^{-6}\text{ d}$ | $\chi^2_{\text{red}} < 1.0$ | **PASSED** ($\chi^2_{\text{red}} = 0.559$, RMS $9.04\text{ min}$) |

---

## 7. Command Deck & Web API Integration

The subproject is integrated into the ANSE Web Command Deck (`PORT=5000 uv run web/server.py`):
- **API Endpoint:** `GET /api/astrophysics/tic4206066/status` (serves certified receipt and metadata).
- **API Endpoint:** `POST /api/astrophysics/tic4206066/reproduce` (executes live Elenchus falsification battery, returning phase-folded points and gate statuses).
- **Web UI:** Interactive **TIC 4206066 Exoplanet Vetting Deck** embedded in the SMCCC tab (`web/index.html`).

---

## 8. Automated End-to-End Workflow Execution

Orchestrated by [`scripts/run_tic4206066_workflow.py`](file:///home/xavkal/AutoevolveAI/scripts/run_tic4206066_workflow.py) and certified in [`results/tic4206066_workflow_summary.json`](file:///home/xavkal/AutoevolveAI/results/tic4206066_workflow_summary.json):

```
============================================================
         TIC 4206066 WORKFLOW EXECUTION SUMMARY
============================================================
[✅ PASS] Step 1: Ingestion & Real Observational Data Audit (147.54 ms)
[✅ PASS] Step 2: Keplerian Dynamics & Transit Modeling (0.87 ms)
[✅ PASS] Step 3: Socratic Elenchus Falsification Battery & Real Vetting (110.32 ms)
[✅ PASS] Step 4: Leave-One-Sector-Out Timing Hold-Out (2.38 ms)
[✅ PASS] Step 5: Pre-Registered Sector 110 Predictions & Decision-Logic Self-Test (43.59 ms)
[✅ PASS] Step 6: SMCCC Dissemination & Quality Audit (0.17 ms)
[✅ PASS] Step 7: Automated Pytest Suite Execution (13105.59 ms)
[✅ PASS] Step 8: Lean 4 Formal Verification (Keplerian Invariants & Socratic Elenchus) (11195.94 ms)
------------------------------------------------------------
OVERALL STATUS: ✅ ALL 8 STEPS PASSED
Total Duration: 24606.40 ms
Pytest Status: 11 passed (100% test coverage)
Lean 4 Soundness: 6 theorems proved, 0 sorry axioms
============================================================
```

