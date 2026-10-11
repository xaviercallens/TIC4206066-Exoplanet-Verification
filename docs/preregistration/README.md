# TIC 4206066: pre-registration of the TESS Sector 110 test

Author: Pavel Rabtsevich, independent researcher, Tarragona, Spain.
This record: doi:10.5281/zenodo.23175179.
Related: preprint doi:10.5281/zenodo.22967456; supporting data and code doi:10.5281/zenodo.22967099 (version 1.0.0,
cited in the preprint; all versions: doi:10.5281/zenodo.22967098).

## What this is, in two minutes

**The signals.** Archival TESS photometry of the nearby K dwarf TIC 4206066 shows two transit-like signals: small,
periodic dips in brightness. *Signal A* repeats every 3.18 days (expected box depth about 440 ppm, i.e. 0.044 %);
*signal B* every 11.13 days (about 590 ppm). Both were found in TESS Sector 98. Signal A is also seen in the older
Sectors 6 and 32; signal B has one counted transit in Sector 6 and none in Sector 32. Neither has been shown to be a planet.

**What happens next.** TESS observes the star again in Sector 110, from 2026-10-31 to 2026-11-26, at 2-min cadence
(Director's Discretionary Target request DDT 100). These are new data that played no part in finding or fitting the
signals.

**Why a pre-registration.** If the analysis of new data is chosen after looking at them, almost any result can be
made to look like a confirmation. This record fixes everything before the Sector 110 observations start: when the
transits must appear, how deep they must be, how the data will be processed, and which results count as *recovered*,
*inconsistent*, *excluded* or *insufficient data*. Frozen code computes outcomes from valid inputs; the registered cases without a usable output require fixed
sentences written from the posted evidence. Rules about
valid inputs and attempts that the code does not check are written out in Appendix G of the document; a result that
breaks any pre-registered rule is invalid.

**What it can show.** Whether the predicted dips reappear at the predicted times (9 times for A, 3 for B, each known
to about 4 min) and depths in independent data, and whether the dip of A comes from the pixels of TIC 4206066.
**What it cannot show.** That either signal is a planet. A recovery is also expected from an eclipsing companion star
too close to resolve; ruling that out needs high-resolution imaging and radial velocities.

**After Sector 110.** The outputs of the frozen scripts and their sha256 will be posted as new versions of this
Zenodo record. The rules in this version cannot be changed afterwards; any later idea is a new, separately labelled
analysis.

## What to read

| you want | file |
|---|---|
| the pre-registration itself | `PREREG_S110_SHORT.pdf` (short overview, start here) and `PREREG_S110_FINAL.pdf` (complete, governs); their LaTeX source is in `protocol_source/` (main file `PREREG_S110_FINAL.tex`) |
| the predicted times | Section 3.2 of the document; `ephemeris_final.json` |
| the expected depths | Section 3.3 of the document; `calib/expected_depths.json` |
| the decision rules | Sections 4 and 5 and Appendices A to I of the document; in code: `s110_final.py`, `s110_localise.py`, `s110_summary.py` |
| where every number comes from | `calib/` (inputs of the ephemeris in `calib/transit_times_A.csv`, `calib/ephemeris_inputs.json`, `calib/ephemeris_sources/`) |
| what every file is | `FILE_GUIDE.md` |
| how the rules evolved | `CHANGELOG.md` and Appendix L of the document (the unpublished drafts v1 and v2 are not part of this deposit) |

Times: the analysis uses only BTJD = BJD_TDB - 2457000, the time stamps of the TESS files. The document also lists,
for convenience, the arrival time of each transit at the Earth in UTC (5.6 to 7.2 min earlier than the BJD_TDB date
written as a calendar date).

## The Zenodo record and the ZIP

The record holds five files: `PREREG_S110_SHORT.pdf` (short overview), `PREREG_S110_FINAL.pdf` (the complete
pre-registration, which governs), this `README.md`, `SHA256SUMS.txt` and `tic4206066_S110_prereg_final.zip`. The ZIP is the complete frozen deposit (a folder `tic4206066_S110_prereg_final/`
with the document and its source, code, predictions, calibrations, validation runs and this README); the two PDFs and
the README outside it are identical copies for reading without unpacking.

Checksums: `SHA256SUMS` inside the ZIP lists the sha256 of every file of the deposit except itself. `SHA256SUMS.txt`
on the record lists the four other uploaded files and every file inside the ZIP, so that after unpacking next to it
`shasum -a 256 -c SHA256SUMS.txt` checks everything. The frozen version is version 1.1 of the protocol: exactly the
files listed in `SHA256SUMS`.

## Reproducing the numbers

Requirements: Python 3.12.13 with the packages of `requirements.txt` (numpy 2.5.3, astropy 8.0.1, scipy 1.18.1,
wotan 1.10, TESS_PRF 0.1.3, numba 0.67.0, llvmlite 0.49.0), on an Apple Silicon Mac (arm64, macOS). Every registered run
refuses to start in another environment, so that the Sector 110 analysis runs exactly as tested. `requirements_lock.txt`
also fixes the indirect dependencies, with distribution hashes. On another platform, edit `requirements.txt` in a copy
of the deposit; the versions actually used are written into every output.

Steps (tested on 2026-10-03 from a fresh copy with freshly downloaded MAST files: every number, trial list and fixed sentence reproduced exactly; 25 to 40 min):

1. `uv venv --python 3.12.13 venv` and `uv pip install --python venv/bin/python --require-hashes -r requirements_lock.txt`
   (or `-r requirements.txt`).
2. `zsh build/fetch_archival_inputs.sh`: downloads the five public archival TESS files (Sectors 6, 32, 98) and the 25
   TESS PRF model files from MAST, puts them next to the deposit where the frozen scripts expect them, and checks their
   sha256. These public files are not included in the deposit.
3. `zsh build/reproduce.sh venv/bin/python`: copies the deposit to `../s110_rebuild/`, rebuilds the expected depths,
   every dry run, the localisation runs, the coverage check and the fixed sentences there, and compares each with the
   deposited file (numbers, trial lists and sentences; run timestamps and output hashes always differ). It never writes
   into the deposit.

Not repeated by these steps:
- The false-alarm trials (hours of CPU): from the deposit folder, run
  `python build/false_alarm_check.py MODE calib/false_alarm_MODE.json NPROC` for MODE = A, B, A_S32, B_S32, A_S32_a3,
  B_S32_a3, then `python build/false_alarm_levels.py`. The counts and bounds are exact; the supplementary confidence
  intervals depend on the order of the trial sets in each file, which follows the order in which parallel jobs finish.
- The archival transit-time measurements and the two posterior fits behind the ephemeris. Their outputs are inputs
  here (`calib/transit_times_A.csv`, `calib/ephemeris_sources/`); `calib/ephemeris_sources/README_sources.md` shows how
  every ephemeris number follows from them.

Note: the field `provenance.expected_depths_sha256` inside `calib/expected_depths.json` is the hash of an earlier
version of that file (it existed when the calibration was rerun); every other output carries the hash of the deposited
file.

## Licences and AI use

Code (`*.py`, `*.sh`): MIT. Everything else: CC BY 4.0. See `LICENSE`; cite as in `CITATION.cff`.
Claude (Anthropic) assisted code and text drafting; Codex (OpenAI) assisted editing and verification.
The author is responsible for the protocol. Internal package checks were separate AI review sessions,
not human peer review. Frozen analysis code is unchanged in v1.1.

## Acknowledgements

This deposit uses data collected by the TESS mission, obtained from the MAST data archive at the Space Telescope Science
Institute (STScI); funding for the TESS mission is provided by NASA's Science Mission Directorate. It uses the TESS-SPOC
(doi:10.17909/t9-wpz1-8s54) and QLP (doi:10.17909/t9-r086-e880) high-level science products and the TESS PRF models
from MAST. This work has made use of data from the European Space Agency (ESA) mission Gaia
(https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC,
https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions,
in particular the institutions participating in the Gaia Multilateral Agreement (Gaia Collaboration et al. 2016,
2023: the Gaia mission; Gaia Data Release 3).
