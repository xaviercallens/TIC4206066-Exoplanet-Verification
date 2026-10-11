# Audit 10_novelty_protocol — REPORT

> Rendered verbatim from the agent's structured result (`_structured/audit__10_novelty_protocol.json`); the sub-agent harness blocked agents from writing .md report files. Scripts, outputs and logs are in this folder.


**Verdict:** PASS_WITH_CHANGES

## Summary

Verdict: PASS_WITH_CHANGES. Every checked number reproduces and the novelty search is complete; the only problems are three small wording fixes (all MINOR).

Deliverables: REPORT.md was not written. The harness blocked it ("Subagents should return findings as text, not write report files"), so the report is this structured output. Scripts and outputs are in /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/.

TASK A (novelty)

Earlier searches:
- The sealed priority re-check (04:46-05:25Z) found nothing, but has two gaps:
  - Its arXiv API calls failed (HTTP 406/429, 0-byte files), so arXiv was covered only by web search on 5 identifiers, with no period strings.
  - ADS was searched through the UI, and only the log lines survive.
- novelty_refresh (10:59Z) re-queried only ExoFOP TOI/CTOI, the ExoFOP target page, NASA toi and exo.MAST. All were empty.

This audit's fresh queries (11:33-11:49Z; every query logged with UTC time, URL, status and sha256 in queries.log, raw responses in raw/):
- ExoFOP TOI list (8,148 rows): no row for TIC 4206066 and nothing within 2'. The nearest TOI is TOI-2451.01 at 43.0'. The 34 period matches (0.1 %, for P/3, P/2, P, 2P, 3P of both signals) are all more than 1 deg away; the closest is TOI-2521.01 (= B/2) at 14.3 deg.
- ExoFOP CTOI list (5,137 rows): no row for the target; the nearest CTOI is at 120'.
- ExoFOP target JSON: tois and ctois empty; imaging, spectroscopy, time series and files all empty.
- exo.MAST DV service: {"TCE":[]}.
- SPOC TCE bulk tables (S6, S32 and S98 single-sector; s0001-s0036/0039/0069/0096 multi-sector): 0 rows for the target and 0 for any of the 49 TIC stars within 2'.
- MAST products: no DV or TCE products for the star, and no TESS-SPOC S98 products.
- NASA Exoplanet Archive: 0 rows in toi, ps and pscomppars, by TIC, by Gaia DR3 id and in a 2' box. The positive control (TOI-700) returns 4 rows.
- SIMBAD: all three identifiers and the 30" cone resolve only to StKM 1-561 (a proper-motion star). Its 10 references are all catalogues or target lists.
- VizieR (166 tables within 15"): no transit, eclipsing-binary or planet table. The only periods are variability entries (8.18, 8.11, 7.98 d, and an ellipsoidal Porb of 4.57 d with Score 0), and none relates to the signal periods.
- T16 (Dataverse v3) and RAVEN (Zenodo, 4 tables): no match. The nearest T16 candidate is 30.2' away; the nearest RAVEN entry is 1.95 deg away.
- arXiv API (run via curl, because python-requests got HTTP 406): 0 hits for 6 identifiers and 15 period strings (each with TESS). A probe confirms numeric phrases are searchable: it finds TOI-908, a different star.
- ADS (UI rendered in the browser pane; no API token was available): 0 hits for every identifier. Positive controls work (1 and 8 hits). The object search returns the 10 SIMBAD catalogue papers. Period searches return only unrelated papers.
- General WebSearch: nothing about the target.

Status sentence: every source it names was actually searched, and its wording matches the strongest permitted form ('no prior public report identified in the searched sources'). One problem: the T16 and RAVEN lists are named without citations.

TASK B (protocol)

What the frozen records and my own recount show:
- Config v1.1 (sha256 80e7ae48..., frozen 20:28:11Z) is the version both production runs started under.
- Input set: 126,246 QLP Sector 98 stars with TIC Tmag <= 11.5. I reproduced this exactly from a fresh copy of the QLP S98 target list plus fresh TIC 8.2 magnitudes for all 628,730 IDs. Of these, 126,070 light curves were searchable (175 corrupt, 1 processing error).
- Search chain: QLP systematics-removed flux, wotan biweight detrending (1 d window), astropy box least squares, then leo-vetter 1.2.0.
- Pool after hard rejects: 4,749 signals. Tier A (all 13 LEO false-alarm tests passed): 141.
- MES ladder on tier A: 125 signals at MES >= 10, 95 at >= 12, 70 at >= 15. The frozen rule takes the lowest rung giving <= 80, so T* = 15, PRIMARY = 70, EXPLORATORY = 71. The ID sets match the manifests exactly.
- The value 15 was not preregistered. The ladder rule was fixed in advance and produced 15 on the final data; amendment 01 says this explicitly. The freeze was a local sha256 record, not a public registration.

TIC 4206066:
- Both signals are in the EXPLORATORY cohort (protocol amendment 01, hash verified), with reason SHORTLIST_MES_BELOW_15.
- MES is 10.597 (3.18 d) and 10.982 (11.13 d). These are LEO-Vetter signal-to-pink-noise values, not SPOC TPS MES.
- Both passed all 13 false-alarm tests and all 4 flux-level false-positive tests (odd_even, vshaped, large, secondary). No LEO test is flagged.
- Other flags (not LEO tests): COMMON_MODE_HOT_FRACTION_0.12 on the 3.18 d signal, and a context flag, GAIA_NEIGHBOR_BRIGHT_ENOUGH_FOR_DEPTH, on both.
- In LEO-Vetter, 'flux-level' covers both the false-alarm and the false-positive tests (compute_flux_metrics feeds both; the paper title is 'LEO-Vetter: Fully Automated Flux- and Pixel-level Vetting'). So 'passed all LEO-vetter flux-level tests' is accurate.

Independent re-run (own detrending and BLS, public leo-vetter 1.2.0, fresh TIC parameters):
- My own BLS recovers both periods.
- In the pipeline's order of processing, MES = 10.602 and 10.975, and all 17 tests pass for both signals.
- At my own BLS ephemeris, the 11.13 d signal fails the 'chases' false-alarm test (mean_chases 0.74-0.76 < 0.78). Its tier-A membership is therefore ephemeris-sensitive.

The Note does not imply the target came from the confirmatory selection.

Environment: installed leo-vetter==1.2.0 into audit_env with uv at 2026-09-25T11:36:52Z. This also installed asteval, brokenaxes, dill and lmfit.

Limitations:
- ADS counts were read from the rendered page; the HTML was not archived.
- The arXiv API searches metadata only, not full text.
- The first arXiv attempt (HTTP 406) and a 'lines=0' logging bug are both annotated in queries.log.

## Numeric rows

| quantity | manuscript | recomputed | status | method / note |
|---|---|---|---|---|
| Number of QLP S98 stars searched with T<=11.5 | 126,246 | 126,246 (fresh QLP S98 target list, 628,730 IDs, and fresh TIC 8.2 Tmag via MAST); 126,070 searchable (175 DATA_CORRUPT, 1 PROCESSING_ERROR) | MATCH | tic_recount.py; data/vetting/targets_processed.csv.gz status counts 126,246 is the input set; 0.14% of light curves were unusable |
| Magnitude limit | T <= 11.5 | TIC 8.2 Tmag <= 11.5 (config v1.1 subset.tmag_max = 11.5; target Tmag 9.6817) | MATCH | config v1.1 copy (sha256 80e7ae48...) and TIC recount  |
| Search data and cadence | Sector 98 QLP light curves, 200 s | hlsp_qlp_tess_ffi_s0098 LLC, TIMEDEL = 200.0 s; sha256 a0341366... identical to the pipeline record | MATCH | fresh MAST download (leo_rerun.py)  |
| Search method | box least squares after biweight detrending | astropy BoxLeastSquares after wotan biweight (1.0 d window); my own re-run recovers 11.1334 d and 3.18273 d | MATCH | config v1.1 search and detrend sections; leo_rerun.py  |
| MES of the 3.18 d signal | 10.6 | pipeline 10.597; independent re-run 10.602 | MATCH | leo_rerun.py (own detrending plus leo-vetter 1.2.0; 11.13 d transits masked, as in the pipeline) LEO-Vetter signal-to-pink-noise MES, not SPOC TPS MES |
| MES of the 11.13 d signal | 11.0 | pipeline 10.982; independent re-run 10.975 | MATCH | leo_rerun.py (unmasked light curve, as in the pipeline)  |
| LEO-Vetter flux-level tests passed | all, for both signals | 13/13 false-alarm plus 4/4 false-positive flux tests passed for both (pipeline and independent re-run at the pipeline ephemerides) | MATCH | candidate.json _tests_failed; signals_all leo_fail_* columns; leo_rerun.py At this audit's own BLS ephemeris the 11.13 d signal fails the 'chases' test (mean_chases 0.74-0.76 < 0.78) |
| Shortlist MES threshold | MES >= 15 cut of the preregistered shortlist | T* = 15, recomputed from the frozen ladder [10,12,15,...] with target_max 80 on the final tier-A pool (125, 95 and 70 signals at 10, 12 and 15) | MINOR_DIFF | protocol_recount.py The value is right, but the rule derived it from the final data; 15 itself was not preregistered, and the freeze was local, not a public registration |
| Cohort of both signals | post hoc exploratory tier | EXPLORATORY (tier A, MES < T* = 15; amendment 01 frozen 2026-09-25T00:19:47Z, hash verified); neither signal is in PRIMARY | MATCH | protocol_recount.py; handoff manifests  |
| Tier-A / PRIMARY / EXPLORATORY sizes | not stated | 141 / 70 / 71 (identical to the manifests, 0 ID differences) | NOT_CHECKED | protocol_recount.py Context only; not a manuscript number |
| Novelty search date | 2026 September 25 (UTC) | sealed re-check 04:46-05:25Z; refresh 10:59Z; this audit 11:33-11:49Z, all on 2026-09-25 UTC | MATCH | queries.log  |
| Named sources searched | ExoFOP, SPOC TCE, NASA EA, T16, RAVEN, SIMBAD, VizieR, ADS, arXiv | all 9 searched (sealed re-check and/or this audit); 0 reports of either signal | MATCH | fresh_queries.py; ads_queries.json; WebSearch The arXiv API gap in the sealed re-check (406/429, no period strings) is closed here |
| ExoFOP TOI/CTOI entries for the target or within 2 arcmin | none (implied) | 0 TOIs (nearest TOI-2451.01 at 43.0'); 0 CTOIs (nearest at 120.1'); 0 period matches within 1 deg | MATCH | fresh_queries.py  |
| SPOC TCEs for the target or TIC stars within 2 arcmin | none (implied) | 0 in the S6/S32/S98 single-sector and s0001-s0036/0039/0069/0096 tables (49 TIC within 2'); exo.MAST {"TCE":[]}; no DV products on MAST | MATCH | fresh_queries.py; MAST product query  |
| NASA Exoplanet Archive rows | none (implied) | 0 in toi, ps and pscomppars by TIC, Gaia id and 2' box (positive control returns 4 rows) | MATCH | fresh_queries.py  |

## Findings

### F1 [MINOR] 'MES >= 15 cut of the preregistered shortlist' overstates what was preregistered

The threshold 15 was not preregistered. The pre-specified rule (tier A, then the lowest rung of the MES ladder [10,12,15,20,...] leaving <= 80 signals) produced 15 on the final data. The rule was frozen as a local sha256 record, not publicly registered, so 'preregistered' suggests a registration readers cannot find. The Note does not imply the target came from the confirmatory selection.

*Evidence:* config/prejev_search.json v1.1 copy (sha256 80e7ae48..., the hash in the search_started events at 20:28:20Z and 21:01:50Z): shortlist.mes_ladder = [10,12,15,...], target_max = 80, no fixed threshold. config/protocol_amendment_01.json statement 1: 'T* is whatever that rule yields on the final data (a preview on ~90% of the data suggested 15; the preview value is not used)'. Recount (protocol_recount.json): tier A = 141; MES>=10: 125, >=12: 95, >=15: 70, so T* = 15.

*Location:* Sect. 1 Target and Data, 2nd sentence of the search description

*Proposed fix:* Replace with: 'missed the MES $\ge15$ threshold set by my pre-specified shortlist rule; I examined them post hoc' (see recommended_wording).

### F2 [MINOR] T16 and RAVEN lists are named without citations

The Status sentence names 'the T16 and RAVEN lists', but neither is cited and neither is in references.bib, so readers cannot identify them. The word budget is tight (texcount 1488/1500), so adding two references means trimming elsewhere.

*Evidence:* manuscript.tex Status section; grep of references.bib finds no Roth or Lafarga entry. ADS: 2026ApJS..284...19R (Roth et al., 'The T16 Planet Hunt'); 2026MNRAS.548ag512L (Lafarga et al., RAVEN; arXiv:2603.22597). Both lists were searched: 0 matches (T16 nearest 30.2'; RAVEN nearest 1.95 deg).

*Location:* Sect. 5 Status, last sentence

*Proposed fix:* Cite Roth et al. 2026 (ApJS 284, 19) and Lafarga et al. 2026 (MNRAS, arXiv:2603.22597) and trim elsewhere to absorb the reference words; otherwise drop 'the T16 and RAVEN lists'.

### F3 [MINOR] 'Both passed' has no noun

Sect. 1 never introduces the two signals before 'Both passed all LEO-vetter flux-level tests'; the reader has to infer them from the abstract.

*Evidence:* manuscript.tex Sect. 1: '...after biweight detrending \citep{Hippke2019}. Both passed all LEO-vetter flux-level tests...'

*Location:* Sect. 1

*Proposed fix:* Write 'Both signals passed ...'.

### F4 [NOTE] 'Passed all LEO-Vetter flux-level tests' verified; one pass is ephemeris-sensitive

In LEO-Vetter, 'flux-level' covers the 13 false-alarm tests and the 4 false-positive tests (odd_even, vshaped, large, secondary); the pixel-level offset test is separate. Both signals pass all 17, which is more than the selection required (13 false-alarm tests only). An independent re-run reproduces this at the pipeline ephemerides. At this audit's own BLS ephemeris for the 11.13 d signal (P = 11.133407 d, T14 = 1.5 h), however, it fails the 'chases' test (mean_chases 0.74-0.76 < 0.78), so its tier-A membership is fragile. The pixel-level test was not run, so neither signal is a LEO 'PC'; the manuscript does not claim it is.

*Evidence:* candidate.json leo_vetter._tests_failed: all 17 false for both signals; data/vetting/signals_all.csv.gz leo_fail_* all False. leo_rerun.json: A (pipeline ephemeris, B masked): MES 10.602, FA_failed [], FP_failed []; B (pipeline ephemeris, unmasked): MES 10.975, [] and []; B at own BLS ephemeris: FA_failed ['chases']. leo-vetter thresholds.check_thresholds adds the 'offset' test only when offset_qual is present. Paper title (2025AJ....170..280K): 'LEO-Vetter: Fully Automated Flux- and Pixel-level Vetting of TESS Planet Candidates...'.

*Location:* Sect. 1

*Proposed fix:* No change required.

### F5 [NOTE] MES is the LEO-Vetter signal-to-pink-noise statistic; values reproduced

The 'multiple-event statistic' is LEO-Vetter's MES (depth over pink noise, Pont et al. 2006), which is how the search config defines it; it is not SPOC TPS MES. The quoted values reproduce independently.

*Evidence:* config v1.1 detection.statistic: 'MES = signal-to-pink-noise ratio, identical definition to leo_vetter get_SES_MES'. Pipeline 10.597 / 10.982; independent 10.602 / 10.975 (leo_rerun.json).

*Location:* Sect. 1

*Proposed fix:* Optional: add '(LEO-Vetter definition)' if the word budget allows.

### F6 [NOTE] 126,246 is the input set; 126,070 light curves were searchable

The count matches TIC 8.2 Tmag <= 11.5 among the 628,730 QLP S98 targets exactly. Of these, 175 light curves were corrupt (< 6,000 good cadences) and 1 failed in processing.

*Evidence:* tic_recount.json: n_tmag_le_11p5 = 126246; data/vetting/vetting_summary.json: OK 126070, DATA_CORRUPT 175, PROCESSING_ERROR 1.

*Location:* Sect. 1

*Proposed fix:* Acceptable as written; optionally write '126,246 stars (126,070 usable)'.

### F7 [NOTE] Novelty claim verified; wording within the permitted strength

Every source named in the Status sentence was searched on 2026-09-25 UTC by the sealed re-check and/or this audit, and no report of either signal exists in them. 'I found no prior public report of either signal in [listed sources]' is equivalent to the permitted 'No prior public report was identified in the searched sources'. Caveats: SPOC 2-min TCEs could only have come from S32; no TESS-SPOC FFI search exists for S6/S32 (those runs start at S36), and there are no S98 TESS-SPOC products; the arXiv API searches metadata only.

*Evidence:* queries.log, fresh_results.json, novelty_summary.json, ads_queries.json (per-source results listed in the summary); MAST product list for the star shows no dvr/dvs/dvt products; the MAST TCE bulk page lists 131 SPOC 2-min tables and no FFI TCE tables.

*Location:* Sect. 5 Status

*Proposed fix:* Keep the sentence. Immediately before submission, re-run ExoFOP TOI/CTOI, NASA toi, exo.MAST and SPOC TCE S107+ and update the date.

### F8 [NOTE] Sealed re-check arXiv gap closed

The sealed priority re-check's arXiv API calls all returned HTTP 406/429 (0-byte files), so it searched only 5 identifiers through the web UI, with no period strings. This audit ran the API via curl for 6 identifiers and 15 period and alias strings (each with TESS): 0 results. A numeric probe (abs:"3.18" AND abs:TESS) returns TOI-908 and TOI-4336 A/TOI-4342, both different stars.

*Evidence:* priority_recheck_tic4206066_SEALED/queries.log lines 'arXiv API ... status=406/429'; this audit's queries.log arXiv lines at 11:40-11:44Z and raw/arxiv_*.xml.

*Location:* 

*Proposed fix:* No change required.

### F9 [NOTE] Spelling of the tool name

The official name is 'LEO-Vetter' (paper title); the manuscript writes 'LEO-vetter', which is the GitHub repository name.

*Evidence:* ADS 2025AJ....170..280K title.

*Location:* Sect. 1

*Proposed fix:* Write 'LEO-Vetter'.

### F10 [NOTE] Outside this audit's scope: Gaia-neighbour depth flag is consistent with the manuscript in the TESS band

The pipeline's context flag GAIA_NEIGHBOR_BRIGHT_ENOUGH_FOR_DEPTH comes from a G-band estimate for Gaia DR3 3220388296975929984 at 11.1" (maximum depth 509 ppm). In the TESS band (TIC 712986249, T = 18.13 against 9.68 for the target) the maximum depth for a total eclipse is ~417 ppm, below the ~500 ppm depth of the 3.18 d signal. This agrees with the manuscript's 'All Gaia DR3 stars able to produce the depth lie >= 43" away'.

*Evidence:* context_decisions.jsonl (tic4206066_p3.182945) gaia_neighbors_bright_enough; raw/mast_tic_cone_2arcmin.csv.

*Location:* 

*Proposed fix:* No change; pointer for the pixel/localization auditor.

## Recommended wording

Sect. 1, sentence after the search description (currently 32 words):

Option A (31 words): "Both signals passed all LEO-Vetter flux-level tests \citep{Kunimoto2025} but, with multiple-event statistics of 10.6 and 11.0, missed the MES $\ge15$ threshold set by my pre-specified shortlist rule; I examined them post hoc."

Option B (24 words, shortest): "Both signals passed all LEO-Vetter flux-level tests \citep{Kunimoto2025}; their multiple-event statistics (10.6, 11.0) fell below the shortlist threshold (15), so I examined them post hoc."

Status sentence: the wording can stay as is; it is within the permitted strength. Either add citations, e.g. "...the T16 \citep{Roth2026} and RAVEN \citep{Lafarga2026} candidate lists..." (2026ApJS..284...19R; 2026MNRAS.548ag512L) and trim elsewhere to stay within 1500 words, or drop "the T16 and RAVEN lists". Strongest permitted form: "No prior public report was identified in the searched sources." Refresh the date after re-running ExoFOP, the NASA toi table, exo.MAST and new SPOC TCE sectors immediately before submission.

## Files

- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/fresh_queries.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/tic_recount.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/protocol_recount.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/leo_rerun.py
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/queries.log
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/fresh_results.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/novelty_summary.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/ads_queries.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/tic_recount.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/protocol_recount.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/leo_rerun.json
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/numeric_rows.csv
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/raw/ (70 raw responses + provenance.jsonl)
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/tic_batches/ (126 TIC Tmag batch CSVs)
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/run_fresh_queries.log
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/run_fresh_queries_arxiv.log
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/tic_recount.log
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/leo_rerun.log
- /Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/10_novelty_protocol/protocol_recount.stdout