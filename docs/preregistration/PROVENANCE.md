# Pre-registration provenance

Zenodo doi:10.5281/zenodo.23175179, "A Pre-registered TESS Sector 110 Test of Two Transit-like Signals in TIC 4206066", v1.1 (2026-10-06).
Verified sha256: `PREREG_S110_FINAL.pdf` = 1d23034d…cb98; `tic4206066_S110_prereg_final.zip` (not copied) = d1ef51d5…f6f21. See `SHA256SUMS.txt`.

## The registered analysis is NOT this repo's runner
- The frozen code (`s110_final.py`, in the Zenodo zip) defines outcomes; it requires a pinned arm64-macOS environment. A run on Linux is only an "additional, separately labelled analysis".
- Primary product: first SPOC 2-min light curve, `PDCSAP_FLUX`, quality mask 17087. Fallbacks (QLP, TESS-SPOC) only if no SPOC 2-min product is public on 2027-02-28.
- No Sector 110 product may be opened before the primary run.

| Item | Registered | `subprojects/tic4206066_sector110_runner.py` |
|---|---|---|
| Signal A transits | 9 (epochs 109-117) | 8 (misses epoch 117) |
| Timing | `ephemeris_final.json` | ~1-2 min offset |
| Expected depth A / B | ~440 / ~590 ppm | 498 ppm single value |
| Decision rules | frozen in `s110_final.py` | own thresholds, not registered |

The runner is an illustrative self-test only. A recovery would not prove a planet.
