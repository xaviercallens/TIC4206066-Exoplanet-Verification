"""Audit 9: light-curve preparation (own code). Inputs: raw MAST products downloaded by the 00_data agent
(sha256 recorded in 00_data/logs/provenance.jsonl and re-hashed here). Output: out/lc_windows.npz with detrended
flux in windows around the 3.18 d (A) and 11.13 d (B) transits.
Choices: QUALITY == 0 for all products (strictest; removes QLP S98 scattered-light bits 29/30 and bit 11),
SPOC/TESS-SPOC PDCSAP_FLUX, QLP SYS_RM_FLUX (the product used in the manuscript for S98; DET_FLUX as variant),
wotan biweight window 0.75 d with both signals' transits masked (2 x T14 wide), empirical per-sector noise from
out-of-transit robust scatter."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from wotan import flatten

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "00_data" / "mast"
FILES = {6: ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR"),
         32: ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR"),
         98: ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX", "DET_FLUX_ERR")}
if len(sys.argv) > 1 and sys.argv[1] == "det":
    FILES[98] = (FILES[98][0], "DET_FLUX", "DET_FLUX_ERR")
PA, T0A, T14A = 3.1827856, 3998.7375, 1.87 / 24     # manuscript ephemeris (BTJD = BJD - 2457000)
PB, T0B, T14B = 11.13274, 4000.865, 1.66 / 24
WIN = 0.30  # half-width of the fitting window around each transit [d]

def phase_dist(t, P, T0):
    return (t - T0 + 0.5 * P) % P - 0.5 * P

out, meta = {}, {}
for sec, (fn, fcol, ecol) in FILES.items():
    p = DATA / fn
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    h = fits.open(p); d = h[1].data
    t = np.array(d["TIME"], float); f = np.array(d[fcol], float)
    q = np.array(d["QUALITY"]).astype(np.int64)
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[ok], f[ok]
    f = f / np.nanmedian(f)
    inA = np.abs(phase_dist(t, PA, T0A)) < T14A          # +-1 T14 (i.e. 2 x T14 wide)
    inB = np.abs(phase_dist(t, PB, T0B)) < T14B
    trend = flatten(t, np.where(inA | inB, np.nan, f), method="biweight", window_length=0.75,
                    return_trend=True, edge_cutoff=0.0, break_tolerance=0.5)[1]
    # wotan returns NaN trend at masked points only if all neighbours are NaN; interpolate masked points
    good = np.isfinite(trend)
    trend = np.interp(t, t[good], trend[good])
    fd = f / trend
    oot = ~(inA | inB)
    mad = 1.4826 * np.median(np.abs(fd[oot] - np.median(fd[oot])))
    exp = float(h[1].header.get("TIMEDEL"))
    # windows
    wA = (np.abs(phase_dist(t, PA, T0A)) < WIN) & ~inB
    wB = (np.abs(phase_dist(t, PB, T0B)) < WIN) & ~inA
    out[f"s{sec}_tA"], out[f"s{sec}_fA"] = t[wA], fd[wA]
    out[f"s{sec}_tB"], out[f"s{sec}_fB"] = t[wB], fd[wB]
    out[f"s{sec}_sig"] = np.array([mad]); out[f"s{sec}_exp"] = np.array([exp])
    nA = len(np.unique(np.round((t[wA] - T0A) / PA)))
    nB = len(np.unique(np.round((t[wB] - T0B) / PB)))
    meta[sec] = dict(file=fn, sha256=sha, flux_col=fcol, n_q0=int(ok.sum()), sigma_oot_ppm=mad * 1e6,
                     exptime_d=exp, nA_windows=nA, nB_windows=nB, npts_A=int(wA.sum()), npts_B=int(wB.sum()),
                     crowdsap=h[1].header.get("CROWDSAP"))
    print(sec, meta[sec])
tag = "" if len(sys.argv) < 2 else "_" + sys.argv[1]
np.savez(HERE / f"out/lc_windows{tag}.npz", **out)
json.dump(meta, open(HERE / f"out/lc_windows{tag}_meta.json", "w"), indent=1)
