"""Shared helpers for the archival false-alarm test of signal A (TIC 4206066).
Times: BTJD = BJD_TDB - 2457000 (mid-exposure) in all three products.
Raw data are read in place (read-only); nothing is copied."""
import os
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))       # .../tess-jev
MAST = os.path.join(ROOT, 'addenda', 'opus55_final_audit_tic4206066', '00_data', 'mast')
FILES = {
    6: 'hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits',      # TESS-SPOC FFI, 30 min, PDCSAP
    32: 'tess2020324010417-s0032-0000000004206066-0200-s_lc.fits',             # SPOC 2 min, PDCSAP
    98: 'hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits',          # QLP 200 s, SYS_RM_FLUX
}
FLUXCOL = {6: 'PDCSAP_FLUX', 32: 'PDCSAP_FLUX', 98: 'SYS_RM_FLUX'}
EXPTIME = {6: 1800 / 86400, 32: 120 / 86400, 98: 200 / 86400}

# --- discovery record (candidates/tic4206066_p3.182945/candidate.json, written 00:09 CEST, before any archival look)
DISC_P = 3.182945191903069
DISC_T0 = 3998.7348383728113
DISC_DUR_H = 1.8
# --- reconciled facts used only for masking / building the null (FINAL_FACTS.csv)
T14A = 1.87 / 24
B_P, B_T0, T14B = 11.13274, 2461000.865 - 2457000.0, 1.66 / 24
FINAL_P, FINAL_T0 = 3.1827850, 2460998.73599 - 2457000.0      # A.P, A.T0 (multi-sector; used ONLY to mask true
                                                                # transits in null trials and for comparison)


def phase_dist(t, P, T0):
    return (t - T0 + 0.5 * P) % P - 0.5 * P


def load(sector):
    with fits.open(os.path.join(MAST, FILES[sector])) as h:
        d = h[1].data
        t = np.array(d['TIME'], float)
        f = np.array(d[FLUXCOL[sector]], float)
        q = np.array(d['QUALITY'], np.int64)
    fin = np.isfinite(t) & np.isfinite(f)
    good = fin & (q == 0)
    dumps = t[np.isfinite(t) & ((q & 32) != 0)]          # momentum-dump cadences (bit 6)
    return t[good], f[good] / np.nanmedian(f[good]), dumps


EDGE = 0.0          # wotan edge_cutoff (per segment); 0 in the headline analysis


def detrend(t, f, mask, window=0.75, edge=None):
    from wotan import flatten
    flat = flatten(t, f, method='biweight', window_length=window, break_tolerance=0.1,
                   edge_cutoff=EDGE if edge is None else edge, return_trend=False, mask=mask)
    return flat


def b_mask(t, half=None):
    half = (0.5 * T14B + 0.25 / 24) if half is None else half      # B in-transit + 15 min margin
    return np.abs(phase_dist(t, B_P, B_T0)) < half


def sigma_from_p(p):
    """one-sided Gaussian equivalent"""
    from scipy.stats import norm
    p = min(max(p, 1e-300), 1.0)
    return float(norm.isf(p))
