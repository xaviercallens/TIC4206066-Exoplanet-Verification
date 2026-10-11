"""PRF helpers for the audit. ONLY the official MAST PRF models (start_s0004, downloaded to data/prf/) are used,
read through TESS_PRF 0.1.3 (Bell & Higgins 2022, ascl:2207.008), which un-interleaves each 117x117 file into 9x9
sub-pixel PRFs of 13x13 pixels (prf[8-i::9, 8-j::9]). No superseded PRF code from any earlier addendum is imported.
Regression tests for the TIC 158955499 PRF mistakes live in tests/test_prf_regressions.py.

Coordinates: stamp pixel indices (0-based, integer = pixel centre), as returned by astropy WCS all_world2pix(.., 0);
TESS_PRF.locate() uses the same convention ("integers refer to the pixel center").
"""
import csv
import os

import numpy as np
from astropy.time import Time

from PRF import TESS_PRF

AUD = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PRFDIR = os.path.join(AUD, "data", "prf")


def get_prf(cam, ccd, sector, col_abs, row_abs):
    return TESS_PRF(cam, ccd, sector, col_abs, row_abs, localdatadir=PRFDIR)


def gaia_sources(radius_arcsec=180.0):
    rows = list(csv.DictReader(open(os.path.join(AUD, "data", "gaia", "gaia_dr3_neighbours_180as.csv"))))
    out = []
    for r in rows:
        if float(r["sep_arcsec"]) > radius_arcsec:
            continue
        out.append(dict(source_id=r["source_id"], ra=float(r["ra"]), dec=float(r["dec"]),
                        pmra=float(r["pmra"]) if r["pmra"] else 0.0, pmdec=float(r["pmdec"]) if r["pmdec"] else 0.0,
                        G=float(r["phot_g_mean_mag"]) if r["phot_g_mean_mag"] else np.nan,
                        bp_rp=float(r["bp_rp"]) if r["bp_rp"] else np.nan, sep=float(r["sep_arcsec"]),
                        ruwe=float(r["ruwe"]) if r["ruwe"] else np.nan))
    return out


def tmag_from_gaia(G, bp_rp):
    """Stassun et al. 2019 (AJ 158, 138) eq. for T from G and BP-RP, valid -0.2 < BP-RP < 3.5 (degrades > 2).
    Missing colour: T = G - 0.430 (their fallback), flagged 'uncertain' with 0.6 mag error (their value)."""
    if not np.isfinite(bp_rp):
        return G - 0.430, 0.6, "no_colour"
    c = bp_rp
    T = G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    flag = "ok" if -0.2 < c < 3.5 else "extrapolated"
    err = 0.006 if c <= 2 else 0.1
    return T, err, flag


def propagate(src, btjd):
    yr = Time(btjd + 2457000.0, format="jd", scale="tdb").jyear
    dt = yr - 2016.0
    ra = src["ra"] + src["pmra"] * dt / 3.6e6 / np.cos(np.radians(src["dec"]))
    de = src["dec"] + src["pmdec"] * dt / 3.6e6
    return ra, de


def scene(sources, wcs, btjd, target_id, shift=(0.0, 0.0), Tlim=19.0):
    """list of (x, y, rel_flux, source_id, T, sep) in stamp pixels; rel_flux relative to the target (from T mags).
    The target is identified by its Gaia source_id only (never by a position match), so the host cannot be duplicated
    into the neighbour list through the TIC(J2000)-vs-Gaia(2016) proper-motion offset (TIC 158955499 lesson)."""
    out = []
    Tt = None
    for s in sources:
        T, eT, fl = tmag_from_gaia(s["G"], s["bp_rp"])
        if s["source_id"] == str(target_id):
            Tt = T
    for s in sources:
        T, eT, fl = tmag_from_gaia(s["G"], s["bp_rp"])
        if not np.isfinite(T) or T > Tlim:
            continue
        ra, de = propagate(s, btjd)
        x, y = wcs.all_world2pix([[ra, de]], 0)[0]
        out.append(dict(x=float(x + shift[0]), y=float(y + shift[1]), rel_flux=10 ** (-0.4 * (T - Tt)), source_id=s["source_id"],
                        T=float(T), T_err=float(eT), T_flag=fl, sep=s["sep"], is_target=(s["source_id"] == str(target_id))))
    return out


def model_image(prf, shape, srcs, fluxes=None):
    img = np.zeros(shape)
    for k, s in enumerate(srcs):
        f = s["rel_flux"] if fluxes is None else fluxes[k]
        if -3 < s["x"] < shape[1] + 2 and -3 < s["y"] < shape[0] + 2:
            img += f * prf.locate(s["x"], s["y"], shape)
    return img
