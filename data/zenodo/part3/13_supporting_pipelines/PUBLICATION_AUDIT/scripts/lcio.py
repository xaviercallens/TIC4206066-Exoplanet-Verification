"""Light-curve loading and detrending used by every analysis step of the audit (one place, one convention).

Flux conventions: every returned flux is a relative flux normalised to its median (dimensionless, 1 = baseline);
depths are reported in ppm = 1e6 x (1 - flux) only at the reporting stage.
Time: BTJD = BJD_TDB - 2457000 (all products declare TIMESYS=TDB, TIMEREF=SOLARSYSTEM, BJDREFI=2457000).
TESScut times carry the CCD-centre barycentric correction; `tesscut_time_fix` adds the target-specific difference
(computed in s06_timing.py; |shift| < 20 s).
"""
import glob
import os

import numpy as np
from astropy.io import fits

AUD = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAST = os.path.join(AUD, "data", "mast")
PHOT = os.path.join(AUD, "out", "phot")

DEFAULT_BITMASK = 17087          # lightkurve TessQualityFlags.DEFAULT_BITMASK
STRAY = 2048                     # SPOC/QLP bit 12 (value 2048) straylight
QLP_FLAG = 1 << 30               # QLP internal poor-quality flag (tess.mit.edu/qlp)
QLP_BIT29 = 1 << 29              # undocumented in QLP docs ("bits 16-30 unused"); set in 44 % of S98 cadences
FROZEN_QLP_MASK = 1610619071     # mask used by the frozen pipeline and the deep follow-up for QLP S98

# preliminary ephemerides (used ONLY to mask transits during detrending; refined values from this audit replace them)
EPH = {
    "A": dict(P=3.1827863, T0=3998.7368, T14=1.86 / 24),
    "B": dict(P=11.132743, T0=4000.8661, T14=1.70 / 24),
}

PRODUCTS = {
    # key: (file, flux column, sector, kind)
    "SPOC2min_S32_PDC": ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX", 32, "spoc"),
    "SPOC2min_S32_SAP": ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "SAP_FLUX", 32, "spoc"),
    "TESSSPOC_S6_PDC": ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX", 6, "spoc"),
    "TESSSPOC_S6_SAP": ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "SAP_FLUX", 6, "spoc"),
    "TESSSPOC_S32_PDC": ("hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits", "PDCSAP_FLUX", 32, "spoc"),
    "QLP_S6": ("hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits", "KSPSAP_FLUX", 6, "qlp_old"),
    "QLP_S6_SAP": ("hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits", "SAP_FLUX", 6, "qlp_old"),
    "QLP_S32": ("hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits", "KSPSAP_FLUX", 32, "qlp_old"),
    "QLP_S32_SAP": ("hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits", "SAP_FLUX", 32, "qlp_old"),
    "QLP_S98_SYS": ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX", 98, "qlp_tglc"),
    "QLP_S98_SAP": ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SAP_FLUX", 98, "qlp_tglc"),
    "QLP_S98_DET": ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "DET_FLUX", 98, "qlp_tglc"),
    "eleanor_S6": ("hlsp_gsfc-eleanor-lite_tess_ffi_s0006-0000000004206066_tess_v1.0_lc.fits", "CORR_FLUX", 6, "eleanor"),
    "TGLC_S6": ("hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits", "cal_aper_flux", 6, "tglc"),
    "TGLC_S6_PSF": ("hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits", "cal_psf_flux", 6, "tglc"),
}

TESSCUT_SHIFT_S = {6: -11.78, 32: -3.28, 98: 0.49}   # median target-minus-CCD-centre barycentric correction [s]


def in_transit(t, key, fac=1.0, eph=None):
    e = (eph or EPH)[key]
    ph = ((t - e["T0"]) / e["P"] + 0.5) % 1 - 0.5
    return np.abs(ph * e["P"]) < fac * e["T14"]


def load_product(key, qmask="default"):
    """qmask for QLP S98: 'frozen' (FROZEN_QLP_MASK), 'default' (hard SPOC bits + QLP 2^30 flag),
    'default+stray' (also drop 2048), 'none_hard' (only hard SPOC bits: keeps 2^29 and 2^30)."""
    fn, col, sec, kind = PRODUCTS[key]
    h = fits.open(os.path.join(MAST, fn))
    d = h[1].data
    t = np.asarray(d["TIME"], float)
    f = np.asarray(d[col], float)
    if kind == "spoc":
        q = np.asarray(d["QUALITY"]).astype(np.int64)
        good = (q & DEFAULT_BITMASK) == 0
    elif kind in ("qlp_old", "qlp_tglc"):
        q = np.asarray(d["QUALITY"]).astype(np.int64)
        if qmask == "frozen":
            good = (q & FROZEN_QLP_MASK) == 0
        elif qmask == "default":
            good = ((q & DEFAULT_BITMASK) == 0) & ((q & QLP_FLAG) == 0)
        elif qmask == "default+stray":
            good = ((q & DEFAULT_BITMASK) == 0) & ((q & QLP_FLAG) == 0) & ((q & STRAY) == 0)
        elif qmask == "none_hard":
            good = (q & DEFAULT_BITMASK) == 0
        else:
            raise ValueError(qmask)
        if kind == "qlp_old" and qmask in ("default", "default+stray"):
            good &= (q & 4096) == 0          # old-format QLP low-precision bit 13 (value 4096)
    elif kind == "eleanor":
        q = np.asarray(d["QUALITY"]).astype(np.int64)
        good = q == 0
    elif kind == "tglc":
        q = np.asarray(d["TESS_flags"]).astype(np.int64) if "TESS_flags" in d.columns.names else np.zeros(len(t), int)
        tg = np.asarray(d["TGLC_flags"]).astype(np.int64) if "TGLC_flags" in d.columns.names else np.zeros(len(t), int)
        good = ((q & DEFAULT_BITMASK) == 0) & (tg == 0)
    else:
        raise ValueError(kind)
    m = np.isfinite(t) & np.isfinite(f) & good & (f > 0)
    t, f, q = t[m], f[m], q[m]
    f = f / np.nanmedian(f)
    cad = float(np.median(np.diff(t)))
    hdr = {**h[0].header, **h[1].header}
    return dict(key=key, t=t, f=f, q=q, sector=sec, cad=cad, kind=kind, file=fn,
                crowdsap=hdr.get("CROWDSAP"), flfrcsap=hdr.get("FLFRCSAP"))


def load_own(name, radius=2.0, drop_stray=False, time_fix=True):
    """Own aperture photometry (s05_photometry.py). name: tesscut_s6/s32/s98, spoc_tpf_s32, tessspoc_tp_s6/_s32."""
    z = np.load(os.path.join(PHOT, f"{name}.npz"))
    t = z["time"].copy()
    sec = int(name.split("_s")[-1])
    if name.startswith("tesscut") and time_fix:
        t = t + TESSCUT_SHIFT_S[sec] / 86400.0
    f = z[f"flux_r{radius}"].astype(float)
    good = np.isfinite(f) & (f > 0)
    if drop_stray:
        good &= ~z["straylight"]
    f = f[good] / np.nanmedian(f[good])
    return dict(key=f"{name}_r{radius}" + ("_nostray" if drop_stray else ""), t=t[good], f=f, q=z["quality"][good],
                sector=sec, cad=float(np.median(np.diff(t[good]))), kind="own", bkg=z["bkg"][good], cx=z["cx"][good],
                cy=z["cy"][good], stray=z["straylight"][good], pos_corr1=z["pos_corr1"][good], pos_corr2=z["pos_corr2"][good])


def detrend(lc, window=0.75, method="biweight", mask_fac=1.5, eph=None, clip_hi=4.0, clip_lo=None, edge_cut=0.0,
            break_tol=0.3):
    """Divide by a wotan trend fitted with BOTH signals masked (+-mask_fac*T14). Upward outliers > clip_hi*MAD removed.
    Returns a copy with 'fd' (detrended), 'trend', 'sig' (robust per-point scatter)."""
    from wotan import flatten
    t, f = lc["t"], lc["f"]
    msk = in_transit(t, "A", mask_fac, eph) | in_transit(t, "B", mask_fac, eph)
    if method == "none":
        tr = np.ones_like(f)
    else:
        kw = dict(method=method, window_length=window, return_trend=True, break_tolerance=break_tol, edge_cutoff=edge_cut)
        _, trm = flatten(t[~msk], f[~msk], **kw)
        ok = np.isfinite(trm)
        tr = np.interp(t, t[~msk][ok], trm[ok])
    fd = f / tr
    r = fd[~msk] - 1
    s = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
    keep = np.isfinite(fd) & (fd < 1 + clip_hi * s)
    if clip_lo is not None:
        keep &= fd > 1 - clip_lo * s
    out = dict(lc)
    for k in list(lc.keys()):
        if isinstance(lc[k], np.ndarray) and len(lc[k]) == len(t):
            out[k] = lc[k][keep]
    out["fd"] = fd[keep]
    out["trend"] = tr[keep]
    out["sig"] = float(s)
    out["detrend"] = dict(window=window, method=method, mask_fac=mask_fac, clip_hi=clip_hi)
    return out


def box_depth(t, f, P, T0, dur, sig=None, err=None):
    """Weighted box depth (ppm) inside |phase|<dur/2 relative to 1 (data already detrended). Returns depth, error, n."""
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    inn = np.abs(ph * P) < dur / 2
    if inn.sum() == 0:
        return np.nan, np.nan, 0
    w = 1.0 / (err[inn] ** 2) if err is not None else np.ones(inn.sum())
    d = 1 - np.sum(w * f[inn]) / np.sum(w)
    e = (1 / np.sqrt(np.sum(w))) if err is not None else sig / np.sqrt(inn.sum())
    return float(d * 1e6), float(e * 1e6), int(inn.sum())
