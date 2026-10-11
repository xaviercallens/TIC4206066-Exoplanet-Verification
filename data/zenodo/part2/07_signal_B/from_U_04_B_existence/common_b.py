"""Common helpers for Audit 4 (Signal B existence), TIC 4206066.

Written from scratch for this audit. Inputs are only the raw public files in ../00_data
(sha256 in ../00_data/DATA_MANIFEST.json). No code or derived tables from earlier analyses are used.

Time system: all products use BTJD = BJD_TDB - 2457000 (TIMESYS=TDB, TIMEREF=SOLARSYSTEM, BJDREFI=2457000).
"""
import os
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "00_data")
MAST = os.path.join(DATA, "mast")
CUT = os.path.join(DATA, "tesscut", "21x21")

# ---- ephemerides used for masking / prediction -------------------------------------------------
# Signal A (manuscript; independently checked by Audit 1).  Only used to mask A.
P_A, T0_A, T14_A = 3.1827856, 3998.7375, 1.87 / 24.0
# Signal B: Sector 98 search value (candidate.json) and manuscript linear ephemeris.
P_B98, T0_B98 = 11.132975671662091, 4000.8681717061445
P_B_MS, T0_B_MS = 11.13274, 4000.865
T14_B = 1.6 / 24.0

# S6 momentum dumps (TESS-SPOC/TESScut QUALITY bit 6; 00_data/logs/momentum_dumps_from_quality.csv)
DUMPS = {
    6: [1468.3912, 1471.5161, 1474.6411, 1481.2868, 1484.4117, 1487.5366],
    32: [2179.4865, 2186.9449, 2193.4867],
    98: [3988.4987, 3995.8208, 4003.0109, 4010.0366, 4017.0992, 4024.5251, 4031.8029, 4038.8444],
}

# SPOC default-like bitmask (lightkurve "default" = 175 | ... ); we write it explicitly from SDPDD Table 32:
# 1 attitude tweak, 2 safe mode, 4 coarse point, 8 earth point, 16 argabrightening, 32 desat, 128 manual exclude,
# 1024 impulsive outlier (2-min) / cosmic ray in collateral (FFI), 2048 straylight (FFI: predicted), 4096 SPOC scattered-light exclude,
# 8192 planet search exclude, 32768 insufficient targets / bad calibration (FFI 15), 65536 ...
SPOC_BAD = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 1024 | 4096 | 8192 | 16384 | 32768 | 65536


def in_A(t, half=None):
    """True within +-half of an A transit centre (default half = T14_A, i.e. a 2*T14 wide mask)."""
    half = T14_A if half is None else half
    ph = ((t - T0_A) / P_A + 0.5) % 1.0 - 0.5
    return np.abs(ph * P_A) < half


def pred_times(P, T0, tmin, tmax):
    n0 = int(np.floor((tmin - T0) / P)) - 1
    n1 = int(np.ceil((tmax - T0) / P)) + 1
    n = np.arange(n0, n1 + 1)
    tc = T0 + n * P
    m = (tc > tmin - 0.2) & (tc < tmax + 0.2)
    return n[m], tc[m]


# ---- loaders ---------------------------------------------------------------------------------------
def load_qlp(sector):
    fn = {6: "hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits",
          32: "hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits",
          98: "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"}[sector]
    d = fits.open(os.path.join(MAST, fn))[1].data
    det = "DET_FLUX" if sector == 98 else "KSPSAP_FLUX"
    out = dict(t=np.asarray(d["TIME"], float), f=np.asarray(d[det], float), sap=np.asarray(d["SAP_FLUX"], float),
               q=np.asarray(d["QUALITY"]).astype(np.int64), bkg=np.asarray(d["SAP_BKG"], float),
               x=np.asarray(d["SAP_X"], float), y=np.asarray(d["SAP_Y"], float))
    if sector == 98:
        out["sysrm"] = np.asarray(d["SYS_RM_FLUX"], float)
    return out


def load_spoc_lc(kind):
    """kind: 's6_tessspoc' (30 min), 's32_2min', 's32_tessspoc' (10 min)."""
    fn = {"s6_tessspoc": "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits",
          "s32_2min": "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits",
          "s32_tessspoc": "hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits"}[kind]
    d = fits.open(os.path.join(MAST, fn))[1].data
    return dict(t=np.asarray(d["TIME"], float), pdc=np.asarray(d["PDCSAP_FLUX"], float), sap=np.asarray(d["SAP_FLUX"], float),
                q=np.asarray(d["QUALITY"]).astype(np.int64), bkg=np.asarray(d["SAP_BKG"], float),
                cx=np.asarray(d["MOM_CENTR1"], float), cy=np.asarray(d["MOM_CENTR2"], float))


def load_eleanor6():
    d = fits.open(os.path.join(MAST, "hlsp_gsfc-eleanor-lite_tess_ffi_s0006-0000000004206066_tess_v1.0_lc.fits"))[1].data
    return dict(t=np.asarray(d["TIME"], float), corr=np.asarray(d["CORR_FLUX"], float), raw=np.asarray(d["RAW_FLUX"], float),
                pca=np.asarray(d["PCA_FLUX"], float), q=np.asarray(d["QUALITY"]).astype(np.int64), bkg=np.asarray(d["FLUX_BKG"], float))


def load_tglc6():
    h = fits.open(os.path.join(MAST, "hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits"))
    d = h[1].data
    cols = d.columns.names
    out = dict(t=np.asarray(d["time"] if "time" in cols else d["TIME"], float), cols=cols)
    for c in cols:
        out[c] = np.asarray(d[c])
    return out


# ---- own TESScut aperture photometry ---------------------------------------------------------------
def tesscut_phot(sector, ap_thresh=0.12, max_r=3.0):
    """Simple aperture photometry on the SPOC-FFI TESScut 21x21 cutout (FLUX not background-subtracted, BACKAPP=F).
    Aperture: pixels within max_r px of the brightest central pixel whose median flux exceeds ap_thresh * peak.
    Background: per-cadence median of the stamp pixels with median flux below the 30th percentile, farther than 5 px from
    the target.  Returns t, flux (aperture sum - n_ap * bkg), bkg per pixel, flux-weighted centroids, quality, aperture."""
    fn = {6: "tess-s0006-1-1_80.859830_-1.326960_21x21_astrocut.fits",
          32: "tess-s0032-1-2_80.859830_-1.326960_21x21_astrocut.fits",
          98: "tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits"}[sector]
    h = fits.open(os.path.join(CUT, fn), memmap=True)
    d = h[1].data
    t = np.asarray(d["TIME"], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    F = np.asarray(d["FLUX"], np.float32)
    good = np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)[:, 200:241]), axis=1)
    med = np.nanmedian(F[good][:: max(1, good.sum() // 400)], axis=0)
    c = np.array(med.shape) // 2
    sub = med[c[0] - 2:c[0] + 3, c[1] - 2:c[1] + 3]
    py, px = np.unravel_index(np.nanargmax(sub), sub.shape)
    py += c[0] - 2
    px += c[1] - 2
    yy, xx = np.mgrid[0:med.shape[0], 0:med.shape[1]]
    rr = np.hypot(yy - py, xx - px)
    bkg_ref = np.nanpercentile(med, 30)
    ap = (rr <= max_r) & ((med - bkg_ref) > ap_thresh * (med[py, px] - bkg_ref))
    bmask = (rr > 5) & (med < np.nanpercentile(med, 30))
    Fb = F[:, bmask]
    bkg = np.nanmedian(Fb, axis=1)
    flux = np.nansum(F[:, ap], axis=1) - ap.sum() * bkg
    w = F[:, ap] - bkg[:, None]
    fx = np.nansum(w * xx[ap][None, :], axis=1) / np.nansum(w, axis=1)
    fy = np.nansum(w * yy[ap][None, :], axis=1) / np.nansum(w, axis=1)
    h.close()
    return dict(t=t, flux=flux, bkg=bkg, cx=fx, cy=fy, q=q, ap=ap, peak=(int(py), int(px)), n_ap=int(ap.sum()))


# ---- detrending ------------------------------------------------------------------------------------
def detrend(t, f, method="biweight", window=0.75, mask=None, break_tol=0.5):
    """Wotan flatten; points in `mask` (True = exclude) are not used to fit the trend but receive the interpolated trend.
    Returns flattened flux (f / trend)."""
    from wotan import flatten
    ok = np.isfinite(t) & np.isfinite(f)
    use = ok & (~mask if mask is not None else True)
    out = np.full_like(f, np.nan, dtype=float)
    if method in ("biweight", "median", "mean", "welsch", "huber"):
        _, tr = flatten(t[use], f[use], method=method, window_length=window, edge_cutoff=0.0, break_tolerance=break_tol,
                        return_trend=True)
        # interpolate the trend to all good points (segment-aware: only within the segment spans)
        trend = np.interp(t[ok], t[use][np.isfinite(tr)], tr[np.isfinite(tr)])
        out[ok] = f[ok] / trend
    elif method == "cosine":
        _, tr = flatten(t[use], f[use], method="cosine", window_length=window, break_tolerance=break_tol, robust=True,
                        return_trend=True)
        trend = np.interp(t[ok], t[use][np.isfinite(tr)], tr[np.isfinite(tr)])
        out[ok] = f[ok] / trend
    else:
        raise ValueError(method)
    return out


# ---- box statistics --------------------------------------------------------------------------------
def box_depth_grid(t, f, grid, dur, flank=None):
    """Box depth (ppm-free, fractional) at each grid centre: mean(flank) - mean(in).  Flanks: dur/2 < |dt| < dur/2 + flank
    (default flank = dur).  t sorted, f finite.  Returns depth, n_in, n_out."""
    flank = dur if flank is None else flank
    o = np.argsort(t)
    t = t[o]
    f = f[o]
    cs = np.r_[0.0, np.cumsum(f)]
    def rng(a, b):
        i0 = np.searchsorted(t, a, "left")
        i1 = np.searchsorted(t, b, "left")
        return cs[i1] - cs[i0], i1 - i0
    s_in, n_in = rng(grid - dur / 2, grid + dur / 2)
    s_l, n_l = rng(grid - dur / 2 - flank, grid - dur / 2)
    s_r, n_r = rng(grid + dur / 2, grid + dur / 2 + flank)
    n_out = n_l + n_r
    with np.errstate(invalid="ignore", divide="ignore"):
        m_in = s_in / n_in
        m_out = (s_l + s_r) / n_out
        depth = m_out - m_in
    return depth, n_in, n_out


def robust_std(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


# ---- transit template and per-event fits ------------------------------------------------------------
def trapezoid(dt, T14, tau):
    """Unit-depth trapezoid: 1 inside flat bottom, linear ingress/egress of length tau, 0 outside |dt| > T14/2."""
    a = np.abs(dt)
    s = np.clip((T14 / 2 - a) / tau, 0.0, 1.0)
    return s


def fit_event(t, f, tc0, T14=None, tau_frac=0.15, win=None, scan_min=60.0, step_min=0.5, return_curve=False):
    """Fit one event: f = a + b*(t-tc0) - D*trap(t-tc).  The data window is FIXED (|t - tc0| < win + scan, default
    win = 3*T14) so that chi2 values at different trial centres are comparable.  tc is scanned over tc0 +- scan_min
    (step step_min) with linear LSQ for (a, b, D) at each tc; the best tc minimises chi2 subject to D > 0.
    Returns dict(tc, dt_min, D, eD_white, n_in, sig_pt, at_edge[, curve])."""
    T14 = T14_B if T14 is None else T14
    win = 3 * T14 if win is None else win
    m = np.isfinite(t) & np.isfinite(f) & (np.abs(t - tc0) < win + scan_min / 1440)
    tt, ff = t[m], f[m]
    if len(tt) < 10:
        return None
    x0 = tt - tc0
    dts = np.arange(-scan_min, scan_min + 1e-9, step_min) / 1440.0
    chis = np.full(len(dts), np.nan)
    Ds = np.full(len(dts), np.nan)
    for i, dtc in enumerate(dts):
        s = trapezoid(x0 - dtc, T14, tau_frac * T14)
        if np.sum(s > 0.5) < 2:
            continue
        A = np.vstack([np.ones_like(x0), x0, -s]).T
        coef, *_ = np.linalg.lstsq(A, ff, rcond=None)
        r = ff - A @ coef
        chis[i] = np.sum(r ** 2)
        Ds[i] = coef[2]
    ok = np.isfinite(chis) & (Ds > 0)
    if not ok.any():
        return None
    i = np.where(ok)[0][np.argmin(chis[ok])]
    dtc = dts[i]
    s = trapezoid(x0 - dtc, T14, tau_frac * T14)
    A = np.vstack([np.ones_like(x0), x0, -s]).T
    coef, *_ = np.linalg.lstsq(A, ff, rcond=None)
    r = ff - A @ coef
    sig = np.std(r[s == 0]) if np.sum(s == 0) > 5 else np.std(r)
    cov = np.linalg.pinv(A.T @ A) * sig ** 2
    out = dict(tc=tc0 + dtc, dt_min=dtc * 1440, D=coef[2], eD_white=np.sqrt(cov[2, 2]), n_in=int(np.sum(s > 0.5)), sig_pt=sig,
               at_edge=bool(abs(dtc) * 1440 > scan_min - 2 * step_min))
    if return_curve:
        out["curve"] = dict(dt_min=(dts * 1440).tolist(), dchi2=((chis - chis[i]) / sig ** 2).tolist(), D=Ds.tolist())
    return out


def fit_event_fixed(t, f, tc, T14=None, tau_frac=0.15, win=None):
    """Depth at a fixed centre with a local linear baseline (same model as fit_event)."""
    T14 = T14_B if T14 is None else T14
    win = 3 * T14 if win is None else win
    m = np.isfinite(t) & np.isfinite(f) & (np.abs(t - tc) < win)
    x = t[m] - tc
    y = f[m]
    s = trapezoid(x, T14, tau_frac * T14)
    if np.sum(s > 0.5) < 3 or np.sum(s == 0) < 6:
        return None
    A = np.vstack([np.ones_like(x), x, -s]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef
    sig = np.std(r[s == 0])
    cov = np.linalg.pinv(A.T @ A) * sig ** 2
    return dict(D=coef[2], eD_white=np.sqrt(cov[2, 2]), n_in=int(np.sum(s > 0.5)), n_out=int(np.sum(s == 0)), sig_pt=sig)
