"""Shared library for Audit 5 (signal B depth) of TIC 4206066.

Written from scratch for this audit. Inputs are only the fresh public files in ../00_data (sha256 in DATA_MANIFEST.json).
No code or derived tables from earlier analyses are used.

Conventions
  time   : BTJD = BJD_TDB - 2457000 (all products carry TIMESYS=TDB, BJDREFI=2457000)
  depth  : "D" = mid-transit depth of a fixed-shape batman template (1 - F(t_c)), in ppm, the same quantity the
           manuscript calls the per-sector depth. Shape (k, b, rho*, u1, u2) is fixed; only its amplitude is fitted.
"""
import os, json, glob, warnings
import numpy as np
from astropy.io import fits
import batman
from wotan import flatten

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "00_data"))
MAST = os.path.join(DATA, "mast")
CUT21 = os.path.join(DATA, "tesscut", "21x21")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
PLOTS = os.path.join(HERE, "plots"); os.makedirs(PLOTS, exist_ok=True)

# ---- ephemerides (manuscript values; used for masking and for the fixed-ephemeris depth fits) ----------------
EPH_A = dict(P=3.1827856, T0=2460998.7375 - 2457000.0, T14=1.87 / 24.0)        # BTJD 3998.7375
EPH_B = dict(P=11.13274, T0=2461000.865 - 2457000.0, T14=1.66 / 24.0)           # BTJD 4000.865
SIG_P_B, SIG_T0_B = 0.00003, 0.003

# ---- stellar / LD inputs --------------------------------------------------------------------------------------
RSTAR, E_RSTAR = 0.623, 0.019           # Rsun (manuscript, Mann et al. relations)
RHO_PRIOR = (2.51, 0.24)                # rho_sun, task-specified Gaussian prior
LD_NOM = (0.32, 0.31)                   # Claret (2017) TESS quadratic, see out/ldc_claret2017.json (ATLAS 0.375/0.278,
LD_SIG = (0.05, 0.05)                   #   PHOENIX-q 0.316/0.313, PHOENIX-r 0.299/0.333 at Teff 4012, logg 4.63)
RSUN_REARTH = 695700.0 / 6378.137       # IAU 2015 nominal R_sun / nominal equatorial R_earth = 109.076
G_CGS = 6.67430e-8
RHO_SUN_CGS = 1.98841e33 / (4.0 / 3.0 * np.pi * 6.957e10 ** 3)   # 1.4097 g/cm3

# ---- quality masks ----------------------------------------------------------------------------------------------
LK_DEFAULT = 17087          # lightkurve 2.6.0 TessQualityFlags.DEFAULT_BITMASK (verified in this environment)
QLP_STRAY = 2048            # bit 12: predicted Earth/Moon straylight (SPOC/TICA-derived)
QLP_B29 = 1 << 29           # undocumented in the QLP pages read here; set on 44% of S98 cadences
QLP_B30 = 1 << 30           # QLP internal poor-quality flag (tess.mit.edu/qlp, 2023-11-21 note)
QLP_SPACECRAFT = 1 | 4 | 8 | 16 | 32 | 128
MASK_S98_BASE = QLP_SPACECRAFT | QLP_STRAY | QLP_B29 | QLP_B30 | 4096   # "scattered-light intervals removed"


def ars_from_rho(rho, P):
    return (rho * RHO_SUN_CGS * G_CGS * (P * 86400.0) ** 2 / (3.0 * np.pi)) ** (1.0 / 3.0)


def phase_dt(t, eph):
    """time from nearest transit centre (d)"""
    ph = ((t - eph["T0"]) / eph["P"] + 0.5) % 1.0 - 0.5
    return ph * eph["P"]


def in_transit(t, eph, half_width):
    return np.abs(phase_dt(t, eph)) < half_width


# ---- template ------------------------------------------------------------------------------------------------------
class Template:
    """batman template with unit mid-transit depth. shape = dict(k, b, rho, u1, u2)."""

    def __init__(self, shape, P=EPH_B["P"]):
        self.shape = dict(shape); self.P = P
        pm = batman.TransitParams()
        pm.t0, pm.per, pm.rp = 0.0, P, shape["k"]
        pm.a = ars_from_rho(shape["rho"], P)
        cosi = shape["b"] / pm.a
        pm.inc = float(np.degrees(np.arccos(cosi)))
        pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, [shape["u1"], shape["u2"]], "quadratic"
        self.pm = pm
        m0 = batman.TransitModel(pm, np.array([0.0]))
        self.mid = 1.0 - float(m0.light_curve(pm)[0])
        k, b, a = shape["k"], shape["b"], pm.a
        arg = np.sqrt(max((1 + k) ** 2 - b ** 2, 0.0)) / (a * np.sin(np.radians(pm.inc)))
        self.t14_d = P / np.pi * np.arcsin(min(arg, 1.0))
        self._cache = {}

    def unit(self, trel, exptime_d, supersample=True):
        """(1 - F)/mid-depth at times trel (d from mid-transit); exposure-integrated if supersample."""
        trel = np.asarray(trel, float)
        ss = int(np.ceil(exptime_d * 1440.0 / 1.0)) if supersample else 1   # <= 1-min sub-samples
        if ss % 2 == 0:
            ss += 1
        if ss > 1:
            m = batman.TransitModel(self.pm, trel, supersample_factor=ss, exp_time=exptime_d)
        else:
            m = batman.TransitModel(self.pm, trel)
        return (1.0 - m.light_curve(self.pm)) / self.mid


def shape_from_T14(b, T14_h, k=0.032, u=LD_NOM, P=EPH_B["P"]):
    """rho giving a requested T14 for given b (circular)"""
    from scipy.optimize import brentq
    def f(lr):
        rho = 10 ** lr; a = ars_from_rho(rho, P); inc = np.arccos(b / a)
        arg = np.sqrt((1 + k) ** 2 - b ** 2) / (a * np.sin(inc))
        return P / np.pi * np.arcsin(min(arg, 1)) * 24 - T14_h
    lr = brentq(f, -1.5, 3.0)
    return dict(k=k, b=b, rho=10 ** lr, u1=u[0], u2=u[1])


# ---- loaders ---------------------------------------------------------------------------------------------------------
def _pack(name, sector, t, f, e, q, cad, meta):
    good = np.isfinite(t) & np.isfinite(f) & (f > 0)
    t, f, e, q = t[good], f[good], e[good], q[good]
    med = np.nanmedian(f)
    return dict(name=name, sector=sector, t=t, f=f / med, e=e / med if e is not None else None, q=q, cad=cad, meta=meta)


def load_spoc(path, name, sector, col="PDCSAP_FLUX", bitmask=LK_DEFAULT, crowd_correct_sap=False):
    h = fits.open(path); d = h[1].data; hd = h[1].header
    t = np.asarray(d["TIME"], float); f = np.asarray(d[col], float)
    e = np.asarray(d[col + "_ERR"], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
    m = (q & bitmask) == 0
    meta = dict(file=os.path.basename(path), col=col, CROWDSAP=hd.get("CROWDSAP"), FLFRCSAP=hd.get("FLFRCSAP"),
                TIMEDEL=hd.get("TIMEDEL"), bitmask=bitmask)
    t, f, e, q = t[m], f[m], e[m], q[m]
    if crowd_correct_sap and col == "SAP_FLUX":
        # SPOC PDC crowding/flux-fraction correction (SDPDD / Kepler DAWG convention):
        # f_corr = (f - (1 - CROWDSAP) * median(f)) / FLFRCSAP ; only the crowding term changes relative depths
        c = hd["CROWDSAP"]; med = np.nanmedian(f)
        f = (f - (1 - c) * med) / hd["FLFRCSAP"]
        meta["crowding_corrected"] = True
    cad = float(hd.get("TIMEDEL"))
    return _pack(name, sector, t, f, e, q, cad, meta)


def load_qlp(path, name, sector, col, bitmask):
    h = fits.open(path); d = h[1].data; hd = h[1].header
    t = np.asarray(d["TIME"], float); f = np.asarray(d[col], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
    ecol = {"SAP_FLUX": None, "KSPSAP_FLUX": "KSPSAP_FLUX_ERR", "DET_FLUX": "DET_FLUX_ERR", "SYS_RM_FLUX": "DET_FLUX_ERR"}.get(col)
    e = np.asarray(d[ecol], float) if ecol and ecol in d.columns.names else np.full_like(f, np.nan)
    m = (q & bitmask) == 0
    meta = dict(file=os.path.basename(path), col=col, bitmask=int(bitmask), TIMEDEL=hd.get("TIMEDEL"),
                BESTAP=hd.get("BESTAP"), CALIB=h[0].header.get("CALIB"))
    return _pack(name, sector, t[m], f[m], e[m], q[m], float(hd["TIMEDEL"]), meta)


def load_tglc_s6(path, col="cal_aper_flux"):
    h = fits.open(path); d = h[1].data; hd = h[1].header
    names = d.columns.names
    t = np.asarray(d["time"], float); f = np.asarray(d[col], float)
    q = np.asarray(d["TESS_flags"]).astype(np.int64) if "TESS_flags" in names else np.zeros(len(t), int)
    tq = np.asarray(d["TGLC_flags"]).astype(np.int64) if "TGLC_flags" in names else np.zeros(len(t), int)
    m = ((q & LK_DEFAULT) == 0) & (tq == 0)
    meta = dict(file=os.path.basename(path), col=col)
    return _pack("S6tglc", 6, t[m], f[m], None, q[m], 1800 / 86400.0, meta) if False else \
        dict(_pack("S6tglc", 6, t[m], f[m], np.full(m.sum(), np.nan), q[m], 1800 / 86400.0, meta))


# ---- own TESScut aperture photometry -----------------------------------------------------------------------------------
_CUBE_CACHE = {}


def load_cube(sector):
    if sector in _CUBE_CACHE:
        return _CUBE_CACHE[sector]
    p = glob.glob(os.path.join(CUT21, f"tess-s{sector:04d}-*_21x21_astrocut.fits"))[0]
    with fits.open(p, memmap=True) as h:
        d = h[1].data
        t = np.array(d["TIME"], float); fl = np.array(d["FLUX"], np.float32); q = np.array(d["QUALITY"]).astype(np.int64)
    _CUBE_CACHE[sector] = (t, fl, q, os.path.basename(p))
    return _CUBE_CACHE[sector]


def own_phot(sector, radius, bitmask=LK_DEFAULT | 2048, extra_bad_times=None, bkg_mode="annulus"):
    """Simple aperture photometry on the SPOC-calibrated FFI cutout (TESScut FLUX; FLUX_BKG is zero, i.e. no
    background subtracted, so the background is estimated here).
    Aperture: pixels whose centre lies within `radius` px of the flux-weighted centroid of the median image.
    Background: per cadence, median of pixels at 6-10 px from the target whose median-image value is below the
    40th percentile of that annulus (star-free sky)."""
    t, fl, q, fn = load_cube(sector)
    good = np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1) & ((q & bitmask) == 0)
    if extra_bad_times is not None:
        good &= ~extra_bad_times(t)
    t, fl, q = t[good], fl[good], q[good]
    med = np.median(fl, axis=0); ny, nx = med.shape; yy, xx = np.mgrid[:ny, :nx]
    c0 = np.hypot(xx - nx // 2, yy - ny // 2) <= 2
    iy, ix = np.unravel_index(np.argmax(np.where(c0, med, -np.inf)), med.shape)
    ann = (np.hypot(xx - ix, yy - iy) >= 6) & (np.hypot(xx - ix, yy - iy) <= 10)
    sky = ann & (med <= np.percentile(med[ann], 40))
    bkg = np.median(fl[:, sky], axis=1)
    sub = med - np.median(med[sky])
    w = np.hypot(xx - ix, yy - iy) <= 3
    cx = np.sum(sub[w] * xx[w]) / np.sum(sub[w]); cy = np.sum(sub[w] * yy[w]) / np.sum(sub[w])
    ap = np.hypot(xx - cx, yy - cy) <= radius
    f = (fl - bkg[:, None, None])[:, ap].sum(axis=1)
    npix = int(ap.sum())
    e = np.sqrt(np.abs(f) + npix * np.abs(bkg)) / 1.0  # rough, unused for fitting (empirical scatter used)
    cad = float(np.median(np.diff(t)))
    meta = dict(file=fn, radius=radius, npix=npix, centroid=(float(cx), float(cy)), bkg_median=float(np.median(bkg)),
                bitmask=int(bitmask))
    return _pack(f"own{sector}_r{radius}", sector, t, f, e, q, cad, meta)


# ---- detrending --------------------------------------------------------------------------------------------------------
def detrend(lc, method="biweight", window=0.75, mask_hw_B=None, mask_hw_A=None, extra_mask=None, break_tol=0.3,
            clip_sigma=4.0, **kw):
    """Global detrending with in-transit A and B cadences excluded from the trend (trend interpolated across them).
    Returns a copy with 'fd' (detrended flux) and 'sig' (robust OOT std)."""
    t, f = lc["t"], lc["f"]
    if mask_hw_B is None:
        mask_hw_B = 1.0 * EPH_B["T14"]   # +-T14 around each B mid-time (twice the half-duration)
    if mask_hw_A is None:
        mask_hw_A = 1.0 * EPH_A["T14"]
    msk = in_transit(t, EPH_B, mask_hw_B) | in_transit(t, EPH_A, mask_hw_A)
    if extra_mask is not None:
        msk |= extra_mask
    if method == "none":
        tr = np.ones_like(f)
    else:
        args = dict(method=method, window_length=window, return_trend=True, break_tolerance=break_tol)
        args.update(kw)
        _, tr0 = flatten(t[~msk], f[~msk], **args)
        ok = np.isfinite(tr0)
        tr = np.interp(t, t[~msk][ok], tr0[ok])
    fd = f / tr
    oot = ~msk & np.isfinite(fd)
    s = 1.4826 * np.median(np.abs(fd[oot] - np.median(fd[oot])))
    keep = np.isfinite(fd) & (fd < 1 + clip_sigma * s)          # clip only upward outliers (flares, cosmics)
    out = dict(lc); out.update(t=t[keep], f=f[keep], fd=fd[keep], q=lc["q"][keep], sig=float(s), detrend=dict(method=method, window=window))
    if lc.get("e") is not None:
        out["e"] = lc["e"][keep]
    return out


# ---- event bookkeeping -------------------------------------------------------------------------------------------------
def events_in(lc, eph=EPH_B, T14=None, half_window=0.3, min_cov=0.5, drop_A_hw=None):
    """list of (epoch, tc, window mask, in-transit coverage fraction) for B events covered by lc"""
    t = lc["t"]; T14 = T14 or eph["T14"]
    if drop_A_hw is None:
        drop_A_hw = 1.0 * EPH_A["T14"]
    notA = ~in_transit(t, EPH_A, drop_A_hw)
    n0 = int(np.floor((t.min() - eph["T0"]) / eph["P"])) - 1; n1 = int(np.ceil((t.max() - eph["T0"]) / eph["P"])) + 1
    evs = []
    cad = lc["cad"]
    for n in range(n0, n1 + 1):
        tc = eph["T0"] + n * eph["P"]
        w = (np.abs(t - tc) < half_window) & notA
        inn = w & (np.abs(t - tc) < T14 / 2)
        cov = inn.sum() * cad / T14
        if cov >= min_cov and w.sum() >= 10:
            evs.append(dict(epoch=n, tc=tc, w=w, cov=float(min(cov, 1.5))))
    return evs


# ---- linear depth fit ---------------------------------------------------------------------------------------------------
def design(lc, evs, tpl, poly=1, use="fd", dt_shift=0.0, supersample=True, per_event_depth=False):
    """Stack windows; columns = per-event polynomial baseline (order `poly`) + depth column(s)."""
    rows_t, rows_f, blocks = [], [], []
    for ev in evs:
        tt = lc["t"][ev["w"]]; ff = lc[use][ev["w"]]
        rows_t.append(tt); rows_f.append(ff); blocks.append((ev, tt))
    N = sum(len(x) for x in rows_t); ne = len(evs); nb = (poly + 1) * ne
    nd = ne if per_event_depth else 1
    X = np.zeros((N, nb + nd)); y = np.concatenate(rows_f); i0 = 0
    for j, (ev, tt) in enumerate(blocks):
        n = len(tt); x = (tt - ev["tc"]) / 0.3
        for p in range(poly + 1):
            X[i0:i0 + n, j * (poly + 1) + p] = x ** p
        s = tpl.unit(tt - ev["tc"] - dt_shift, lc["cad"], supersample=supersample)
        X[i0:i0 + n, nb + (j if per_event_depth else 0)] = -s
        i0 += n
    return X, y, nb


def fit_depth(lc, tpl, evs=None, poly=1, use="fd", dt_shift=0.0, supersample=True, per_event_depth=False, sigma=None):
    """Weighted linear LSQ. sigma per point = robust OOT std of the window residuals (sector-level) unless given.
    Returns depth(s) [ppm], white error(s) [ppm], chi2, n."""
    evs = evs if evs is not None else events_in(lc)
    if not evs:
        return None
    X, y, nb = design(lc, evs, tpl, poly, use, dt_shift, supersample, per_event_depth)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c
    s = sigma if sigma is not None else lc.get("sig", 1.4826 * np.median(np.abs(r - np.median(r))))
    cov = np.linalg.pinv(X.T @ X) * s ** 2
    D = c[nb:] * 1e6; eD = np.sqrt(np.diag(cov)[nb:]) * 1e6
    chi2 = float(np.sum((r / s) ** 2))
    return dict(D=D.tolist(), eD=eD.tolist(), chi2=chi2, n=int(len(y)), sig_ppm=float(s * 1e6),
                rms_resid_ppm=float(np.std(r) * 1e6), epochs=[ev["epoch"] for ev in evs], cov=[ev["cov"] for ev in evs])


def wmean(d, e):
    d = np.asarray(d, float); e = np.asarray(e, float); w = 1 / e ** 2
    return float(np.sum(w * d) / w.sum()), float(1 / np.sqrt(w.sum()))


def s98_vs_archival(dep):
    """dep: {6:(D,e),32:(D,e),98:(D,e)} -> dict with archival mean, difference, sigma, common depth, chi2s"""
    a, ea = wmean([dep[6][0], dep[32][0]], [dep[6][1], dep[32][1]])
    d98, e98 = dep[98]
    diff = d98 - a; sig = diff / np.hypot(e98, ea)
    allD = [dep[s][0] for s in (6, 32, 98)]; allE = [dep[s][1] for s in (6, 32, 98)]
    mu, emu = wmean(allD, allE)
    chi2_common = float(np.sum(((np.array(allD) - mu) / np.array(allE)) ** 2))
    return dict(archival=(a, ea), S98=(d98, e98), diff=diff, nsigma=float(sig), common=(mu, emu), chi2_common_3sectors=chi2_common)


def jdump(obj, name):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(str(type(o)))
    with open(os.path.join(OUT, name), "w") as fh:
        json.dump(obj, fh, indent=1, default=conv)


# ---- standard data sets -------------------------------------------------------------------------------------------------
F_S6_SPOC = os.path.join(MAST, "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits")
F_S32_2M = os.path.join(MAST, "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits")
F_S32_SPOC10 = os.path.join(MAST, "hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits")
F_S6_QLP = os.path.join(MAST, "hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits")
F_S32_QLP = os.path.join(MAST, "hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits")
F_S98_QLP = os.path.join(MAST, "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")
F_S6_TGLC = os.path.join(MAST, "hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits")
F_S6_ELEANOR = os.path.join(MAST, "hlsp_gsfc-eleanor-lite_tess_ffi_s0006-0000000004206066_tess_v1.0_lc.fits")


def baseline_sets():
    """the manuscript's product choice: TESS-SPOC S6 PDCSAP 30 min, SPOC S32 2-min PDCSAP, QLP S98 SYS_RM_FLUX 200 s"""
    return {6: load_spoc(F_S6_SPOC, "S6_tessspoc_pdc", 6),
            32: load_spoc(F_S32_2M, "S32_spoc2m_pdc", 32),
            98: load_qlp(F_S98_QLP, "S98_qlp_sysrm", 98, "SYS_RM_FLUX", MASK_S98_BASE)}


def qlp_s98_bad_time_fn(bitmask=MASK_S98_BASE):
    """returns a function flagging times (e.g. TESScut cadences) that fall on QLP-S98 cadences with bits in bitmask"""
    h = fits.open(F_S98_QLP); d = h[1].data
    tq = np.asarray(d["TIME"], float); qq = np.asarray(d["QUALITY"]).astype(np.int64)
    bad = (qq & bitmask) != 0
    order = np.argsort(tq); tq, bad = tq[order], bad[order]

    def fn(t):
        j = np.clip(np.searchsorted(tq, t), 1, len(tq) - 1)
        jj = np.where(np.abs(tq[j - 1] - t) < np.abs(tq[j] - t), j - 1, j)
        near = np.abs(tq[jj] - t) < 150 / 86400.0
        return near & bad[jj]
    return fn
