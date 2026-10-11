"""R-S6 referee check: shared helpers (written from scratch for this check; no audit code is imported).

Times: BTJD = BJD_TDB - 2457000 (all three products: TIMESYS=TDB, BJDREFI=2457000).
Data (fresh 00_data/mast, sha256 recorded in every output):
  S6  TESS-SPOC 30-min PDCSAP_FLUX   (FFI, TIMEDEL 1800 s)
  S32 SPOC 2-min PDCSAP_FLUX          (TIMEDEL 120 s)
  S98 QLP 200-s SYS_RM_FLUX           (TIMEDEL 200 s)
Quality: QUALITY == 0 in all three. For QLP S98 this removes bit 12 (2048, scattered light), the QLP internal bit 2^30
and the undocumented 2^29 bit that brackets the scattered-light intervals (the manuscript's "scattered-light intervals
removed"); for S6 it removes the momentum-dump cadences (bit 6, value 32) and the manual-exclude cadences next to them.
"""
import os, hashlib
import numpy as np
from astropy.io import fits
import batman

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "00_data", "mast")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

# manuscript ephemerides (used only to place windows and masks; everything is refitted)
P_A, T0_A = 3.1827856, 2460998.7375 - 2457000.0
T14_A = 1.87 / 24.0
P_B, T0_B = 11.13274, 2461000.865 - 2457000.0
T14_B = 1.65 / 24.0

# stellar parameters (manuscript Sec. 1)
R_STAR, E_R_STAR = 0.623, 0.019
RHO_MED, RHO_EP, RHO_EM = 2.51, 0.25, 0.22
RSUN_OVER_REARTH = 6.957e8 / 6.3781e6          # IAU 2015 nominal
G_CGS = 6.67430e-8
RHO_SUN = 1.98841e33 / (4.0 / 3.0 * np.pi * 6.957e10 ** 3)

FILES = {
    6: ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX"),
    32: ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX"),
    98: ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX"),
}
# integration time span per cadence (s) = TIMEDEL; supersampling factor (sub-steps <= 2 min)
CAD_S = {6: 1800.0, 32: 120.0, 98: 200.0}
SS = {6: 15, 32: 2, 98: 4}

# S6 momentum dumps: cadences with QUALITY bit 6 (value 32) in the TESS-SPOC S6 file (see load_dumps_s6)


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def load(sector, qmode="q0"):
    fn, col = FILES[sector]
    path = os.path.join(DATA, fn)
    h = fits.open(path)
    d = h[1].data
    t = np.asarray(d["TIME"], float)
    f = np.asarray(d[col], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    good = np.isfinite(t) & np.isfinite(f) & (f > 0)
    if qmode == "q0":
        good &= q == 0
    elif qmode == "lkdefault":          # lightkurve 'default' bitmask 175 would be looser; 17087 = the package's choice
        good &= (q & 17087) == 0
    tt, ff = t[good], f[good]
    ff = ff / np.nanmedian(ff)
    meta = dict(file=fn, sha256=sha256(path), column=col, n_rows=int(len(t)), n_kept=int(good.sum()), qmode=qmode,
                TIMEDEL_s=float(h[1].header["TIMEDEL"]) * 86400.0)
    return tt, ff, meta


def load_dumps_s6():
    fn, _ = FILES[6]
    d = fits.open(os.path.join(DATA, fn))[1].data
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    return np.asarray(d["TIME"], float)[(q & 32) != 0]


def phase_dt(t, P, T0):
    ep = np.round((t - T0) / P)
    return t - (T0 + ep * P), ep.astype(int)


def ars_from_rho(rho, P):
    return (G_CGS * rho * RHO_SUN * (P * 86400.0) ** 2 / (3.0 * np.pi)) ** (1.0 / 3.0)


def t14_hours(P, k, b, ars):
    inc = np.arccos(np.clip(b / ars, 0, 1))
    arg = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (ars * np.sin(inc))
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1)) * 24.0


def claret_u():
    """Claret (2017) TESS quadratic LD, PHOENIX-COND quasi-spherical (Type q), LSM, Z=0, bilinear in Teff/logg,
    read from the VizieR table saved by audit 07 (data file only).  Teff=4012 K, logg=4.632."""
    from astropy.table import Table
    tab = Table.read(os.path.join(ROOT, "07_A_shape_radii", "claret", "tableab.ecsv"))
    tab = tab[(tab["Type"] == "q") & (tab["Z"] == 0.0)]
    teff, logg = 4012.0, 4.438 + np.log10(0.606) - 2 * np.log10(0.623)

    def v(T, g):
        r = tab[(tab["Teff"] == T) & np.isclose(tab["logg"], g)]
        return np.array([float(r["aLSM"][0]), float(r["bLSM"][0])])
    gw = (logg - 4.5) / 0.5
    tw = (teff - 4000.0) / 100.0
    lo = v(4000.0, 4.5) * (1 - gw) + v(4000.0, 5.0) * gw
    hi = v(4100.0, 4.5) * (1 - gw) + v(4100.0, 5.0) * gw
    return lo * (1 - tw) + hi * tw


class Model:
    """batman quadratic-LD circular model on fixed times, supersampled over the cadence (exp_time = TIMEDEL)."""

    def __init__(self, t, sector, ss=None):
        self.t = np.asarray(t, float)
        pm = batman.TransitParams()
        pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w = 0.0, P_A, 0.02, 12.0, 89.5, 0.0, 90.0
        pm.u, pm.limb_dark = [0.3, 0.3], "quadratic"
        self.pm = pm
        ss = SS[sector] if ss is None else ss
        self.ss = ss
        if ss > 1:
            self.m = batman.TransitModel(pm, self.t, supersample_factor=ss, exp_time=CAD_S[sector] / 86400.0)
        else:
            self.m = batman.TransitModel(pm, self.t)

    def flux(self, T0, P, k, b, ars, u1, u2):
        pm = self.pm
        pm.t0, pm.per, pm.rp, pm.a = T0, P, k, ars
        pm.inc = np.degrees(np.arccos(min(b / ars, 1.0)))
        pm.u = [u1, u2]
        return self.m.light_curve(pm)


def windows(sector, half=0.30, qmode="q0", t_extra_mask=None):
    """Raw normalised flux within +-half d of every predicted A transit; B transits (|dt_B| < T14_B) removed."""
    t, f, meta = load(sector, qmode)
    dA, ep = phase_dt(t, P_A, T0_A)
    dB, _ = phase_dt(t, P_B, T0_B)
    m = (np.abs(dA) < half) & (np.abs(dB) > T14_B)
    if t_extra_mask is not None:
        m &= ~t_extra_mask(t)
    return t[m], f[m], ep[m], meta


def coverage(t, tc, cad_s):
    nexp = T14_A / (cad_s / 86400.0)
    inn = np.abs(t - tc) < T14_A / 2
    pre = (t < tc - T14_A / 2) & (t > tc - 0.3)
    post = (t > tc + T14_A / 2) & (t < tc + 0.3)
    return float(inn.sum() / nexp), int(pre.sum()), int(post.sum())


def pct(x, q=(2.275, 15.865, 50.0, 84.135, 97.725)):
    return np.percentile(np.asarray(x), q).tolist()
