"""Verifier common helpers (independent of the auditor's common.py). Times in BTJD (TIMESYS=TDB, BJDREFI=2457000)."""
import os
import numpy as np
import batman
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "00_data", "mast")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
F = {6: ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX", 1800.0, 15),
     32: ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX", 120.0, 3),
     98: ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX", 200.0, 5)}
PA, T0A, T14A = 3.1827856, 2460998.7375 - 2457000, 1.87 / 24
PB, T0B, T14B = 11.13274, 2461000.865 - 2457000, 1.65 / 24
RS, ERS, MS, EMS = 0.623, 0.019, 0.606, 0.015
RHO0, RHOP, RHOM = 2.51, 0.25, 0.22
LD_PHX = (0.3156, 0.3131)   # own interpolation of VizieR J/A+A/600/A30 tableab (Type q, LSM) at 4012 K, logg 4.63
LD_ATL = (0.3750, 0.2780)
REARTH_PER_RSUN = 6.957e8 / 6.3781e6          # IAU 2015 nominal
G = 6.67430e-8; RHOSUN = 1.98841e33 / (4 / 3 * np.pi * 6.957e10 ** 3)


def ars(rho, P):
    return (G * rho * RHOSUN * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


def t14(P, k, b, a):
    inc = np.arccos(b / a)
    return P / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (a * np.sin(inc)), 0, 1))


def load(sector, col=None):
    fn, c, cad, ss = F[sector]
    d = fits.open(os.path.join(D, fn))[1].data
    t = np.asarray(d["TIME"], float); f = np.asarray(d[col or c], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    g = (q == 0) & np.isfinite(t) & np.isfinite(f)
    return t[g], f[g] / np.nanmedian(f[g]), cad, ss


def segment_of(t, tc, gap=0.1):
    """boolean mask of points in the same contiguous segment (gaps > gap d split) as time tc."""
    br = np.where(np.diff(t) > gap)[0]
    starts = np.r_[0, br + 1]; ends = np.r_[br, len(t) - 1]
    m = np.zeros(len(t), bool)
    for a, b in zip(starts, ends):
        if t[a] - 0.05 <= tc <= t[b] + 0.05:
            m[a:b + 1] = True
    return m


class Model:
    def __init__(self, t, cad, ss):
        p = batman.TransitParams()
        p.t0, p.per, p.rp, p.a, p.inc, p.ecc, p.w, p.u, p.limb_dark = 0, PA, 0.02, 12, 89.5, 0, 90, [0.3, 0.3], "quadratic"
        self.p = p
        self.m = batman.TransitModel(p, t, supersample_factor=ss, exp_time=cad / 86400) if ss > 1 else batman.TransitModel(p, t)

    def __call__(self, T0, P, k, b, a, u1, u2):
        p = self.p
        p.t0, p.per, p.rp, p.a, p.inc, p.u = T0, P, k, a, np.degrees(np.arccos(b / a)), [u1, u2]
        return self.m.light_curve(p)
