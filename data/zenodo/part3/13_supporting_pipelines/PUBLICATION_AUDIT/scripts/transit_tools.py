"""Transit model + data-window helpers for the audit (batman quadratic LD, circular orbit).

Conventions: rho_star in solar units (rho_sun = 1.40978 g/cm^3 for M_sun = 1.98841e33 g, R_sun = 6.957e10 cm);
a/R* = (G rho P^2 / 3 pi)^(1/3); depth in relative flux; times in BTJD (BJD_TDB - 2457000).
"""
import numpy as np
import batman

import lcio

G_CGS = 6.6743e-8
RHO_SUN = 1.98841e33 / (4.0 / 3.0 * np.pi * 6.957e10 ** 3)

# Claret (2017) J/A+A/600/A30 table15 (quasi-spherical PHOENIX-COND, TESS, LSM), bilinear in Teff/logg
_LD = {(4000, 4.5): (0.3390, 0.2994), (4000, 5.0): (0.2260, 0.3679), (4100, 4.5): (0.3955, 0.2618),
       (4100, 5.0): (0.2673, 0.3371), (3900, 4.5): (0.2946, 0.3316), (3900, 5.0): (0.1991, 0.3902)}


def claret_ld(teff, logg):
    T0 = 3900 if teff < 4000 else 4000
    T1 = T0 + 100
    tw = (teff - T0) / 100.0
    gw = (logg - 4.5) / 0.5
    v = lambda T: np.array(_LD[(T, 4.5)]) * (1 - gw) + np.array(_LD[(T, 5.0)]) * gw
    return tuple(v(T0) * (1 - tw) + v(T1) * tw)


def ars_from_rho(rho, P):
    return (rho * RHO_SUN * G_CGS * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1.0 / 3.0)


def rho_from_ars(ars, P):
    return 3 * np.pi * ars ** 3 / (G_CGS * (P * 86400.0) ** 2) / RHO_SUN


def t14_hours(P, k, b, ars):
    arg = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (ars * np.sqrt(np.clip(1 - (b / ars) ** 2, 1e-12, None)))
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1)) * 24.0


class Model:
    """Caches one batman TransitModel per dataset (supersampling set from the cadence)."""

    def __init__(self, t, cad_days):
        self.t = t
        self.cad = cad_days
        self.ss = max(1, int(np.ceil(cad_days * 1440 / 2.0)))
        pm = batman.TransitParams()
        pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w = 0.0, 3.0, 0.02, 12.0, 89.0, 0.0, 90.0
        pm.u, pm.limb_dark = [0.3, 0.3], "quadratic"
        self.pm = pm
        if self.ss > 1:
            self.m = batman.TransitModel(pm, t, supersample_factor=self.ss, exp_time=cad_days)
        else:
            self.m = batman.TransitModel(pm, t)

    def lc(self, t0, P, k, b, ars, u1, u2):
        if b >= ars:
            return np.ones_like(self.t)
        pm = self.pm
        pm.t0, pm.per, pm.rp, pm.a = t0, P, k, ars
        pm.inc = np.degrees(np.arccos(b / ars))
        pm.u = [u1, u2]
        return self.m.light_curve(pm)


def windows(lcs, sig, eph, half=0.30, other_fac=1.0):
    """Cut +-half d windows around every transit of `sig`, dropping points within other_fac*T14 of the other signal.
    lcs: list of detrended light-curve dicts (lcio.detrend output)."""
    other = "B" if sig == "A" else "A"
    e = eph[sig]
    W = []
    for lc in lcs:
        t, f = lc["t"], lc["fd"]
        ph = ((t - e["T0"]) / e["P"] + 0.5) % 1 - 0.5
        m = (np.abs(ph * e["P"]) < half) & ~lcio.in_transit(t, other, other_fac, eph)
        if m.sum() < 10:
            continue
        W.append(dict(key=lc["key"], sector=lc["sector"], t=t[m], f=f[m], sig=lc["sig"], cad=lc["cad"],
                      epoch=np.round((t[m] - e["T0"]) / e["P"]).astype(int)))
    return W


def beta_red(res, cad, scales_min=(20, 40, 60, 90, 120)):
    """Time-averaging red-noise factor (Pont+2006 / Winn+2008): beta = sigma_N / (sigma_1/sqrt(N) sqrt(M/(M-1))),
    averaged over bin sizes spanning ingress..transit timescales. Returns max(1, median beta)."""
    s1 = np.std(res)
    betas = []
    for sm in scales_min:
        n = max(1, int(round(sm / (cad * 1440))))
        M = len(res) // n
        if M < 5 or n < 2:
            continue
        rb = res[: M * n].reshape(M, n).mean(axis=1)
        exp = s1 / np.sqrt(n) * np.sqrt(M / (M - 1.0))
        betas.append(np.std(rb) / exp)
    return float(max(1.0, np.median(betas))) if betas else 1.0
