"""Shared constants and helpers for audit 07 (Signal A shape tests and radii). Written from scratch for this audit.

Times are BTJD = BJD_TDB - 2457000 (all products: TIMESYS=TDB, BJDREFI=2457000; see 00_data/README.md sec. 6).
"""
import os
import numpy as np
import batman

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "00_data", "mast")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

# Ephemerides quoted in the manuscript (used only to place windows / masks; P, T0 of A are refitted)
P_A, T0_A = 3.1827856, 2460998.7375 - 2457000.0      # Sec. 2
T14_A_H = 1.87
P_B, T0_B = 11.13274, 2461000.865 - 2457000.0        # Sec. 3
T14_B_H = 1.65

# Stellar parameters (manuscript Sec. 1)
R_STAR, E_R_STAR = 0.623, 0.019
M_STAR, E_M_STAR = 0.606, 0.015
RHO_MED, RHO_EP, RHO_EM = 2.51, 0.25, 0.22
TEFF, LOGG = 4012.0, 4.438 + np.log10(0.606) - 2 * np.log10(0.623)   # logg_sun = 4.438 (IAU nominal) -> 4.632

# IAU 2015 nominal radii (Resolution B3; also astropy.constants R_sun, R_earth)
R_SUN_M, R_EARTH_M = 6.957e8, 6.3781e6
RSUN_OVER_REARTH = R_SUN_M / R_EARTH_M                 # 109.076
G_CGS = 6.67430e-8
RHO_SUN_CGS = 1.98841e33 / (4.0 / 3.0 * np.pi * (6.957e10) ** 3)   # 1.4097 g cm^-3 (nominal M_sun from GM_sun/G)

FILES = {
    6: ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX", 1800.0),
    32: ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX", 120.0),
    98: ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX", 200.0),
}
# supersampling factor per cadence (sub-exposure ~ 2 min for 30-min, 1 min for 2-min, 50 s for 200 s)
SUPERSAMPLE = {6: 15, 32: 2, 98: 4}


def ars_from_rho(rho, P):
    """a/R* for a circular orbit from rho* (solar units) and P (d): a/R* = (G rho P^2 / 3pi)^(1/3)."""
    return (G_CGS * rho * RHO_SUN_CGS * (P * 86400.0) ** 2 / (3.0 * np.pi)) ** (1.0 / 3.0)


def t14_hours(P, k, b, ars):
    """Total duration, circular orbit (Winn 2010, eq. 14)."""
    inc = np.arccos(b / ars)
    arg = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (ars * np.sin(inc))
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1)) * 24.0


def phase_dt(t, P, T0):
    """Time from nearest mid-transit (d) and epoch number."""
    ep = np.round((t - T0) / P)
    return t - (T0 + ep * P), ep.astype(int)


class BatmanSet:
    """One batman model per sector, supersampled over the cadence."""

    def __init__(self, t, cadence_s, ss):
        self.t = np.asarray(t, float)
        pm = batman.TransitParams()
        pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w = 0.0, P_A, 0.02, 12.0, 89.5, 0.0, 90.0
        pm.u, pm.limb_dark = [0.3, 0.3], "quadratic"
        self.pm = pm
        exp = cadence_s / 86400.0
        if ss > 1:
            self.m = batman.TransitModel(pm, self.t, supersample_factor=ss, exp_time=exp)
        else:
            self.m = batman.TransitModel(pm, self.t)

    def flux(self, T0, P, k, b, ars, u1, u2):
        pm = self.pm
        pm.t0, pm.per, pm.rp, pm.a = T0, P, k, ars
        pm.inc = np.degrees(np.arccos(min(b / ars, 1.0)))
        pm.u = [u1, u2]
        return self.m.light_curve(pm)


def claret_ld():
    """Bilinear interpolation of Claret (2017) TESS quadratic LDCs (VizieR J/A+A/600/A30).
    Returns dict with PHOENIX-COND quasi-spherical (Type q, table 15 in the paper; VizieR 'tableab', LSM) and
    ATLAS (table 25, Z=0, xi=2 km/s, LSM)."""
    from astropy.table import Table
    res = {}
    t = Table.read(os.path.join(HERE, "claret", "tableab.ecsv"))
    t = t[(t["Type"] == "q") & (t["Z"] == 0.0)]

    def interp(tab, Tlo, Thi, alab, blab):
        def val(T, g):
            r = tab[(tab["Teff"] == T) & (np.isclose(tab["logg"], g))]
            assert len(r) == 1, (T, g, len(r))
            return np.array([float(r[alab][0]), float(r[blab][0])])
        gw = (LOGG - 4.5) / 0.5
        tw = (TEFF - Tlo) / (Thi - Tlo)
        vlo = val(Tlo, 4.5) * (1 - gw) + val(Tlo, 5.0) * gw
        vhi = val(Thi, 4.5) * (1 - gw) + val(Thi, 5.0) * gw
        return vlo * (1 - tw) + vhi * tw
    res["phoenix_q_LSM"] = interp(t, 4000.0, 4100.0, "aLSM", "bLSM").tolist()
    t2 = Table.read(os.path.join(HERE, "claret", "table25.ecsv"))
    t2 = t2[(t2["Z"] == 0.0) & (t2["xi"] == 2.0)]
    res["atlas_Z0_xi2_LSM"] = interp(t2, 4000.0, 4250.0, "aLSM", "bLSM").tolist()
    return res


def pct(x, w=None, q=(2.275, 15.865, 50.0, 84.135, 97.725)):
    x = np.asarray(x)
    if w is None:
        return np.percentile(x, q).tolist()
    o = np.argsort(x)
    cw = np.cumsum(w[o]) / np.sum(w)
    return [float(np.interp(qq / 100.0, cw, x[o])) for qq in q]
