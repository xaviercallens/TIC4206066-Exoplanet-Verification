"""Audit 9: transit-model machinery (own code; batman, Kreidberg 2015). Phase-folded windows from c02_lc_prep.py.
Model for a planet transiting a host that supplies a fraction fh of the total TESS flux:
    F = 1 - fh * (1 - batman(k, a/R_host, b, u1, u2))
(fh = 1: host = target; fh < 1: host = unresolved companion or target diluted by a companion).
Circular orbit parametrised by the host density: a/R = (G rho P^2 / 3pi)^(1/3). For eccentric orbits the transit
shape at fixed (k, b) is, to first order, that of a circular orbit with an 'apparent' density
rho_app = rho_true * [(1 + e sin w)/sqrt(1 - e^2)]^3 (Kipping 2010; Dawson & Johnson 2012, Eq. 39)."""
from pathlib import Path
import numpy as np
import batman

HERE = Path(__file__).resolve().parent
G = 6.67430e-8; RHO_SUN = 1.98841e33 / (4 / 3 * np.pi * 6.957e10 ** 3)   # g cm^-3 (IAU nominal)
SIGNALS = {"A": dict(P=3.1827856, T0=3998.7375), "B": dict(P=11.13274, T0=4000.865)}
SS = {6: 15, 32: 1, 98: 5}   # supersampling factors (30 min, 2 min, 200 s cadences)

def aR_from_rho(rho_solar, P):
    return (G * rho_solar * RHO_SUN * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)

class Folded:
    def __init__(self, sig="A", npz="lc_windows.npz", sectors=(6, 32, 98), sigma_scale=None):
        d = np.load(HERE / "out" / npz)
        P, T0 = SIGNALS[sig]["P"], SIGNALS[sig]["T0"]
        self.P = P; self.sectors = sectors
        self.t, self.f, self.s, self.models, self.params = {}, {}, {}, {}, {}
        for sec in sectors:
            t = d[f"s{sec}_t{sig}"]; f = d[f"s{sec}_f{sig}"]
            ph = (t - T0 + 0.5 * P) % P - 0.5 * P
            self.t[sec] = ph; self.f[sec] = f
            sc = 1.0 if sigma_scale is None else sigma_scale.get(sec, 1.0)
            self.s[sec] = float(d[f"s{sec}_sig"][0]) * sc
            pm = batman.TransitParams()
            pm.t0 = 0.0; pm.per = P; pm.rp = 0.02; pm.a = 12.; pm.inc = 89.; pm.ecc = 0.; pm.w = 90.
            pm.limb_dark = "quadratic"; pm.u = [0.4, 0.2]
            exp = float(d[f"s{sec}_exp"][0])
            self.models[sec] = batman.TransitModel(pm, ph, supersample_factor=SS[sec], exp_time=exp)
            self.params[sec] = pm
        self.n = sum(len(self.f[s]) for s in sectors)

    def model(self, sec, k, aR, b, u1, u2, fh=1.0, dt=0.0):
        pm = self.params[sec]
        cosi = b / aR
        if not (0 <= cosi < 1):
            return None
        pm.rp = k; pm.a = aR; pm.inc = np.degrees(np.arccos(cosi)); pm.u = [u1, u2]; pm.t0 = dt
        m = self.models[sec].light_curve(pm)
        return 1 - fh * (1 - m)

    def chi2(self, k, aR, b, u1, u2, fh=1.0, dt=0.0, sectors=None):
        c = 0.0
        for sec in (sectors or self.sectors):
            m = self.model(sec, k, aR, b, u1, u2, fh, dt)
            if m is None:
                return 1e30
            c += np.sum(((self.f[sec] - m) / self.s[sec]) ** 2)
        return c

    def resid(self, k, aR, b, u1, u2, fh=1.0, dt=0.0):
        return {sec: self.f[sec] - self.model(sec, k, aR, b, u1, u2, fh, dt) for sec in self.sectors}

def t14_hours(P, aR, k, b):
    """Total duration, circular orbit (Winn 2010 Eq. 14)."""
    inc = np.arccos(b / aR)
    arg = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(inc))
    return P * 24 / np.pi * np.arcsin(np.clip(arg, 0, 1))

# Claret (2017) TESS quadratic LD (LSM, ATLAS) as tabulated in the TRICERATOPS package data file (read as data only)
import pandas as _pd
_LD = _pd.read_csv("/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/lib/python3.11/site-packages/triceratops/data/ldc_tess.csv")
def claret_ld(teff, logg, Z=0.0):
    t = _LD[_LD.Z == _LD.Z.values[np.abs(_LD.Z.values - Z).argmin()]]
    tt = t.Teff.values[np.abs(t.Teff.values - np.clip(teff, 3500, 50000)).argmin()]
    t = t[t.Teff == tt]
    gg = t.logg.values[np.abs(t.logg.values - np.clip(logg, 3.5, 5.0)).argmin()]
    r = t[t.logg == gg].iloc[0]
    return float(r.aLSM), float(r.bLSM), int(tt), float(gg)
