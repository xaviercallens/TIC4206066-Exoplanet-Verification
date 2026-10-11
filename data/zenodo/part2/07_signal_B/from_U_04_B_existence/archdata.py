"""Archival (Sector 6 and 32) light-curve variants, A removed, detrended; box-depth grids for the link test."""
import numpy as np
import common_b as c

NOMINAL_CAD = {"s6": 1800 / 86400., "s32_2min": 120 / 86400., "s32_10min": 600 / 86400.}


def get(sector, phot, window=0.75, method="biweight"):
    """Return t, f (normalised, detrended, A removed).  phot in {'pipe','qlp','own','own05','eleanor','tglc_psf','tglc_aper','sap'}."""
    if sector == 6:
        if phot in ("pipe", "sap"):
            d = c.load_spoc_lc("s6_tessspoc")
            t, f = d["t"], (d["pdc"] if phot == "pipe" else d["sap"])
            good = np.isfinite(t) & np.isfinite(f) & ((d["q"] & c.SPOC_BAD) == 0)
        elif phot == "qlp":
            d = c.load_qlp(6)
            t, f = d["t"], d["f"]
            good = np.isfinite(t) & np.isfinite(f) & (d["q"] == 0)
        elif phot in ("own", "own05"):
            p = np.load(c.HERE + "/out/ownphot_s6_%s.npz" % ("ap12" if phot == "own" else "ap05"))
            t, f, q = p["t"], p["flux"], p["q"]
            good = np.isfinite(t) & np.isfinite(f) & ((q & (c.SPOC_BAD | 2048)) == 0)
        elif phot == "eleanor":
            d = c.load_eleanor6()
            t, f = d["t"], d["corr"]
            good = np.isfinite(t) & np.isfinite(f) & (d["q"] == 0)
        elif phot in ("tglc_psf", "tglc_aper"):
            d = c.load_tglc6()
            t = d["t"]
            f = np.asarray(d["cal_psf_flux"] if phot == "tglc_psf" else d["cal_aper_flux"], float)
            good = np.isfinite(t) & np.isfinite(f) & (np.asarray(d["TESS_flags"]) == 0) & (np.asarray(d["TGLC_flags"]) == 0)
        else:
            raise ValueError(phot)
    elif sector == 32:
        if phot in ("pipe", "sap"):
            d = c.load_spoc_lc("s32_2min")
            t, f = d["t"], (d["pdc"] if phot == "pipe" else d["sap"])
            good = np.isfinite(t) & np.isfinite(f) & ((d["q"] & c.SPOC_BAD) == 0)
        elif phot == "pipe10":
            d = c.load_spoc_lc("s32_tessspoc")
            t, f = d["t"], d["pdc"]
            good = np.isfinite(t) & np.isfinite(f) & ((d["q"] & c.SPOC_BAD) == 0)
        elif phot == "qlp":
            d = c.load_qlp(32)
            t, f = d["t"], d["f"]
            good = np.isfinite(t) & np.isfinite(f) & (d["q"] == 0)
        elif phot in ("own", "own05"):
            p = np.load(c.HERE + "/out/ownphot_s32_%s.npz" % ("ap12" if phot == "own" else "ap05"))
            t, f, q = p["t"], p["flux"], p["q"]
            good = np.isfinite(t) & np.isfinite(f) & ((q & (c.SPOC_BAD | 2048)) == 0)
        else:
            raise ValueError(phot)
    t, f = t[good], f[good] / np.nanmedian(f[good])
    inA = c.in_A(t)
    fd = c.detrend(t, f, method, window, mask=inA)
    r = c.robust_std(fd)
    keep = np.isfinite(fd) & ~inA & (fd < 1 + 5 * r)   # remove A and upward outliers only
    return t[keep], fd[keep]


class BoxGrid:
    """Box depth and calibrated uncertainty on a dense grid of trial centres."""

    def __init__(self, t, f, dur=c.T14_B, step=0.5 / 1440, min_cov=0.5):
        self.t, self.f = t, f
        self.dur = dur
        self.g0 = t.min() - 0.2
        self.grid = np.arange(self.g0, t.max() + 0.2, step)
        self.step = step
        d, nin, nout = c.box_depth_grid(t, f, self.grid, dur)
        cad = np.median(np.diff(np.sort(t)))
        n_nom = dur / cad
        sig_pt = c.robust_std(f - 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            sw = sig_pt * np.sqrt(1 / nin + 1 / nout)
        ok = (nin >= min_cov * n_nom) & (nout >= min_cov * 2 * n_nom)
        z = d[ok] / sw[ok]
        self.beta = max(1.0, c.robust_std(z))
        self.depth = np.where(ok, d, np.nan)
        self.sig = np.where(ok, sw * self.beta, np.nan)
        self.sig_pt, self.cad, self.n_nom = sig_pt, cad, n_nom
        self.tmin, self.tmax = t.min(), t.max()

    def lookup(self, tc):
        i = np.rint((tc - self.g0) / self.step).astype(np.int64)
        inside = (i >= 0) & (i < len(self.grid))
        i = np.clip(i, 0, len(self.grid) - 1)
        d = np.where(inside, self.depth[i], np.nan)
        s = np.where(inside, self.sig[i], np.nan)
        return d, s

    def snr_series(self):
        return self.grid, self.depth / self.sig


def combined_snr(grids, Tref, Pvec):
    """grids: list of BoxGrid; Tref scalar; Pvec array of periods.  Returns combined S/N per period (weighted mean depth /
    its error over all predicted events with coverage), and number of events used."""
    num = np.zeros(len(Pvec)); den = np.zeros(len(Pvec)); nev = np.zeros(len(Pvec), int)
    for g in grids:
        a = (g.tmin - 0.2 - Tref) / Pvec
        b = (g.tmax + 0.2 - Tref) / Pvec
        nlo = int(np.floor(min(a.min(), b.min()))) - 1
        nhi = int(np.ceil(max(a.max(), b.max()))) + 1
        for n in range(nlo, nhi + 1):
            tc = Tref + n * Pvec
            m = (tc > g.tmin - 0.1) & (tc < g.tmax + 0.1)
            if not m.any():
                continue
            d, s = g.lookup(tc)
            ok = m & np.isfinite(d) & np.isfinite(s)
            num[ok] += d[ok] / s[ok] ** 2
            den[ok] += 1 / s[ok] ** 2
            nev[ok] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = num / np.sqrt(den)
    snr[den == 0] = 0.0
    return snr, nev
