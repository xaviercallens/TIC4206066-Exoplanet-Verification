"""Sector 98 light-curve variants used by the S98-only search and null tests (A removed)."""
import numpy as np
import common_b as c

B_SEGS = [(3995.80, 4002.80), (4010.02, 4016.89), (4024.51, 4031.60), (4038.83, 4045.91)]   # b-half orbits (DRN134)


def in_bsegs(t):
    m = np.zeros_like(t, bool)
    for a, b in B_SEGS:
        m |= (t >= a) & (t <= b)
    return m


def inB(t, P=c.P_B98, T0=c.T0_B98, half=None):
    half = c.T14_B if half is None else half
    ph = ((t - T0) / P + 0.5) % 1.0 - 0.5
    return np.abs(ph * P) < half


def get(variant):
    """Returns t, f (normalised, A removed, finite only)."""
    if variant.startswith("QLP"):
        q = c.load_qlp(98)
        t, qq = q["t"], q["q"]
        clean = (qq == 0)
        if variant == "QLP_DET_raw":
            f = q["f"]
        elif variant == "QLP_DET_bw075":
            m = clean & np.isfinite(q["f"])
            f = np.full_like(t, np.nan)
            f[m] = c.detrend(t[m], q["f"][m], "biweight", 0.75, mask=c.in_A(t[m]))
        elif variant == "QLP_SYSRM_bw075":
            m = clean & np.isfinite(q["sysrm"])
            f = np.full_like(t, np.nan)
            f[m] = c.detrend(t[m], q["sysrm"][m], "biweight", 0.75, mask=c.in_A(t[m]))
        else:
            raise ValueError(variant)
        m = clean & np.isfinite(f) & ~c.in_A(t)
    elif variant.startswith("OWN"):
        p = np.load(c.HERE + "/out/ownphot_s98_ap12.npz")
        t, fl, qq, bkg = p["t"], p["flux"], p["q"], p["bkg"]
        good = np.isfinite(fl) & np.isfinite(t) & ((qq & (c.SPOC_BAD | 2048)) == 0)
        if variant == "OWN_bhalves_bw075":
            sel = good & in_bsegs(t)
        elif variant == "OWN_bkgcut_bw075":
            ref = np.median(bkg[good & in_bsegs(t)])
            sel = good & (bkg < 1.5 * ref)
        else:
            raise ValueError(variant)
        f = np.full_like(t, np.nan)
        f[sel] = c.detrend(t[sel], fl[sel] / np.nanmedian(fl[sel]), "biweight", 0.75, mask=c.in_A(t[sel]))
        # 5-sigma upper outlier clip only (keeps dips)
        r = c.robust_std(f[sel])
        f[sel & (f > 1 + 5 * r)] = np.nan
        m = np.isfinite(f) & ~c.in_A(t)
    else:
        raise ValueError(variant)
    return t[m], f[m]


def period_grid(T, pmin=1.0, pmax=20.0, dmin=1.0 / 24, os_=5):
    P = [pmin]
    while P[-1] < pmax:
        P.append(P[-1] + dmin * P[-1] ** 2 / (T * os_))
    return np.array(P)


DURS = np.array([1.0, 1.3, 1.6, 2.0, 2.5, 3.0]) / 24.0


def bls_spectrum(t, f, periods, dy=None):
    from astropy.timeseries import BoxLeastSquares
    dy = c.robust_std(f) if dy is None else dy
    b = BoxLeastSquares(t, f, dy=np.full_like(f, dy))
    r = b.power(periods, DURS, objective="likelihood")
    pw = np.asarray(r.power)
    snr = np.asarray(r.depth) / np.asarray(r.depth_err)
    return r, pw, snr


def sde(pw):
    return (pw - np.mean(pw)) / np.std(pw)
