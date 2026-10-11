"""Build detrended light curves for S6 (TESS-SPOC 30-min PDCSAP), S32 (SPOC 2-min PDCSAP), S98 (QLP 200-s SYS_RM_FLUX).

Quality: SPOC/TESS-SPOC keep QUALITY == 0 (SDPDD Rev F Table 32). QLP S98 baseline keeps QUALITY == 0, i.e. drops
bit 12 (2048, predicted Earth/Moon straylight), 2^30 (QLP internal poor-quality flag; tess.mit.edu/qlp), 32/4
(desaturation / coarse point) and the undocumented 2^29 bit, which in this file brackets the straylight intervals at the
start/end of each orbit.  Variant 'b29': additionally keep cadences whose ONLY flag is 2^29.
Detrending: wotan 1.10 biweight (Hippke et al. 2019), window 0.75 d, break_tolerance 0.5 d, both signals masked
(+-1.0 T14 around mid-transit, i.e. a 2*T14 wide mask).  Upward 4-sigma (MAD) outliers removed after detrending.
Variant 'mask05': additionally mask phase 0.5 of A (for the secondary-eclipse test).
"""
import sys, json, hashlib
import numpy as np
from astropy.io import fits
from wotan import flatten, transit_mask
from common import *

VARIANT = sys.argv[1] if len(sys.argv) > 1 else "base"   # base | b29 | mask05 | win05 | det98 | raw
WINDOW = 0.5 if VARIANT == "win05" else 0.75


def load(sector):
    fn, col, cad = FILES[sector]
    path = os.path.join(DATA, fn)
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    d = fits.open(path)[1].data
    t = np.asarray(d["TIME"], float)
    if sector == 98 and VARIANT == "det98":
        col = "DET_FLUX"
    f = np.asarray(d[col], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    if sector == 98 and VARIANT == "b29":
        good = (q == 0) | (q == (1 << 29))
    else:
        good = q == 0
    good &= np.isfinite(t) & np.isfinite(f)
    return t[good], f[good] / np.nanmedian(f[good]), cad, dict(file=fn, sha256=sha, column=col,
                                                              n_total=int(len(t)), n_kept=int(good.sum()))


def main():
    info = {}
    for s in (6, 32, 98):
        t, f, cad, meta = load(s)
        mA = transit_mask(time=t, period=P_A, duration=2 * T14_A_H / 24.0, T0=T0_A)
        mB = transit_mask(time=t, period=P_B, duration=2 * T14_B_H / 24.0, T0=T0_B)
        m = mA | mB
        if VARIANT == "mask05":
            m |= transit_mask(time=t, period=P_A, duration=2 * T14_A_H / 24.0, T0=T0_A + 0.5 * P_A)
        if VARIANT == "raw":      # no detrending: per-window linear baselines are profiled in the joint fit
            fl, tr = f.copy(), np.ones_like(f)
        else:
            fl, tr = flatten(t, f, method="biweight", window_length=WINDOW, break_tolerance=0.5, return_trend=True,
                             mask=m)
        ok = np.isfinite(fl)
        r = fl[ok & ~m] - 1.0
        mad = 1.4826 * np.median(np.abs(r - np.median(r)))
        if VARIANT != "raw":
            ok &= ~((fl - 1.0) > 4 * mad)   # upward outliers only
        np.savez(os.path.join(OUT, f"lc_{VARIANT}_S{s}.npz"), t=t[ok], f=fl[ok], trend=tr[ok], raw=f[ok], cad=cad)
        meta.update(n_after_detrend=int(ok.sum()), mad_ppm=float(mad * 1e6), window_d=WINDOW)
        info[s] = meta
        print(s, meta)
    json.dump(info, open(os.path.join(OUT, f"prep_{VARIANT}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
