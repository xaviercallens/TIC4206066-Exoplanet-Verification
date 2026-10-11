"""Verifier check: Lomb-Scargle of the Sector 32 SPOC 2-min light curve (PDCSAP and SAP), read-only from
deep_followup_tic4206066/lc, to test whether the ~8 d, ~1 mmag catalogue signal is visible in the project data.
Transit windows (A: P=3.182785, T0=2460998.736-2457000; B: P=11.13274, T0=2461000.865-2457000, +-1.5 h) masked."""
import json, numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.timeseries import LombScargle
H = Path(__file__).resolve().parent
f = Path("/Users/paulsmith/space <3/tess-jev/addenda/deep_followup_tic4206066/lc/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits")
d = fits.getdata(f, 1)
t = d["TIME"]; q = d["QUALITY"]
out = {"file": str(f)}
def mask_tr(t):
    m = np.ones_like(t, bool)
    for P, T0 in ((3.182785, 1998.736), (11.13274, 2000.865)):
        ph = ((t - T0 + 0.5 * P) % P) - 0.5 * P
        m &= np.abs(ph) > 1.5 / 24
    return m
for col in ("PDCSAP_FLUX", "SAP_FLUX"):
    y = d[col]; ok = np.isfinite(t) & np.isfinite(y) & (q == 0)
    ok &= mask_tr(np.where(np.isfinite(t), t, 0))
    tt, yy = t[ok], y[ok] / np.nanmedian(y[ok]) - 1
    # bin to 30 min for speed
    b = np.floor((tt - tt[0]) / (30 / 1440)).astype(int)
    tb = np.bincount(b, tt) / np.bincount(b); yb = np.bincount(b, yy) / np.bincount(b); nb = np.bincount(b)
    g = nb > 5; tb, yb = tb[g], yb[g]
    fr = np.linspace(1 / 20, 1 / 0.5, 20000)
    ls = LombScargle(tb, yb); p = ls.power(fr)
    i = np.argmax(p); per = 1 / fr[i]
    f8 = 1 / 8.0; j = np.argmin(np.abs(fr - f8))
    mdl = ls.model(tb, fr[i])
    # power / FAP at 7.9-8.3 d
    sel = (1 / fr > 7.8) & (1 / fr < 8.4)
    k = np.argmax(np.where(sel, p, 0))
    out[col] = dict(n=int(tb.size), span_d=float(np.ptp(tb)), rms_ppm=float(np.std(yb) * 1e6),
                    top_period_d=float(per), top_power=float(p[i]), top_fap_baluev=float(ls.false_alarm_probability(p[i])),
                    top_semi_amp_ppm=float((mdl.max() - mdl.min()) / 2 * 1e6),
                    best_7p8_8p4_period=float(1 / fr[k]), best_7p8_8p4_power=float(p[k]),
                    best_7p8_8p4_fap=float(ls.false_alarm_probability(p[k])),
                    semi_amp_at_that_ppm=float((ls.model(tb, fr[k]).max() - ls.model(tb, fr[k]).min()) / 2 * 1e6))
json.dump(out, open(H / "v02_s32_rotation_ls.json", "w"), indent=1); print(json.dumps(out, indent=1))
