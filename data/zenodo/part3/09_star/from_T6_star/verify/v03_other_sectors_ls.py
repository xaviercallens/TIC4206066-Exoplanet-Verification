"""Verifier check: same LS test on Sector 6 (TESS-SPOC 30-min PDCSAP/SAP) and Sector 98 (QLP SAP_FLUX / KSPSAP_FLUX),
and on the S32 PDCSAP split in halves (orbit 1 / orbit 2), to see whether the 7.97 d S32 signal recurs.
Read-only inputs from deep_followup_tic4206066/lc."""
import json, numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.timeseries import LombScargle
H = Path(__file__).resolve().parent
L = Path("/Users/paulsmith/space <3/tess-jev/addenda/deep_followup_tic4206066/lc")
def mask_tr(t):
    m = np.ones_like(t, bool)
    for P, T0 in ((3.182785, 1998.736), (11.13274, 2000.865)):
        ph = ((t - T0 + 0.5 * P) % P) - 0.5 * P
        m &= np.abs(ph) > 1.5 / 24
    return m
def ls_test(t, y, lab):
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]; m = mask_tr(t); t, y = t[m], y[m] / np.nanmedian(y[m]) - 1
    b = np.floor((t - t[0]) / (30 / 1440)).astype(int); nb = np.bincount(b); g = nb > 0
    tb = (np.bincount(b, t)[g] / nb[g]); yb = (np.bincount(b, y)[g] / nb[g])
    fr = np.linspace(1 / 20, 1 / 0.5, 20000); ls = LombScargle(tb, yb); p = ls.power(fr); i = np.argmax(p)
    sel = (1 / fr > 7.5) & (1 / fr < 8.7); k = np.argmax(np.where(sel, p, 0))
    amp = lambda fq: float((ls.model(tb, fq).max() - ls.model(tb, fq).min()) / 2 * 1e6)
    return dict(label=lab, n=int(tb.size), span_d=float(np.ptp(tb)), rms_ppm=float(np.std(yb) * 1e6), top_period=float(1 / fr[i]),
                top_power=float(p[i]), top_semi_amp_ppm=amp(fr[i]), p7p5_8p7_period=float(1 / fr[k]), p7p5_8p7_power=float(p[k]),
                p7p5_8p7_semi_amp_ppm=amp(fr[k]), rank_of_8d_peak=int(np.sum(p > p[k])))
out = []
d = fits.getdata(L / "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", 1)
q = d["QUALITY"] == 0
for c in ("PDCSAP_FLUX", "SAP_FLUX"):
    out.append(ls_test(d["TIME"][q], d[c][q], f"S6 TESS-SPOC {c}"))
d = fits.getdata(L / "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", 1)
q = d["QUALITY"] == 0
for c in ("SAP_FLUX", "KSPSAP_FLUX"):
    if c in d.columns.names: out.append(ls_test(d["TIME"][q], d[c][q], f"S98 QLP {c}"))
d = fits.getdata(L / "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", 1)
q = d["QUALITY"] == 0; t = d["TIME"]; mid = np.nanmedian(t[q])
gap = t[q][np.argmax(np.diff(t[q]))]
for lab, s in (("S32 PDCSAP orbit1", q & (t <= gap)), ("S32 PDCSAP orbit2", q & (t > gap))):
    out.append(ls_test(t[s], d["PDCSAP_FLUX"][s], lab))
json.dump(out, open(H / "v03_other_sectors_ls.json", "w"), indent=1)
for o in out: print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items()})
