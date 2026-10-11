"""Verifier (audit 08): independent S32 (SPOC 2-min PDCSAP) + S98 (QLP SAP) folds for TRICERATOPS.
Own code; choices deliberately differ from the auditor's f01 (window 1.0 d, mask +-0.75*T14, 5-min bins,
S98 quality: drop any nonzero QUALITY). Run with audit_env.
"""
import json, os, hashlib
import numpy as np
from astropy.io import fits
from wotan import flatten

H = os.path.dirname(os.path.abspath(__file__))
MAST = os.path.join(H, "..", "..", "00_data", "mast")
OUT = os.path.join(H, "folds"); os.makedirs(OUT, exist_ok=True)
EPH = {"A": (3.1827856, 3998.7375, 1.87 / 24), "B": (11.13274, 4000.865, 1.65 / 24)}
SRC = {"S32": ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX"),
       "S98": ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SAP_FLUX")}


def ph(t, s, dT0=0.0):
    P, T0, _ = EPH[s]
    return (t - (T0 + dT0) + 0.5 * P) % P - 0.5 * P


meta = {"eph": EPH, "inputs": {}, "folds": {}}
LC = {}
for sec, (fn, col) in SRC.items():
    p = os.path.join(MAST, fn)
    meta["inputs"][sec] = dict(file=fn, col=col, sha256=hashlib.sha256(open(p, "rb").read()).hexdigest())
    d = fits.getdata(p, 1)
    t = np.array(d["TIME"], float); f = np.array(d[col], float); q = np.array(d["QUALITY"])
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    intr = np.zeros(t.size, bool)
    for s in EPH:
        intr |= np.abs(ph(t, s)) < 0.75 * EPH[s][2]
    ff = f.copy(); ff[intr] = np.nan          # wotan ignores NaN in the trend fit -> transits masked
    _, tr = flatten(t, ff, method="biweight", window_length=1.0, edge_cutoff=0.0, break_tolerance=0.5, return_trend=True)
    good = np.isfinite(tr)
    tr = np.interp(t, t[good], tr[good])
    LC[sec] = (t, f / tr)
    meta["inputs"][sec]["n_good"] = int(t.size)


def fold(sig, secs, binmin=5.0, dT0=0.0):
    other = "B" if sig == "A" else "A"
    xs, ys = [], []
    for s in secs:
        t, f = LC[s]
        m = (np.abs(ph(t, other)) > 0.75 * EPH[other][2]) & (np.abs(ph(t, sig, dT0)) < 0.25)
        xs.append(ph(t, sig, dT0)[m]); ys.append(f[m])
    x = np.concatenate(xs); y = np.concatenate(ys)
    e = np.arange(-0.25, 0.25 + 1e-9, binmin / 1440.0)
    k = np.digitize(x, e)
    rows = np.array([(x[k == i].mean(), np.median(y[k == i]) if False else y[k == i].mean(), y[k == i].size)
                     for i in range(1, e.size) if (k == i).sum() >= 3])
    T14 = EPH[sig][2]
    oot = np.abs(rows[:, 0]) > 0.8 * T14
    sig_oot = float(np.std(rows[oot, 1], ddof=1))
    inn = np.abs(rows[:, 0]) < 0.3 * T14
    depth = float(1 - rows[inn, 1].mean())
    name = f"v_{sig}_{'+'.join(secs)}_{int(binmin)}min" + ("" if dT0 == 0 else f"_dT{dT0:+.4f}")
    np.savetxt(os.path.join(OUT, name + ".csv"), rows[:, :2], delimiter=",", header="dt_d,flux")
    meta["folds"][name] = dict(sig=sig, secs=secs, binmin=binmin, n_pts=int(x.size), n_bins=int(len(rows)),
                               sigma_oot=sig_oot, depth_ppm=depth * 1e6)
    print(name, x.size, len(rows), "sigma_oot %.0f ppm depth %.0f ppm" % (sig_oot * 1e6, depth * 1e6))


for s in ("A", "B"):
    fold(s, ["S32", "S98"])
    fold(s, ["S32", "S98"], binmin=10.0)
json.dump(meta, open(os.path.join(OUT, "folds_meta.json"), "w"), indent=1)
