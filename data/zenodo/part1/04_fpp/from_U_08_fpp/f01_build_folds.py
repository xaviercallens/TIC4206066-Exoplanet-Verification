"""Audit 8: build my own detrended, folded S32 (SPOC 2-min PDCSAP) + S98 (QLP) light curves for TRICERATOPS.
Run with audit_env. Inputs: fresh 00_data/mast files only. Writes folds/*.csv and folds/folds_meta.json.
"""
import json, os, hashlib
import numpy as np
from astropy.io import fits
from wotan import flatten

H = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(H, "..", "00_data", "mast")
OUT = os.path.join(H, "folds"); os.makedirs(OUT, exist_ok=True)
EPH = {"A": dict(P=3.1827856, T0=3998.7375, T14=1.87/24),
       "B": dict(P=11.13274, T0=4000.865, T14=1.65/24)}
FILES = {"S32": "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits",
         "S98": "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"}

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

def dphase(t, sig, T0=None):
    e = EPH[sig]; T0 = e["T0"] if T0 is None else T0
    return ((t - T0) / e["P"] + 0.5) % 1 * e["P"] - 0.5 * e["P"]

def load(name):
    p = os.path.join(DATA, FILES[name]); d = fits.open(p)[1].data
    if name == "S32":
        t, f, q = d["TIME"], d["PDCSAP_FLUX"], d["QUALITY"]
        cad = 120.0
    else:
        t, f, q = d["TIME"], d["SAP_FLUX"], d["QUALITY"]   # QLP raw best-aperture flux; own detrending
        cad = 200.0
    m = np.isfinite(t) & np.isfinite(f) & (q == 0)        # any flagged cadence dropped (removes S98 straylight bit 2048 and QLP bits 29/30)
    t, f = np.asarray(t[m], float), np.asarray(f[m], float)
    f = f / np.nanmedian(f)
    intr = np.zeros(len(t), bool)
    for s in EPH:
        intr |= np.abs(dphase(t, s)) < EPH[s]["T14"]      # +-1 T14 around each centre
    trend_oot = flatten(t[~intr], f[~intr], method="biweight", window_length=0.75, break_tolerance=0.5, return_trend=True)[1]
    ok = np.isfinite(trend_oot)
    trend = np.interp(t, t[~intr][ok], trend_oot[ok])
    fd = f / trend
    sig = 1.4826 * np.nanmedian(np.abs(fd[~intr] - np.nanmedian(fd[~intr])))
    return dict(t=t, fd=fd, sig=sig, cad=cad, sha=sha(p), n=len(t))

def build(sig, sectors, binmin=6.0, T0shift=0.0, tag=""):
    other = "B" if sig == "A" else "A"
    T, F, E = [], [], []
    for s in sectors:
        L = LC[s]
        keep = np.abs(dphase(L["t"], other)) >= EPH[other]["T14"]
        ph = dphase(L["t"], sig, EPH[sig]["T0"] + T0shift)
        w = keep & (np.abs(ph) < 0.25)
        T.append(ph[w]); F.append(L["fd"][w]); E.append(np.full(w.sum(), L["sig"]))
    t, f, e = map(np.concatenate, (T, F, E))
    o = np.argsort(t); t, f, e = t[o], f[o], e[o]
    if binmin is None:
        rows = np.c_[t, f, e]
    else:
        edges = np.arange(-0.25, 0.25 + 1e-9, binmin / 1440); idx = np.digitize(t, edges); rows = []
        for i in range(1, len(edges)):
            m = idx == i
            if m.sum() < 3: continue
            ww = 1 / e[m] ** 2
            rows.append((t[m].mean(), np.sum(ww * f[m]) / ww.sum(), 1 / np.sqrt(ww.sum())))
        rows = np.array(rows)
    T14 = EPH[sig]["T14"]
    oot = np.abs(rows[:, 0]) > 0.75 * T14
    core = np.abs(rows[:, 0]) < 0.25 * T14
    meta = dict(signal=sig, sectors=sectors, bin_min=binmin, T0shift_d=T0shift, n_raw=int(len(t)), n_rows=int(len(rows)),
                sigma_empirical_oot=float(np.std(rows[oot, 1], ddof=1)), sigma_formal_mean=float(np.mean(rows[:, 2])),
                depth_core_ppm=float((1 - np.average(rows[core, 1], weights=1 / rows[core, 2] ** 2)) * 1e6))
    name = f"fold_{sig}_{'+'.join(sectors)}_{'unb' if binmin is None else str(int(binmin))+'min'}{tag}"
    np.savetxt(os.path.join(OUT, name + ".csv"), rows, delimiter=",", header="t_from_mid_d,flux,flux_err")
    meta["file"] = name + ".csv"
    return name, meta

LC = {s: load(s) for s in FILES}
META = dict(ephemerides=EPH, inputs={s: dict(file=FILES[s], sha256=LC[s]["sha"], n_good=LC[s]["n"], sigma_point=LC[s]["sig"]) for s in LC}, folds={})
for sig in ("A", "B"):
    for sec in (["S32", "S98"], ["S32"], ["S98"]):
        n, m = build(sig, sec); META["folds"][n] = m
    for b in (2.0, 10.0):
        n, m = build(sig, ["S32", "S98"], binmin=b); META["folds"][n] = m
    n, m = build(sig, ["S32", "S98"], binmin=None); META["folds"][n] = m
    # epoch sensitivity: shift T0 by +-1 sigma-ish (A: 0.0005 d ~ 0.7 min; B: 0.003 d)
    dT = 0.0005 if sig == "A" else 0.003
    for s_, tg in ((+dT, "_T0p"), (-dT, "_T0m")):
        n, m = build(sig, ["S32", "S98"], T0shift=s_, tag=tg); META["folds"][n] = m
json.dump(META, open(os.path.join(OUT, "folds_meta.json"), "w"), indent=1)
for k, v in META["folds"].items():
    print(k, v["n_raw"], v["n_rows"], "sig_emp=%.0f ppm sig_form=%.0f ppm depth=%.0f ppm" % (v["sigma_empirical_oot"]*1e6, v["sigma_formal_mean"]*1e6, v["depth_core_ppm"]))
