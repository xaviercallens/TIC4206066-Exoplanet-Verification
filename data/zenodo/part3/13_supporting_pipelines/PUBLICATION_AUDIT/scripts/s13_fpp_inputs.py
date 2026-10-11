"""Step 13a: clean TRICERATOPS inputs built from this audit's own products (audit_env).
Folds: pipeline light curves (S32 SPOC 2-min PDCSAP, S98 QLP SYS_RM_FLUX with the frozen mask = clean intervals,
S6 TESS-SPOC PDCSAP separately), detrended with BOTH signals masked, the OTHER signal's transits removed (+-1.0 T14),
folded on the audited ephemeris (out/ephem_current.json), +-0.25 d, binned to 6 min (inverse-variance).
Apertures (absolute CCD col,row as TRICERATOPS expects): SPOC optimal-aperture bit (value 2) of the S6 TESS-SPOC TP
and the S32 SPOC 2-min TPF; S98: 3x3 pixels centred on the target pixel (QLP BESTAP = '3x3 square').
"""
import glob
import json
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
from audit_common import AUD, save_json  # noqa

EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
OUT = os.path.join(AUD, "out", "fpp", "inputs")
os.makedirs(OUT, exist_ok=True)
LCS = {
    "S6": lcio.detrend(lcio.load_product("TESSSPOC_S6_PDC"), window=0.75, eph=EPH),
    "S32": lcio.detrend(lcio.load_product("SPOC2min_S32_PDC"), window=0.75, eph=EPH),
    "S98": lcio.detrend(lcio.load_product("QLP_S98_SYS", "frozen"), window=0.75, eph=EPH),
}


def fold(sig, sets, binmin=6.0, half=0.25):
    e = EPH[sig]
    other = "B" if sig == "A" else "A"
    T, F, E = [], [], []
    for s in sets:
        L = LCS[s]
        keep = ~lcio.in_transit(L["t"], other, 1.0, EPH)
        ph = (((L["t"] - e["T0"]) / e["P"] + 0.5) % 1 - 0.5) * e["P"]
        w = keep & (np.abs(ph) < half)
        T.append(ph[w])
        F.append(L["fd"][w])
        E.append(np.full(w.sum(), L["sig"]))
    t, f, er = np.concatenate(T), np.concatenate(F), np.concatenate(E)
    o = np.argsort(t)
    t, f, er = t[o], f[o], er[o]
    edges = np.arange(-half, half + 1e-9, binmin / 1440)
    idx = np.digitize(t, edges)
    bt, bf, be = [], [], []
    for i in range(1, len(edges)):
        m = idx == i
        if m.sum() < 3:
            continue
        ww = 1 / er[m] ** 2
        bt.append(np.mean(t[m]))
        bf.append(np.sum(ww * f[m]) / ww.sum())
        be.append(1 / np.sqrt(ww.sum()))
    return np.array(bt), np.array(bf), np.array(be), len(t)


meta = dict(ephemeris=EPH, folds={}, notes="6-min bins; S98 = QLP SYS_RM_FLUX frozen mask (clean intervals)")
for sig in ("A", "B"):
    for name, sets in {"S32+S98": ["S32", "S98"], "S32only": ["S32"], "S98only": ["S98"], "S6only": ["S6"]}.items():
        binm = 30.0 if name == "S6only" else 6.0
        bt, bf, be, n = fold(sig, sets, binmin=binm)
        fn = os.path.join(OUT, f"{sig}_{name}_binned.csv")
        np.savetxt(fn, np.c_[bt, bf, be], delimiter=",", header="t_from_mid_d,flux,flux_err")
        intr = np.abs(bt) < 0.35 * EPH[sig]["T14"]
        meta["folds"][f"{sig}_{name}"] = dict(file=os.path.relpath(fn, AUD), n_raw=int(n), n_bins=int(len(bt)), bin_min=binm,
                                              mean_err_ppm=float(np.mean(be) * 1e6), in_transit_mean_depth_ppm=float((1 - np.mean(bf[intr])) * 1e6))
        print(sig, name, meta["folds"][f"{sig}_{name}"])

aps = {}
for sec, fn in ((6, "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits"),
                (32, "tess2020324010417-s0032-0000000004206066-0200-s_tp.fits")):
    h = fits.open(os.path.join(lcio.MAST, fn))
    ap = h[2].data
    c0, r0 = h[1].header["1CRV4P"], h[1].header["2CRV4P"]
    rows, cols = np.where((ap & 2) > 0)
    aps[str(sec)] = dict(pixels=[[int(c0 + c), int(r0 + r)] for r, c in zip(rows, cols)],
                         source=f"{fn} APERTURE ext, bit value 2 (optimal aperture); 1CRV4P={c0}, 2CRV4P={r0}")
h = fits.open(glob.glob(os.path.join(AUD, "data/tesscut/tess-s0098-*_astrocut.fits"))[0])
w = WCS(h[2].header)
c0, r0 = h[1].header["1CRV4P"], h[1].header["2CRV4P"]
x, y = w.all_world2pix(80.85982666, -1.32696401, 0)
cx, cy = int(np.round(c0 + x)), int(np.round(r0 + y))
aps["98"] = dict(pixels=[[cx + i, cy + j] for i in (-1, 0, 1) for j in (-1, 0, 1)],
                 source=f"QLP BESTAP '3x3 square' centred on the target pixel ({cx},{cy}) from TESScut WCS")
json.dump(aps, open(os.path.join(OUT, "apertures.json"), "w"), indent=1)
meta["apertures"] = {k: dict(n_pix=len(v["pixels"]), source=v["source"]) for k, v in aps.items()}
# compare with the earlier apertures file
old = json.load(open(os.path.join(AUD, "..", "deep_followup_tic4206066", "fpp", "inputs", "apertures.json")))
meta["apertures_identical_to_deep_followup"] = {k: sorted(map(tuple, aps[k]["pixels"])) == sorted(map(tuple, old[k]["pixels"])) for k in aps}
save_json(meta, os.path.join(OUT, "inputs_meta.json"))
print(json.dumps(meta["apertures"], indent=1), meta["apertures_identical_to_deep_followup"])
