"""Step 7: independent recovery of both signals with astropy BoxLeastSquares (a different BLS implementation from the
frozen pipeline's), on (i) S98 alone with three quality-mask choices and two independent photometries, (ii) each
archival sector alone, (iii) all sectors combined. Iterative: after the strongest peak, mask it and search again.
SNR here = BLS depth / depth_err (white) ; SDE = (peak - mean)/std of the power spectrum (log-likelihood objective).
Harmonic/alias bookkeeping: every top peak is compared with P_A, P_B and k*P, P/k (k = 2, 3, 4) within 0.5 %.
"""
import json
import os
import sys

import numpy as np
from astropy.timeseries import BoxLeastSquares

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
from audit_common import save_json  # noqa

PA, PB = 3.1827863, 11.132743
DUR = np.array([0.04, 0.055, 0.07, 0.085, 0.10, 0.13])


def relation(p):
    rels = []
    for name, P0 in (("A", PA), ("B", PB)):
        for lab, fac in (("P", 1), ("2P", 2), ("3P", 3), ("4P", 4), ("P/2", 1 / 2), ("P/3", 1 / 3), ("P/4", 1 / 4),
                         ("2P/3", 2 / 3), ("3P/2", 3 / 2)):
            if abs(p / (P0 * fac) - 1) < 0.005:
                rels.append(f"{name}:{lab}")
    return rels


def bin_lc(t, f, width):
    if width is None:
        return t, f
    edges = np.arange(t.min(), t.max() + width, width)
    idx = np.digitize(t, edges)
    bt, bf = [], []
    for i in np.unique(idx):
        m = idx == i
        if m.sum() >= 1:
            bt.append(t[m].mean())
            bf.append(f[m].mean())
    return np.array(bt), np.array(bf)


def search(t, f, pmin=0.5, pmax=None, n_peaks=3, nit=3, os_=5, mask_found=True):
    base = t.max() - t.min()
    pmax = pmax or min(30.0, base / 2)
    # log-uniform period grid: d ln P = dur_min / (baseline * oversampling)  (phase drift < dur_min/os over baseline)
    dlnp = DUR.min() / (base * os_)
    periods = np.exp(np.arange(np.log(pmin), np.log(pmax), dlnp))
    res_all = []
    tt, ff = t.copy(), f.copy()
    for it in range(nit):
        err = np.full(len(tt), 1.4826 * np.median(np.abs(ff - np.median(ff))))
        bls = BoxLeastSquares(tt, ff, err)
        r = bls.power(periods, DUR, objective="likelihood")
        pw = r.power
        sde = (pw - np.mean(pw)) / np.std(pw)
        i = int(np.argmax(pw))
        st = bls.compute_stats(r.period[i], r.duration[i], r.transit_time[i])
        snr = float(r.depth[i] / r.depth_err[i])
        # top distinct peaks
        order = np.argsort(pw)[::-1]
        tops = []
        for j in order:
            if all(abs(r.period[j] / x["P"] - 1) > 0.01 for x in tops):
                tops.append(dict(P=float(r.period[j]), SDE=float(sde[j]), depth_ppm=float(r.depth[j] * 1e6),
                                 snr=float(r.depth[j] / r.depth_err[j]), dur_h=float(r.duration[j] * 24),
                                 T0=float(r.transit_time[j]), rel=relation(r.period[j])))
            if len(tops) >= n_peaks + 3:
                break
        found = dict(iteration=it, P=float(r.period[i]), T0=float(r.transit_time[i]), dur_h=float(r.duration[i] * 24),
                     depth_ppm=float(r.depth[i] * 1e6), depth_err_ppm=float(r.depth_err[i] * 1e6), snr=snr, SDE=float(sde[i]),
                     rel=relation(r.period[i]), n_transits_with_data=int(np.sum(st["per_transit_count"] > 0)),
                     depth_odd_ppm=float(st["depth_odd"][0] * 1e6), depth_even_ppm=float(st["depth_even"][0] * 1e6),
                     depth_half_ppm=float(st["depth_half"][0] * 1e6), depth_phased_ppm=float(st["depth_phased"][0] * 1e6),
                     top_peaks=tops, n_periods=len(periods))
        res_all.append(found)
        if not mask_found:
            break
        ph = ((tt - r.transit_time[i]) / r.period[i] + 0.5) % 1 - 0.5
        keep = np.abs(ph * r.period[i]) > 1.5 * r.duration[i]
        tt, ff = tt[keep], ff[keep]
    return res_all


def prep(lc, window=0.75):
    d = lcio.detrend(lc, window=window, mask_fac=0.0)   # NO transit mask: blind search must not know the signals
    return d["t"], d["fd"]


out = {}
cases = {
    "S98_QLP_SYS_frozenmask": lambda: lcio.load_product("QLP_S98_SYS", "frozen"),
    "S98_QLP_SYS_defaultmask": lambda: lcio.load_product("QLP_S98_SYS", "default"),
    "S98_QLP_SYS_defaultmask_nostray": lambda: lcio.load_product("QLP_S98_SYS", "default+stray"),
    "S98_TESScut_own_r2": lambda: lcio.load_own("tesscut_s98", 2.0),
    "S98_TESScut_own_r2_nostray": lambda: lcio.load_own("tesscut_s98", 2.0, drop_stray=True),
    "S32_SPOC2min_PDC": lambda: lcio.load_product("SPOC2min_S32_PDC"),
    "S32_TESScut_own_r2": lambda: lcio.load_own("tesscut_s32", 2.0),
    "S6_TESSSPOC_PDC": lambda: lcio.load_product("TESSSPOC_S6_PDC"),
    "S6_TESScut_own_r2.5": lambda: lcio.load_own("tesscut_s6", 2.5),
}
for name, fn in cases.items():
    lc = fn()
    t, f = prep(lc)
    if lc["cad"] < 0.0015:     # 2-min data: bin to 6 min for speed (transit 1.6-1.9 h)
        t, f = bin_lc(t, f, 6 / 1440)
    res = search(t, f)
    out[name] = res
    print(f"== {name}: n={len(t)}")
    for r in res:
        print(f"   it{r['iteration']}: P={r['P']:.5f} dur={r['dur_h']:.2f}h depth={r['depth_ppm']:.0f}±{r['depth_err_ppm']:.0f} "
              f"SNR={r['snr']:.1f} SDE={r['SDE']:.1f} rel={r['rel']} ntr={r['n_transits_with_data']} half={r['depth_half_ppm']:.0f}")
        print("        tops:", [(round(x['P'], 4), round(x['SDE'], 1), x['rel']) for x in r["top_peaks"][:5]])
    save_json(out, "out/s07_search.json")

# combined blind search: all three sectors (pipeline set and own set), 30-min bins
for name, keys in {"ALL_pipeline": ["TESSSPOC_S6_PDC", "SPOC2min_S32_PDC", ("QLP_S98_SYS", "default")],
                   "ALL_own": [("own", "tesscut_s6", 2.5), ("own", "tesscut_s32", 2.0), ("own", "tesscut_s98", 2.0)]}.items():
    T, F = [], []
    for k in keys:
        if isinstance(k, tuple) and k[0] == "own":
            lc = lcio.load_own(k[1], k[2])
        elif isinstance(k, tuple):
            lc = lcio.load_product(k[0], k[1])
        else:
            lc = lcio.load_product(k)
        t, f = prep(lc)
        t, f = bin_lc(t, f, 30 / 1440)
        T.append(t)
        F.append(f / np.median(f))
    t, f = np.concatenate(T), np.concatenate(F)
    res = search(t, f, pmin=1.0, pmax=15.0, os_=3)
    out[name] = res
    print(f"== {name}: n={len(t)}")
    for r in res:
        print(f"   it{r['iteration']}: P={r['P']:.6f} dur={r['dur_h']:.2f}h depth={r['depth_ppm']:.0f}±{r['depth_err_ppm']:.0f} "
              f"SNR={r['snr']:.1f} SDE={r['SDE']:.1f} rel={r['rel']} ntr={r['n_transits_with_data']}")
        print("        tops:", [(round(x['P'], 5), round(x['SDE'], 1), x['rel']) for x in r["top_peaks"][:5]])
    save_json(out, "out/s07_search.json")
