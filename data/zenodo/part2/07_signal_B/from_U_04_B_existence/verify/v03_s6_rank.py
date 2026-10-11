"""V03: Sector 6 -- are the two predicted B dips the two strongest non-overlapping 1.6-h boxes outside A?
Own box estimator (in-box mean vs 1 after detrending, beta-calibrated), 4 detrendings x 4 photometries."""
import json, numpy as np
from math import comb
import vcommon as v
from v02_link import Grid

TB = [v.T0_MS + n * v.P_MS for n in (-227, -226)]
DUMPS6 = []


def get(phot):
    if phot in ("pdc", "sap"):
        d, t, q = v.load("ts6")
        f = np.asarray(d["PDCSAP_FLUX" if phot == "pdc" else "SAP_FLUX"], float)
        ok = np.isfinite(t) & np.isfinite(f) & ((q & v.SPOC_BAD) == 0)
        dumps = t[(q & 32) > 0]
    elif phot == "qlp":
        d, t, q = v.load("qlp6")
        f = np.asarray(d["KSPSAP_FLUX"], float)
        ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
        dumps = np.array([])
    elif phot == "tglc_aper":
        from astropy.io import fits
        d = fits.open(v.os.path.join(v.MAST, v.F["tglc6"]))[1].data
        t = np.asarray(d["time"], float); f = np.asarray(d["cal_aper_flux"], float)
        ok = np.isfinite(t) & np.isfinite(f) & (np.asarray(d["TESS_flags"]) == 0) & (np.asarray(d["TGLC_flags"]) == 0)
        dumps = np.array([])
    return t[ok], f[ok] / np.nanmedian(f[ok]), dumps


def detr(t, f, meth, w):
    from wotan import flatten
    A = v.inA(t)
    _, tr = flatten(t[~A], f[~A], method=meth, window_length=w, return_trend=True, break_tolerance=0.4)
    ok = np.isfinite(tr)
    fd = f / np.interp(t, t[~A][ok], tr[ok])
    s = v.rstd(fd[~A])
    keep = np.isfinite(fd) & ~A & (fd < 1 + 4 * s)
    return t[keep], fd[keep]


res = {}
_, _, dumps = get("pdc")
for phot in ("pdc", "sap", "qlp", "tglc_aper"):
    t, f, _ = get(phot)
    for meth, w in (("biweight", 0.5), ("biweight", 0.75), ("biweight", 1.0), ("median", 0.75)):
        tt, ff = detr(t, f, meth, w)
        g = Grid(tt, ff)
        z = g.dep / g.err
        ok = np.isfinite(z)
        x, s = g.g[ok], z[ok]
        order = np.argsort(-s)
        boxes = []
        blocked = np.zeros(len(x), bool)
        for i in order:
            if blocked[i]:
                continue
            boxes.append((float(x[i]), float(s[i])))
            blocked |= np.abs(x - x[i]) < v.DUR
            if len(boxes) > 400:
                break
        ranks = []
        for tb in TB:
            c = [(k + 1, b) for k, b in enumerate(boxes) if abs(b[0] - tb) < v.DUR / 2]
            ranks.append((c[0][0], round(c[0][1][1], 2), round(c[0][1][0], 3)) if c else None)
        others = [(round(b[0], 3), round(b[1], 2), round(float(np.min(np.abs(dumps - b[0])) * 24), 1) if len(dumps) else None)
                  for b in boxes if all(abs(b[0] - tb) > v.DUR for tb in TB)][:4]
        N = int(round(np.sum(ok) * (g.g[1] - g.g[0]) / v.DUR))
        key = "%s|%s%.2f" % (phot, meth, w)
        res[key] = dict(B=ranks, top_others=others, N_indep=N, top2_are_B=bool(ranks[0] and ranks[1] and sorted([ranks[0][0], ranks[1][0]]) == [1, 2]))
        print(key, res[key], flush=True)
N = res["pdc|biweight0.75"]["N_indep"]
K = 2 * 0.004 * 226.5 / v.DUR
res["p_rank"] = dict(N=N, K=K, p_3spacings=K * 3 / comb(N, 2), p_1spacing=K / comb(N, 2))
print(res["p_rank"])
res["dumps_s6"] = [float(x) for x in dumps]
json.dump(res, open(v.HERE + "/out/v03_s6_rank.json", "w"), indent=1)
