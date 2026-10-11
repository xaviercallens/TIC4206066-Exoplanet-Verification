"""Time-averaging beta factor (Pont et al. 2006; Winn et al. 2008) of detrended out-of-transit residuals per sector
(BASE products), for bin sizes 0.5-2 x T14_B. Output out/s06_beta.json"""
from b_run import *
from s_configs import CONFIGS
cfg = CONFIGS["BASE"]; res = {}
for sec in (6, 32, 98):
    lc = get_lc(cfg["products"][sec], spoc_mask=17087); d = detrend(lc)
    oot = ~(in_transit(d["t"], EPH_B, EPH_B["T14"]) | in_transit(d["t"], EPH_A, EPH_A["T14"]))
    t, r = d["t"][oot], d["fd"][oot] - 1
    s1 = np.std(r); out = {}
    for fac in (0.5, 0.75, 1.0, 1.5, 2.0):
        L = fac * EPH_B["T14"]; edges = np.arange(t.min(), t.max() + L, L); idx = np.digitize(t, edges)
        nper = np.bincount(idx); full = [i for i in np.unique(idx) if nper[i] >= 0.8 * L / d["cad"]]
        means = np.array([r[idx == i].mean() for i in full]); N = np.median([nper[i] for i in full]); M = len(full)
        sN = np.std(means); exp = s1 / np.sqrt(N) * np.sqrt(M / (M - 1))
        out[str(fac)] = dict(beta=float(sN / exp), N=float(N), M=int(M))
    res[str(sec)] = dict(product=lc["name"], sigma1_ppm=float(s1 * 1e6), bins=out, beta_median=float(np.median([v["beta"] for v in out.values()])))
    print(sec, lc["name"], {k: round(v["beta"], 2) for k, v in out.items()}, flush=True)
jdump(res, "s06_beta.json")
