"""T6_star rerun 1: re-read the ARIADNE BMA pickle of companion_check_tic4206066 (read-only) and settle which
Teff / R / L values and 1-sigma errors the run actually supports.

Checks:
 (a) percentiles of ARIADNE's 'weighted_average' arrays (sum_i w_i * independent draw from grid i)
 (b) percentiles of ARIADNE's 'weighted_samples' arrays (weighted KDE resample of the concatenated traces)
 (c) a proper mixture: N*w_i draws from each grid posterior ('originals')
 (d) ARIADNE 1.5.0's own out_filler HDR (utils.estimate_pdf 300-pt KDE grid + estimate_cdf 499-bin histogram),
     copied verbatim below, applied to (a) and (b) -> must reproduce best_fit_average.dat / best_fit_sample.dat
 (e) the same HDR with a consistent grid (pdf and cdf on the same 300-pt KDE grid) applied to (c)
Source of copied functions: astroARIADNE 1.5.0 utils.py (sha256 f2468c57...,
/private/tmp/.../dd0bfc95.../ariadne_env/lib/python3.11/site-packages/astroARIADNE/utils.py L25-55, L93-139).
Output: out/t6_01_ariadne_bma.json
"""
import json, pickle
from pathlib import Path
import numpy as np
from scipy.special import erf
from scipy.stats import gaussian_kde

ROOT = Path("/Users/paulsmith/space <3/tess-jev")
AR = ROOT / "addenda/companion_check_tic4206066/sed_out"
OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(exist_ok=True)


class _Stub:
    def __init__(self, *a, **k): pass
    def __setstate__(self, st): self.__dict__["_state"] = st


class _U(pickle.Unpickler):
    def find_class(self, mod, name):
        try:
            return super().find_class(mod, name)
        except Exception:
            return type(name, (_Stub,), {})


# ---- verbatim from astroARIADNE 1.5.0 utils.py ----
def estimate_pdf(distribution):
    kde = gaussian_kde(distribution)
    xmin, xmax = distribution.min(), distribution.max()
    xx = np.linspace(xmin, xmax, 300)
    pdf = kde(xx)
    return xx, pdf


def estimate_cdf(distribution, hdr=False):
    h, hx = np.histogram(distribution, density=True, bins=499)
    cdf = np.zeros(500)
    if hdr:
        idx = np.argsort(h)[::-1]
        cdf[1:] = np.cumsum(h[idx]) * np.diff(hx)
    else:
        cdf[1:] = np.cumsum(h) * np.diff(hx)
    return cdf


def credibility_interval_hdr(xx, pdf, cdf, sigma=1.):
    best = xx[np.argmax(pdf)]
    z = erf(sigma / np.sqrt(2))
    idx = np.argsort(pdf)[::-1]
    idx_hdr = np.where(cdf >= z)[0][0]
    hdr = pdf[idx][:idx_hdr]
    hdr_min = hdr.min()
    low = xx[pdf > hdr_min].min()
    high = xx[pdf > hdr_min].max()
    return best, low, high
# ---- end verbatim ----


def ariadne_hdr(s):
    s = np.asarray(s, float)
    xx, pdf = estimate_pdf(s)
    cdf = estimate_cdf(s, hdr=True)
    best, lo, up = credibility_interval_hdr(xx, pdf, cdf, 1)
    return dict(best=float(best), plus=float(up - best), minus=float(best - lo))


def consistent_hdr(s, level=0.682689492):
    s = np.asarray(s, float)
    xx, pdf = estimate_pdf(s)
    p = pdf / pdf.sum()
    order = np.argsort(pdf)[::-1]
    k = np.searchsorted(np.cumsum(p[order]), level) + 1
    sel = xx[order[:k]]
    best = xx[np.argmax(pdf)]
    return dict(best=float(best), plus=float(sel.max() - best), minus=float(best - sel.min()))


def pct(s):
    s = np.asarray(s, float)
    lo, med, hi = np.percentile(s, [15.8655, 50, 84.1345])
    return dict(median=float(med), plus=float(hi - med), minus=float(med - lo), std=float(s.std()), n=int(s.size))


bma = _U(open(AR / "BMA.pkl", "rb")).load()
w = {k: float(v) for k, v in bma["weights"].items()}
rng = np.random.default_rng(20260925)
res = {"input": str(AR / "BMA.pkl"), "weights": w, "sqrt_sum_w2": float(np.sqrt(sum(v * v for v in w.values()))), "par": {}}
for par in ("teff", "rad", "lum", "z", "dist"):
    if par not in bma["weighted_average"]:
        continue
    wa = np.asarray(bma["weighted_average"][par], float)
    ws = np.asarray(bma["weighted_samples"][par], float)
    pool = [rng.choice(np.asarray(o[par], float), int(round(w[g] * 400000))) for g, o in bma["originals"].items() if par in o]
    mix = np.concatenate(pool)
    d = {"weighted_average_pct": pct(wa), "weighted_samples_pct": pct(ws), "mixture_pct": pct(mix),
         "per_grid_median": {g: float(np.median(o[par])) for g, o in bma["originals"].items() if par in o},
         "per_grid_std": {g: float(np.std(o[par])) for g, o in bma["originals"].items() if par in o}}
    # full arrays as in the ARIADNE run
    sub = lambda a: a  # full arrays, as ARIADNE did
    d["ariadne_hdr_on_weighted_average"] = ariadne_hdr(sub(wa))
    d["ariadne_hdr_on_weighted_samples"] = ariadne_hdr(sub(ws))
    d["ariadne_hdr_on_mixture"] = ariadne_hdr(sub(mix))
    d["consistent_hdr_on_weighted_average"] = consistent_hdr(sub(wa))
    d["consistent_hdr_on_mixture"] = consistent_hdr(sub(mix))
    if par in bma.get("best_fit_averaged", {}):
        d["pickle_best_fit_averaged"] = [float(bma["best_fit_averaged"][par])] + [float(x) for x in bma["uncertainties_averaged"][par]]
        d["pickle_best_fit_samples"] = [float(bma["best_fit_samples"][par])] + [float(x) for x in bma["uncertainties_samples"][par]]
    res["par"][par] = d
json.dump(res, open(OUT / "t6_01_ariadne_bma.json", "w"), indent=1)
for p, d in res["par"].items():
    print(p, "WA", {k: round(v, 4) for k, v in d["weighted_average_pct"].items()},
          "\n   MIX", {k: round(v, 4) for k, v in d["mixture_pct"].items()},
          "\n   WS", {k: round(v, 4) for k, v in d["weighted_samples_pct"].items()},
          "\n   ARIADNE-HDR(WA)", d["ariadne_hdr_on_weighted_average"], "ARIADNE-HDR(WS)", d["ariadne_hdr_on_weighted_samples"],
          "\n   ARIADNE-HDR(MIX)", d["ariadne_hdr_on_mixture"], "consistentHDR(WA)", d["consistent_hdr_on_weighted_average"],
          "consistentHDR(MIX)", d["consistent_hdr_on_mixture"],
          "\n   pickle best_fit_averaged", d.get("pickle_best_fit_averaged"), "samples", d.get("pickle_best_fit_samples"))
