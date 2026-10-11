"""V02: archival link test (Sectors 6 and 32), own implementation.
Statistic: max over a +-0.004 d period scan (centred on the S98 search ephemeris used by the manuscript) of the combined
1.6-h box S/N at the predicted archival times.  Box S/N per event = (1 - mean in-box flux of the detrended light curve) /
(sigma_pt / sqrt(n_in) * beta), beta calibrated so that the robust std of box S/N over the whole sector is 1.
Nulls (data kept as observed, i.e. B's archival dips remain in the data):
  EXCL : P0 ~ U(9, 13.5) excluding |P0 - 11.133| < 0.1, T0 uniform over one period (what both the manuscript and the
         auditor did, apart from the manuscript removing B's dips);
  NEAR : P0 ~ U(11.033, 11.233) (the excluded strip) -> fraction of the prior strip whose scan reproduces the observed
         statistic; the full-prior 'self-hit' FAP = f_near * 0.2/4.5 + f_excl * 4.3/4.5.
  REM  : EXCL with B's archival windows removed (manuscript's null).
"""
import sys, json, time
import numpy as np
from astropy.io import fits
import vcommon as v

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
PC, TREF = 11.132975671662091, 4000.8681717061445      # manuscript scan centre (S98 search values)
W = 0.004
NSTEP = 1601
BARCH = [1473.7328, 1484.8655, 2175.0957, 2186.2284, 2197.3612]


def prep(kind):
    d, t, q = v.load(kind)
    f = np.asarray(d["PDCSAP_FLUX"], float)
    ok = np.isfinite(t) & np.isfinite(f) & ((q & v.SPOC_BAD) == 0)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    A = v.inA(t)
    fd = v.biweight(t, f, A, 0.75)
    s = v.rstd(fd[~A])
    keep = np.isfinite(fd) & ~A & (fd < 1 + 4 * s)
    return t[keep], fd[keep]


class Grid:
    def __init__(self, t, f, step=1.0 / 1440):
        o = np.argsort(t); t, f = t[o], f[o]
        self.t0 = t.min() - 0.3; self.step = step
        g = np.arange(self.t0, t.max() + 0.3, step)
        cs = np.r_[0, np.cumsum(f)]
        i0 = np.searchsorted(t, g - v.DUR / 2); i1 = np.searchsorted(t, g + v.DUR / 2)
        n = i1 - i0
        cad = np.median(np.diff(t))
        nnom = v.DUR / cad
        sig = v.rstd(f)
        with np.errstate(invalid="ignore", divide="ignore"):
            dep = 1 - (cs[i1] - cs[i0]) / n
            e = sig / np.sqrt(n)
        ok = n >= 0.5 * nnom
        z = dep[ok] / e[ok]
        self.beta = max(1.0, v.rstd(z))
        self.dep = np.where(ok, dep, np.nan); self.err = np.where(ok, e * self.beta, np.nan)
        self.g = g; self.tmin, self.tmax = t.min(), t.max()

    def look(self, tc):
        i = np.rint((tc - self.t0) / self.step).astype(np.int64)
        inside = (i >= 0) & (i < len(self.g))
        i = np.clip(i, 0, len(self.g) - 1)
        return np.where(inside, self.dep[i], np.nan), np.where(inside, self.err[i], np.nan)


def remove_B(t, f):
    m = np.ones_like(t, bool)
    for tc in BARCH:
        m &= np.abs(t - tc) > 1.5 * v.DUR
    return t[m], f[m]


def stat(grids, P0, T0, w):
    P = P0 + np.linspace(-w, w, NSTEP)
    num = np.zeros(NSTEP); den = np.zeros(NSTEP)
    for g in grids:
        nlo = int(np.floor((g.tmin - T0) / P.min())) - 1; nhi = int(np.ceil((g.tmax - T0) / P.max())) + 1
        for n in range(min(nlo, nhi), max(nlo, nhi) + 1):
            tc = T0 + n * P
            if tc.max() < g.tmin - 0.1 or tc.min() > g.tmax + 0.1:
                continue
            d, e = g.look(tc)
            ok = np.isfinite(d)
            num[ok] += d[ok] / e[ok] ** 2; den[ok] += 1 / e[ok] ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(den > 0, num / np.sqrt(den), 0.0)
    j = int(np.argmax(s))
    return s[j], P[j]


if __name__ == "__main__":
    t6, f6 = prep("ts6"); t32, f32 = prep("s32")
    G = {"6": Grid(t6, f6), "32": Grid(t32, f32)}
    GR = {"6": Grid(*remove_B(t6, f6)), "32": Grid(*remove_B(t32, f32))}
    out = dict(beta={k: g.beta for k, g in G.items()}, runs=[])
    for secs in (["6", "32"], ["6"], ["32"]):
        grids = [G[s] for s in secs]; gridsR = [GR[s] for s in secs]
        obs, Pobs = stat(grids, PC, TREF, W)
        per = []
        for s in secs:
            for n in range(-300, 0):
                tc = TREF + n * Pobs
                if G[s].tmin - 0.05 < tc < G[s].tmax + 0.05:
                    d, e = G[s].look(np.array([tc]))
                    per.append((s, n, round(float(tc), 4), None if not np.isfinite(d[0]) else round(float(d[0] * 1e6)),
                                None if not np.isfinite(e[0]) else round(float(e[0] * 1e6))))
        res = dict(sectors=secs, obs=float(obs), P_obs=float(Pobs), per_event=per)
        t0 = time.time()
        for mode in ("EXCL", "NEAR", "REM"):
            vals = np.empty(NTR)
            for k in range(NTR):
                if mode == "NEAR":
                    P0 = rng.uniform(11.033, 11.233)
                else:
                    while True:
                        P0 = rng.uniform(9.0, 13.5)
                        if abs(P0 - 11.133) > 0.1:
                            break
                T0 = TREF + rng.uniform(-0.5, 0.5) * P0
                vals[k] = stat(gridsR if mode == "REM" else grids, P0, T0, W * P0 / PC)[0]
            nge = int(np.sum(vals >= obs))
            res[mode] = dict(n=NTR, n_ge=nge, frac=nge / NTR, q=[float(x) for x in np.percentile(vals, [50, 99, 99.9, 100])],
                             n_ge_90pct_obs=int(np.sum(vals >= 0.9 * obs)))
            print(secs, mode, "obs %.2f P %.6f" % (obs, Pobs), res[mode], "%.0fs" % (time.time() - t0), flush=True)
        fe, fn = res["EXCL"]["frac"], res["NEAR"]["frac"]
        res["FAP_full_prior_9_13.5"] = fn * 0.2 / 4.5 + fe * 4.3 / 4.5
        res["FAP_full_prior_upper95"] = ((res["NEAR"]["n_ge"] + 3) / NTR) * 0.2 / 4.5 + ((res["EXCL"]["n_ge"] + 3) / NTR) * 4.3 / 4.5
        out["runs"].append(res)
        print(secs, "per-event", per, "FAP(full prior)", res["FAP_full_prior_9_13.5"], flush=True)
        json.dump(out, open(v.HERE + "/out/v02_link.json", "w"), indent=1)
