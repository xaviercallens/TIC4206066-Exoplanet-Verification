"""Step 8: joint transit fits (batman + emcee), independent of the deep follow-up code.

usage: s08_globalfit.py SIGNAL SET PRIOR [S98MASK] [WINDOW] [TAG]
  SIGNAL A|B ; SET pipe|own|own2min ; PRIOR free|mann ; S98MASK frozen|default ; WINDOW biweight length (d)
Data sets
  pipe   : S6 TESS-SPOC PDCSAP (30 min), S32 SPOC PDCSAP (2 min), S98 QLP SYS_RM_FLUX (200 s)
  own    : TESScut aperture photometry of this audit, S6 r=2.5 px, S32 r=2.0 px (10 min), S98 r=2.0 px (200 s,
           straylight-flagged cadences removed); ownall: same but S98 straylight cadences kept
  own2min: as own but S32 from the SPOC 2-min TPF (own aperture r=2.0 px, no PDC)
Model: circular orbit, quadratic LD (Gaussian priors sigma 0.1 around Claret 2017 table15 LSM at the adopted
Teff/logg), per-data-set flux offset and jitter. Free-rho prior: uniform in ln(rho) on [0.1, 50] rho_sun
(posterior under a uniform-in-rho prior is obtained by importance weights w = rho and reported too).
The OTHER signal's transits are excluded from the fitted windows (+-1.0 T14) and both are masked in detrending.
Convergence: emcee autocorrelation time tau checked every 2000 steps; stop when N > 60 max(tau) and tau changes < 2 %.
"""
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np
import emcee

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import transit_tools as tt  # noqa
from audit_common import AUD, save_json  # noqa

SIGNAL = sys.argv[1]
SET = sys.argv[2]
PRIOR = sys.argv[3]
S98MASK = sys.argv[4] if len(sys.argv) > 4 else "frozen"
WIN = float(sys.argv[5]) if len(sys.argv) > 5 else 0.75
TAG = sys.argv[6] if len(sys.argv) > 6 else f"{SIGNAL}_{SET}_{PRIOR}_{S98MASK}_w{WIN}"

ST = json.load(open(os.path.join(AUD, "out", "s04_stellar.json")))
RHO_MANN = ST["density_rho_sun"]["median"]
E_RHO_MANN = 0.5 * (ST["density_rho_sun"]["plus1"] + ST["density_rho_sun"]["minus1"])
TEFF, LOGG = 4030.0, ST["logg_cgs"]["median"]
U1, U2 = tt.claret_ld(TEFF, LOGG)

EPH = {k: dict(v) for k, v in lcio.EPH.items()}
ephf = os.path.join(AUD, "out", "ephem_current.json")
if os.path.exists(ephf):
    EPH.update(json.load(open(ephf)))


def load_set():
    if SET == "pipe":
        lcs = [lcio.load_product("TESSSPOC_S6_PDC"), lcio.load_product("SPOC2min_S32_PDC"),
               lcio.load_product("QLP_S98_SYS", S98MASK)]
    elif SET == "own":
        lcs = [lcio.load_own("tesscut_s6", 2.5), lcio.load_own("tesscut_s32", 2.0),
               lcio.load_own("tesscut_s98", 2.0, drop_stray=True)]
    elif SET == "ownall":
        lcs = [lcio.load_own("tesscut_s6", 2.5), lcio.load_own("tesscut_s32", 2.0), lcio.load_own("tesscut_s98", 2.0)]
    elif SET == "own2min":
        lcs = [lcio.load_own("tesscut_s6", 2.5), lcio.load_own("spoc_tpf_s32", 2.0),
               lcio.load_own("tesscut_s98", 2.0, drop_stray=True)]
    else:
        raise ValueError(SET)
    return [lcio.detrend(lc, window=WIN, eph=EPH) for lc in lcs]


LCS = load_set()
W = tt.windows(LCS, SIGNAL, EPH, half=0.30)
MODELS = [tt.Model(w["t"], w["cad"]) for w in W]
NW = len(W)
E = EPH[SIGNAL]
NAMES = ["T0", "P", "lnk", "b", "lnrho", "u1", "u2"] + [f"c_{w['key']}" for w in W] + [f"lnjit_{w['key']}" for w in W]


def unpack(th):
    T0, P, lnk, b, lnrho, u1, u2 = th[:7]
    c = th[7:7 + NW]
    lj = th[7 + NW:7 + 2 * NW]
    return T0, P, np.exp(lnk), b, np.exp(lnrho), u1, u2, c, lj


def lnprob(th):
    T0, P, k, b, rho, u1, u2, c, lj = unpack(th)
    if not (abs(T0 - E["T0"]) < 0.1 and abs(P - E["P"]) < 0.01 and 0.005 < k < 0.2 and 0 <= b < 1 + k
            and 0.1 < rho < 50 and u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0):
        return -np.inf
    if np.any(np.abs(c - 1) > 0.01) or np.any(lj < np.log(1e-6)) or np.any(lj > np.log(1e-2)):
        return -np.inf
    lp = -0.5 * ((u1 - U1) / 0.1) ** 2 - 0.5 * ((u2 - U2) / 0.1) ** 2
    if PRIOR == "mann":
        lp += -0.5 * ((rho - RHO_MANN) / E_RHO_MANN) ** 2 + np.log(rho)   # sampling in ln(rho): p(ln rho) = N(rho) * rho
    ars = tt.ars_from_rho(rho, P)
    if b >= ars:
        return -np.inf
    for w, m, ci, lji in zip(W, MODELS, c, lj):
        mod = ci * m.lc(T0, P, k, b, ars, u1, u2)
        s2 = w["sig"] ** 2 + np.exp(2 * lji)
        lp += -0.5 * np.sum((w["f"] - mod) ** 2 / s2 + np.log(2 * np.pi * s2))
    return lp


def neg(th):
    v = lnprob(th)
    return 1e25 if not np.isfinite(v) else -v


def main():
    from scipy.optimize import minimize
    rng = np.random.default_rng(12345)
    k0 = 0.021 if SIGNAL == "A" else 0.028
    best = None
    for b0 in (0.1, 0.4, 0.7, 0.85):
        for r0 in (1.0, 2.5, 5.0):
            x0 = np.r_[E["T0"], E["P"], np.log(k0), b0, np.log(r0), U1, U2, np.ones(NW), np.full(NW, np.log(1e-4))]
            r = minimize(neg, x0, method="Nelder-Mead", options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-4))
            if best is None or r.fun < best.fun:
                best = r
    xb = best.x
    ndim = len(xb)
    nwalk = max(48, 2 * ndim + 2)
    p0 = xb + 1e-4 * rng.normal(size=(nwalk, ndim)) * np.maximum(np.abs(xb), 1e-2) * 0.01
    p0[:, 3] = np.clip(xb[3] + 0.1 * rng.normal(size=nwalk), 0.0, 0.95)
    p0[:, 4] = xb[4] + 0.1 * rng.normal(size=nwalk)
    p0[:, 0] = xb[0] + 5e-4 * rng.normal(size=nwalk)
    p0[:, 1] = xb[1] + 2e-6 * rng.normal(size=nwalk)
    for i in range(nwalk):
        while not np.isfinite(lnprob(p0[i])):
            p0[i] = xb + 1e-5 * rng.normal(size=ndim)
            p0[i, 3] = abs(p0[i, 3])
    t0 = time.time()
    # serial sampling: one likelihood call is ~0.2 ms, so process-pool overhead would dominate
    s = emcee.EnsembleSampler(nwalk, ndim, lnprob)
    old_tau = np.inf
    tau_hist = []
    for _ in s.sample(p0, iterations=80000, progress=False):
        if s.iteration % 2000:
            continue
        tau = s.get_autocorr_time(tol=0)
        tau_hist.append((s.iteration, float(np.max(tau))))
        conv = np.all(tau * 60 < s.iteration) and np.all(np.abs(old_tau - tau) / tau < 0.02)
        old_tau = tau
        if conv:
            break
    tau = s.get_autocorr_time(tol=0)
    burn = int(5 * np.max(tau))
    thin = max(1, int(0.5 * np.min(tau)))
    ch = s.get_chain(discard=burn, thin=thin, flat=True)
    lnp = s.get_log_prob(discard=burn, thin=thin, flat=True)
    np.savez_compressed(os.path.join(AUD, "out", "fits", f"chain_{TAG}.npz"), chain=ch, lnp=lnp, names=np.array(NAMES))
    T0, P, k, b, rho = ch[:, 0], ch[:, 1], np.exp(ch[:, 2]), ch[:, 3], np.exp(ch[:, 4])
    ars = tt.ars_from_rho(rho, P)
    inc = np.degrees(np.arccos(b / ars))
    T14 = tt.t14_hours(P, k, b, ars)
    Rs = np.random.default_rng(1).normal(ST["radius_mann2015_MK"]["median"], ST["radius_mann2015_MK"]["std"], len(k))
    Ms = np.random.default_rng(2).normal(ST["mass_mann2019_table6_n5"]["median"], ST["mass_mann2019_table6_n5"]["std"], len(k))
    Rp = k * Rs * 109.076  # R_sun / R_earth (IAU nominal: 695700 / 6378.1)
    a_au = (Ms * (P / 365.25) ** 2) ** (1 / 3)
    pct = lambda v, w=None: [float(x) for x in np.percentile(v, [2.275, 15.865, 50, 84.135, 97.725])]

    def wpct(v, w):
        o = np.argsort(v)
        cw = np.cumsum(w[o]) / np.sum(w)
        return [float(np.interp(q, cw, v[o])) for q in (0.02275, 0.15865, 0.5, 0.84135, 0.97725)]

    res = dict(tag=TAG, signal=SIGNAL, set=SET, prior=PRIOR, s98mask=S98MASK, window=WIN, ld_center=[U1, U2],
               rho_mann=[RHO_MANN, E_RHO_MANN], datasets=[dict(key=w["key"], sector=w["sector"], n=len(w["t"]),
                                                                 sig_ppm=w["sig"] * 1e6, cad_min=w["cad"] * 1440,
                                                                 epochs=sorted(set(w["epoch"].tolist()))) for w in W],
               n_walkers=nwalk, n_steps=int(s.iteration), tau=tau.tolist(), tau_max=float(np.max(tau)),
               n_over_tau=float(s.iteration / np.max(tau)), burn=burn, thin=thin, n_samples=len(ch),
               acceptance=float(np.mean(s.acceptance_fraction)), tau_history=tau_hist, runtime_s=time.time() - t0,
               map_lnprob=float(-best.fun), map=[float(x) for x in xb],
               T0=pct(T0), P=pct(P), k=pct(k), depth_k2_ppm=pct(k ** 2 * 1e6), b=pct(b), rho=pct(rho), a_over_Rs=pct(ars),
               inc_deg=pct(inc), T14_h=pct(T14), Rp_Rearth=pct(Rp), a_AU=pct(a_au), u1=pct(ch[:, 5]), u2=pct(ch[:, 6]),
               rho_uniform_prior_reweighted=wpct(rho, rho) if PRIOR == "free" else None,
               b_uniform_rho_reweighted=wpct(b, rho) if PRIOR == "free" else None,
               frac_samples_b_gt_0p9=float(np.mean(b > 0.9)), frac_rho_at_upper_edge=float(np.mean(rho > 40)),
               frac_rho_at_lower_edge=float(np.mean(rho < 0.15)),
               corr_T0_P=float(np.corrcoef(T0, P)[0, 1]))
    save_json(res, os.path.join(AUD, "out", "fits", f"fit_{TAG}.json"))
    print(json.dumps({k_: res[k_] for k_ in ["tag", "n_steps", "tau_max", "n_over_tau", "acceptance", "T0", "P", "k",
                                              "depth_k2_ppm", "b", "rho", "T14_h", "Rp_Rearth",
                                              "rho_uniform_prior_reweighted"]}, indent=0))


if __name__ == "__main__":
    os.makedirs(os.path.join(AUD, "out", "fits"), exist_ok=True)
    main()
