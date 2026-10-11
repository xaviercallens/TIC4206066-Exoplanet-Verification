"""Joint batman + emcee fit of all Signal-A transits (S6 30-min, S32 2-min, S98 200-s), circular orbit.

usage: s03_jointfit.py RHO LD NOISE [DATAVARIANT] [NSTEPS]
  RHO   : mann (split-normal 2.51 +0.25/-0.22 rho_sun) | uni (uniform in rho, 0.05-30) | logu (uniform in ln rho, 0.05-30)
  LD    : claret (Gaussian sigma 0.1 on u1,u2 around Claret 2017 PHOENIX-COND q, LSM) | atlas (same, ATLAS Z=0 xi=2) |
          free (Kipping 2013 q1,q2 ~ U(0,1))
  NOISE : jit (free white-noise sigma per sector) | beta (sigma fixed to beta_s * sigma_s from the baseline MAP; red noise)
Parameters: T0, P, k, b, ln rho, LD(2), per-sector offset (3), per-sector ln sigma (3; not for 'beta').
Model: batman quadratic LD, supersampled per cadence (S6 x15 over 1800 s, S32 x2 over 120 s, S98 x4 over 200 s).
Data: +-0.30 d around each A transit; points within +-1.0 T14_B (1.65 h) of a B transit removed.
"""
import sys, json, time
import numpy as np
import emcee
from scipy.optimize import minimize
from common import *

RHO, LD, NOISE = sys.argv[1], sys.argv[2], sys.argv[3]
DV = sys.argv[4] if len(sys.argv) > 4 else "base"
NSTEPS = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
EXCL = [int(x) for x in os.environ.get("EXCL_S98_EPOCHS", "").split(",") if x.strip()]   # optional S98 epochs to drop
TAG = f"{RHO}_{LD}_{NOISE}_{DV}" + ("_excl" + "_".join(map(str, EXCL)) if EXCL else "")
LDC = claret_ld()
U0 = np.array(LDC["atlas_Z0_xi2_LSM"] if LD == "atlas" else LDC["phoenix_q_LSM"])
SIG_U = 0.1
SECT = (6, 32, 98)

# ---------------- data
D = []
for s in SECT:
    z = np.load(os.path.join(OUT, f"lc_{DV}_S{s}.npz"))
    t, f = z["t"], z["f"]
    dA, _ = phase_dt(t, P_A, T0_A)
    dB, _ = phase_dt(t, P_B, T0_B)
    m = (np.abs(dA) < 0.30) & (np.abs(dB) > T14_B_H / 24.0)
    if s == 98 and EXCL:
        m &= ~np.isin(np.round((t - T0_A) / P_A).astype(int), EXCL)
    dd = dict(s=s, t=t[m], f=f[m], mod=BatmanSet(t[m], float(z["cad"]), SUPERSAMPLE[s]))
    if DV == "raw":   # per-window linear baseline, profiled analytically (weighted equally within a window)
        ep = np.round((t[m] - T0_A) / P_A).astype(int)
        _, wid = np.unique(ep, return_inverse=True)
        x = t[m] - (T0_A + ep * P_A)
        nw_ = wid.max() + 1
        S0 = np.bincount(wid, minlength=nw_); S1 = np.bincount(wid, x, nw_); S2 = np.bincount(wid, x * x, nw_)
        dd.update(wid=wid, x=x, S0=S0, S1=S1, S2=S2, det=S0 * S2 - S1 ** 2)
    D.append(dd)
NPTS = sum(len(d["t"]) for d in D)

beta_fix = None
if NOISE == "beta":
    bj = json.load(open(os.path.join(OUT, f"beta_{DV}.json")))
    beta_fix = np.array([bj[str(s)]["sigma_white"] * bj[str(s)]["beta"] for s in SECT])

NAMES = ["T0", "P", "k", "b", "lnrho", "l1", "l2"] + [f"c{s}" for s in SECT] + \
        ([f"lns{s}" for s in SECT] if NOISE == "jit" else [])


def ld_from(l1, l2):
    if LD == "free":
        q1, q2 = l1, l2
        if not (0 < q1 < 1 and 0 < q2 < 1):
            return None
        sq = np.sqrt(q1)
        return 2 * sq * q2, sq * (1 - 2 * q2)
    u1, u2 = l1, l2
    if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0):
        return None
    return u1, u2


def log_prior(th):
    T0, P, k, b, lnrho, l1, l2 = th[:7]
    rho = np.exp(lnrho)
    if not (0.004 < k < 0.1 and 0 <= b < 1 + k and 0.05 < rho < 30 and abs(T0 - T0_A) < 0.05 and abs(P - P_A) < 0.001):
        return -np.inf
    u = ld_from(l1, l2)
    if u is None:
        return -np.inf
    lp = 0.0
    if RHO == "mann":
        sg = RHO_EP if rho > RHO_MED else RHO_EM
        lp += -0.5 * ((rho - RHO_MED) / sg) ** 2 + lnrho          # density prior in rho, sampled in ln rho
    elif RHO == "uni":
        lp += lnrho                                               # p(ln rho) = p(rho) * rho, p(rho)=const
    elif RHO == "logu":
        lp += 0.0
    if DV == "raw":
        lp += -0.5 * np.sum((np.array(th[7:10]) / 1e-4) ** 2)   # offsets are absorbed by the profiled baselines
    if LD in ("claret", "atlas"):
        lp += -0.5 * np.sum(((np.array(u) - U0) / SIG_U) ** 2)
    if NOISE == "jit":
        ls = th[10:13]
        if np.any(ls < np.log(1e-5)) or np.any(ls > np.log(1e-2)):
            return -np.inf
    return lp


def models(th):
    T0, P, k, b, lnrho, l1, l2 = th[:7]
    u1, u2 = ld_from(l1, l2)
    ars = ars_from_rho(np.exp(lnrho), P)
    if b >= ars:
        return None
    return [d["mod"].flux(T0, P, k, b, ars, u1, u2) + th[7 + i] for i, d in enumerate(D)]


def log_like(th, return_parts=False):
    ms = models(th)
    if ms is None:
        return -np.inf
    ll = 0.0
    for i, d in enumerate(D):
        sig = np.exp(th[10 + i]) if NOISE == "jit" else beta_fix[i]
        r = d["f"] - ms[i]
        if DV == "raw":
            nw_ = len(d["S0"])
            Sy = np.bincount(d["wid"], r, nw_); Sxy = np.bincount(d["wid"], r * d["x"], nw_)
            a0 = (d["S2"] * Sy - d["S1"] * Sxy) / d["det"]; a1 = (d["S0"] * Sxy - d["S1"] * Sy) / d["det"]
            r = r - a0[d["wid"]] - a1[d["wid"]] * d["x"]
        ll += -0.5 * np.sum(r ** 2) / sig ** 2 - len(r) * np.log(sig)
    return ll


def log_prob(th):
    lp = log_prior(th)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_like(th)
    return lp + ll if np.isfinite(ll) else -np.inf


def start():
    l12 = [0.36, 0.25] if LD == "free" else list(U0)   # q1=(u1+u2)^2, q2=u1/(2(u1+u2))
    if LD == "free":
        uu = U0.sum()
        l12 = [uu ** 2, U0[0] / (2 * uu)]
    th = [T0_A, P_A, 0.0212, 0.35, np.log(2.5)] + l12 + [0.0, 0.0, 0.0]
    if NOISE == "jit":
        th += [np.log(2.0e-4), np.log(8.5e-4), np.log(6.2e-4)]
    return np.array(th)


def main():
    t0 = time.time()
    x0 = start()
    tt = time.time(); [log_prob(x0) for _ in range(20)]; per_call = (time.time() - tt) / 20
    print(TAG, "npts", NPTS, "per-call ms", per_call * 1e3, flush=True)
    neg = lambda x: -log_prob(x) if np.isfinite(log_prob(x)) else 1e25
    best = x0
    for it in range(3):
        r = minimize(neg, best, method="Nelder-Mead", options=dict(maxiter=20000, maxfev=20000, xatol=1e-9, fatol=1e-3))
        best = r.x
    print("MAP", dict(zip(NAMES, np.round(best, 6))), -r.fun, flush=True)
    ndim = len(best)
    nw = max(40, 3 * ndim)
    rng = np.random.default_rng(42)
    scale = np.array([5e-4, 2e-6, 5e-4, 0.05, 0.03, 0.03, 0.03] + [2e-5] * 3 + ([0.02] * 3 if NOISE == "jit" else []))
    p0 = []
    while len(p0) < nw:
        c = best + scale * rng.standard_normal(ndim)
        if c[3] < 0:
            c[3] = -c[3]
        if np.isfinite(log_prob(c)):
            p0.append(c)
    sampler = emcee.EnsembleSampler(nw, ndim, log_prob,
                                    moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    taus = []
    for i, _ in enumerate(sampler.sample(np.array(p0), iterations=NSTEPS, progress=False)):
        if (i + 1) % 2000 == 0:
            tau = sampler.get_autocorr_time(tol=0)
            taus.append([i + 1, float(np.max(tau))])
            print(f"step {i+1} tau_max {np.max(tau):.1f} acc {np.mean(sampler.acceptance_fraction):.3f} "
                  f"elapsed {time.time()-t0:.0f}s", flush=True)
            if (i + 1) > 60 * np.max(tau) and (i + 1) >= 6000:
                break
    tau = sampler.get_autocorr_time(tol=0)
    burn = int(3 * np.max(tau))
    thin = max(1, int(np.max(tau) / 4))
    ch = sampler.get_chain(discard=burn, thin=thin, flat=True)
    lpch = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
    np.savez_compressed(os.path.join(OUT, f"chain_{TAG}.npz"), chain=ch, logp=lpch, names=NAMES)
    # derived
    T0, P, k, b, rho = ch[:, 0], ch[:, 1], ch[:, 2], ch[:, 3], np.exp(ch[:, 4])
    u = np.array([ld_from(a, c) for a, c in zip(ch[:, 5], ch[:, 6])])
    ars = ars_from_rho(rho, P)
    t14 = t14_hours(P, k, b, ars)
    rs = rng.normal(R_STAR, E_R_STAR, len(k))
    rp = k * rs * RSUN_OVER_REARTH
    rp_fixedR = k * R_STAR * RSUN_OVER_REARTH
    # mid-transit model depth (unbinned) for a random subset
    import batman
    pm = batman.TransitParams(); pm.t0 = 0; pm.per = P_A; pm.ecc = 0; pm.w = 90; pm.limb_dark = "quadratic"
    idx = rng.choice(len(k), size=min(3000, len(k)), replace=False)
    dmid = []
    for j in idx:
        pm.rp, pm.a = k[j], ars[j]; pm.inc = np.degrees(np.arccos(b[j] / ars[j])); pm.u = list(u[j])
        dmid.append(1 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0])
    dmid = np.array(dmid)
    res = dict(tag=TAG, rho_prior=RHO, ld=LD, noise=NOISE, data_variant=DV, npts=NPTS,
               n_per_sector={str(d["s"]): int(len(d["t"])) for d in D}, ld_center=U0.tolist(), sig_u=SIG_U,
               supersample=SUPERSAMPLE, n_walkers=nw, n_steps=int(sampler.iteration), tau=tau.tolist(),
               tau_max=float(np.max(tau)), n_over_tau=float(sampler.iteration / np.max(tau)), burn=burn, thin=thin,
               n_samples=int(len(ch)), acceptance=float(np.mean(sampler.acceptance_fraction)), tau_history=taus,
               map=dict(zip(NAMES, best.tolist())), map_logp=float(-r.fun), runtime_s=time.time() - t0,
               quantiles="2.275,15.865,50,84.135,97.725 percent",
               T0=pct(T0), P=pct(P), k=pct(k), k_mean_std=[float(k.mean()), float(k.std())], depth_k2_ppm=pct(k ** 2 * 1e6),
               depth_mid_model_ppm=pct(dmid * 1e6), b=pct(b), rho=pct(rho), ars=pct(ars), T14_h=pct(t14),
               u1=pct(u[:, 0]), u2=pct(u[:, 1]), Rp_Re=pct(rp), Rp_Re_mean_std=[float(rp.mean()), float(rp.std())],
               Rp_Re_fixedRstar=pct(rp_fixedR),
               rho_upper=dict(q95=float(np.percentile(rho, 95)), q97_7=float(np.percentile(rho, 97.725)),
                              q99=float(np.percentile(rho, 99)), q99_87=float(np.percentile(rho, 99.865))),
               frac_b_gt_0p9=float(np.mean(b > 0.9)), frac_rho_lt_0p1=float(np.mean(rho < 0.1)))
    if NOISE == "jit":
        res["sigma_ppm"] = {str(s): pct(np.exp(ch[:, 10 + i]) * 1e6) for i, s in enumerate(SECT)}
    json.dump(res, open(os.path.join(OUT, f"fit_{TAG}.json"), "w"), indent=1)
    print(json.dumps({kk: res[kk] for kk in ["k", "b", "rho", "T14_h", "Rp_Re", "rho_upper", "tau_max", "n_over_tau"]}))


if __name__ == "__main__":
    main()
