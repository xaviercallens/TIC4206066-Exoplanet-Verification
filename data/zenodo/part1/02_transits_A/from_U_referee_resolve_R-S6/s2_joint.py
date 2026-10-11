"""R-S6 (b), (c)-depth, (d): joint batman + emcee fit of Signal A with and without Sector 6.

usage: s2_joint.py SECTORS BASE NOISE [options as KEY=VAL]
  SECTORS : 6,32,98 | 32,98
  BASE    : raw  (raw normalised flux, one linear baseline per transit window, profiled analytically)
            bw   (wotan biweight 0.75 d, A and B masked over 2*T14, then the same per-window linear baseline)
  NOISE   : jit  (free white sigma per sector) | beta (sigma = beta * sigma_white per sector, from MAP residuals;
            time-averaging method, Pont et al. 2006 / Winn et al. 2008, bins 30-120 min inside windows)
  options : SS6=15 (S6 supersampling; 1 = none)  SPLIT6=1 (separate k6 and a timing offset dt6 [min] for S6)
            EXCL=98:-1 (drop sector:epoch windows, comma separated)  NSTEPS=20000  TAG=...
Model: circular orbit, quadratic LD with Gaussian prior sigma 0.1 around Claret (2017), Mann split-normal density prior
2.51 +0.25/-0.22 rho_sun.  Windows +-0.30 d around each A transit, B cadences (|dt_B| < 1.65 h) removed.
"""
import sys, json, time
import numpy as np
import emcee
from scipy.optimize import minimize
from rs6_common import *

SECT = tuple(int(s) for s in sys.argv[1].split(","))
BASE, NOISE = sys.argv[2], sys.argv[3]
OPT = dict(a.split("=", 1) for a in sys.argv[4:])
SS6 = int(OPT.get("SS6", 15))
SPLIT6 = OPT.get("SPLIT6", "0") == "1"
EXCL = [tuple(int(y) for y in x.split(":")) for x in OPT.get("EXCL", "").split(",") if x.strip()]
NSTEPS = int(OPT.get("NSTEPS", 20000))
RHO = OPT.get("RHO", "mann")     # mann | uni (uniform in rho, 0.05-30 rho_sun: density free)
TAG = OPT.get("TAG", f"S{'_'.join(map(str, SECT))}_{BASE}_{NOISE}_ss6{SS6}" + ("_split6" if SPLIT6 else "") + ("" if RHO == "mann" else f"_rho{RHO}") +
              ("_excl" + "_".join(f"{a}e{b}" for a, b in EXCL) if EXCL else ""))
U0 = claret_u()
SIG_U = 0.1
rng = np.random.default_rng(7)


def get_data(sec):
    t, f, meta = load(sec)
    dA, ep = phase_dt(t, P_A, T0_A)
    dB, _ = phase_dt(t, P_B, T0_B)
    if BASE == "bw":
        from wotan import flatten, transit_mask
        m = transit_mask(time=t, period=P_A, duration=2 * T14_A, T0=T0_A) | \
            transit_mask(time=t, period=P_B, duration=2 * T14_B, T0=T0_B)
        fl = flatten(t, f, method="biweight", window_length=0.75, break_tolerance=0.5, mask=m)
        ok = np.isfinite(fl)
        t, f, dA, ep, dB = t[ok], fl[ok], dA[ok], ep[ok], dB[ok]
    sel = (np.abs(dA) < 0.30) & (np.abs(dB) > T14_B)
    for s_, e_ in EXCL:
        if s_ == sec:
            sel &= ep != e_
    t, f, ep = t[sel], f[sel], ep[sel]
    x = t - (T0_A + ep * P_A)
    _, wid = np.unique(ep, return_inverse=True)
    nw = wid.max() + 1
    S0 = np.bincount(wid, minlength=nw).astype(float); S1 = np.bincount(wid, x, nw); S2 = np.bincount(wid, x * x, nw)
    ss = SS6 if sec == 6 else None
    return dict(s=sec, t=t, f=f, ep=ep, x=x, wid=wid, nw=nw, S0=S0, S1=S1, S2=S2, det=S0 * S2 - S1 ** 2,
                mod=Model(t, sec, ss), meta=meta, epochs=sorted(set(ep.tolist())))


D = [get_data(s) for s in SECT]
NS = len(D)
NAMES = ["T0", "P", "k", "b", "lnrho", "u1", "u2"] + [f"lns{d['s']}" for d in D] + (["k6", "dt6"] if SPLIT6 else [])
I_S = 7
I_K6 = 7 + NS
BETA = None
LS_PIN = None


def resid(th, d):
    T0, P, k, b, lnrho, u1, u2 = th[:7]
    ars = ars_from_rho(np.exp(lnrho), P)
    if d["s"] == 6 and SPLIT6:
        k_, T0_ = th[I_K6], T0 + th[I_K6 + 1] / 1440.0
    else:
        k_, T0_ = k, T0
    m = d["mod"].flux(T0_, P, k_, b, ars, u1, u2)
    r = d["f"] - m
    Sy = np.bincount(d["wid"], r, d["nw"]); Sxy = np.bincount(d["wid"], r * d["x"], d["nw"])
    a0 = (d["S2"] * Sy - d["S1"] * Sxy) / d["det"]; a1 = (d["S0"] * Sxy - d["S1"] * Sy) / d["det"]
    return r - a0[d["wid"]] - a1[d["wid"]] * d["x"]


def log_prior(th):
    T0, P, k, b, lnrho, u1, u2 = th[:7]
    rho = np.exp(lnrho)
    if not (0.004 < k < 0.1 and 0 <= b < 1 + k and 0.05 < rho < 30 and abs(T0 - T0_A) < 0.02 and abs(P - P_A) < 4e-5):
        return -np.inf
    if not (u1 > 0 and u1 + u2 < 1 and u1 + 2 * u2 > 0):
        return -np.inf
    if b >= ars_from_rho(rho, P):
        return -np.inf
    if RHO == "mann":
        sg = RHO_EP if rho > RHO_MED else RHO_EM
        lp = -0.5 * ((rho - RHO_MED) / sg) ** 2 + lnrho
    else:
        lp = lnrho          # uniform in rho, sampled in ln rho
    lp += -0.5 * np.sum(((np.array([u1, u2]) - U0) / SIG_U) ** 2)
    if NOISE == "jit":
        ls = th[I_S:I_S + NS]
        if np.any(ls < np.log(1e-5)) or np.any(ls > np.log(1e-2)):
            return -np.inf
    else:   # beta mode: the jitter slots are unused; pin them to keep the sampler well defined
        lp += -0.5 * np.sum(((np.asarray(th[I_S:I_S + NS]) - LS_PIN) / 0.01) ** 2)
    if SPLIT6:
        if not (0.004 < th[I_K6] < 0.1 and abs(th[I_K6 + 1]) < 120):
            return -np.inf
    return lp


def log_like(th):
    ll = 0.0
    for i, d in enumerate(D):
        r = resid(th, d)
        sig = np.exp(th[I_S + i]) if NOISE == "jit" else BETA[i]
        ll += -0.5 * np.sum(r * r) / sig ** 2 - len(r) * np.log(sig)
    return ll


def log_prob(th):
    lp = log_prior(th)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_like(th)
    return lp + ll if np.isfinite(ll) else -np.inf


def beta_factor(r, t, wid, cad_s, bins_min=(30, 60, 90, 120)):
    s1 = np.std(r)
    out = {}
    for bm in bins_min:
        n = int(round(bm * 60.0 / cad_s))
        if n < 2:
            continue
        means = []
        for w in np.unique(wid):
            rw = r[wid == w]
            m = len(rw) // n
            if m < 1:
                continue
            means.extend(rw[:m * n].reshape(m, n).mean(1).tolist())
        means = np.array(means)
        M = len(means)
        if M < 5:
            continue
        exp = s1 / np.sqrt(n) * np.sqrt(M / (M - 1.0))
        out[bm] = float(np.std(means) / exp)
    b = max(1.0, float(np.median(list(out.values())))) if out else 1.0
    return b, out, float(s1)


def main():
    global BETA
    t0 = time.time()
    x0 = [3998.7360, 3.1827830, 0.0210, 0.35, np.log(2.5), U0[0], U0[1]]   # start near the s1 timing solution
    x0 += [np.log(np.std(d["f"])) for d in D]
    if SPLIT6:
        x0 += [0.0210, 0.0]
    x0 = np.array(x0)
    # MAP with free jitter first (also used to define beta)
    NOISE_SAVE = NOISE
    globals()["NOISE"] = "jit"
    neg = lambda x: -log_prob(x) if np.isfinite(log_prob(x)) else 1e25
    best = x0
    for _ in range(4):
        r = minimize(neg, best, method="Nelder-Mead", options=dict(maxiter=30000, maxfev=30000, xatol=1e-10, fatol=1e-4))
        best = r.x
    globals()["NOISE"] = NOISE_SAVE
    betas = {}
    bl = []
    for i, d in enumerate(D):
        res = resid(best, d)
        b_, per, s1 = beta_factor(res, d["t"], d["wid"], CAD_S[d["s"]])
        betas[str(d["s"])] = dict(beta=b_, by_bin=per, sigma_white=s1, n=int(len(res)))
        bl.append(b_ * s1)
    BETA = np.array(bl)
    globals()["LS_PIN"] = np.array(best[I_S:I_S + NS])
    print(TAG, "npts", {d["s"]: len(d["t"]) for d in D}, "epochs", {d["s"]: d["epochs"] for d in D}, flush=True)
    print("MAP", dict(zip(NAMES, np.round(best, 7))), "betas", betas, flush=True)
    ndim = len(best)
    nw = max(40, 3 * ndim)
    scale = np.array([5e-4, 2e-6, 5e-4, 0.05, 0.03, 0.03, 0.03] + [0.02] * NS + ([5e-4, 2.0] if SPLIT6 else []))
    p0 = []
    while len(p0) < nw:
        c = best + scale * rng.standard_normal(ndim)
        c[3] = abs(c[3])
        if np.isfinite(log_prob(c)):
            p0.append(c)
    sampler = emcee.EnsembleSampler(nw, ndim, log_prob, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    for i, _ in enumerate(sampler.sample(np.array(p0), iterations=NSTEPS, progress=False)):
        if (i + 1) % 2000 == 0:
            tau = sampler.get_autocorr_time(tol=0)
            print(f"step {i+1} tau_max {np.max(tau):.1f} acc {np.mean(sampler.acceptance_fraction):.3f} {time.time()-t0:.0f}s", flush=True)
            if (i + 1) >= 8000 and (i + 1) > 60 * np.max(tau):
                break
    tau = sampler.get_autocorr_time(tol=0)
    burn = int(3 * np.max(tau)); thin = max(1, int(np.max(tau) / 4))
    ch = sampler.get_chain(discard=burn, thin=thin, flat=True)
    np.savez_compressed(os.path.join(OUT, f"chain_{TAG}.npz"), chain=ch, names=NAMES)
    T0, P, k, b, rho = ch[:, 0], ch[:, 1], ch[:, 2], ch[:, 3], np.exp(ch[:, 4])
    ars = ars_from_rho(rho, P)
    t14 = t14_hours(P, k, b, ars)
    rs = rng.normal(R_STAR, E_R_STAR, len(k))
    rp = k * rs * RSUN_OVER_REARTH
    res = dict(tag=TAG, rho_prior=RHO, sectors=SECT, base=BASE, noise=NOISE, ss6=SS6, split6=SPLIT6, excl=EXCL,
               npts={str(d["s"]): int(len(d["t"])) for d in D}, epochs={str(d["s"]): d["epochs"] for d in D},
               n_windows={str(d["s"]): int(d["nw"]) for d in D}, supersample={str(d["s"]): int(d["mod"].ss) for d in D},
               data={str(d["s"]): d["meta"] for d in D}, ld_center=U0.tolist(), betas=betas,
               n_walkers=nw, n_steps=int(sampler.iteration), tau=tau.tolist(), n_over_tau=float(sampler.iteration / np.max(tau)),
               n_samples=int(len(ch)), acceptance=float(np.mean(sampler.acceptance_fraction)),
               map=dict(zip(NAMES, best.tolist())), runtime_s=time.time() - t0,
               quantiles="2.275,15.865,50,84.135,97.725",
               T0=pct(T0), P=pct(P), k=pct(k), k_mean_std=[float(k.mean()), float(k.std())], b=pct(b), rho=pct(rho),
               T14_h=pct(t14), depth_k2_ppm=pct(k ** 2 * 1e6), Rp_Re=pct(rp), Rp_mean_std=[float(rp.mean()), float(rp.std())],
               rho_q97_7=float(np.percentile(rho, 97.725)), corr_T0_P=float(np.corrcoef(T0, P)[0, 1]))
    if SPLIT6:
        k6, dt6 = ch[:, I_K6], ch[:, I_K6 + 1]
        res.update(k6=pct(k6), dt6_min=pct(dt6), k6_minus_k=pct(k6 - k), k6_over_k=pct(k6 / k),
                   z_k6=float(np.mean(k6 - k) / np.std(k6 - k)), z_dt6=float(np.mean(dt6) / np.std(dt6)),
                   depth_ratio_k6sq_over_ksq=pct((k6 / k) ** 2))
    if NOISE == "jit":
        res["sigma_ppm"] = {str(d["s"]): pct(np.exp(ch[:, I_S + i]) * 1e6) for i, d in enumerate(D)}
    json.dump(res, open(os.path.join(OUT, f"fit_{TAG}.json"), "w"), indent=1, default=float)
    print(json.dumps({kk: res[kk] for kk in ["k", "b", "T14_h", "Rp_Re", "P", "T0", "n_over_tau"] +
                      (["k6", "dt6_min", "k6_minus_k", "z_k6", "z_dt6"] if SPLIT6 else [])}))


if __name__ == "__main__":
    main()
