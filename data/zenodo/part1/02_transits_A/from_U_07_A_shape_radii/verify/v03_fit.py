"""Verifier joint fit (own code): batman quadratic LD, supersampled (S6 x15/1800 s, S32 x3/120 s, S98 x5/200 s),
circular orbit, on locally detrended light curves from v02_prep.py (|dt| < 0.25 d windows).
usage: v03_fit.py LCTAG RHO LD EXCL NSTEPS
  RHO: mann (split normal 2.51 +0.25/-0.22 in rho) | uni (flat in rho, 0.05-30) | logu (flat in ln rho)
  LD : claret (N(0.316,0.1), N(0.313,0.1)) | free (Kipping q1,q2 ~ U(0,1))
  EXCL: comma list of 'sector:epoch' to drop, or 'none'
Params: T0, P, k, b, ln rho, l1, l2, ln sigma_6, ln sigma_32, ln sigma_98 (white jitter per sector)."""
import sys, json, time
import numpy as np, emcee
from scipy.optimize import minimize
from vcommon import *
LCTAG, RHO, LD, EXCL, NST = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
TAG = f"{LCTAG}_{RHO}_{LD}_{EXCL.replace(':', 'e').replace(',', '+')}"
z = np.load(os.path.join(OUT, f"lc_{LCTAG}.npz"))
keep = np.abs(z["t"] - (T0A + z["e"] * PA)) < 0.25
if EXCL != "none":
    for it in EXCL.split(","):
        s_, e_ = map(int, it.split(":"))
        keep &= ~((z["s"] == s_) & (z["e"] == e_))
SEC = (6, 32, 98)
DS = []
for s in SEC:
    m = keep & (z["s"] == s)
    DS.append(dict(s=s, t=z["t"][m], f=z["f"][m], mod=Model(z["t"][m], F[s][2], F[s][3]), ne=len(np.unique(z["e"][m]))))
NEP = sum(d["ne"] for d in DS)
U0 = np.array(LD_PHX)


def uu(l1, l2):
    if LD == "free":
        if not (0 < l1 < 1 and 0 < l2 < 1):
            return None
        s = np.sqrt(l1); return 2 * s * l2, s * (1 - 2 * l2)
    if not (l1 > 0 and l1 + l2 < 1 and l1 + 2 * l2 > 0):
        return None
    return l1, l2


def lnp(th):
    T0, P, k, b, lr, l1, l2 = th[:7]
    ls = th[7:]
    if not (0.005 < k < 0.08 and 0 <= b < 1 + k and np.log(0.05) < lr < np.log(30) and np.all(np.abs(ls + 7.5) < 3.5)):
        return -np.inf
    u = uu(l1, l2)
    if u is None:
        return -np.inf
    rho = np.exp(lr)
    lp = 0.0
    if RHO == "mann":
        lp += -0.5 * ((rho - RHO0) / (RHOP if rho > RHO0 else RHOM)) ** 2 + lr
    elif RHO == "uni":
        lp += lr
    if LD == "claret":
        lp += -0.5 * np.sum(((np.array(u) - U0) / 0.1) ** 2)
    a = ars(rho, P)
    if b >= a:
        return -np.inf
    ll = 0.0
    for i, d in enumerate(DS):
        r = d["f"] - d["mod"](T0, P, k, b, a, u[0], u[1])
        s2 = np.exp(2 * ls[i])
        ll += -0.5 * np.sum(r * r) / s2 - 0.5 * len(r) * np.log(s2)
    return lp + ll


def main():
    t0 = time.time()
    l12 = list(U0) if LD == "claret" else [(U0.sum()) ** 2, U0[0] / (2 * U0.sum())]
    x0 = np.array([T0A, PA, 0.021, 0.3, np.log(2.5)] + l12 + [np.log(2e-4), np.log(8.5e-4), np.log(6.2e-4)])
    r = minimize(lambda x: -lnp(x) if np.isfinite(lnp(x)) else 1e25, x0, method="Nelder-Mead",
                 options=dict(maxiter=15000, maxfev=15000, xatol=1e-9, fatol=1e-4))
    r = minimize(lambda x: -lnp(x) if np.isfinite(lnp(x)) else 1e25, r.x, method="Nelder-Mead",
                 options=dict(maxiter=15000, maxfev=15000, xatol=1e-9, fatol=1e-4))
    best = r.x
    print(TAG, "MAP", np.round(best, 6), "lnp", -r.fun, "epochs", NEP, flush=True)
    nw, nd = 40, len(best)
    rng = np.random.default_rng(2026)
    sc = np.array([1e-3, 3e-6, 5e-4, 0.05, 0.03, 0.03, 0.03, 0.01, 0.01, 0.01])
    p0 = []
    while len(p0) < nw:
        c = best + sc * rng.standard_normal(nd); c[3] = abs(c[3])
        if np.isfinite(lnp(c)):
            p0.append(c)
    sam = emcee.EnsembleSampler(nw, nd, lnp)
    for i, _ in enumerate(sam.sample(np.array(p0), iterations=NST)):
        if (i + 1) % 1000 == 0:
            print("step", i + 1, "tau_max %.0f" % np.max(sam.get_autocorr_time(tol=0)), "%.0fs" % (time.time() - t0), flush=True)
    tau = sam.get_autocorr_time(tol=0); tm = float(np.max(tau))
    ch = sam.get_chain(discard=int(3 * tm), thin=max(1, int(tm / 5)), flat=True)
    np.savez_compressed(os.path.join(OUT, f"chain_{TAG}.npz"), chain=ch)
    T0, P, k, b, rho = ch[:, 0], ch[:, 1], ch[:, 2], ch[:, 3], np.exp(ch[:, 4])
    u = np.array([uu(a, c) for a, c in zip(ch[:, 5], ch[:, 6])])
    a = ars(rho, P); T14 = t14(P, k, b, a) * 24
    Rs = np.random.default_rng(1).normal(RS, ERS, len(k)); Rp = k * Rs * REARTH_PER_RSUN
    q = lambda x: np.percentile(x, [2.275, 15.865, 50, 84.135, 97.725]).tolist()
    res = dict(tag=TAG, lc=LCTAG, rho_prior=RHO, ld=LD, excl=EXCL, n_epochs=NEP, npts=int(sum(len(d["t"]) for d in DS)),
               nsteps=NST, tau_max=tm, n_over_tau=NST / tm, nsamp=len(ch), acc=float(np.mean(sam.acceptance_fraction)),
               k=q(k), b=q(b), rho=q(rho), T14_h=q(T14), Rp=q(Rp), Rp_mean_sd=[float(Rp.mean()), float(Rp.std())],
               u1=q(u[:, 0]), u2=q(u[:, 1]), rho_q977=float(np.percentile(rho, 97.725)), sig_ppm=[q(np.exp(ch[:, 7 + i]) * 1e6)[2] for i in range(3)],
               map=best.tolist(), runtime_s=time.time() - t0)
    json.dump(res, open(os.path.join(OUT, f"fit_{TAG}.json"), "w"), indent=1)
    print(json.dumps({kk: res[kk] for kk in ("tag", "n_epochs", "k", "b", "T14_h", "Rp", "rho_q977", "n_over_tau")}), flush=True)


if __name__ == "__main__":
    main()
