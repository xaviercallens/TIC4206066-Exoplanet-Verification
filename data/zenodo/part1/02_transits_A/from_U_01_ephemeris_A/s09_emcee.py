"""Route 2b (part 2): emcee posterior of the global transit fit (P, T0 free; per-window linear baselines
profiled; uniform priors: T0, P, k in (0.005,0.08), b in [0,1+k), aRs in (3,40), q1,q2 in (0,1)).
Per-sector sigma inflated by the red-noise beta factor (time-averaging method, bins 30-120 min) of the
out-of-transit detrended flux."""
import numpy as np, json, sys, time
import emcee
from a_common import *
from a_global import build_windows, window_chi2, TAGS_BASE

def beta_factor(tag):
    z = np.load(f'lc/{tag}.npz'); t = z['t']; y = z['flat']; cad = float(z['exptime'])
    oot = ~(transit_mask(t, MS_PA, MS_T0A, T14A) | transit_mask(t, MS_PB, MS_T0B, T14B))
    t = t[oot]; y = y[oot] - 1; s1 = np.std(y)
    betas = []
    for binmin in [30, 45, 60, 90, 120]:
        nb = binmin / 1440. / cad
        if nb < 2: continue
        edges = np.arange(t.min(), t.max() + binmin / 1440., binmin / 1440.)
        idx = np.digitize(t, edges)
        means = []; counts = []
        for k in np.unique(idx):
            w = idx == k
            if w.sum() >= 0.8 * nb: means.append(y[w].mean()); counts.append(w.sum())
        if len(means) < 10: continue
        M = np.mean(counts); expected = s1 / np.sqrt(M) * np.sqrt(len(means) / (len(means) - 1))
        betas.append(np.std(means) / expected)
    return float(np.median(betas)) if betas else 1.0, betas

variant = sys.argv[1] if len(sys.argv) > 1 else 'base'
tags = dict(TAGS_BASE)
W = build_windows(tags=tags, order=1)
betas = {}
for sec, tag in tags.items():
    b, bl = beta_factor(tag); betas[sec] = dict(beta=max(1.0, b), per_bin=bl)
for w in W:
    w['sig'] = w['sig'] * betas[w['sector']]['beta']
print('betas', {s: round(v['beta'], 3) for s, v in betas.items()})
lsq = json.load(open('out_s05_globalfit_lsq.json'))['raw_lin']['params']
p0 = np.array([lsq[n][0] for n in ['T0', 'P', 'k', 'b', 'aRs', 'q1', 'q2']])

def lnprob(p):
    T0, P, k, b, aRs, q1, q2 = p
    if not (0.005 < k < 0.08 and 0 <= b < 1 + k and 3 < aRs < 40 and 0 < q1 < 1 and 0 < q2 < 1):
        return -np.inf
    if b >= aRs: return -np.inf
    return -0.5 * window_chi2(W, T0, P, k, b, aRs, q1, q2)

ndim, nw = 7, 32
nsteps = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
rng = np.random.default_rng(42)
scale = np.array([5e-4, 1e-6, 1e-3, 0.05, 0.3, 0.05, 0.05])
start = p0 + scale * rng.normal(size=(nw, ndim))
start[:, 3] = np.abs(start[:, 3]); start[:, 5:] = np.clip(start[:, 5:], 0.02, 0.98)
t0 = time.time()
sampler = emcee.EnsembleSampler(nw, ndim, lnprob)
sampler.run_mcmc(start, nsteps, progress=False)
chain = sampler.get_chain()
lp = sampler.get_log_prob()
np.savez_compressed(f'out_s09_emcee_chain_{variant}.npz', chain=chain[::5], lnp=lp[::5])
burn = nsteps // 3
flat = chain[burn:].reshape(-1, ndim)
try:
    tau = sampler.get_autocorr_time(discard=burn, quiet=True).tolist()
except Exception as ex:
    tau = str(ex)
names = ['T0', 'P', 'k', 'b', 'aRs', 'q1', 'q2']
res = dict(variant=variant, betas=betas, nwalkers=nw, nsteps=nsteps, burn=burn, autocorr=tau,
           acceptance=float(np.mean(sampler.acceptance_fraction)), runtime_s=time.time() - t0,
           n_windows=len(W), n_per_sector={s: sum(w['sector'] == s for w in W) for s in (6, 32, 98)})
for i, n in enumerate(names):
    q = np.percentile(flat[:, i], [2.275, 15.865, 50, 84.135, 97.725])
    res[n] = dict(median=float(q[2]), lo=float(q[2] - q[1]), hi=float(q[3] - q[2]), p2=float(q[0]), p98=float(q[4]),
                  mean=float(flat[:, i].mean()), std=float(flat[:, i].std()))
# T0 at weighted centre and at manuscript reference
C = np.cov(flat[:, 0], flat[:, 1]); nc = -C[0, 1] / C[1, 1]
tc = flat[:, 0] + nc * flat[:, 1]
res['n_center'] = float(nc); res['T0_center'] = [float(np.median(tc)), float(np.std(tc))]
res['corr_T0_P'] = float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]))
res['P_minus_ms_sigma'] = float((res['P']['median'] - MS_PA) / np.hypot(res['P']['std'], MS_ePA))
res['T0_minus_ms_min'] = float((res['T0']['median'] - MS_T0A) * 1440)
print(json.dumps({k: v for k, v in res.items() if k not in ('betas',)}, indent=1))
json.dump(res, open(f'out_s09_emcee_{variant}.json', 'w'), indent=1)
