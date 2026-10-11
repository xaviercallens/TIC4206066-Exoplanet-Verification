"""Joint emcee fit of signal B shape: k = Rp/R*, b, rho*, u1, u2, T0, P over all 7 archival+S98 events.

Usage: python s02_joint_fit.py <prior: prior|free|freeU> <mask: default|keep512> [nsteps]
  prior  : rho* ~ N(2.51, 0.24) rho_sun ; free: log-uniform rho in [0.3, 100] ; freeU: uniform rho in [0.3, 100]
  mask   : SPOC quality mask = lightkurve default 17087 (manuscript) or 16575 (keeps the PDC impulsive-outlier bit 512)
Data: S6 TESS-SPOC PDCSAP 30 min, S32 SPOC 2-min PDCSAP, S98 QLP SYS_RM_FLUX 200 s (scattered-light/QLP-bad bits removed),
biweight 0.75 d detrending with A and B in-transit cadences excluded from the trend, signal A subtracted with a fixed
template (k=0.0212, b=0.35, mid-depth 522 ppm), windows +-0.3 d around each B event, per-event linear baseline profiled
analytically (linear least squares inside the likelihood), per-sector white noise fixed at the robust OOT scatter.
Circular orbit. Supersampling to <= 1-min sub-exposures. Output: out/s02_joint_<prior>_<mask>.json (+ chain .npz)
"""
import sys, time
import numpy as np
import emcee
from b_run import *

prior = sys.argv[1] if len(sys.argv) > 1 else "prior"
maskname = sys.argv[2] if len(sys.argv) > 2 else "default"
nsteps = int(sys.argv[3]) if len(sys.argv) > 3 else 6000
spoc_mask = 17087 if maskname == "default" else 16575
rng = np.random.default_rng(42)

P = {6: "S6_tessspoc_pdc", 32: "S32_spoc2m_pdc", 98: "S98_qlp_sysrm"}
segs = []
for sec in (6, 32, 98):
    lc = get_lc(P[sec], spoc_mask=spoc_mask)
    lcd = detrend(lc)
    fdA = lcd["fd"] + 522.0e-6 * tplA().unit(phase_dt(lcd["t"], EPH_A), lcd["cad"], True)
    evs = events_in(lcd, EPH_B, T14=EPH_B["T14"], half_window=0.3, min_cov=0.25, drop_A_hw=1e-9)
    for ev in evs:
        w = ev["w"]
        ss = int(np.ceil(lcd["cad"] * 1440)); ss += (ss % 2 == 0)
        segs.append(dict(sec=sec, epoch=ev["epoch"], t=lcd["t"][w], y=fdA[w], sig=lcd["sig"], cad=lcd["cad"], ss=ss,
                         X0=np.vstack([np.ones(w.sum()), (lcd["t"][w] - ev["tc"]) / 0.3]).T))
print("events:", [(s["sec"], s["epoch"], len(s["t"])) for s in segs], flush=True)

# one batman model object per segment (params updated in place)
pm = batman.TransitParams()
pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w, pm.u, pm.limb_dark = EPH_B["T0"], EPH_B["P"], 0.03, 30.0, 89.0, 0.0, 90.0, [0.32, 0.31], "quadratic"
for s in segs:
    s["m"] = batman.TransitModel(pm, s["t"], supersample_factor=s["ss"], exp_time=s["cad"]) if s["ss"] > 1 else batman.TransitModel(pm, s["t"])
    XtX = s["X0"].T @ s["X0"]; s["proj"] = s["X0"] @ np.linalg.solve(XtX, s["X0"].T)

names = ["k", "b", "log10rho", "u1", "u2", "T0", "P"]


def lnprob(th):
    k, b, lr, u1, u2, T0, Pp = th
    rho = 10 ** lr
    if not (0.005 < k < 0.2 and 0 <= b < 1 + k and -0.53 < lr < 2.0):
        return -np.inf
    if not (u1 + u2 < 1 and u1 > 0 and u1 + 2 * u2 > 0):
        return -np.inf
    if abs(T0 - EPH_B["T0"]) > 0.05 or abs(Pp - EPH_B["P"]) > 0.002:
        return -np.inf
    lp = -0.5 * ((u1 - LD_NOM[0]) / LD_SIG[0]) ** 2 - 0.5 * ((u2 - LD_NOM[1]) / LD_SIG[1]) ** 2
    if prior == "prior":
        lp += -0.5 * ((rho - RHO_PRIOR[0]) / RHO_PRIOR[1]) ** 2
    elif prior == "freeU":
        lp += lr * np.log(10)          # uniform in rho  <=> p(log rho) ~ rho
    a = ars_from_rho(rho, Pp)
    if b >= a:
        return -np.inf
    pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.u = T0, Pp, k, a, np.degrees(np.arccos(b / a)), [u1, u2]
    chi2 = 0.0
    for s in segs:
        mod = s["m"].light_curve(pm)
        # profile per-event linear baseline (multiplicative on the model ~ additive at these depths)
        r = s["y"] - mod
        r = r - s["proj"] @ r
        chi2 += np.sum(r * r) / s["sig"] ** 2
    return lp - 0.5 * chi2


ndim, nw = len(names), 40
p0 = np.array([0.033 if prior == "prior" else 0.030, 0.85 if prior == "prior" else 0.5, np.log10(2.51 if prior == "prior" else 10.0),
               0.32, 0.31, EPH_B["T0"], EPH_B["P"]])
pos = p0 + np.array([0.001, 0.02, 0.02, 0.01, 0.01, 0.0005, 0.00001]) * rng.standard_normal((nw, ndim))
pos[:, 1] = np.abs(pos[:, 1])
t0 = time.time()
sampler = emcee.EnsembleSampler(nw, ndim, lnprob)
sampler.run_mcmc(pos, nsteps, progress=False)
chain = sampler.get_chain(discard=nsteps // 3, thin=5, flat=True)
try:
    tau = sampler.get_autocorr_time(discard=nsteps // 3, quiet=True).tolist()
except Exception as e:
    tau = str(e)
k, b, lr, u1, u2, T0, Pp = chain.T
rho = 10 ** lr
Rs = RSTAR + E_RSTAR * rng.standard_normal(len(k))
Rp = k * Rs * RSUN_REARTH
a = ars_from_rho(rho, Pp); inc = np.arccos(b / a)
T14 = Pp / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (a * np.sin(inc)), 0, 1)) * 24
# mid-transit depth of each posterior sample (batman at t=T0)
idx = rng.choice(len(k), size=min(3000, len(k)), replace=False)
mid = []
for i in idx:
    pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.u = 0.0, Pp[i], k[i], a[i], np.degrees(inc[i]), [u1[i], u2[i]]
    mid.append(1 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0])
mid = np.array(mid) * 1e6


def q(x):
    p = np.percentile(x, [2.275, 15.865, 50, 84.135, 97.725])
    return dict(med=float(p[2]), lo=float(p[2] - p[1]), hi=float(p[3] - p[2]), p2=float(p[0]), p98=float(p[4]), mean=float(np.mean(x)), sd=float(np.std(x)))


best = chain[np.argmax(sampler.get_log_prob(discard=nsteps // 3, thin=5, flat=True))]
res = dict(prior=prior, mask=maskname, spoc_mask=spoc_mask, nsteps=nsteps, nwalkers=nw, acceptance=float(np.mean(sampler.acceptance_fraction)),
           autocorr=tau, runtime_s=time.time() - t0, n_events=len(segs), events=[(s["sec"], s["epoch"]) for s in segs],
           k=q(k), b=q(b), rho=q(rho), u1=q(u1), u2=q(u2), T0=q(T0), P=q(Pp), Rp_earth=q(Rp), T14_h=q(T14), mid_depth_ppm=q(mid),
           k2_ppm=q(k ** 2 * 1e6), Rp_from_k_median=float(np.median(k) * RSTAR * RSUN_REARTH),
           best=dict(zip(names, map(float, best))), rho_prior=RHO_PRIOR if prior == "prior" else None,
           frac_b_gt_0p9=float(np.mean(b > 0.9)), frac_grazing=float(np.mean(b > 1 - k)))
jdump(res, f"s02_joint_{prior}_{maskname}.json")
np.savez_compressed(os.path.join(OUT, f"s02_chain_{prior}_{maskname}.npz"), chain=chain, names=names)
print(json.dumps({kk: res[kk] for kk in ("acceptance", "autocorr", "k", "b", "rho", "Rp_earth", "T14_h", "mid_depth_ppm", "frac_grazing")}, indent=1))
