"""Audit 2, steps 4-5: global (joint) fixed-shape fit for P and T0 with formal and beta-scaled errors; TTV tests on
the per-transit O-C: chi2, sinusoid periodogram with permutation and empirical-noise FAP, 95% upper limits on the
sinusoidal semi-amplitude (Neyman construction with the empirical injection error distribution)."""
import json
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy import stats
import t2_common as C

rng = np.random.default_rng(5)
ev = pd.read_csv(f"{C.OUT}/t3_events_base.csv")
meta = json.load(open(f"{C.OUT}/t3_meta_base.json"))
t5 = json.load(open(f"{C.OUT}/t5_results.json"))
beta = {int(k): v for k, v in t5["betas"].items()}
P0, T00 = meta["ephem_window"]["P"], meta["ephem_window"]["T0"]
tpl = C.Template(); T14 = tpl.t14
lcs = {s: C.load_lc(s) for s in (6, 32, 98)}
out = {}

# ------------------------------------------------------------------ global fit
def build(sel):
    W = []
    for _, r in ev[sel].iterrows():
        lc = lcs[r.sector]; _, tcs = C.a_epochs(lc["t"], P0, T00)
        tt, yy, cov = C.prepare_window(lc, r.tc_pred, tcs, T14)
        W.append(dict(t=tt, y=yy, n=int(r.epoch), s=int(r.sector), cad=lc["cad"], rms=float(np.std(yy)), ctr=r.tc_pred))
    return W

def gfit(W, use_beta):
    secs = [6, 32, 98]
    def unpack(p):
        return p[0], p[1], dict(zip(secs, p[2:5])), p[5:].reshape(len(W), 2)
    def resid(p):
        T0, P, D, B = unpack(p); rr = []
        for i, w in enumerate(W):
            x = w["t"] - w["ctr"]
            m = B[i, 0] + B[i, 1] * x - D[w["s"]] * tpl(w["t"] - (T0 + w["n"] * P), w["cad"])
            sig = w["rms"] * (beta[w["s"]] if use_beta else 1.0)
            rr.append((w["y"] - m) / sig)
        return np.concatenate(rr)
    p0 = np.r_[T00, P0, 5e-4, 5e-4, 5e-4, np.tile([1.0, 0.0], len(W))]
    sc = np.r_[1e-3, 1e-6, 1e-4, 1e-4, 1e-4, np.tile([1e-4, 1e-4], len(W))]
    r = least_squares(resid, p0, x_scale=sc, diff_step=None, method="trf")
    J = r.jac; cov = np.linalg.inv(J.T @ J)
    chi2 = float(np.sum(r.fun ** 2)); dof = len(r.fun) - len(p0)
    return dict(T0=float(r.x[0]), P=float(r.x[1]), eT0=float(np.sqrt(cov[0, 0])), eP=float(np.sqrt(cov[1, 1])),
                D6=float(r.x[2]), D32=float(r.x[3]), D98=float(r.x[4]), chi2=chi2, dof=int(dof), use_beta=use_beta,
                n_windows=len(W))

W20 = build(ev.sel)
out["global_fit_20"] = {"white": gfit(W20, False), "beta": gfit(W20, True)}
W23 = build(np.ones(len(ev), bool))
out["global_fit_23"] = {"white": gfit(W23, False), "beta": gfit(W23, True)}
for k, v in out.items():
    for kk, vv in v.items():
        print(k, kk, f"P={vv['P']:.8f}+-{vv['eP']:.8f} T0={vv['T0']:.5f}+-{vv['eT0']:.5f} chi2={vv['chi2']:.0f}/{vv['dof']} D={vv['D6']*1e6:.0f},{vv['D32']*1e6:.0f},{vv['D98']*1e6:.0f}")

# ------------------------------------------------------------------ TTV tests
inj = {s: pd.read_csv(f"{C.OUT}/t4_inj_S{s}.csv") for s in (6, 32, 98)}
for s, d in inj.items():
    d["sel"] = (d.coverage >= 0.7) & (d.snr >= 2) & (~d.edge)
    d["sel60"] = (d.coverage >= 0.7) & (d.snr_alt >= 2) & (~d.edge_alt)

def design(n, f=None):
    cols = [np.ones_like(n), n]
    if f is not None:
        cols += [np.sin(2 * np.pi * f * n), np.cos(2 * np.pi * f * n)]
    return np.vstack(cols).T

def wls_chi2(A, y, e):
    Aw = A / e[:, None]; yw = y / e
    c, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    return float(np.sum((yw - Aw @ c) ** 2)), c

ttv = {}
for label, rows, tcol, poolcol, dtcol in (("R144", ev[ev.sel], "tc_obs", "sel", "dt"),
                                          ("R60", ev[(ev.coverage >= 0.7) & (ev.snr_alt >= 2) & (~ev.edge_alt)], "tc_obs_alt", "sel60", "dt_alt")):
    n = rows.epoch.values.astype(float); y = rows[tcol].values * 1440.0  # minutes
    sd = {s: float(inj[s][dtcol][inj[s][poolcol]].std(ddof=1)) for s in inj}
    pools = {s: inj[s][dtcol][inj[s][poolcol]].values for s in inj}
    e = rows.sector.map(sd).values
    y0 = y - (T00 + P0 * n) * 1440.0   # remove a reference ephemeris for numerical conditioning
    chi_lin, c_lin = wls_chi2(design(n), y0, e)
    # frequency grid in cycles per epoch: from 1/(2*baseline) to 0.5 (Nyquist for integer epochs)
    base = n.max() - n.min()
    freqs = np.arange(1.0 / (2 * base), 0.5, 1.0 / (10 * base))
    def periodogram(yy, ee):
        cl, _ = wls_chi2(design(n), yy, ee)
        return np.array([cl - wls_chi2(design(n, f), yy, ee)[0] for f in freqs]), cl
    dchi, _ = periodogram(y0, e)
    jbest = int(np.argmax(dchi))
    # FAP 1: permutation of (residual, sigma) pairs among epochs
    oc = y0 - design(n) @ c_lin
    Nperm = 2000
    mx_perm = []
    for k in range(Nperm):
        pi = rng.permutation(len(n))
        mx_perm.append(periodogram(oc[pi], e[pi])[0].max())
    mx_perm = np.array(mx_perm)
    # FAP 2: noise drawn from the empirical injection distribution (per sector), no TTV
    mx_emp = []
    for k in range(Nperm):
        yy = np.array([rng.choice(pools[s]) for s in rows.sector])
        mx_emp.append(periodogram(yy, e)[0].max())
    mx_emp = np.array(mx_emp)
    # 95% upper limit on the semi-amplitude (Neyman): at each frequency the amplitude A such that 95% of
    # simulated data sets (sinusoid of amplitude A, random phase, empirical noise) give a fitted amplitude larger
    # than the observed fitted amplitude.  Coarse frequency grid to keep the cost down.
    fgrid = freqs[::5]
    ul = []
    Nsim = 400
    noise = np.array([[rng.choice(pools[s]) for s in rows.sector] for _ in range(Nsim)])  # (Nsim, nev)
    for f in fgrid:
        A = design(n, f); Aw = A / e[:, None]
        L = np.linalg.pinv(Aw)                          # (4, nev), applied to y/e
        a_obs = np.hypot(*(L @ (y0 / e))[2:4])
        phases = rng.uniform(0, 2 * np.pi, Nsim)
        def frac_above(Aamp):
            sig = Aamp * np.sin(2 * np.pi * f * n[None, :] + phases[:, None])
            Y = (sig + noise) / e[None, :]
            cc = Y @ L.T
            return np.mean(np.hypot(cc[:, 2], cc[:, 3]) > a_obs)
        lo, hi = 0.0, 5.0
        while frac_above(hi) < 0.95 and hi < 5000:
            hi *= 2
        for _ in range(25):
            mid = 0.5 * (lo + hi)
            if frac_above(mid) >= 0.95: hi = mid
            else: lo = mid
        ul.append(dict(f_cyc_per_epoch=float(f), P_ttv_epochs=float(1 / f), P_ttv_days=float(P0 / f), A_obs_min=float(a_obs), A95_min=float(hi)))
    ul = pd.DataFrame(ul)
    ul.to_csv(f"{C.OUT}/t6_ttv_upper_limits_{label}.csv", index=False)
    # summarise ranges
    def summ(m):
        u = ul[m]
        return dict(n=int(len(u)), median_A95_min=float(u.A95_min.median()), p90_A95_min=float(u.A95_min.quantile(0.9)),
                    min_A95_min=float(u.A95_min.min()), max_A95_min=float(u.A95_min.max()))
    sum_ul = {"P_ttv < 20 epochs (< 64 d, within-sector sampling)": summ(ul.P_ttv_epochs < 20),
              "20-200 epochs (64-640 d)": summ((ul.P_ttv_epochs >= 20) & (ul.P_ttv_epochs < 200)),
              ">= 200 epochs (> 640 d; partly degenerate with the linear ephemeris)": summ(ul.P_ttv_epochs >= 200)}
    chi_mc = []
    for k in range(20000):
        yy = np.array([rng.choice(pools[s]) for s in rows.sector])
        chi_mc.append(wls_chi2(design(n), yy, e)[0])
    chi_mc = np.array(chi_mc)
    ttv[label] = dict(n_events=int(len(n)), sigmas_min=sd, chi2_linear=chi_lin, dof=int(len(n) - 2),
                      p_chi2_gauss=float(stats.chi2.sf(chi_lin, len(n) - 2)),
                      p_chi2_empirical=float(np.mean(chi_mc >= chi_lin)),
                      best_freq=float(freqs[jbest]), best_P_ttv_epochs=float(1 / freqs[jbest]), best_P_ttv_days=float(P0 / freqs[jbest]),
                      best_dchi2=float(dchi[jbest]),
                      FAP_permutation=float(np.mean(mx_perm >= dchi[jbest])),
                      FAP_empirical_noise=float(np.mean(mx_emp >= dchi[jbest])),
                      n_freqs=int(len(freqs)), freq_step=float(freqs[1] - freqs[0]),
                      upper_limits_summary=sum_ul)
    np.savez(f"{C.OUT}/t6_periodogram_{label}.npz", freqs=freqs, dchi=dchi, mx_perm=mx_perm, mx_emp=mx_emp, n=n, oc=oc, e=e)
    print(label, json.dumps(ttv[label], indent=1))
out["ttv"] = ttv
json.dump(out, open(f"{C.OUT}/t6_results.json", "w"), indent=1)
