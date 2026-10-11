"""T2_A_shape r02: own V-vs-U shape test for signal A with three noise treatments, on two stored detrendings from U:
  LC1 = U/07_A_shape_radii/verify/out/lc_o2_hw0.45_SYS_RM_FLUX.npz (per-transit local quadratic baselines, segment-aware, B removed)
  LC2 = U/07_A_shape_radii/out/lc_base_S{6,32,98}.npz (wotan biweight 0.75 d; B cadences removed here, |dt_B| < 1.65 h)
Models (independent code): exposure-integrated symmetric trapezoid (depth, T14, T12 free) vs triangle (T12 = T14/2),
plus one flux offset per sector; data |t - tc| < 0.25 d. Ephemeris fixed (U07 raw-flux joint-fit medians).
Noise treatments: (1) white: per-sector sigma = std of out-of-transit residuals; (2) beta: sigma*beta, beta = time-averaging
factor (Pont+2006 / Winn+2008) of trapezoid residuals, median over bins 30,60,90,120 min, computed per sector inside the
transit windows, floor 1; (3) event bootstrap: resample transits (sector, epoch) with replacement (white sigma), 400 draws,
report percentiles of Delta chi2 = chi2(V) - chi2(U); plus leave-one-transit-out minimum.
Delta chi2 > 0 means the V (triangle) fits worse."""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
ROOT = "/Users/paulsmith/space <3/tess-jev"
U7 = os.path.join(ROOT, "addenda/opus55_final_audit_tic4206066/07_A_shape_radii")
W = os.path.join(ROOT, "addenda/tic4206066_FINAL_RECONCILIATION/_work/T2_A_shape")
PA, T0A = 3.1827845765743845, 3998.736563362111          # U07 out/tests_results.json median_params (biweight chain); raw chain differs < 1e-5 d
PB, T0B, T14B = 11.13274, 4000.865, 1.65 / 24
CAD = {6: (1800.0, 15), 32: (120.0, 3), 98: (200.0, 5)}
T14G = 1.87 / 24
rng = np.random.default_rng(20260925)

def load_lc1():
    z = np.load(os.path.join(U7, "verify/out/lc_o2_hw0.45_SYS_RM_FLUX.npz"))
    return z["t"], z["f"], z["s"].astype(int), z["e"].astype(int)

def load_lc2():
    T, Fl, S, E = [], [], [], []
    for s in (6, 32, 98):
        z = np.load(os.path.join(U7, f"out/lc_base_S{s}.npz"))
        t, f = z["t"], z["f"]
        e = np.round((t - T0A) / PA).astype(int)
        dB = t - (T0B + np.round((t - T0B) / PB) * PB)
        m = np.isfinite(f) & (np.abs(dB) > T14B)
        T.append(t[m]); Fl.append(f[m]); S.append(np.full(m.sum(), s)); E.append(e[m])
    return np.concatenate(T), np.concatenate(Fl), np.concatenate(S), np.concatenate(E)

def trap(x, T14, T12, d):
    x = np.abs(x); y = np.zeros_like(x); y[x <= T14 / 2 - T12] = 1.0
    ing = (x > T14 / 2 - T12) & (x < T14 / 2); y[ing] = (T14 / 2 - x[ing]) / T12
    return d * y

class Data:
    def __init__(self, t, f, s, e, excl=()):
        dt = t - (T0A + e * PA); keep = np.abs(dt) < 0.25
        for (ss, ee) in excl: keep &= ~((s == ss) & (e == ee))
        self.dt, self.f, self.s, self.e, self.t = dt[keep], f[keep], s[keep], e[keep], t[keep]
        self.sub = {}
        for sec in (6, 32, 98):
            cad, n = CAD[sec]; off = ((np.arange(n) + 0.5) / n - 0.5) * cad / 86400; self.sub[sec] = (off, n)
        self.events = sorted(set(zip(self.s.tolist(), self.e.tolist())))
    def select(self, evlist):
        idx = []
        for (ss, ee) in evlist: idx.append(np.where((self.s == ss) & (self.e == ee))[0])
        return np.concatenate(idx)

def model(p, mode, dt, s, sub):
    d, T14 = p[0], p[1]; T12 = T14 / 2 if mode == "tri" else p[2]; offs = p[2:] if mode == "tri" else p[3:]
    out = np.empty_like(dt)
    for i, sec in enumerate((6, 32, 98)):
        m = s == sec
        if not m.any(): continue
        off, n = sub[sec]; x = (dt[m][:, None] + off[None, :]).ravel()
        out[m] = 1 + offs[i] - trap(x, T14, max(T12, 1e-5), d).reshape(-1, n).mean(axis=1)
    return out

def fit(mode, dt, f, s, sig_pt, sub):
    best = None
    for T14g in (0.072, 0.078, 0.09, 0.105):
        for T12g in ((0.004, 0.012, 0.025) if mode == "trap" else (None,)):
            p0 = [5e-4, T14g] + ([] if mode == "tri" else [T12g]) + [0.0, 0.0, 0.0]
            lb = [0, 0.03] + ([] if mode == "tri" else [1e-4]) + [-2e-3] * 3
            ub = [3e-3, 0.2] + ([] if mode == "tri" else [0.1]) + [2e-3] * 3
            r = least_squares(lambda p: (f - model(p, mode, dt, s, sub)) / sig_pt, p0, bounds=(lb, ub), x_scale="jac")
            if mode == "trap" and r.x[2] > r.x[1] / 2 + 1e-9: continue
            c2 = float(np.sum(r.fun ** 2))
            if best is None or c2 < best[0]: best = (c2, r.x)
    return best

def beta_factor(res, t, e, s_arr, sec, cad):
    """time-averaging beta for one sector, bins inside each transit window only"""
    m = s_arr == sec; r, tt, ee = res[m], t[m], e[m]; sig1 = np.std(r); betas = {}
    for bmin in (30, 60, 90, 120):
        nb = max(int(round(bmin * 60 / cad)), 1)
        if nb < 2 and sec == 6 and bmin < 60: continue
        means = []
        for ev in np.unique(ee):
            rr = r[ee == ev]; order = np.argsort(tt[ee == ev]); rr = rr[order]
            for k in range(len(rr) // nb): means.append(rr[k * nb:(k + 1) * nb].mean())
        means = np.array(means); M = len(means)
        if M < 5 or nb < 2: continue
        exp = sig1 / np.sqrt(nb) * np.sqrt(M / (M - 1)); betas[bmin] = float(np.std(means) / exp)
    b = float(np.median(list(betas.values()))) if betas else 1.0
    return max(b, 1.0), betas

def analyse(D, label, nboot=400):
    oot = np.abs(D.dt) > 0.75 * T14G
    sig = {sec: float(np.std(D.f[(D.s == sec) & oot])) for sec in (6, 32, 98)}
    sig_pt = np.array([sig[x] for x in D.s])
    U = fit("trap", D.dt, D.f, D.s, sig_pt, D.sub); V = fit("tri", D.dt, D.f, D.s, sig_pt, D.sub)
    res_u = D.f - model(U[1], "trap", D.dt, D.s, D.sub)
    betas = {sec: beta_factor(res_u, D.t, D.e, D.s, sec, CAD[sec][0]) for sec in (6, 32, 98)}
    sig_b = sig_pt * np.array([betas[x][0] for x in D.s])
    Ub = fit("trap", D.dt, D.f, D.s, sig_b, D.sub); Vb = fit("tri", D.dt, D.f, D.s, sig_b, D.sub)
    npts = len(D.f); chi2red_u = U[0] / (npts - 6)
    # leave-one-transit-out
    loo = []
    for ev in D.events:
        idx = D.select([x for x in D.events if x != ev])
        u = fit("trap", D.dt[idx], D.f[idx], D.s[idx], sig_pt[idx], D.sub); v = fit("tri", D.dt[idx], D.f[idx], D.s[idx], sig_pt[idx], D.sub)
        loo.append((ev, v[0] - u[0]))
    # event bootstrap
    bs = []
    for _ in range(nboot):
        pick = [D.events[i] for i in rng.integers(0, len(D.events), len(D.events))]
        idx = D.select(pick)
        u = fit("trap", D.dt[idx], D.f[idx], D.s[idx], sig_pt[idx], D.sub); v = fit("tri", D.dt[idx], D.f[idx], D.s[idx], sig_pt[idx], D.sub)
        bs.append(v[0] - u[0])
    bs = np.array(bs)
    r = dict(label=label, n_points=npts, n_transits=len(D.events), sigma_white_ppm={k: v * 1e6 for k, v in sig.items()},
             chi2red_U_white=chi2red_u,
             white=dict(dchi2=V[0] - U[0], dchi2_rescaled_to_chi2red1=(V[0] - U[0]) / chi2red_u,
                        U=dict(depth_ppm=U[1][0] * 1e6, T14_h=U[1][1] * 24, T12_h=U[1][2] * 24), V=dict(depth_ppm=V[1][0] * 1e6, T14_h=V[1][1] * 24)),
             beta=dict(betas={k: v[0] for k, v in betas.items()}, betas_by_bin={k: v[1] for k, v in betas.items()}, dchi2=Vb[0] - Ub[0]),
             loo_min=dict(event=min(loo, key=lambda x: x[1])[0], dchi2=min(x[1] for x in loo)),
             loo_all=[[int(a[0]), int(a[1]), float(b)] for a, b in loo],
             bootstrap=dict(n=int(nboot), median=float(np.median(bs)), p2_5=float(np.percentile(bs, 2.5)), p16=float(np.percentile(bs, 16)),
                            frac_le_0=float(np.mean(bs <= 0)), frac_lt_9=float(np.mean(bs < 9))))
    print(json.dumps({k: v for k, v in r.items() if k != "loo_all"}), flush=True)
    return r

if __name__ == "__main__":
    nboot = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    out = []
    t, f, s, e = load_lc1()
    out.append(analyse(Data(t, f, s, e), "LC1_localquad_all", nboot))
    out.append(analyse(Data(t, f, s, e, excl=[(98, -1)]), "LC1_localquad_excl98e-1", nboot))
    t, f, s, e = load_lc2()
    out.append(analyse(Data(t, f, s, e), "LC2_biweight_all", nboot))
    out.append(analyse(Data(t, f, s, e, excl=[(98, -1)]), "LC2_biweight_excl98e-1", nboot))
    json.dump(out, open(os.path.join(W, "out/r02_vshape.json"), "w"), indent=1)
