"""Shape / EB tests for Signal A, using the baseline joint-fit posterior (fit_mann_claret_jit_base).

1) red-noise beta per sector (time-averaging, Pont+2006 / Winn+2008) on joint-fit residuals -> beta_base.json
2) depth scaling of the fixed best-fit shape (supersampled per cadence, per-transit constant offset):
   per transit, per sector, odd vs even.  Realistic errors: the same fit repeated with the whole ephemeris shifted by
   random phase offsets (injection-free 'fake transit' sets, windows avoiding real A/B transits) -> empirical sigma.
3) secondary eclipse: box of width T14 at phase 0.5 (light curves detrended with phase 0.5 masked), empirical sigma
   from random same-width boxes; plus a scan over all phases (widths 0.5, 1, 2 x T14) on the base light curves.
4) V-shape: supersampled trapezoid (free T12) vs triangle (T12 = T14/2); delta chi2 with white and beta-scaled errors.
5) sqrt(depth) sanity check of Rp/R*.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from common import *

rng = np.random.default_rng(7)
FIT = json.load(open(os.path.join(OUT, "fit_mann_claret_jit_base.json")))
ch = np.load(os.path.join(OUT, "chain_mann_claret_jit_base.npz"))["chain"]
med = np.median(ch, axis=0)
T0, P, k, b, lnrho, u1, u2 = med[:7]
ars = float(ars_from_rho(np.exp(lnrho), P))
T14 = float(t14_hours(P, k, b, ars)) / 24.0
SECT = (6, 32, 98)
R = {}


def load(dv, s):
    z = np.load(os.path.join(OUT, f"lc_{dv}_S{s}.npz"))
    return z["t"], z["f"], float(z["cad"])


def shape(t, cad, s, tc_shift=0.0):
    """unit-depth template: (1 - F)/delta_mid for the median parameters, supersampled."""
    m = BatmanSet(t, cad, SUPERSAMPLE[s])
    F = m.flux(T0 + tc_shift, P, k, b, ars, u1, u2)
    return (1 - F)


# unsupersampled mid-transit depth of the median model
import batman
pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = k; pm.a = ars; pm.inc = np.degrees(np.arccos(b / ars))
pm.ecc = 0; pm.w = 90; pm.u = [u1, u2]; pm.limb_dark = "quadratic"
DMID = 1 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0]
R["median_params"] = dict(T0=T0, P=P, k=k, b=b, rho=float(np.exp(lnrho)), u1=u1, u2=u2, ars=ars, T14_h=T14 * 24,
                          depth_mid_ppm=DMID * 1e6)

# ---------------------------------------------------------------- 1) beta
beta = {}
for i, s in enumerate(SECT):
    t, f, cad = load("base", s)
    dA, _ = phase_dt(t, P, T0); dB, _ = phase_dt(t, P_B, T0_B)
    m = (np.abs(dA) < 0.3) & (np.abs(dB) > T14_B_H / 24)
    t, f = t[m], f[m]
    off = FIT["map"][f"c{s}"]
    res = f - (1 - shape(t, cad, s)) - off
    s1 = np.std(res)
    # per window time-averaging
    ep = np.round((t - T0) / P).astype(int)
    betas = {}
    for bm in (30, 60, 90, 120):
        n = int(round(bm * 60 / cad))
        if n < 2:
            continue
        binned = []
        for e in np.unique(ep):
            r = res[ep == e]
            M = len(r) // n
            if M < 1:
                continue
            binned += list(r[:M * n].reshape(M, n).mean(axis=1))
        binned = np.array(binned); M = len(binned)
        if M < 8:
            continue
        betas[bm] = float(np.std(binned) / (s1 / np.sqrt(n) * np.sqrt(M / (M - 1.0))))
    bval = max(1.0, float(np.median(list(betas.values()))))
    beta[str(s)] = dict(sigma_white=float(s1), beta=bval, betas_by_bin_min=betas, n=int(len(res)))
json.dump(beta, open(os.path.join(OUT, "beta_base.json"), "w"), indent=1)
R["beta"] = beta
print("beta", beta)


# ---------------------------------------------------------------- 2) depth scaling
def depth_fit(sets, shift=0.0):
    """sets: list of (t, f, cad, s) already restricted. One common depth, per-transit offsets.
    returns depth (units of DMID-template, i.e. mid-transit depth in ppm), formal sigma (white, per-sector sigma)."""
    A = []; y = []; w = []; groups = []
    gid = 0
    rows = []
    for (t, f, cad, s) in sets:
        tmp = shape(t, cad, s, shift)
        ep = np.round((t - T0 - shift) / P).astype(int)
        sig = beta[str(s)]["sigma_white"]
        for e in np.unique(ep):
            mm = ep == e
            rows.append((tmp[mm] / DMID, f[mm], sig))
    # analytic: for each window, remove offset by projecting out the mean (weighted)
    num = 0.0; den = 0.0
    for tmp, f, sig in rows:
        x = tmp - tmp.mean(); yy = (1 - f) - (1 - f).mean()
        num += np.sum(x * yy) / sig ** 2; den += np.sum(x * x) / sig ** 2
    if den == 0:
        return np.nan, np.nan
    return num / den, 1 / np.sqrt(den)


def select(dv, s, parity=None, shift=0.0, sector_epochs=None):
    t, f, cad = load(dv, s)
    dA, epA = phase_dt(t, P, T0)
    dB, _ = phase_dt(t, P_B, T0_B)
    dS, epS = phase_dt(t, P, T0 + shift)
    m = (np.abs(dS) < 0.3) & (np.abs(dB) > T14_B_H / 24)
    if shift != 0.0:
        m &= np.abs(dA) > T14 * 1.0          # keep real transits out of fake windows
    if parity is not None:
        m &= (np.abs(epS) % 2) == parity
    return (t[m], f[m], cad, s)


def empirical(fn, nfake=400):
    vals = []
    shifts = rng.uniform(0.12, 0.88, nfake) * P
    for sh in shifts:
        if abs(sh - 0.5 * P) < 0.06 * P:
            continue
        v = fn(sh)
        if np.isfinite(v):
            vals.append(v)
    vals = np.array(vals)
    return float(1.4826 * np.median(np.abs(vals - np.median(vals)))), float(np.std(vals)), len(vals)


out2 = {}
# per sector
for s in SECT:
    d, e = depth_fit([select("base", s)])
    emp_mad, emp_std, n = empirical(lambda sh: depth_fit([select("base", s, shift=sh)], sh)[0])
    out2[f"S{s}"] = dict(depth_ppm=d * 1e6, sig_formal_ppm=e * 1e6, sig_emp_mad_ppm=emp_mad * 1e6,
                         sig_emp_std_ppm=emp_std * 1e6, n_fake=n)
# all
d, e = depth_fit([select("base", s) for s in SECT])
emp = empirical(lambda sh: depth_fit([select("base", s, shift=sh) for s in SECT], sh)[0])
out2["all"] = dict(depth_ppm=d * 1e6, sig_formal_ppm=e * 1e6, sig_emp_mad_ppm=emp[0] * 1e6,
                   sig_emp_std_ppm=emp[1] * 1e6)
# odd / even
for par, lab in ((1, "odd"), (0, "even")):
    d, e = depth_fit([select("base", s, parity=par) for s in SECT])
    emp = empirical(lambda sh: depth_fit([select("base", s, parity=par, shift=sh) for s in SECT], sh)[0])
    out2[lab] = dict(depth_ppm=d * 1e6, sig_formal_ppm=e * 1e6, sig_emp_mad_ppm=emp[0] * 1e6,
                     sig_emp_std_ppm=emp[1] * 1e6)


def sect_chi2(key):
    ds = np.array([out2[f"S{s}"]["depth_ppm"] for s in SECT])
    ss = np.array([max(out2[f"S{s}"]["sig_formal_ppm"], out2[f"S{s}"][key]) for s in SECT])
    w = 1 / ss ** 2; mu = np.sum(w * ds) / np.sum(w)
    return dict(mean=float(mu), chi2=float(np.sum((ds - mu) ** 2 * w)), dof=2)


out2["sector_consistency_emp_std"] = sect_chi2("sig_emp_std_ppm")
for key in ("sig_formal_ppm", "sig_emp_std_ppm"):
    so = max(out2["odd"]["sig_formal_ppm"], out2["odd"][key]); se = max(out2["even"]["sig_formal_ppm"], out2["even"][key])
    out2[f"odd_even_diff_nsigma_{key}"] = float((out2["odd"]["depth_ppm"] - out2["even"]["depth_ppm"]) / np.hypot(so, se))
# per transit
pt = []
for s in SECT:
    t, f, cad, _ = select("base", s)
    ep = np.round((t - T0) / P).astype(int)
    for e in np.unique(ep):
        mm = ep == e
        inn = np.abs(t[mm] - (T0 + e * P)) < T14 / 2
        if inn.sum() == 0:
            continue
        d, er = depth_fit([(t[mm], f[mm], cad, s)])
        pt.append(dict(sector=s, epoch=int(e), n_in=int(inn.sum()), depth_ppm=d * 1e6, sig_white_ppm=er * 1e6))
out2["per_transit"] = pt
dd = np.array([p["depth_ppm"] for p in pt]); ee = np.array([p["sig_white_ppm"] for p in pt])
w = 1 / ee ** 2; mu = np.sum(w * dd) / np.sum(w)
out2["per_transit_chi2"] = dict(n=len(pt), wmean=float(mu), chi2=float(np.sum((dd - mu) ** 2 * w)), dof=len(pt) - 1)
R["depths"] = out2
print(json.dumps({kk: v for kk, v in out2.items() if kk != "per_transit"}, indent=1))


# ---------------------------------------------------------------- 3) secondary
def box_depth(dv, phase, width_d, sectors=SECT, excl_real=True):
    """flux in box (mean, relative to local +-0.3 d baseline) centred on the given phase; sector-weighted."""
    vals = []; ws = []
    for s in sectors:
        t, f, cad = load(dv, s)
        dA, _ = phase_dt(t, P, T0); dB, _ = phase_dt(t, P_B, T0_B)
        ok = (np.abs(dA) > T14) & (np.abs(dB) > T14_B_H / 24)
        dP, ep = phase_dt(t, P, T0 + phase * P)
        inn = ok & (np.abs(dP) < width_d / 2)
        base = ok & (np.abs(dP) >= width_d / 2) & (np.abs(dP) < 0.3)
        if inn.sum() < 3 or base.sum() < 10:
            continue
        sig = beta[str(s)]["sigma_white"]
        # per-epoch local baseline
        num = 0.0; nin = 0
        for e in np.unique(ep[inn]):
            bi = base & (ep == e); ii = inn & (ep == e)
            if bi.sum() < 5:
                continue
            num += np.sum(f[ii] - np.mean(f[bi])); nin += ii.sum()
        if nin == 0:
            continue
        vals.append(num / nin); ws.append(nin / sig ** 2)
    vals, ws = np.array(vals), np.array(ws)
    return float(np.sum(vals * ws) / np.sum(ws)), float(1 / np.sqrt(np.sum(ws)))


out3 = {}
for dv in ("mask05", "base"):
    v, e = box_depth(dv, 0.5, T14)
    rnd = []
    for ph in rng.uniform(0.1, 0.9, 600):
        if abs(ph - 0.5) < 0.04:
            continue
        rnd.append(box_depth(dv, ph, T14)[0])
    rnd = np.array(rnd)
    out3[dv] = dict(flux_ppm=v * 1e6, sig_formal_ppm=e * 1e6, sig_emp_std_ppm=float(np.std(rnd) * 1e6),
                    sig_emp_mad_ppm=float(1.4826 * np.median(np.abs(rnd - np.median(rnd))) * 1e6),
                    mean_random_ppm=float(np.mean(rnd) * 1e6), n_random=len(rnd))
# scan across phases (eccentric orbits)
scan = {}
for wf in (0.5, 1.0, 2.0):
    phs = np.arange(0.06, 0.94, 0.002)
    vals = np.array([box_depth("base", ph, wf * T14)[0] for ph in phs])
    sd = 1.4826 * np.median(np.abs(vals - np.median(vals)))
    j = np.argmin(vals)
    scan[str(wf)] = dict(min_flux_ppm=float(vals[j] * 1e6), phase=float(phs[j]), sig_emp_ppm=float(sd * 1e6),
                         snr_min=float(vals[j] / sd), max_flux_ppm=float(vals.max() * 1e6))
    np.savetxt(os.path.join(OUT, f"secondary_scan_w{wf}.txt"), np.c_[phs, vals * 1e6], header="phase flux_ppm")
out3["scan"] = scan
R["secondary"] = out3
print(json.dumps(out3, indent=1))


# ---------------------------------------------------------------- 4) V-shape
def trap_unit(dt, T14d, T12d):
    x = np.abs(dt)
    y = np.zeros_like(x)
    full = x <= T14d / 2 - T12d
    y[full] = 1.0
    ing = (x > T14d / 2 - T12d) & (x < T14d / 2)
    y[ing] = (T14d / 2 - x[ing]) / T12d
    return y


DATA = []
for s in SECT:
    t, f, cad, _ = select("base", s)
    ss = SUPERSAMPLE[s]
    off = (np.arange(ss) + 0.5) / ss - 0.5
    tsub = (t[:, None] + off[None, :] * cad / 86400.0).ravel()
    ep = np.round((t - T0) / P).astype(int)
    DATA.append(dict(s=s, t=t, f=f, tsub=tsub, ss=ss, ep=ep))


def resid(p, mode, scale):
    depth, T14d, T12d = p[0], p[1], (p[1] / 2 if mode == "tri" else p[2])
    offs = p[2 if mode == "tri" else 3:]
    out = []
    for i, d in enumerate(DATA):
        dt, _ = phase_dt(d["tsub"], P, T0)
        m = 1 - depth * trap_unit(dt, T14d, max(T12d, 1e-5)).reshape(-1, d["ss"]).mean(axis=1)
        sig = beta[str(d["s"])]["sigma_white"] * (beta[str(d["s"])]["beta"] if scale else 1.0)
        out.append((d["f"] - m - offs[i]) / sig)
    return np.concatenate(out)


out4 = {}
for scale in (False, True):
    best = {}
    for mode in ("trap", "tri"):
        bestc = None
        for T14g in (0.07, 0.078, 0.09):
            for T12g in (0.003, 0.006, 0.012, 0.02):
                p0 = [5e-4, T14g] + ([] if mode == "tri" else [T12g]) + [0, 0, 0]
                lb = [0, 0.02] + ([] if mode == "tri" else [1e-4]) + [-1e-3] * 3
                ub = [5e-3, 0.2] + ([] if mode == "tri" else [0.1]) + [1e-3] * 3
                r = least_squares(resid, p0, args=(mode, scale), bounds=(lb, ub))
                c2 = float(np.sum(r.fun ** 2))
                if mode == "trap" and r.x[2] > r.x[1] / 2:
                    continue
                if bestc is None or c2 < bestc[0]:
                    bestc = (c2, r.x.tolist())
        best[mode] = dict(chi2=bestc[0], params=bestc[1])
    out4["beta_scaled" if scale else "white"] = dict(trap=best["trap"], tri=best["tri"],
                                                    dchi2=best["tri"]["chi2"] - best["trap"]["chi2"])
R["vshape"] = out4
print(json.dumps(out4, indent=1))

# ---------------------------------------------------------------- 5) sqrt(depth) sanity
mu = np.sqrt(1 - b ** 2)
Imu = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2
Iavg = 1 - u1 / 3 - u2 / 6
fac = Imu / Iavg
d_all = out2["all"]["depth_ppm"] * 1e-6  # mid-transit depth of the scaled template
R["sqrt_depth"] = dict(ld_factor=float(fac), k_from_depth=float(np.sqrt(d_all / fac)), k_naive=float(np.sqrt(d_all)),
                       Rp_from_depth=float(np.sqrt(d_all / fac) * R_STAR * RSUN_OVER_REARTH),
                       Rp_naive_500ppm=float(np.sqrt(500e-6) * R_STAR * RSUN_OVER_REARTH))
print(R["sqrt_depth"])
json.dump(R, open(os.path.join(OUT, "tests_results.json"), "w"), indent=1)
