"""Step 10: odd/even, secondary eclipse, phase scan, harmonic (P/2, 2P, P/3, 3P) and V/U-shape tests, both signals.

Errors: every box depth gets two errors, (i) white = robust scatter / sqrt(N_in), (ii) EMPIRICAL = standard deviation
of box depths of the same width placed at random phases that avoid both signals' transits (this absorbs red noise and
detrending residuals). Significances quoted with the empirical error. The phase scan reports the largest |z| over all
phases and the fraction of random phases with |z| at least as large (look-elsewhere).
Odd/even uses the per-event depths of s09 (red-noise-inflated) plus a direct box comparison with empirical errors.
Shape: trapezoid fit (total duration T14, ingress fraction x = T12/T14 in [0, 0.5], depth) on the combined fold;
Delta chi2 between the best U-shape (x free) and a pure V (x = 0.5) with chi2 rescaled to chi2_red = 1.
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
from audit_common import AUD, save_json  # noqa

EPH = {k: dict(v) for k, v in lcio.EPH.items()}
ephf = os.path.join(AUD, "out", "ephem_current.json")
if os.path.exists(ephf):
    EPH.update(json.load(open(ephf)))
rng = np.random.default_rng(99)


def load(setname):
    if setname == "pipe":
        lcs = [lcio.load_product("TESSSPOC_S6_PDC"), lcio.load_product("SPOC2min_S32_PDC"), lcio.load_product("QLP_S98_SYS", "frozen")]
    else:
        lcs = [lcio.load_own("tesscut_s6", 2.5), lcio.load_own("tesscut_s32", 2.0), lcio.load_own("tesscut_s98", 2.0, drop_stray=True)]
    return [lcio.detrend(lc, window=0.75, eph=EPH) for lc in lcs]


def box(t, f, w, P, T0, phase, dur):
    ph = ((t - T0) / P - phase + 0.5) % 1 - 0.5
    m = np.abs(ph * P) < dur / 2
    if m.sum() < 3:
        return np.nan, 0
    return 1 - np.sum(w[m] * f[m]) / np.sum(w[m]), int(m.sum())


def empirical_sigma(t, f, w, P, T0, dur, ntrial=3000, avoid=None):
    """box depths at random phases avoiding both signals' transits (+-1.5 T14) -> robust std."""
    vals = []
    for _ in range(ntrial):
        ph = rng.uniform(0.1, 0.9)
        d, n = box(t, f, w, P, T0, ph, dur)
        if n < 3:
            continue
        # reject if the box overlaps either signal's transits
        ph_all = ((t - T0) / P - ph + 0.5) % 1 - 0.5
        m = np.abs(ph_all * P) < dur / 2
        if avoid is not None and np.any(avoid[m]):
            continue
        vals.append(d)
    vals = np.array(vals)
    return float(1.4826 * np.median(np.abs(vals - np.median(vals)))), float(np.std(vals)), len(vals), vals


def trapezoid(tr, depth, T14, x):
    """symmetric trapezoid: total duration T14, ingress duration x*T14 (x in (0, 0.5])."""
    a = np.abs(tr)
    ti = x * T14
    out = np.zeros_like(tr)
    flat = a <= T14 / 2 - ti
    ing = (a > T14 / 2 - ti) & (a < T14 / 2)
    out[flat] = depth
    out[ing] = depth * (T14 / 2 - a[ing]) / ti
    return out


def run(sig, setname):
    lcs = load(setname)
    e = EPH[sig]
    other = "B" if sig == "A" else "A"
    P, T0, T14 = e["P"], e["T0"], e["T14"]
    res = {}
    # combine all sectors, weights 1/sig^2 per data set
    t = np.concatenate([lc["t"] for lc in lcs])
    f = np.concatenate([lc["fd"] for lc in lcs])
    w = np.concatenate([np.full(len(lc["t"]), 1 / lc["sig"] ** 2) for lc in lcs])
    keep = ~lcio.in_transit(t, other, 1.0, EPH)
    t, f, w = t[keep], f[keep], w[keep]
    near_sig = lcio.in_transit(t, sig, 1.5, EPH)
    s_rob, s_std, ntr, vals = empirical_sigma(t, f, w, P, T0, T14, avoid=near_sig)
    d0, n0 = box(t, f, w, P, T0, 0.0, T14)
    white = 1 / np.sqrt(np.sum(w[np.abs((((t - T0) / P + 0.5) % 1 - 0.5) * P) < T14 / 2]))
    res["primary"] = dict(depth_ppm=d0 * 1e6, n_in=n0, err_white_ppm=white * 1e6, err_emp_robust_ppm=s_rob * 1e6,
                          err_emp_std_ppm=s_std * 1e6, n_random=ntr, snr_emp=d0 / s_rob, snr_white=d0 / white,
                          frac_random_ge_primary=float(np.mean(vals >= d0)))
    # secondary at phase 0.5 (circular) and a sweep 0.35-0.65 for eccentric cases
    d5, n5 = box(t, f, w, P, T0, 0.5, T14)
    res["secondary_0p5"] = dict(depth_ppm=d5 * 1e6, n_in=n5, err_emp_ppm=s_rob * 1e6, z=d5 / s_rob,
                                upper_3sigma_ppm=(max(d5, 0) + 3 * s_rob) * 1e6)
    # full phase scan
    phases = np.arange(0.05, 0.9501, 0.005)
    sc = []
    for ph in phases:
        d, n = box(t, f, w, P, T0, ph, T14)
        sc.append((ph, d * 1e6, n))
    sc = np.array(sc)
    z = sc[:, 1] / (s_rob * 1e6)
    i = int(np.nanargmax(np.abs(z)))
    res["phase_scan"] = dict(max_abs_z=float(np.abs(z[i])), phase_at_max=float(sc[i, 0]), depth_at_max_ppm=float(sc[i, 1]),
                             frac_random_abs_ge=float(np.mean(np.abs(vals) >= np.abs(sc[i, 1]) * 1e-6)),
                             n_independent_phases=int(round(0.9 * P / T14)),
                             scan=[[float(a), float(b)] for a, b, _ in sc])
    # P/3 check (dips at phases 1/3, 2/3)
    res["P_over_3_check"] = {str(ph): dict(depth_ppm=box(t, f, w, P, T0, ph, T14)[0] * 1e6,
                                           z=box(t, f, w, P, T0, ph, T14)[0] / s_rob) for ph in (1 / 3, 2 / 3)}
    # odd/even by direct box (empirical error from random phases at 2P)
    n_ep = np.round((t - T0) / P).astype(int)
    oe = {}
    for par in (0, 1):
        m = (n_ep % 2) == par
        d, n = box(t[m], f[m], w[m], P, T0, 0.0, T14)
        s_r, _, _, _ = empirical_sigma(t[m], f[m], w[m], P, T0, T14, avoid=near_sig[m])
        oe["even" if par == 0 else "odd"] = dict(depth_ppm=d * 1e6, err_emp_ppm=s_r * 1e6, n_in=n)
    dz = (oe["odd"]["depth_ppm"] - oe["even"]["depth_ppm"]) / np.hypot(oe["odd"]["err_emp_ppm"], oe["even"]["err_emp_ppm"])
    res["odd_even_box"] = dict(**oe, diff_sigma=float(dz))
    # n mod 3 grouping (3P test)
    g3 = {}
    for r3 in (0, 1, 2):
        m = (n_ep % 3) == r3
        d, n = box(t[m], f[m], w[m], P, T0, 0.0, T14)
        s_r, _, _, _ = empirical_sigma(t[m], f[m], w[m], P, T0, T14, avoid=near_sig[m])
        g3[str(r3)] = dict(depth_ppm=d * 1e6, err_emp_ppm=s_r * 1e6, n_in=n)
    ds = np.array([g3[k]["depth_ppm"] for k in g3])
    es = np.array([g3[k]["err_emp_ppm"] for k in g3])
    okg = np.isfinite(ds) & np.isfinite(es) & (es > 0)
    mu = np.sum(ds[okg] / es[okg] ** 2) / np.sum(1 / es[okg] ** 2)
    res["mod3"] = dict(groups=g3, chi2_const=float(np.sum(((ds[okg] - mu) / es[okg]) ** 2)), dof=int(okg.sum() - 1))
    # shape: trapezoid on the fold within +-3 T14, errors rescaled
    ph = (((t - T0) / P + 0.5) % 1 - 0.5) * P
    m = np.abs(ph) < 3 * T14
    tr, ff, ww = ph[m], f[m], w[m]

    def chi2(p, xfix=None):
        dep, dur, x, c = p if xfix is None else (p[0], p[1], xfix, p[2])
        if not (0 < dep < 0.01 and 0.02 < dur < 0.3 and 0.001 <= x <= 0.5):
            return 1e30
        mod = c - trapezoid(tr, dep, dur, x)
        return np.sum(ww * (ff - mod) ** 2)

    best_u = min((minimize(chi2, [5e-4 if sig == "A" else 8e-4, T14, x0, 1.0], method="Nelder-Mead",
                           options=dict(xatol=1e-7, fatol=1e-6, maxiter=20000)) for x0 in (0.05, 0.15, 0.3)), key=lambda r: r.fun)
    best_v = minimize(lambda p: chi2(p, 0.5), [5e-4 if sig == "A" else 8e-4, T14, 1.0], method="Nelder-Mead",
                      options=dict(xatol=1e-7, fatol=1e-6, maxiter=20000))
    dof = m.sum() - 4
    scale = best_u.fun / dof
    # profile over x
    xs = np.linspace(0.02, 0.5, 25)
    prof = [minimize(lambda p: chi2(p, xx), [best_u.x[0], best_u.x[1], best_u.x[3]], method="Nelder-Mead",
                     options=dict(xatol=1e-7, fatol=1e-6, maxiter=20000)).fun for xx in xs]
    prof = (np.array(prof) - best_u.fun) / scale
    res["shape"] = dict(U=dict(depth_ppm=best_u.x[0] * 1e6, T14_h=best_u.x[1] * 24, ingress_frac=best_u.x[2]),
                        V=dict(depth_ppm=best_v.x[0] * 1e6, T14_h=best_v.x[1] * 24),
                        dchi2_V_minus_U=float((best_v.fun - best_u.fun) / scale), chi2_scale=float(scale),
                        x_profile=[[float(a), float(b)] for a, b in zip(xs, prof)],
                        x_1sigma_upper=float(xs[np.where(prof <= 1)[0][-1]]) if np.any(prof <= 1) else None,
                        x_2sigma_upper=float(xs[np.where(prof <= 4)[0][-1]]) if np.any(prof <= 4) else None)
    return res


out = {"ephemeris": EPH}
for sig in ("A", "B"):
    for setname in ("pipe", "own"):
        r = run(sig, setname)
        out[f"{sig}_{setname}"] = r
        p = r["primary"]
        print(f"{sig} {setname}: primary {p['depth_ppm']:.0f} ppm, err white {p['err_white_ppm']:.0f} emp {p['err_emp_robust_ppm']:.0f} "
              f"-> S/N emp {p['snr_emp']:.1f} (white {p['snr_white']:.1f}); random >= primary: {p['frac_random_ge_primary']:.4f}")
        s = r["secondary_0p5"]
        print(f"    secondary(0.5): {s['depth_ppm']:.0f} ± {s['err_emp_ppm']:.0f} ppm (z={s['z']:.2f}), 3sigma upper {s['upper_3sigma_ppm']:.0f}")
        ps = r["phase_scan"]
        print(f"    phase scan: max|z| {ps['max_abs_z']:.2f} at {ps['phase_at_max']:.3f} ({ps['depth_at_max_ppm']:.0f} ppm); "
              f"~{ps['n_independent_phases']} independent phases")
        print(f"    P/3 phases: {[(k, round(v['depth_ppm']), round(v['z'], 2)) for k, v in r['P_over_3_check'].items()]}")
        oe = r["odd_even_box"]
        print(f"    odd {oe['odd']['depth_ppm']:.0f}±{oe['odd']['err_emp_ppm']:.0f} even {oe['even']['depth_ppm']:.0f}±{oe['even']['err_emp_ppm']:.0f} "
              f"-> {oe['diff_sigma']:.2f} sigma")
        print(f"    mod3 chi2 {r['mod3']['chi2_const']:.2f}/{r['mod3']['dof']}: {[(k, v['depth_ppm']) for k, v in r['mod3']['groups'].items()]}")
        sh = r["shape"]
        print(f"    shape: U {sh['U']} ; V {sh['V']} ; dchi2(V-U) = {sh['dchi2_V_minus_U']:.2f}; ingress frac <= {sh['x_2sigma_upper']} (2 sigma)")
save_json(out, "out/s10_eb_tests.json")
