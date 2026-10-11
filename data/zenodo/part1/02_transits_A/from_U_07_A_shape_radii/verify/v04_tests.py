"""Verifier shape/EB tests on raw normalized flux with LOCAL quadratic baselines fitted simultaneously per window
(same contiguous segment only; windows +-0.45 d). Template shape fixed at given parameters.
 - per-sector / odd / even / all depths (mid-transit depth of the template) with formal and EMPIRICAL sigma
   (same estimator at randomly shifted ephemerides that avoid real A and B transits)
 - phase-0.5 box flux (width T14) relative to local baseline, empirical sigma from random phases
 - V-shape: supersampled trapezoid vs triangle on the locally detrended light curve (from v02_prep), white and
   beta-rescaled errors.
usage: v04_tests.py K B RHO U1 U2 COL98 LABEL"""
import sys, json
import numpy as np
from scipy.optimize import least_squares
from vcommon import *
K, B, RHO, U1, U2 = map(float, sys.argv[1:6]); COL98 = sys.argv[6]; LABEL = sys.argv[7]
A = float(ars(RHO, PA)); T14 = float(t14(PA, K, B, A))
rng = np.random.default_rng(99)
LC = {s: load(s, COL98 if s == 98 else None) for s in (6, 32, 98)}
import batman
pm = batman.TransitParams(); pm.t0 = 0; pm.per = PA; pm.rp = K; pm.a = A; pm.inc = np.degrees(np.arccos(B / A))
pm.ecc = 0; pm.w = 90; pm.u = [U1, U2]; pm.limb_dark = "quadratic"
DMID = 1 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0]


def windows(s, shift=0.0, parity=None, excl=()):
    """yield (t, f, x, epoch) for each window centred on T0A+shift+e*PA with >=50% coverage; A/B real transits removed
    when shift != 0."""
    t, f, cad, ss = LC[s]
    dA = t - (T0A + np.round((t - T0A) / PA) * PA); dB = t - (T0B + np.round((t - T0B) / PB) * PB)
    c0 = T0A + shift
    for e in range(int(np.floor((t.min() - c0) / PA)), int(np.ceil((t.max() - c0) / PA)) + 1):
        if parity is not None and abs(e) % 2 != parity:
            continue
        if (s, e) in excl:
            continue
        tc = c0 + e * PA
        w = (np.abs(t - tc) < 0.45) & segment_of(t, tc) & (np.abs(dB) > T14B)
        if shift != 0.0:
            w &= np.abs(dA) > 0.75 * T14A
        inn = w & (np.abs(t - tc) < T14 / 2)
        if inn.sum() * cad / 86400 < 0.5 * T14:
            continue
        oot = w & (np.abs(t - tc) > 0.75 * T14A)
        if oot.sum() < 6:
            continue
        yield t[w], f[w], t[w] - tc, e, cad, ss


def depth(sets, shift=0.0, parity=None, excl=(), per_window=False):
    num = den = 0.0; rows = []
    for s in sets:
        W = list(windows(s, shift, parity, excl))
        if not W:
            continue
        parts = []
        for t, f, x, e, cad, ss in W:
            tm = Model(t, cad, ss)
            tmp = (1 - tm(T0A + shift + e * PA - 0.0, PA, K, B, A, U1, U2)) / DMID
            X = np.c_[np.ones_like(x), x, x * x]
            H = X @ np.linalg.pinv(X)
            parts.append((tmp - H @ tmp, (1 - f) - H @ (1 - f), e))
        sig = np.std(np.concatenate([p[1] - (np.sum(p[0] * p[1]) / max(np.sum(p[0] ** 2), 1e-30)) * p[0] for p in parts]))
        for xt, yy, e in parts:
            num += np.sum(xt * yy) / sig ** 2; den += np.sum(xt * xt) / sig ** 2
            if per_window and np.sum(xt * xt) > 0:
                rows.append(dict(sector=s, epoch=int(e), depth=float(np.sum(xt * yy) / np.sum(xt * xt) * 1e6),
                                 sig=float(sig / np.sqrt(np.sum(xt * xt)) * 1e6)))
    if den == 0:
        return np.nan, np.nan, rows
    return num / den * 1e6, 1 / np.sqrt(den) * 1e6, rows


def emp(fn, n=250):
    v = []
    for sh in rng.uniform(0.1, 0.9, n) * PA:
        if abs(sh / PA - 0.5) < 0.05:
            continue
        x = fn(sh)
        if np.isfinite(x):
            v.append(x)
    v = np.array(v)
    return float(np.std(v)), float(1.4826 * np.median(np.abs(v - np.median(v)))), float(np.mean(v)), len(v)


R = dict(label=LABEL, template=dict(k=K, b=B, rho=RHO, u1=U1, u2=U2, ars=A, T14_h=T14 * 24, depth_mid_ppm=DMID * 1e6))
EX = ((98, -1),)
for key, sets, par, ex in [("S6", (6,), None, ()), ("S32", (32,), None, ()), ("S98", (98,), None, ()), ("S98_excl-1", (98,), None, EX),
                           ("all", (6, 32, 98), None, ()), ("all_excl98-1", (6, 32, 98), None, EX),
                           ("odd", (6, 32, 98), 1, ()), ("even", (6, 32, 98), 0, ()),
                           ("odd_excl98-1", (6, 32, 98), 1, EX), ("even_excl98-1", (6, 32, 98), 0, EX)]:
    d, e, _ = depth(sets, 0.0, par, ex)
    es = emp(lambda sh: depth(sets, sh, par, ex)[0])
    R[key] = dict(depth_ppm=d, sig_formal_ppm=e, sig_emp_std_ppm=es[0], sig_emp_mad_ppm=es[1], emp_mean_ppm=es[2], n_emp=es[3])
    print(key, {k: round(v, 1) for k, v in R[key].items()}, flush=True)
_, _, rows = depth((6, 32, 98), per_window=True)
R["per_transit"] = rows
R["n_transits_measured"] = len(rows)
for a_, b_ in (("odd", "even"), ("odd_excl98-1", "even_excl98-1")):
    so = max(R[a_]["sig_formal_ppm"], R[a_]["sig_emp_std_ppm"]); se = max(R[b_]["sig_formal_ppm"], R[b_]["sig_emp_std_ppm"])
    R[f"{a_}_minus_{b_}_nsig"] = (R[a_]["depth_ppm"] - R[b_]["depth_ppm"]) / np.hypot(so, se)
ds = np.array([R[k]["depth_ppm"] for k in ("S6", "S32", "S98")]); ss = np.array([max(R[k]["sig_formal_ppm"], R[k]["sig_emp_std_ppm"]) for k in ("S6", "S32", "S98")])
w = 1 / ss ** 2; mu = np.sum(w * ds) / np.sum(w); R["sector_chi2"] = dict(mean=mu, chi2=float(np.sum(w * (ds - mu) ** 2)), dof=2)
print("sector chi2", R["sector_chi2"], "odd-even nsig", R["odd_minus_even_nsig"], R["odd_excl98-1_minus_even_excl98-1_nsig"], flush=True)


# ------------------------------------------------ phase 0.5 box (flux relative to local quadratic baseline)
TWO_SIDED = False   # v07 sets True (baseline >= 0.15 d on both sides)


def boxflux(phase, width):
    num = den = 0.0
    for s in (6, 32, 98):
        t, f, cad, ss = LC[s]
        dA = t - (T0A + np.round((t - T0A) / PA) * PA); dB = t - (T0B + np.round((t - T0B) / PB) * PB)
        ok = (np.abs(dA) > 0.75 * T14A) & (np.abs(dB) > T14B)
        c0 = T0A + phase * PA
        vals = []
        for e in range(int(np.floor((t.min() - c0) / PA)), int(np.ceil((t.max() - c0) / PA)) + 1):
            tc = c0 + e * PA
            w = ok & (np.abs(t - tc) < 0.45) & segment_of(t, tc)
            inn = w & (np.abs(t - tc) < width / 2)
            oot = w & (np.abs(t - tc) > 0.75 * T14A)
            if inn.sum() * cad / 86400 < 0.5 * width or oot.sum() < 6:
                continue
            if TWO_SIDED and not ((t[oot] < tc - 0.15).any() and (t[oot] > tc + 0.15).any()):
                continue
            x = t - tc
            c = np.polyfit(x[oot], f[oot], 2)
            r = f[w] / np.polyval(c, x[w]) - 1
            vals.append((r[np.abs(x[w]) < width / 2], np.std(r[np.abs(x[w]) >= 0.75 * T14A])))
        if not vals:
            continue
        sig = np.median([v[1] for v in vals])
        for r, _ in vals:
            num += np.sum(r) / sig ** 2; den += len(r) / sig ** 2
    return (num / den * 1e6, 1 / np.sqrt(den) * 1e6) if den > 0 else (np.nan, np.nan)


v05, e05 = boxflux(0.5, T14)
rnd = np.array([boxflux(ph, T14)[0] for ph in rng.uniform(0.1, 0.9, 300) if abs(ph - 0.5) > 0.04])
rnd = rnd[np.isfinite(rnd)]
R["phase0p5"] = dict(flux_ppm=v05, sig_formal_ppm=e05, sig_emp_std_ppm=float(np.std(rnd)), emp_mean_ppm=float(np.mean(rnd)),
                     n_random=len(rnd), note="flux relative to baseline: negative = dip")
print("phase0.5", R["phase0p5"], flush=True)
json.dump(R, open(os.path.join(OUT, f"tests_{LABEL}.json"), "w"), indent=1, default=float)
