"""Per-sector, odd/even and all-transit depths with a detrending-independent method: raw normalized flux
(PDCSAP S6/S32, QLP SYS_RM_FLUX S98, QUALITY==0), fixed transit shape from a joint-fit posterior median, a per-window
linear baseline fitted simultaneously (linear least squares).  Empirical errors: identical fits with the whole
ephemeris shifted by random phase offsets (windows avoiding the real A and B transits).
usage: s06_depths_localbaseline.py [CHAIN_TAG]   (default mann_claret_jit_raw)
Depths are the mid-transit depth of the scaled template (unbinned model), in ppm."""
import sys, json
import numpy as np
import batman
from common import *

TAG = sys.argv[1] if len(sys.argv) > 1 else "mann_claret_jit_raw"
ch = np.load(os.path.join(OUT, f"chain_{TAG}.npz"))["chain"]
T0, P, k, b, lnrho, u1, u2 = np.median(ch, axis=0)[:7]
ars = float(ars_from_rho(np.exp(lnrho), P)); T14 = float(t14_hours(P, k, b, ars)) / 24
pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = k; pm.a = ars; pm.inc = np.degrees(np.arccos(b / ars))
pm.ecc = 0; pm.w = 90; pm.u = [u1, u2]; pm.limb_dark = "quadratic"
DMID = 1 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0]
rng = np.random.default_rng(11)
SECT = (6, 32, 98)
RAW = {s: np.load(os.path.join(OUT, f"lc_raw_S{s}.npz")) for s in SECT}


def sel(s, shift=0.0, parity=None):
    z = RAW[s]; t, f, cad = z["t"], z["f"], float(z["cad"])
    dA, _ = phase_dt(t, P, T0); dB, _ = phase_dt(t, P_B, T0_B); dS, eS = phase_dt(t, P, T0 + shift)
    m = (np.abs(dS) < 0.3) & (np.abs(dB) > T14_B_H / 24)
    if shift != 0.0:
        m &= np.abs(dA) > T14
    if parity is not None:
        m &= (np.abs(eS) % 2) == parity
    return t[m], f[m], cad, s


def fit(sets, shift=0.0):
    num = 0.0; den = 0.0
    for t, f, cad, s in sets:
        if len(t) == 0:
            continue
        tmp = (1 - BatmanSet(t, cad, SUPERSAMPLE[s]).flux(T0 + shift, P, k, b, ars, u1, u2)) / DMID
        ep = np.round((t - T0 - shift) / P).astype(int)
        # residual scatter per sector for weights: from a baseline-only fit
        resid_all = []; parts = []
        for e in np.unique(ep):
            mm = ep == e
            x = t[mm] - (T0 + shift + e * P)
            X = np.c_[np.ones(mm.sum()), x]
            H = X @ np.linalg.pinv(X)
            xt = tmp[mm] - H @ tmp[mm]; yy = (1 - f[mm]) - H @ (1 - f[mm])
            parts.append((xt, yy))
        sig = np.std(np.concatenate([p[1] for p in parts]))
        for xt, yy in parts:
            num += np.sum(xt * yy) / sig ** 2; den += np.sum(xt * xt) / sig ** 2
    return (num / den, 1 / np.sqrt(den)) if den > 0 else (np.nan, np.nan)


def emp(fn, n=300):
    v = []
    for sh in rng.uniform(0.12, 0.88, n) * P:
        if abs(sh - 0.5 * P) < 0.06 * P:
            continue
        x = fn(sh)
        if np.isfinite(x):
            v.append(x)
    return float(np.std(v))


R = dict(chain=TAG, depth_mid_template_ppm=DMID * 1e6, T14_h=T14 * 24)
for s in SECT:
    d, e = fit([sel(s)])
    R[f"S{s}"] = dict(depth=d * 1e6, sig_formal=e * 1e6, sig_emp=emp(lambda sh: fit([sel(s, sh)], sh)[0]) * 1e6)
d, e = fit([sel(s) for s in SECT])
R["all"] = dict(depth=d * 1e6, sig_formal=e * 1e6, sig_emp=emp(lambda sh: fit([sel(s, sh) for s in SECT], sh)[0]) * 1e6)
for par, lab in ((1, "odd"), (0, "even")):
    d, e = fit([sel(s, parity=par) for s in SECT])
    R[lab] = dict(depth=d * 1e6, sig_formal=e * 1e6,
                  sig_emp=emp(lambda sh: fit([sel(s, sh, par) for s in SECT], sh)[0]) * 1e6)
so = max(R["odd"]["sig_formal"], R["odd"]["sig_emp"]); se = max(R["even"]["sig_formal"], R["even"]["sig_emp"])
R["odd_minus_even_nsigma"] = (R["odd"]["depth"] - R["even"]["depth"]) / np.hypot(so, se)
ds = np.array([R[f"S{s}"]["depth"] for s in SECT]); ss = np.array([max(R[f"S{s}"]["sig_formal"], R[f"S{s}"]["sig_emp"]) for s in SECT])
w = 1 / ss ** 2; mu = np.sum(w * ds) / np.sum(w)
R["sector_chi2"] = dict(mean=float(mu), chi2=float(np.sum(w * (ds - mu) ** 2)), dof=2)
# S98 without the edge transit (epoch -1) and the 3-point partial (epoch 8)
t, f, cad, s = sel(98)
keep = ~np.isin(np.round((t - T0) / P).astype(int), [-1, 8])
d, e = fit([(t[keep], f[keep], cad, 98)])
R["S98_excl_-1_8"] = dict(depth=d * 1e6, sig_formal=e * 1e6)
print(json.dumps(R, indent=1, default=float))
json.dump(R, open(os.path.join(OUT, f"depths_localbaseline_{TAG}.json"), "w"), indent=1, default=float)
