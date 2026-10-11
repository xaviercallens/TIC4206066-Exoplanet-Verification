"""Verifier V-shape test: supersampled symmetric trapezoid (depth, T14, T12 free) vs triangle (T12 = T14/2) on the
locally detrended light curve (v02_prep output), |dt| < 0.25 d, per-sector offsets; chi2 with per-sector white sigma
(OOT std) and with sigma * beta (time-averaging beta from trapezoid residuals, 30-120 min bins within windows).
usage: v05_vshape.py LCTAG EXCL"""
import sys, json
import numpy as np
from scipy.optimize import least_squares
from vcommon import *
LCTAG, EXCL = sys.argv[1], sys.argv[2]
z = np.load(os.path.join(OUT, f"lc_{LCTAG}.npz"))
dt = z["t"] - (T0A + z["e"] * PA)
keep = np.abs(dt) < 0.25
if EXCL != "none":
    for it in EXCL.split(","):
        s_, e_ = map(int, it.split(":")); keep &= ~((z["s"] == s_) & (z["e"] == e_))
DS = []
for s in (6, 32, 98):
    m = keep & (z["s"] == s)
    cad, ss = F[s][2], F[s][3]
    off = ((np.arange(ss) + 0.5) / ss - 0.5) * cad / 86400
    DS.append(dict(s=s, dt=dt[m], f=z["f"][m], e=z["e"][m], dsub=(dt[m][:, None] + off[None, :]).ravel(), ss=ss, cad=cad,
                   sig=float(np.std(z["f"][m][np.abs(dt[m]) > 0.75 * T14A]))))


def trap(x, T14, T12, d):
    x = np.abs(x); y = np.zeros_like(x)
    y[x <= T14 / 2 - T12] = 1
    ing = (x > T14 / 2 - T12) & (x < T14 / 2)
    y[ing] = (T14 / 2 - x[ing]) / T12
    return d * y


def model(p, mode, i):
    d, T14 = p[0], p[1]; T12 = T14 / 2 if mode == "tri" else p[2]
    o = p[2 if mode == "tri" else 3:][i]
    D = DS[i]
    return 1 + o - trap(D["dsub"], T14, max(T12, 1e-5), d).reshape(-1, D["ss"]).mean(axis=1)


def resid(p, mode, beta):
    return np.concatenate([(D["f"] - model(p, mode, i)) / (D["sig"] * beta[i]) for i, D in enumerate(DS)])


def best(mode, beta):
    out = None
    for T14g in (0.07, 0.078, 0.09, 0.105):
        for T12g in (0.004, 0.01, 0.02, 0.035):
            p0 = [5e-4, T14g] + ([] if mode == "tri" else [min(T12g, T14g / 2 - 1e-3)]) + [0, 0, 0]
            lb = [0, 0.03] + ([] if mode == "tri" else [1e-4]) + [-1e-3] * 3
            ub = [3e-3, 0.2] + ([] if mode == "tri" else [0.1]) + [1e-3] * 3
            r = least_squares(resid, p0, args=(mode, beta), bounds=(lb, ub), x_scale="jac")
            if mode == "trap" and r.x[2] > r.x[1] / 2 + 1e-9:
                continue
            c2 = float(np.sum(r.fun ** 2))
            if out is None or c2 < out[0]:
                out = (c2, r.x)
    return out


one = np.ones(3)
tr = best("trap", one)
# beta from trapezoid residuals
betas = []
for i, D in enumerate(DS):
    res = D["f"] - model(tr[1], "trap", i)
    s1 = np.std(res); bl = []
    for bm in (30, 60, 90, 120):
        n = int(round(bm * 60 / D["cad"]))
        if n < 2:
            continue
        bb = []
        for e in np.unique(D["e"]):
            r = res[D["e"] == e]; M = len(r) // n
            if M >= 1:
                bb += list(r[:M * n].reshape(M, n).mean(axis=1))
        bb = np.array(bb); M = len(bb)
        if M >= 8:
            bl.append(np.std(bb) / (s1 / np.sqrt(n) * np.sqrt(M / (M - 1.0))))
    betas.append(max(1.0, float(np.median(bl))) if bl else 1.0)
betas = np.array(betas)
R = dict(lc=LCTAG, excl=EXCL, sig_white_ppm=[D["sig"] * 1e6 for D in DS], beta_values=betas.tolist(), npts=int(sum(len(D["f"]) for D in DS)))
for lab, bb in (("white", one), ("beta", betas)):
    t_ = best("trap", bb); v_ = best("tri", bb)
    R[lab] = dict(chi2_trap=t_[0], chi2_tri=v_[0], dchi2=v_[0] - t_[0], trap=dict(depth_ppm=t_[1][0] * 1e6, T14_h=t_[1][1] * 24, T12_h=t_[1][2] * 24),
                  tri=dict(depth_ppm=v_[1][0] * 1e6, T14_h=v_[1][1] * 24))
print(json.dumps(R, indent=1))
json.dump(R, open(os.path.join(OUT, f"vshape_{LCTAG}_{EXCL.replace(':', 'e').replace(',', '+')}.json"), "w"), indent=1)
