"""R-S6 (a): is the S32 <-> S98 cycle count unique without Sector 6?  Light-curve level test.

For each cycle-count alias j (P_j = (T98 - T32)/(N0 + j) from s1), search T0 (+-40 min) and P (+-4e-5 d) around the alias
and record the best detection statistic Delta-chi2 = chi2(no transit) - chi2(transit) of the fixed-shape, cadence-integrated
template on S32 (SPOC 2-min PDCSAP) + S98 (QLP SYS_RM_FLUX, scattered light removed), B cadences removed, one linear
baseline per +-0.3 d window (profiled), common depth > 0, per-sector sigma fixed (point-to-point MAD / sqrt 2).
"""
import json
import numpy as np
from rs6_common import *
import batman

S1 = json.load(open(os.path.join(OUT, "s1_timing_R60.json")))
CC = S1["cycle_count"]
U = claret_u()


def tpl_binned():
    tau = np.arange(-0.40, 0.40 + 1e-9, 2.0 / 86400.0)
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp = 0.0, P_A, 0.0212
    pm.a = ars_from_rho(2.51, P_A); pm.inc = np.degrees(np.arccos(0.35 / pm.a)); pm.ecc, pm.w = 0.0, 90.0
    pm.u, pm.limb_dark = list(U), "quadratic"
    fl = batman.TransitModel(pm, tau).light_curve(pm)
    s = (1 - fl) / (1 - fl[np.argmin(np.abs(tau))])
    out = {}
    for sec in (32, 98):
        n = int(round(CAD_S[sec] / 2.0))
        out[sec] = np.convolve(s, np.ones(n) / n, mode="same")
    return tau, out


TAU, TPL = tpl_binned()
LC = {}
for sec in (32, 98):
    t, f, meta = load(sec)
    dB, _ = phase_dt(t, P_B, T0_B)
    k = np.abs(dB) > T14_B
    t, f = t[k], f[k]
    d = np.diff(f)
    sig = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0)   # point-to-point scatter (robust to trends)
    LC[sec] = (t, f, sig)


def stat(T0, P):
    A = 0.0; B = 0.0
    for sec in (32, 98):
        t, f, sig = LC[sec]
        x, ep = phase_dt(t, P, T0)
        m = np.abs(x) < 0.30
        x, ep, ff = x[m], ep[m], f[m]
        s = np.interp(x, TAU, TPL[sec], left=0, right=0)
        _, wid = np.unique(ep, return_inverse=True)
        nw = wid.max() + 1
        S0 = np.bincount(wid, minlength=nw).astype(float); S1_ = np.bincount(wid, x, nw); S2 = np.bincount(wid, x * x, nw)
        det = S0 * S2 - S1_ ** 2
        det[det == 0] = np.inf

        def proj(y):
            Sy = np.bincount(wid, y, nw); Sxy = np.bincount(wid, y * x, nw)
            a0 = (S2 * Sy - S1_ * Sxy) / det; a1 = (S0 * Sxy - S1_ * Sy) / det
            return y - a0[wid] - a1[wid] * x
        sp, fp = proj(s), proj(ff)
        A += sp @ sp / sig ** 2; B += sp @ fp / sig ** 2
    # depth D = -B/A with model f = base - D s  -> use B' = sum s_p (-f_p)
    return A, B


def dchi2(T0, P):
    A, B = stat(T0, P)
    # model f = base - D s: D = -(s.f)/(s.s); improvement = (s.f)^2/(s.s) if D > 0
    D = -B / A
    return (B * B / A) if D > 0 else 0.0, D


def main():
    res = []
    T98, E98 = CC["T98_BTJD"], CC["E98"]
    for a in CC["aliases"]:
        Pj = a["P"]
        T0j = T98 - E98 * Pj
        best = (-1, None, None, None)
        for dT in np.arange(-40, 40.01, 1.0) / 1440.0:
            for dP in np.arange(-4e-5, 4.0001e-5, 1e-6):
                v, D = dchi2(T0j + dT, Pj + dP)
                if v > best[0]:
                    best = (v, T0j + dT, Pj + dP, D)
        res.append(dict(j=a["j"], N=a["N"], P_alias=Pj, best_dchi2=float(best[0]), best_T0=float(best[1]), best_P=float(best[2]),
                        depth_ppm=float(best[3] * 1e6)))
        print(res[-1], flush=True)
    b0 = [r for r in res if r["j"] == 0][0]["best_dchi2"]
    for r in res:
        r["dchi2_below_true"] = b0 - r["best_dchi2"]
    json.dump(dict(method=__doc__, sigma_ppm={str(k): LC[k][2] * 1e6 for k in LC}, aliases=res), open(os.path.join(OUT, "s1b_alias.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
