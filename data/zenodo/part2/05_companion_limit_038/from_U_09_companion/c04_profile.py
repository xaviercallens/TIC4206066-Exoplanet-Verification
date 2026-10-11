"""Audit 9, step 6/8: can a planet transiting an unresolved bound companion reproduce the transit shape?
For each companion (host) mass M_c: R_c (Mann 2015 via Mann 2019 M_Ks), rho_c, TESS flux fraction f_c
(Mamajek sequence + Stassun 2019 Eq. 1, relative to the observed total T=9.68 at d = 1000/28.075 pc), Claret LD.
 (a) circular orbit: rho fixed to rho_c, dilution fixed to f_c; free k, b  -> chi2_circ(M_c)
 (b) eccentric orbit: rho_app free (photo-eccentric), dilution f_c -> chi2_ecc(M_c) and rho_app upper limit
     -> minimum e at the most favourable omega (transit at apoapsis): e_min = (1-x^2)/(1+x^2), x=(rho_max/rho_c)^(1/3).
 Reference: target as host (f=1), rho free and rho fixed to the Mann value.
Usage: c04_profile.py A|B [npz] [ldfree]"""
import sys, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
from c03_transit_model import Folded, aR_from_rho, claret_ld, t14_hours
from c00_relations import companion_props

HERE = Path(__file__).resolve().parent
SIG = sys.argv[1] if len(sys.argv) > 1 else "A"
NPZ = sys.argv[2] if len(sys.argv) > 2 else "lc_windows.npz"
LDFREE = len(sys.argv) > 3 and sys.argv[3] == "ldfree"
PLX = 28.075018; T_OBS = 9.6817          # Gaia DR3 parallax; TESS mag from the S98 QLP header (TIC)
DM = 5 * np.log10(1000 / PLX / 10); MT_TOT = T_OBS - DM
fo = Folded(SIG, NPZ)
P = fo.P
u_t = claret_ld(4012, 4.63)[:2]
print("target LD", u_t, "M_T,tot", MT_TOT)

def q_to_u(q1, q2):
    s = np.sqrt(q1); return 2 * s * q2, s * (1 - 2 * q2)

KMAX = 1.0
def best_fit(aR=None, fh=1.0, u=u_t, starts_b=(0.1, 0.4, 0.7, 0.9), k0=None, ldfree=False):
    """minimise chi2 over k, b (and log rho if aR is None; q1,q2 if ldfree). Returns (chi2, dict)."""
    kk0 = k0 if k0 is not None else np.sqrt(5e-4 / fh) if SIG == "A" else np.sqrt(1e-3 / fh)
    best = (np.inf, None)
    for b0 in starts_b:
        x0 = [np.log(kk0), b0]
        if aR is None: x0.append(np.log(2.5))
        if ldfree: x0 += [0.36, 0.3]
        def fun(x):
            k = np.exp(x[0]); b = x[1]
            if b < 0 or b > 1 + k or k > KMAX: return 1e30
            j = 2
            if aR is None:
                a_ = aR_from_rho(np.exp(x[j]), P); j += 1
            else:
                a_ = aR
            if ldfree:
                q1, q2 = x[j], x[j + 1]
                if not (0 < q1 < 1 and 0 < q2 < 1): return 1e30
                uu = q_to_u(q1, q2)
            else:
                uu = u
            return fo.chi2(k, a_, b, uu[0], uu[1], fh)
        r = minimize(fun, x0, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-4, maxiter=4000, maxfev=8000))
        if r.fun < best[0]:
            k = np.exp(r.x[0]); b = r.x[1]
            d = dict(k=k, b=b)
            j = 2
            if aR is None:
                d["rho"] = float(np.exp(r.x[j])); d["aR"] = aR_from_rho(d["rho"], P); j += 1
            else:
                d["aR"] = aR
            if ldfree:
                d["u"] = q_to_u(r.x[j], r.x[j + 1])
            d["T14_h"] = float(t14_hours(P, d["aR"], k, b))
            d["depth_ppm_obs"] = float(fh * k ** 2 * 1e6)
            best = (r.fun, d)
    return best

t0 = time.time()
res = {"signal": SIG, "npz": NPZ, "ldfree": LDFREE, "n_points": fo.n, "sigma_ppm": {s: fo.s[s] * 1e6 for s in fo.sectors}}
c_free, d_free = best_fit(None, ldfree=LDFREE)
res["target_free_rho"] = dict(chi2=c_free, **{k: (float(v) if np.ndim(v) == 0 else list(map(float, v))) for k, v in d_free.items()})
print("target free:", c_free, d_free, time.time() - t0)
aR_m = aR_from_rho(2.51, P)
c_mann, d_mann = best_fit(aR_m, ldfree=LDFREE)
res["target_mann_rho"] = dict(chi2=c_mann, dchi2=c_mann - c_free, **{k: (float(v) if np.ndim(v) == 0 else list(map(float, v))) for k, v in d_mann.items()})
print("target Mann rho:", c_mann - c_free, d_mann)
# null (no transit)
c_null = sum(np.sum(((fo.f[s] - 1) / fo.s[s]) ** 2) for s in fo.sectors)
res["chi2_null"] = c_null; res["dchi2_transit"] = c_null - c_free
print("transit detection dchi2", c_null - c_free)

rows = []

for M in np.round(np.arange(0.08, 0.6061, 0.02), 3).tolist() + [0.38, 0.54, 0.606]:
    p = companion_props(M)
    MT1 = float(companion_props(0.606)["MT"]); fc = 1 / (1 + 10 ** (-0.4 * (MT1 - float(p["MT"]))))  # F_c/(F_1+F_c), primary = 0.606 Msun sequence point
    rho = float(M / p["R"] ** 3)
    logg = float(np.log10(6.674e-8 * M * 1.989e33 / (p["R"] * 6.957e10) ** 2))
    u = claret_ld(float(p["Teff"]), logg)[:2]
    aR = aR_from_rho(rho, P)
    KMAX = min(1.0, 20.0 / (float(p["R"]) * 109.08))   # R_p <= 20 R_E (TRICERATOPS sample_rp maximum)
    c_circ, d_circ = best_fit(aR, fh=fc, u=u, ldfree=LDFREE)
    # eccentric: rho_app free
    c_ecc, d_ecc = best_fit(None, fh=fc, u=u, ldfree=LDFREE)
    # rho_app upper limit (Delta chi2 = 4 w.r.t. c_ecc, profiling k, b)
    def prof(lr):
        c, _ = best_fit(aR_from_rho(np.exp(lr), P), fh=fc, u=u, starts_b=(0.1, 0.5, 0.85), ldfree=LDFREE)
        return c - c_ecc - 4.0
    lr0 = np.log(d_ecc["rho"])
    try:
        hi = lr0 + 0.1
        while prof(hi) < 0 and hi < lr0 + 6: hi += 0.4
        rho_max = float(np.exp(brentq(prof, lr0, hi, xtol=2e-3))) if hi < lr0 + 6 else np.inf
    except Exception as ex:
        rho_max = np.nan
    x = (rho_max / rho) ** (1 / 3)
    e_min = float(max(0.0, (1 - x ** 2) / (1 + x ** 2))) if np.isfinite(x) else np.nan
    row = dict(M=M, R=float(p["R"]), Teff=float(p["Teff"]), MT=float(p["MT"]), f_c=fc, rho=rho, aR=aR, u=u,
               chi2_circ=c_circ, dchi2_circ=c_circ - c_free, k_circ=d_circ["k"], b_circ=d_circ["b"],
               T14_circ_h=d_circ["T14_h"], Rp_circ_Re=d_circ["k"] * float(p["R"]) * 109.08,
               chi2_ecc=c_ecc, dchi2_ecc=c_ecc - c_free, rho_app_best=d_ecc["rho"], k_ecc=d_ecc["k"], b_ecc=d_ecc["b"],
               rho_app_max_2sig=rho_max, e_min_w270=e_min, Rp_ecc_Re=d_ecc["k"] * float(p["R"]) * 109.08)
    rows.append(row)
    print(f"M={M:.3f} R={row['R']:.3f} f={fc:.4f} rho={rho:.2f} | circ dchi2={row['dchi2_circ']:7.1f} k={row['k_circ']:.3f} b={row['b_circ']:.2f} T14={row['T14_circ_h']:.2f}h | ecc dchi2={row['dchi2_ecc']:6.1f} rho_app={row['rho_app_best']:.2f} max={rho_max:.2f} e_min={e_min:.3f}  ({time.time()-t0:.0f}s)", flush=True)
res["companions"] = rows
tag = f"{SIG}{'_' + NPZ.replace('.npz','') if NPZ != 'lc_windows.npz' else ''}{'_ldfree' if LDFREE else ''}"
json.dump(res, open(HERE / f"out/profile_{tag}.json", "w"), indent=1, default=float)
