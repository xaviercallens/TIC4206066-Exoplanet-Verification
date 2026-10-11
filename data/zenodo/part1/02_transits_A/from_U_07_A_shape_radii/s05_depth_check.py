"""Detrending-independent depth check: raw normalized flux (PDCSAP / QLP SYS_RM / DET), per-transit local quadratic
baseline fitted jointly with a box-free batman template (shape fixed at the baseline posterior median); compares with
the wotan-detrended depths. Also compares the wotan trend inside vs outside the mask."""
import json
import numpy as np
from astropy.io import fits
from common import *
ch = np.load(os.path.join(OUT, "chain_mann_claret_jit_base.npz"))["chain"]
T0, P, k, b, lnrho, u1, u2 = np.median(ch, axis=0)[:7]
ars = float(ars_from_rho(np.exp(lnrho), P)); T14 = float(t14_hours(P, k, b, ars)) / 24
import batman
pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = k; pm.a = ars; pm.inc = np.degrees(np.arccos(b / ars))
pm.ecc = 0; pm.w = 90; pm.u = [u1, u2]; pm.limb_dark = "quadratic"
DMID = 1 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0]
out = {}
def fit_local(t, f, cad, s, half=0.3, order=2):
    tmp = 1 - BatmanSet(t, cad, SUPERSAMPLE[s]).flux(T0, P, k, b, ars, u1, u2)
    ep = np.round((t - T0) / P).astype(int)
    cols = [tmp / DMID]
    for e in np.unique(ep):
        mm = (ep == e).astype(float); x = (t - (T0 + e * P)) * mm
        for o in range(order + 1):
            cols.append(mm * x ** o)
    A = np.array(cols).T
    c, *_ = np.linalg.lstsq(A, 1 - f, rcond=None)
    r = (1 - f) - A @ c
    sig = np.std(r)
    cov = np.linalg.inv(A.T @ A) * sig ** 2
    return c[0] * 1e6, np.sqrt(cov[0, 0]) * 1e6, len(np.unique(ep))
for s, cols in ((6, ["PDCSAP_FLUX", "SAP_FLUX"]), (32, ["PDCSAP_FLUX", "SAP_FLUX"]), (98, ["SYS_RM_FLUX", "DET_FLUX", "SAP_FLUX"])):
    fn, _, cad = FILES[s]
    d = fits.open(os.path.join(DATA, fn))[1].data
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    t0 = np.asarray(d["TIME"], float)
    for col in cols:
        f0 = np.asarray(d[col], float)
        g = (q == 0) & np.isfinite(f0) & np.isfinite(t0)
        t, f = t0[g], f0[g] / np.nanmedian(f0[g])
        dA, _ = phase_dt(t, P, T0); dB, _ = phase_dt(t, P_B, T0_B)
        m = (np.abs(dA) < 0.3) & (np.abs(dB) > T14_B_H / 24)
        for order in (1, 2):
            out[f"S{s}_{col}_poly{order}"] = fit_local(t[m], f[m], cad, s, order=order)
    z = np.load(os.path.join(OUT, f"lc_base_S{s}.npz"))
    t, f = z["t"], z["f"]
    dA, _ = phase_dt(t, P, T0); dB, _ = phase_dt(t, P_B, T0_B)
    m = (np.abs(dA) < 0.3) & (np.abs(dB) > T14_B_H / 24)
    out[f"S{s}_wotan_base_poly0"] = fit_local(t[m], f[m], cad, s, order=0)
for kk, v in out.items():
    print(f"{kk:30s} depth_mid {v[0]:7.1f} +- {v[1]:5.1f} ppm  (n_transit windows {v[2]})")
json.dump(dict(depth_mid_template_ppm=DMID * 1e6, results={kk: list(map(float, v)) for kk, v in out.items()}),
          open(os.path.join(OUT, "depth_check_local_baseline.json"), "w"), indent=1)
