"""V04: B linear ephemeris from per-event mid-times (S6 x2 TESS-SPOC PDCSAP; S98 E0, E1 QLP q==0 own-detrended SYS_RM),
trapezoid (T14 1.6 h, ingress 15%), cadence-integrated, local linear baseline, tc scanned on a 0.25-min grid.
Timing errors from injection-recovery of the same template at random epochs in the same light curve."""
import json, numpy as np
import vcommon as v

rng = np.random.default_rng(404)


def trap(dt, T14=v.DUR, tau=0.15 * v.DUR):
    return np.clip((T14 / 2 - np.abs(dt)) / tau, 0, 1)


def model_unit(x, cad, nss):
    off = (np.arange(nss) + 0.5) / nss - 0.5
    return trap(x[:, None] + off[None, :] * cad).mean(1)


def fit_tc(t, f, tc0, cad, nss, scan=60.0, step=0.25, win=3 * v.DUR):
    m = np.isfinite(f) & (np.abs(t - tc0) < win + scan / 1440)
    tt, ff = t[m], f[m]
    best = None
    for dtm in np.arange(-scan, scan + 1e-9, step):
        x = tt - tc0 - dtm / 1440
        s = model_unit(x, cad, nss)
        if np.sum(s > 0.5) < 2:
            continue
        X = np.vstack([np.ones_like(x), x, -s]).T
        c, *_ = np.linalg.lstsq(X, ff, rcond=None)
        chi = np.sum((ff - X @ c) ** 2)
        if c[2] > 0 and (best is None or chi < best[0]):
            best = (chi, dtm, c[2])
    return None if best is None else dict(dt_min=best[1], D=best[2])


def prep_s6():
    d, t, q = v.load("ts6")
    f = np.asarray(d["PDCSAP_FLUX"], float)
    ok = np.isfinite(t) & np.isfinite(f) & ((q & v.SPOC_BAD) == 0)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    A = v.inA(t)
    B = np.zeros_like(t, bool)
    for n in (-227, -226):
        B |= np.abs(t - (v.T0_MS + n * v.P_MS)) < v.DUR
    fd = v.biweight(t, f, A | B, 0.75)
    return t[~A], fd[~A], 1800 / 86400, 15


def prep_s98():
    d, t, q = v.load("qlp98")
    f = np.asarray(d["SYS_RM_FLUX"], float)
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    A = v.inA(t)
    B = np.zeros_like(t, bool)
    for n in range(-1, 5):
        B |= np.abs(t - (v.T0_MS + n * v.P_MS)) < v.DUR
    fd = v.biweight(t, f, A | B, 0.75)
    return t[~A], fd[~A], 200 / 86400, 3


def inj_err(t, f, cad, nss, depth, avoid, n=250):
    dts = []
    k = 0
    while len(dts) < n and k < 50 * n:
        k += 1
        tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
        if np.min(np.abs(tc - np.array(avoid))) < 4 * v.DUR:
            continue
        w = np.abs(t - tc) < v.DUR / 2
        if w.sum() < 0.8 * v.DUR / cad:
            continue
        fi = f - depth * model_unit(t - tc, cad, nss)
        r = fit_tc(t, fi, tc, cad, nss, scan=90)
        if r is not None:
            dts.append(r["dt_min"])
    dts = np.array(dts)
    return float(np.std(dts)), float(v.rstd(dts)), len(dts)


ev = []
for name, prep, epochs in (("S6", prep_s6, (-227, -226)), ("S98", prep_s98, (0, 1))):
    t, f, cad, nss = prep()
    avoid = [v.T0_MS + n * v.P_MS for n in range(-230, 6)]
    for n in epochs:
        tc0 = v.T0_MS + n * v.P_MS
        r = fit_tc(t, f, tc0, cad, nss)
        e_sd, e_r, ne = inj_err(t, f, cad, nss, r["D"], avoid)
        ev.append(dict(sector=name, epoch=n, tc=tc0 + r["dt_min"] / 1440, dt_min=r["dt_min"], D_ppm=r["D"] * 1e6, e_inj_sd_min=e_sd, e_inj_rstd_min=e_r, n_inj=ne))
        print(ev[-1], flush=True)


def linfit(E, T, S):
    A = np.vstack([np.ones_like(E), E]).T
    Wm = np.diag(1 / S ** 2)
    cov = np.linalg.inv(A.T @ Wm @ A)
    p = cov @ A.T @ Wm @ T
    chi2 = float(np.sum(((T - A @ p) / S) ** 2))
    return dict(T0=float(p[0]), P=float(p[1]), eT0=float(np.sqrt(cov[0, 0])), eP=float(np.sqrt(cov[1, 1])), chi2=chi2, dof=len(T) - 2)


E = np.array([e["epoch"] for e in ev], float); T = np.array([e["tc"] for e in ev])
out = {}
for k in ("e_inj_sd_min", "e_inj_rstd_min"):
    S = np.array([e[k] for e in ev]) / 1440
    out[k] = linfit(E, T, S)
    r = out[k]
    print(k, "P = %.6f +- %.6f  T0 = %.4f +- %.4f  chi2 %.2f/%d" % (r["P"], r["eP"], r["T0"], r["eT0"], r["chi2"], r["dof"]))
json.dump(dict(events=ev, fits=out, manuscript=dict(P=11.13274, eP=0.00003, T0=4000.865, eT0=0.003)), open(v.HERE + "/out/v04_ephem.json", "w"), indent=1)
