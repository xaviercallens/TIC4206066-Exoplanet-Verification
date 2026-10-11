"""Verifier: independent quick re-derivation of S32/S98 3.18 d (and 11.13 d) difference-image positions.
Own code. Differences from auditor: POS_CORR handled as extra per-event regressors (not a global pre-decorrelation);
target reference = PRF fit of the median image (single star, free bkg); errors from event bootstrap (not injections).
"""
import os, sys, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import minimize
import PRF

D = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data"
OUT = os.path.dirname(os.path.abspath(__file__))
F = {"S32_tpf": D + "/mast/tess2020324010417-s0032-0000000004206066-0200-s_tp.fits",
     "S98_ffi": D + "/tesscut/21x21/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits",
     "S6_ffi": D + "/tesscut/21x21/tess-s0006-1-1_80.859830_-1.326960_21x21_astrocut.fits"}
QLP98 = D + "/mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"
EPH = {"A": (3.1827856, 2460998.7375 - 2457000, 1.87 / 24), "B": (11.13274, 2461000.865 - 2457000, 1.65 / 24)}
BITS = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 512 | 2048 | 4096 | 16384


def load(key):
    h = fits.open(F[key]); d = h[1].data
    t = np.array(d["TIME"], float); fl = np.array(d["FLUX"], float); q = np.array(d["QUALITY"])
    pc1 = np.array(d["POS_CORR1"], float); pc2 = np.array(d["POS_CORR2"], float)
    ok = np.isfinite(t) & ((q & BITS) == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), 1)
    if key == "S98_ffi":
        ql = fits.getdata(QLP98, 1); tq = np.array(ql["TIME"], float); qq = np.array(ql["QUALITY"])
        badq = tq[(qq & (2048 | (1 << 29))) > 0]
        j = np.searchsorted(np.sort(badq), t)
        sb = np.sort(badq)
        near = np.zeros(len(t), bool)
        for jj in (j - 1, j):
            jj = np.clip(jj, 0, len(sb) - 1); near |= np.abs(sb[jj] - t) < 150 / 86400.
        ok &= ~near
    h1 = h[1].header
    return dict(t=t[ok], F=fl[ok], pc1=pc1[ok], pc2=pc2[ok], wcs=WCS(h[2].header), cam=h[0].header["CAMERA"],
                ccd=h[0].header["CCD"], sector=h[0].header["SECTOR"], c0=h1["1CRV4P"], r0=h1["2CRV4P"],
                dt=h1["TIMEDEL"])


def phasemask(t, eph, pad):
    P, T0, T14 = eph
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    return np.abs(ph * P) < T14 / 2 + pad


def events(ds, sig, W=4 / 24):
    P, T0, T14 = EPH[sig]
    other = "B" if sig == "A" else "A"
    om = phasemask(ds["t"], EPH[other], 0.5 / 24)
    n = np.arange(np.ceil((ds["t"].min() - T0) / P), np.floor((ds["t"].max() - T0) / P) + 1)
    out = []
    for tc in T0 + n * P:
        idx = np.where((np.abs(ds["t"] - tc) <= T14 / 2 + W) & ~om)[0]
        if len(idx) < 5: continue
        tt = ds["t"][idx]
        nin = np.sum(np.abs(tt - tc) < T14 / 2); npre = np.sum(tt < tc - T14 / 2); npost = np.sum(tt > tc + T14 / 2)
        need = 0.3 * W / ds["dt"]
        if nin >= 0.6 * T14 / ds["dt"] and npre >= need and npost >= need:
            out.append((tc, idx))
    return out


def ediff(ds, tc, idx, sig, usepc):
    t = ds["t"][idx]; x = t - tc
    T14 = EPH[sig][2]; dt = ds["dt"]
    lo = np.maximum(t - dt / 2, tc - T14 / 2); hi = np.minimum(t + dt / 2, tc + T14 / 2)
    box = np.clip(hi - lo, 0, None) / dt
    cols = [np.ones_like(x), x]
    if usepc:
        for p in (ds["pc1"][idx], ds["pc2"][idx]):
            if np.all(np.isfinite(p)): cols.append(p - p.mean())
    cols.append(box)
    X = np.vstack(cols).T
    Y = ds["F"][idx].reshape(len(idx), -1)
    Ci = np.linalg.pinv(X.T @ X); b = Ci @ X.T @ Y
    r = Y - X @ b; s2 = (r ** 2).sum(0) / max(len(idx) - X.shape[1], 1)
    sh = ds["F"].shape[1:]
    return (-b[-1]).reshape(sh), (s2 * Ci[-1, -1]).reshape(sh)


def combine(ds, evs, sig, usepc):
    Ds, Vs = zip(*[ediff(ds, tc, idx, sig, usepc) for tc, idx in evs])
    Ds, Vs = np.array(Ds), np.array(Vs); w = 1 / Vs
    return (w * Ds).sum(0) / w.sum(0), 1 / w.sum(0)


class Fitter:
    def __init__(self, ds):
        self.ds = ds
        img = np.nanmedian(ds["F"], 0); self.img = img
        ny, nx = img.shape
        # target from Gaia DR3 propagated
        import pandas as pd
        g = pd.read_csv(D + "/catalogs/gaiadr3_gaia_source_cone240as.csv")
        r = g[g.source_id == 3220388198192519424].iloc[0]
        ep = 2000 + (np.median(ds["t"]) + 2457000 - 2451545) / 365.25
        dec = r.dec + r.pmdec * (ep - 2016) / 3.6e6; ra = r.ra + r.pmra * (ep - 2016) / 3.6e6 / np.cos(np.radians(dec))
        x0, y0 = ds["wcs"].all_world2pix([ra], [dec], 0); self.xw, self.yw = float(x0[0]), float(y0[0])
        self.ra, self.dec = ra, dec
        self.prf = PRF.TESS_PRF(ds["cam"], ds["ccd"], ds["sector"], ds["c0"] + self.xw, ds["r0"] + self.yw,
                                localdatadir=D + "/prf/start_s0004")
        xc, yc = int(round(self.xw)), int(round(self.yw)); h = 4
        self.sl = (slice(max(0, yc - h), min(ny, yc + h + 1)), slice(max(0, xc - h), min(nx, xc + h + 1)))
        self.shape = img.shape
        # OOT reference position: fit median image with single PRF
        self.xt, self.yt = self.fit(img, np.clip(img, 1, None) + 100.0)[:2]

    def chi2(self, d, s, x, y):
        P = self.prf.locate(x, y, self.shape)[self.sl].ravel()
        A = np.vstack([P, np.ones_like(P)]).T / s[:, None]; bb = d / s
        c, *_ = np.linalg.lstsq(A, bb, rcond=None); r = bb - A @ c
        return r @ r, c

    def fit(self, D, V):
        d = D[self.sl].ravel(); s = np.sqrt(V[self.sl].ravel())
        best = min(((self.chi2(d, s, self.xw + a, self.yw + b)[0], self.xw + a, self.yw + b)
                    for a in np.arange(-2, 2.01, .5) for b in np.arange(-2, 2.01, .5)))
        r = minimize(lambda p: self.chi2(d, s, p[0], p[1])[0], best[1:], method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-4))
        c2, c = self.chi2(d, s, *r.x)
        return r.x[0], r.x[1], c[0]

    def to_arcsec(self, x, y):
        w = self.ds["wcs"]
        ra0, de0 = w.all_pix2world([self.xt], [self.yt], 0); ra1, de1 = w.all_pix2world([x], [y], 0)
        return float((ra1[0] - ra0[0]) * np.cos(np.radians(de0[0])) * 3600), float((de1[0] - de0[0]) * 3600)


def run(key, sig, usepc, nboot=200, drop=None, seed=0):
    ds = load(key); f = Fitter(ds)
    evs = events(ds, sig)
    if drop is not None:
        evs = [e for e in evs if abs(e[0] - drop) > 0.05]
    Dm, Vm = combine(ds, evs, sig, usepc)
    x, y, a = f.fit(Dm, Vm); off = f.to_arcsec(x, y)
    rng = np.random.default_rng(seed); bo = []
    per = [ediff(ds, tc, idx, sig, usepc) for tc, idx in evs]
    for k in range(nboot):
        ii = rng.integers(0, len(evs), len(evs))
        Ds = np.array([per[i][0] for i in ii]); Vs = np.array([per[i][1] for i in ii]); w = 1 / Vs
        xb, yb, _ = f.fit((w * Ds).sum(0) / w.sum(0), 1 / w.sum(0)); bo.append(f.to_arcsec(xb, yb))
    bo = np.array(bo)
    jk = []
    for k in (range(len(evs)) if len(evs) >= 3 else []):
        Ds = np.array([per[i][0] for i in range(len(evs)) if i != k]); Vs = np.array([per[i][1] for i in range(len(evs)) if i != k]); w = 1 / Vs
        xb, yb, _ = f.fit((w * Ds).sum(0) / w.sum(0), 1 / w.sum(0)); jk.append((float(evs[k][0]),) + f.to_arcsec(xb, yb))
    return dict(key=key, sig=sig, usepc=usepc, drop=drop, n=len(evs), tcs=[float(e[0]) for e in evs], off=off,
                boot_sd=bo.std(0, ddof=1).tolist(), boot_med=np.median(bo, 0).tolist(), jackknife=jk,
                oot_ref_minus_wcs_px=[f.xt - f.xw, f.yt - f.yw], amp=float(a))


if __name__ == "__main__":
    res = []
    for args in [("S32_tpf", "A", False), ("S32_tpf", "A", True), ("S32_tpf", "B", False), ("S32_tpf", "B", True),
                 ("S98_ffi", "A", False), ("S98_ffi", "A", True), ("S98_ffi", "B", False)]:
        r = run(*args, nboot=100); res.append(r)
        print(r["key"], r["sig"], "pc" if r["usepc"] else "raw", r["n"], np.round(r["off"], 1), "bootSD", np.round(r["boot_sd"], 1),
              "ref-wcs px", np.round(r["oot_ref_minus_wcs_px"], 3), flush=True)
    r = run("S98_ffi", "A", False, nboot=100, drop=3995.555); res.append(r)
    print("S98 A drop3995.555", r["n"], np.round(r["off"], 1), np.round(r["boot_sd"], 1))
    json.dump(res, open(OUT + "/v1_loc.json", "w"), indent=1)
