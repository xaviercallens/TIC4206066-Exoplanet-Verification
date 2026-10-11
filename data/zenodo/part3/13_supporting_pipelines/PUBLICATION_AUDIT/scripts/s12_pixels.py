"""Step 12: pixel-level localisation of both signals with the OFFICIAL TESS PRF (TESS_PRF 0.1.3, MAST start_s0004
models), in every usable pixel data set: TESScut FFI cubes S6 (30 min), S32 (10 min), S98 (200 s) and the SPOC 2-min
TPF (S32).

Per data set and signal:
 1 cadence mask: lightkurve DEFAULT_BITMASK; S98 straylight cadences dropped (their background is not modelled well);
   FFI cubes background-subtracted per cadence (median of faint pixels > 5 px from target); TPF FLUX already is.
 2 astrometric/pointing solution: the median out-of-transit image is fitted with the PRF scene model of all Gaia DR3
   sources (T < 19, positions propagated to the sector epoch) with a free global shift (dx, dy), target flux and
   background -> PRF-based target position; neighbour positions carry the same shift.
 3 difference image per transit (linear per-pixel baseline from 0.75-2.5 T14, in-transit |dt| < 0.4 T14), averaged.
   Cadences within 1.0 T14 of the OTHER signal are excluded.
 4 noise map: per-pixel std of mean difference images built the same way at 400 random fake epoch sets (same number
   of epochs, away from both signals) -> includes red noise and systematics.
 5 fits: D = A * PRF(x, y) + c on pixels within 4 px of the target, weighted by the noise map: free position (grid
   0.05 px over +-3 px) and fixed positions (target, each Gaia neighbour).
 6 injection-recovery IN THE REAL PIXELS: (a) model-free on-target dips (real background-subtracted pixels within 3 px
   of the target multiplied by 1 - depth*box) at 300 fake epoch sets, fitted exactly like the data; (b) dips on each
   neighbour able to host the signal (neighbour PRF x flux, total flux loss matched) -> distributions of the fitted
   offset and of chi2(target) - chi2(neighbour).
"""
import glob
import json
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import prf_tools as pt  # noqa
from audit_common import AUD, save_json  # noqa

TARGET_ID = 3220388198192519424
PIX = 21.0  # arcsec per TESS pixel (nominal)
EPH = {k: dict(v) for k, v in lcio.EPH.items()}
ephf = os.path.join(AUD, "out", "ephem_current.json")
if os.path.exists(ephf):
    EPH.update(json.load(open(ephf)))

CUBES = {
    "S6_ffi": dict(glob="data/tesscut/tess-s0006-*_astrocut.fits", kind="ffi", sector=6, cam=1, ccd=1),
    "S32_ffi": dict(glob="data/tesscut/tess-s0032-*_astrocut.fits", kind="ffi", sector=32, cam=1, ccd=2),
    "S32_tpf2min": dict(glob="data/mast/tess2020324010417-s0032-*_tp.fits", kind="tpf", sector=32, cam=1, ccd=2),
    "S98_ffi": dict(glob="data/tesscut/tess-s0098-*_astrocut.fits", kind="ffi", sector=98, cam=1, ccd=2),
}


# ------------------------------------------------------------------------------------------------ utilities
def inject_on_target(cube, bkg, t, tc, dur, depth, x0, y0, radius=3.0):
    """MODEL-FREE injection into real pixels: background-subtracted pixel values within `radius` px of (x0, y0) are
    multiplied by (1 - depth) during |t - tc| < dur/2; the background is added back. No PRF model is involved."""
    out = cube.copy()
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    inside = np.hypot(xx - x0, yy - y0) <= radius
    tcs = np.atleast_1d(tc)
    box = np.zeros(len(t), bool)
    for c in tcs:
        box |= np.abs(t - c) < dur / 2
    for k in np.where(box)[0]:
        out[k][inside] = (cube[k][inside] - bkg[k]) * (1 - depth) + bkg[k]
    return out


def required_depth(depth_obs, fT_ap, fN_ap):
    """eclipse depth a neighbour needs so that the aperture light curve drops by depth_obs."""
    return depth_obs * (fT_ap + fN_ap) / fN_ap


def required_depths_example(sector_cube="S98_ffi", radius=2.0, depth=5e-4):
    """used by tests/test_prf_regressions.py (R4): for the real Gaia neighbours within 60'' in the S98 stamp, compare the
    eclipse depth needed with fluxes INSIDE a radius-2 px target aperture (PRF-integrated) with the full-inclusion
    value. For an aperture centred on the target the aperture requirement can only be larger."""
    C = load_cube(sector_cube)
    shape = C["cube"].shape[1:]
    sc = pt.scene(C["sources"], C["wcs"], C["tmid"], TARGET_ID)
    tg = [s for s in sc if s["is_target"]][0]
    prf = pt.get_prf(C["cam"], C["ccd"], C["sector"], C["col0"] + tg["x"], C["row0"] + tg["y"])
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    ap = np.hypot(xx - tg["x"], yy - tg["y"]) <= radius
    mT = prf.locate(tg["x"], tg["y"], shape)
    ex = []
    for s in sc:
        if s["is_target"] or s["sep"] > 60:
            continue
        mN = prf.locate(s["x"], s["y"], shape) * s["rel_flux"]
        fT_ap, fN_ap = mT[ap].sum(), mN[ap].sum()
        if fN_ap <= 0:
            continue
        ex.append(dict(source_id=s["source_id"], sep=s["sep"], req_aperture=float(required_depth(depth, fT_ap, fN_ap)),
                       req_full_inclusion=float(required_depth(depth, mT.sum(), mN.sum())),
                       frac_target_in_ap=float(fT_ap / mT.sum()), frac_neighbour_in_ap=float(fN_ap / mN.sum())))
    return ex


def load_cube(name, drop_stray=True):
    c = CUBES[name]
    path = glob.glob(os.path.join(AUD, c["glob"]))[0]
    h = fits.open(path)
    d = h[1].data
    t = np.asarray(d["TIME"], float)
    fl = np.asarray(d["FLUX"], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    good = np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1) & ((q & lcio.DEFAULT_BITMASK) == 0)
    if drop_stray and c["kind"] == "ffi" and c["sector"] == 98:
        good &= (q & lcio.STRAY) == 0
    t, fl, q = t[good], fl[good], q[good]
    if c["kind"] == "ffi":
        t = t + lcio.TESSCUT_SHIFT_S[c["sector"]] / 86400.0
    w = WCS(h[2].header)
    ny, nx = fl.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    src = pt.gaia_sources(180)
    tmid = float(np.median(t))
    sc0 = pt.scene(src, w, tmid, TARGET_ID)
    tg = [s for s in sc0 if s["is_target"]][0]
    if c["kind"] == "ffi":
        med = np.nanmedian(fl, axis=0)
        far = np.hypot(xx - tg["x"], yy - tg["y"]) > 5
        bpix = far & (med < np.nanpercentile(med[far], 30))
        bkg = np.nanmedian(fl[:, bpix], axis=1)
        fl = fl - bkg[:, None, None]
    else:
        bkg = np.zeros(len(t))
    pc1 = np.asarray(d["POS_CORR1"], float)[good] if "POS_CORR1" in d.columns.names else np.zeros(len(t))
    pc2 = np.asarray(d["POS_CORR2"], float)[good] if "POS_CORR2" in d.columns.names else np.zeros(len(t))
    col0, row0 = h[1].header.get("1CRV4P"), h[1].header.get("2CRV4P")
    return dict(name=name, t=t, cube=fl, bkg=bkg, wcs=w, col0=col0, row0=row0, tmid=tmid, sources=src, pc1=pc1, pc2=pc2,
                cad=float(np.median(np.diff(t))), **c)


def astrometric_solution(C, prf):
    """fit the median image with the PRF scene model: global shift + target flux + background (+ neighbour fluxes fixed
    to their T-mag ratios)."""
    med = np.nanmedian(C["cube"], axis=0)
    ny, nx = med.shape
    sc = pt.scene(C["sources"], C["wcs"], C["tmid"], TARGET_ID, Tlim=19)
    tg = [s for s in sc if s["is_target"]][0]
    yy, xx = np.mgrid[:ny, :nx]
    fitpix = np.hypot(xx - tg["x"], yy - tg["y"]) <= 5

    def resid(p):
        dx, dy = p
        srcs = [dict(s, x=s["x"] + dx, y=s["y"] + dy) for s in sc]
        m = pt.model_image(prf, med.shape, srcs)
        A = np.vstack([m[fitpix], np.ones(fitpix.sum())]).T
        coef, *_ = np.linalg.lstsq(A, med[fitpix], rcond=None)
        return np.sum((med[fitpix] - A @ coef) ** 2), coef

    best = None
    for dx0 in (-0.2, 0.0, 0.2):
        for dy0 in (-0.2, 0.0, 0.2):
            r = minimize(lambda p: resid(p)[0], [dx0, dy0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-6))
            if best is None or r.fun < best.fun:
                best = r
    chi, coef = resid(best.x)
    srcs = [dict(s, x=s["x"] + best.x[0], y=s["y"] + best.x[1]) for s in sc]
    model = pt.model_image(prf, med.shape, srcs) * coef[0] + coef[1]
    frac_resid = float(np.sqrt(np.sum((med[fitpix] - model[fitpix]) ** 2) / np.sum(med[fitpix] ** 2)))
    return dict(shift_px=[float(best.x[0]), float(best.x[1])], target_flux=float(coef[0]), bkg=float(coef[1]),
                frac_rms_resid=frac_resid), srcs, model


def epochs_for(C, sig):
    e = EPH[sig]
    other = "B" if sig == "A" else "A"
    t = C["t"]
    T14 = e["T14"]
    okc = ~lcio.in_transit(t, other, 1.0, EPH)
    n0 = int(np.floor((t.min() - e["T0"]) / e["P"])) - 1
    n1 = int(np.ceil((t.max() - e["T0"]) / e["P"])) + 1
    eps = []
    for n in range(n0, n1 + 1):
        tc = e["T0"] + n * e["P"]
        i_in = okc & (np.abs(t - tc) < 0.4 * T14)
        i_l = okc & (t - tc < -0.75 * T14) & (t - tc > -2.5 * T14)
        i_r = okc & (t - tc > 0.75 * T14) & (t - tc < 2.5 * T14)
        if i_in.sum() >= 2 and i_l.sum() >= 2 and i_r.sum() >= 2:
            eps.append((n, tc))
    return eps


def diff_image(C, cube, tcs, T14, other_mask):
    t = C["t"]
    ny, nx = cube.shape[1:]
    D = []
    for tc in tcs:
        i_in = other_mask & (np.abs(t - tc) < 0.4 * T14)
        i_o = other_mask & (np.abs(t - tc) > 0.75 * T14) & (np.abs(t - tc) < 2.5 * T14)
        if i_in.sum() < 2 or i_o.sum() < 4:
            continue
        X = np.vstack([np.ones(i_o.sum()), t[i_o] - tc]).T
        coef, *_ = np.linalg.lstsq(X, cube[i_o].reshape(i_o.sum(), -1), rcond=None)
        base = coef[0].reshape(ny, nx)
        D.append(base - cube[i_in].mean(axis=0))
    return np.array(D)


class GridFitter:
    """precomputed PRF models on a position grid -> vectorised linear fits D = A*PRF + c."""

    def __init__(self, prf, shape, x0, y0, half=3.0, step=0.05, rfit=4.0):
        self.gx = np.arange(x0 - half, x0 + half + 1e-9, step)
        self.gy = np.arange(y0 - half, y0 + half + 1e-9, step)
        ny, nx = shape
        yy, xx = np.mgrid[:ny, :nx]
        self.pix = np.hypot(xx - x0, yy - y0) <= rfit
        self.models = np.array([[prf.locate(x, y, shape)[self.pix] for x in self.gx] for y in self.gy])  # (ny_g,nx_g,npix)
        self.prf, self.shape = prf, shape

    def fit(self, D, sig):
        d = D[self.pix]
        w = 1.0 / sig[self.pix] ** 2
        M = self.models
        Sw, Sd, Sm = w.sum(), (w * d).sum(), (M * w).sum(-1)
        Smm, Smd = (M * M * w).sum(-1), (M * (w * d)).sum(-1)
        det = Smm * Sw - Sm ** 2
        A = (Smd * Sw - Sm * Sd) / det
        c = (Sd - A * Sm) / Sw
        chi2 = (w * d * d).sum() - 2 * A * Smd - 2 * c * Sd + A * A * Smm + 2 * A * c * Sm + c * c * Sw
        j = np.unravel_index(np.argmin(chi2), chi2.shape)
        return dict(x=float(self.gx[j[1]]), y=float(self.gy[j[0]]), A=float(A[j]), c=float(c[j]), chi2=float(chi2[j]),
                    chi2map=chi2)

    def fit_fixed(self, D, sig, x, y):
        d = D[self.pix]
        w = 1.0 / sig[self.pix] ** 2
        m = self.prf.locate(x, y, self.shape)[self.pix]
        X = np.vstack([m, np.ones_like(m)]).T
        Wm = np.diag(w)
        coef = np.linalg.solve(X.T @ Wm @ X, X.T @ Wm @ d)
        r = d - X @ coef
        return dict(A=float(coef[0]), c=float(coef[1]), chi2=float(np.sum(w * r * r)))


def random_epoch_sets(C, sig, n_ep, nsets, rng):
    e = EPH[sig]
    t = C["t"]
    T14 = e["T14"]
    bad = lcio.in_transit(t, "A", 3.0, EPH) | lcio.in_transit(t, "B", 3.0, EPH)
    cand = t[~bad]
    sets = []
    tries = 0
    while len(sets) < nsets and tries < nsets * 200:
        tries += 1
        cs = []
        for _ in range(n_ep * 20):
            c = rng.choice(cand)
            if (np.abs(t - c) < 0.4 * T14).sum() < 2:
                continue
            l = ((t - c < -0.75 * T14) & (t - c > -2.5 * T14)).sum()
            r_ = ((t - c > 0.75 * T14) & (t - c < 2.5 * T14)).sum()
            win = np.abs(t - c) < 2.6 * T14
            if l < 2 or r_ < 2 or bad[win].any() or any(abs(c - x) < 6 * T14 for x in cs):
                continue
            cs.append(c)
            if len(cs) == n_ep:
                break
        if len(cs) == n_ep:
            sets.append(cs)
    return sets


def analyse(name, sig, n_noise=400, n_inj=300, seed=7):
    rng = np.random.default_rng(seed)
    C = load_cube(name)
    shape = C["cube"].shape[1:]
    sc0 = pt.scene(C["sources"], C["wcs"], C["tmid"], TARGET_ID)
    tg0 = [s for s in sc0 if s["is_target"]][0]
    prf = pt.get_prf(C["cam"], C["ccd"], C["sector"], C["col0"] + tg0["x"], C["row0"] + tg0["y"])
    astro, srcs, oot_model = astrometric_solution(C, prf)
    tg = [s for s in srcs if s["is_target"]][0]
    e = EPH[sig]
    other = "B" if sig == "A" else "A"
    okc = ~lcio.in_transit(C["t"], other, 1.0, EPH)
    eps = epochs_for(C, sig)
    if len(eps) == 0:
        return dict(name=name, signal=sig, n_epochs=0)
    tcs = [x[1] for x in eps]
    Dreal = diff_image(C, C["cube"], tcs, e["T14"], okc)
    Dm = Dreal.mean(axis=0)
    # noise map from fake epoch sets
    fsets = random_epoch_sets(C, sig, len(Dreal), n_noise, rng)
    fakes = np.array([diff_image(C, C["cube"], s, e["T14"], np.ones(len(C["t"]), bool)).mean(axis=0) for s in fsets])
    sig_map = np.std(fakes, axis=0)
    sig_map = np.where(sig_map > 0, sig_map, np.nanmedian(sig_map))
    G = GridFitter(prf, shape, tg["x"], tg["y"])
    free = G.fit(Dm, sig_map)
    fix_t = G.fit_fixed(Dm, sig_map, tg["x"], tg["y"])
    depth_ppm = fix_t["A"] / astro["target_flux"] * 1e6
    # neighbours that could host the signal (T-mag flux ratio with full inclusion; conservative T for missing colour)
    neigh = []
    for s in srcs:
        if s["is_target"]:
            continue
        fr = s["rel_flux"] * (10 ** (0.4 * s["T_err"]) if s["T_flag"] != "ok" else 1.0)   # brighter by 1 sigma if uncertain
        req_full = depth_ppm * 1e-6 * (1 + fr) / fr
        dist = np.hypot(s["x"] - tg["x"], s["y"] - tg["y"])
        if dist * PIX > 150:
            continue
        ff = G.fit_fixed(Dm, sig_map, s["x"], s["y"])
        neigh.append(dict(source_id=s["source_id"], T=s["T"], T_flag=s["T_flag"], sep_arcsec=float(dist * PIX),
                          x=s["x"], y=s["y"], rel_flux=s["rel_flux"], req_eclipse_full_inclusion=float(req_full),
                          chi2=ff["chi2"], dchi2_vs_target=float(ff["chi2"] - fix_t["chi2"]), A=ff["A"]))
    neigh.sort(key=lambda r: r["sep_arcsec"])
    # injections on the target (model-free, real pixels)
    isets = random_epoch_sets(C, sig, len(Dreal), n_inj, np.random.default_rng(seed + 1))
    inj = []
    for s in isets:
        cube_i = inject_on_target(C["cube"] + C["bkg"][:, None, None], C["bkg"], C["t"], s, e["T14"], depth_ppm * 1e-6,
                                  tg["x"], tg["y"], 3.0) - C["bkg"][:, None, None]
        Di = diff_image(C, cube_i, s, e["T14"], np.ones(len(C["t"]), bool)).mean(axis=0)
        fr_ = G.fit(Di, sig_map)
        ft = G.fit_fixed(Di, sig_map, tg["x"], tg["y"])
        rec = dict(dx=fr_["x"] - tg["x"], dy=fr_["y"] - tg["y"], chi2_t=ft["chi2"], A=ft["A"])
        for nb in neigh[:12]:
            rec[f"dchi2_{nb['source_id']}"] = G.fit_fixed(Di, sig_map, nb["x"], nb["y"])["chi2"] - ft["chi2"]
        inj.append(rec)
    idx_ = np.array([[r["dx"], r["dy"]] for r in inj])
    cov = np.cov(idx_.T)
    mu = idx_.mean(axis=0)
    off = np.array([free["x"] - tg["x"], free["y"] - tg["y"]])
    maha = float(np.sqrt((off - mu) @ np.linalg.inv(cov) @ (off - mu)))
    frac_further = float(np.mean(np.hypot(idx_[:, 0], idx_[:, 1]) >= np.hypot(*off)))
    # injections on capable neighbours (flux loss matched to the observed aperture loss)
    nb_inj = {}
    for nb in [n for n in neigh if n["req_eclipse_full_inclusion"] <= 1.0][:6]:
        recs = []
        m_nb = prf.locate(nb["x"], nb["y"], shape) * nb["rel_flux"] * astro["target_flux"]
        m_tg = prf.locate(tg["x"], tg["y"], shape) * astro["target_flux"]
        ap = np.hypot(*np.mgrid[:shape[0], :shape[1]][::-1] - np.array([tg["x"], tg["y"]])[:, None, None]) <= 2.0
        loss = depth_ppm * 1e-6 * (m_tg[ap].sum() + m_nb[ap].sum())
        ecl = loss / m_nb[ap].sum()
        for s in isets[:150]:
            cube_i = C["cube"].copy()
            box = np.zeros(len(C["t"]), bool)
            for c in s:
                box |= np.abs(C["t"] - c) < e["T14"] / 2
            cube_i[box] -= ecl * m_nb
            Di = diff_image(C, cube_i, s, e["T14"], np.ones(len(C["t"]), bool)).mean(axis=0)
            fr_ = G.fit(Di, sig_map)
            ft = G.fit_fixed(Di, sig_map, tg["x"], tg["y"])
            fn = G.fit_fixed(Di, sig_map, nb["x"], nb["y"])
            recs.append(dict(dx=fr_["x"] - tg["x"], dy=fr_["y"] - tg["y"], dchi2=fn["chi2"] - ft["chi2"]))
        dchi_real = nb["dchi2_vs_target"]
        dn = np.array([r["dchi2"] for r in recs])
        dt_ = np.array([r[f"dchi2_{nb['source_id']}"] for r in inj if f"dchi2_{nb['source_id']}" in r])
        nb_inj[nb["source_id"]] = dict(sep_arcsec=nb["sep_arcsec"], eclipse_needed_aperture=float(ecl),
                                       real_dchi2=dchi_real, inj_on_neighbour_dchi2_median=float(np.median(dn)),
                                       inj_on_neighbour_dchi2_p95=float(np.percentile(dn, 95)),
                                       frac_neighbour_inj_with_dchi2_ge_real=float(np.mean(dn >= dchi_real)),
                                       inj_on_target_dchi2_median=float(np.median(dt_)) if len(dt_) else None,
                                       frac_target_inj_with_dchi2_le_real=float(np.mean(dt_ <= dchi_real)) if len(dt_) else None,
                                       inj_on_neighbour_offset_median_arcsec=float(np.median(np.hypot([r["dx"] for r in recs], [r["dy"] for r in recs])) * PIX))
    # sky offsets
    t_sky = C["wcs"].all_pix2world([[tg["x"] - astro["shift_px"][0], tg["y"] - astro["shift_px"][1]]], 0)[0]
    f_sky = C["wcs"].all_pix2world([[free["x"] - astro["shift_px"][0], free["y"] - astro["shift_px"][1]]], 0)[0]
    dra = (f_sky[0] - t_sky[0]) * np.cos(np.radians(t_sky[1])) * 3600
    dde = (f_sky[1] - t_sky[1]) * 3600
    res = dict(name=name, signal=sig, n_epochs=len(Dreal), epochs=[int(x[0]) for x in eps], cadence_min=C["cad"] * 1440,
               astrometry=astro, target_px=[tg["x"], tg["y"]], depth_ppm_from_prf_fit=depth_ppm,
               free_fit=dict(x=free["x"], y=free["y"], dx_px=float(off[0]), dy_px=float(off[1]), offset_arcsec=float(np.hypot(*off) * PIX),
                             dRA_arcsec=float(dra), dDec_arcsec=float(dde), chi2=free["chi2"]),
               target_fixed=fix_t, n_fit_pixels=int(G.pix.sum()),
               injection_on_target=dict(n=len(inj), mean_dx_px=float(mu[0]), mean_dy_px=float(mu[1]),
                                        sd_dx_px=float(np.sqrt(cov[0, 0])), sd_dy_px=float(np.sqrt(cov[1, 1])),
                                        sd_arcsec=float(np.sqrt(np.trace(cov) / 2) * PIX),
                                        mahalanobis_real=maha, frac_injections_further_than_real=frac_further,
                                        median_recovered_A=float(np.median([r["A"] for r in inj]))),
               neighbours=neigh[:15], neighbour_injections=nb_inj, n_noise_sets=len(fakes))
    return res


if __name__ == "__main__":
    which = sys.argv[1:] or list(CUBES)
    out = {}
    for name in which:
        for sig in ("A", "B"):
            try:
                r = analyse(name, sig)
            except Exception as ex:  # keep going, record the failure
                r = dict(name=name, signal=sig, error=repr(ex))
            out[f"{name}_{sig}"] = r
            if "free_fit" in r:
                ff, ij = r["free_fit"], r["injection_on_target"]
                print(f"{name:12s} {sig}: n_ep={r['n_epochs']} depth(PRF)={r['depth_ppm_from_prf_fit']:.0f} ppm; "
                      f"offset {ff['offset_arcsec']:.1f}\" (dRA {ff['dRA_arcsec']:+.1f}, dDec {ff['dDec_arcsec']:+.1f}); "
                      f"on-target injections sd {ij['sd_arcsec']:.1f}\", Mahalanobis {ij['mahalanobis_real']:.2f}, "
                      f"frac further {ij['frac_injections_further_than_real']:.2f}; astro shift {r['astrometry']['shift_px']}", flush=True)
                for k, v in r["neighbour_injections"].items():
                    print(f"      neighbour {k} {v['sep_arcsec']:.1f}\": needs {v['eclipse_needed_aperture']*100:.0f}% eclipse; real dchi2(N-T) "
                          f"{v['real_dchi2']:.1f}; neighbour-inj median {v['inj_on_neighbour_dchi2_median']:.1f} (frac>=real "
                          f"{v['frac_neighbour_inj_with_dchi2_ge_real']:.3f}); target-inj median {v['inj_on_target_dchi2_median']}", flush=True)
            else:
                print(name, sig, r.get("error", "no epochs"), flush=True)
            save_json(out, f"out/s12_pixels_{'_'.join(which)}.json")
