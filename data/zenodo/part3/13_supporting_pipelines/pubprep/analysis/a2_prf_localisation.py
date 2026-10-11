"""Independent pixel-level localisation of both signals with the official TESS PRF (Step 3).

Written from scratch for this check (no code from the deep follow-up, the TIC 158955499 audit or the publication
audit is imported). Inputs are raw data only: TESScut cubes (S6, S32, S98), the SPOC S32 2-min TPF, the MAST
start_s0004 PRF files (read with TESS_PRF 0.1.3 = PyPI 'tess-prf', which un-interleaves the 9x9 sub-pixel PRFs) and the
Gaia DR3 census (proper motions propagated to each sector's epoch).

Method (per pixel set and signal):
 1 cadence selection: SPOC quality DEFAULT_BITMASK (17087) removed; straylight (2048) removed in S98; cadences within
   1 T14 of the other signal removed; FFI cubes background-subtracted (median of faint pixels > 5 px from the target).
 2 scene: the median image is fitted with a PRF scene of all Gaia sources (T < 19) with a free global shift; this
   fixes the target's stamp position (x_t, y_t) at the sector epoch.
 3 difference image: for each transit, per-pixel linear baseline from cadences 0.75-2.5 T14 from mid-transit minus
   the mean in-transit image (|dt| < 0.4 T14); averaged over the transits.
 4 noise: per-pixel scatter of 300 difference images built from random fake epoch sets away from both signals.
 5 fits: D = A PRF(x, y) + c over pixels within 4 px of the target; free position (0.1-px grid, +-3 px, parabolic
   refinement) and fixed positions (target, each Gaia neighbour).
 6 calibration: 300 model-free injections (real background-subtracted pixels within 3 px of the target multiplied by
   1 - depth during fake transits), fitted identically -> mean and covariance of the recovered offset -> z of the real
   offset; plus PRF-shaped injections on each neighbour able to host the signal (aperture flux loss matched).
 7 S6: momentum-dump times (quality bit 32) are listed against the transit times.
Outputs: out/a2_localisation.json, out/a2_diffimages.png
"""
import glob, json, os, sys, warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PRF import TESS_PRF
warnings.filterwarnings("ignore")

H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "out"); os.makedirs(OUT, exist_ok=True)
ADD = os.path.abspath(os.path.join(H, "..", ".."))
AUDD = os.path.join(ADD, "tic4206066_PUBLICATION_AUDIT", "data")     # raw downloads only (TESScut, PRF, Gaia census)
LC = os.path.join(ADD, "deep_followup_tic4206066", "lc")
TARGET = "3220388198192519424"; PIX = 21.0
LK_DEFAULT = 17087; STRAY = 2048
TESSCUT_SHIFT_S = {6: -11.78, 32: -3.28, 98: 0.49}
EPH = {"A": dict(P=3.1827853, T0=3998.73679, T14=1.874 / 24), "B": dict(P=11.132745, T0=4000.86587, T14=1.66 / 24)}
SETS = {"S6_ffi": dict(path=glob.glob(os.path.join(AUDD, "tesscut", "tess-s0006-*_astrocut.fits"))[0], kind="ffi", sector=6, cam=1, ccd=1),
        "S32_tpf": dict(path=os.path.join(LC, "tess2020324010417-s0032-0000000004206066-0200-s_tp.fits"), kind="tpf", sector=32, cam=1, ccd=2),
        "S32_ffi": dict(path=glob.glob(os.path.join(AUDD, "tesscut", "tess-s0032-*_astrocut.fits"))[0], kind="ffi", sector=32, cam=1, ccd=2),
        "S98_ffi": dict(path=glob.glob(os.path.join(AUDD, "tesscut", "tess-s0098-*_astrocut.fits"))[0], kind="ffi", sector=98, cam=1, ccd=2)}
rng = np.random.default_rng(2026)
CLEAN98 = os.environ.get("CLEAN98", "1") == "1"   # default: drop QLP-2^29 (scattered-light) cadences in S98

def in_tr(t, key, fac=1.0):
    e = EPH[key]; ph = ((t - e["T0"]) / e["P"] + 0.5) % 1 - 0.5
    return np.abs(ph * e["P"]) < fac * e["T14"]

def tmag(G, bprp):
    if not np.isfinite(bprp): return G - 0.43
    c = bprp; return G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473   # Stassun+2019 eq. 1

def gaia():
    import csv
    rows = list(csv.DictReader(open(os.path.join(AUDD, "gaia", "gaia_dr3_neighbours_180as.csv"))))
    out = []
    for r in rows:
        g = float(r["phot_g_mean_mag"]) if r["phot_g_mean_mag"] else np.nan
        if not np.isfinite(g): continue
        bprp = float(r["bp_rp"]) if r["bp_rp"] else np.nan
        out.append(dict(id=r["source_id"], ra=float(r["ra"]), dec=float(r["dec"]), pmra=float(r["pmra"] or 0), pmdec=float(r["pmdec"] or 0),
                        T=tmag(g, bprp), G=g, sep=float(r["sep_arcsec"])))
    return out

def load(name):
    c = SETS[name]; h = fits.open(c["path"]); d = h[1].data
    t = np.asarray(d["TIME"], float); fl = np.asarray(d["FLUX"], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
    dumps = t[np.isfinite(t) & ((q & 32) != 0)]
    good = np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1) & ((q & LK_DEFAULT) == 0)
    if c["sector"] == 98: good &= (q & STRAY) == 0
    if c["sector"] == 98 and CLEAN98:
        # exclude the Earth-scattered-light intervals that QLP flags with bit 2^29 (not in the TESScut QUALITY column):
        # a cube cadence is dropped when the QLP cadence nearest in time (< 60 s) carries the flag
        dq = fits.open(os.path.join(LC, "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"))[1].data
        tq = np.asarray(dq["TIME"], float); qq = np.asarray(dq["QUALITY"]).astype(np.int64); m = np.isfinite(tq); tq, qq = tq[m], qq[m]
        j = np.clip(np.searchsorted(tq, t), 1, len(tq) - 1); j = np.where(np.abs(tq[j - 1] - t) < np.abs(tq[j] - t), j - 1, j)
        flagged = (np.abs(tq[j] - t) < 60 / 86400) & ((qq[j] & (1 << 29)) != 0)
        good &= ~flagged
    t, fl, q = t[good], fl[good], q[good]
    if c["kind"] == "ffi": t = t + TESSCUT_SHIFT_S[c["sector"]] / 86400
    w = WCS(h[2].header); ny, nx = fl.shape[1:]; yy, xx = np.mgrid[:ny, :nx]
    # background for FFI cubes (TPF FLUX is already background-subtracted)
    if c["kind"] == "ffi":
        med = np.median(fl, axis=0); cen = np.hypot(xx - nx / 2, yy - ny / 2) < 3
        iy, ix = np.unravel_index(np.argmax(np.where(cen, med, -np.inf)), med.shape)
        far = np.hypot(xx - ix, yy - iy) > 5; bp = far & (med < np.percentile(med[far], 30))
        bkg = np.median(fl[:, bp], axis=1); fl = fl - bkg[:, None, None]
    else:
        bkg = np.zeros(len(t))
    col0, row0 = h[1].header["1CRV4P"], h[1].header["2CRV4P"]
    return dict(name=name, t=t, cube=fl, bkg=bkg, wcs=w, col0=col0, row0=row0, cad=float(np.median(np.diff(t))), dumps=dumps, **c)

def scene(C, srcs, shift=(0, 0), Tlim=19.0):
    """stamp positions of Gaia sources at the sector epoch (proper motion propagated from 2016.0)"""
    yr = Time(np.median(C["t"]) + 2457000, format="jd", scale="tdb").jyear; dt = yr - 2016.0
    out = []
    for s in srcs:
        if s["T"] > Tlim: continue
        ra = s["ra"] + s["pmra"] * dt / 3.6e6 / np.cos(np.radians(s["dec"])); de = s["dec"] + s["pmdec"] * dt / 3.6e6
        x, y = C["wcs"].all_world2pix([[ra, de]], 0)[0]
        out.append(dict(s, x=float(x + shift[0]), y=float(y + shift[1]), is_target=s["id"] == TARGET))
    return out

def model(prf, shape, srcs, T_ref):
    img = np.zeros(shape)
    for s in srcs:
        if -3 < s["x"] < shape[1] + 2 and -3 < s["y"] < shape[0] + 2:
            img += 10 ** (-0.4 * (s["T"] - T_ref)) * prf.locate(s["x"], s["y"], shape)
    return img

def astrometry(C, prf, srcs):
    med = np.median(C["cube"], axis=0); shape = med.shape; ny, nx = shape; yy, xx = np.mgrid[:ny, :nx]
    sc = scene(C, srcs); tg = [s for s in sc if s["is_target"]][0]; fit = np.hypot(xx - tg["x"], yy - tg["y"]) <= 5
    def chi(p):
        m = model(prf, shape, scene(C, srcs, p), tg["T"]); A = np.vstack([m[fit], np.ones(fit.sum())]).T
        coef, *_ = np.linalg.lstsq(A, med[fit], rcond=None); return np.sum((med[fit] - A @ coef) ** 2), coef
    best = None
    for x0 in (-0.3, 0, 0.3):
        for y0 in (-0.3, 0, 0.3):
            r = minimize(lambda p: chi(p)[0], [x0, y0], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-7))
            if best is None or r.fun < best.fun: best = r
    _, coef = chi(best.x); sc = scene(C, srcs, best.x)
    return dict(shift=[float(v) for v in best.x], flux=float(coef[0])), sc

def diff_images(C, cube, tcs, T14, allowed):
    t = C["t"]; D = []
    for tc in tcs:
        i_in = allowed & (np.abs(t - tc) < 0.4 * T14); i_o = allowed & (np.abs(t - tc) > 0.75 * T14) & (np.abs(t - tc) < 2.5 * T14)
        if i_in.sum() < 2 or i_o.sum() < 4: continue
        if (t[i_o] < tc).sum() < 2 or (t[i_o] > tc).sum() < 2:   # one-sided baseline: constant instead of linear
            D.append(cube[i_o].mean(axis=0) - cube[i_in].mean(axis=0)); continue
        X = np.vstack([np.ones(i_o.sum()), t[i_o] - tc]).T
        coef, *_ = np.linalg.lstsq(X, cube[i_o].reshape(i_o.sum(), -1), rcond=None)
        D.append(coef[0].reshape(cube.shape[1:]) - cube[i_in].mean(axis=0))
    return np.array(D)

class Fitter:
    def __init__(self, prf, shape, x0, y0, half=3.0, step=0.1, rfit=4.0):
        self.gx = np.arange(x0 - half, x0 + half + 1e-9, step); self.gy = np.arange(y0 - half, y0 + half + 1e-9, step)
        ny, nx = shape; yy, xx = np.mgrid[:ny, :nx]; self.pix = np.hypot(xx - x0, yy - y0) <= rfit; self.shape = shape; self.prf = prf
        self.M = np.array([[prf.locate(x, y, shape)[self.pix] for x in self.gx] for y in self.gy])
    def _lin(self, m, d, w):
        Sw, Sd, Sm = w.sum(), (w * d).sum(), (m * w).sum(-1); Smm, Smd = (m * m * w).sum(-1), (m * (w * d)).sum(-1)
        det = Smm * Sw - Sm ** 2; A = (Smd * Sw - Sm * Sd) / det; c = (Sd - A * Sm) / Sw
        chi2 = (w * d * d).sum() - 2 * A * Smd - 2 * c * Sd + A * A * Smm + 2 * A * c * Sm + c * c * Sw
        return A, c, chi2
    def free(self, D, sig):
        d = D[self.pix]; w = 1 / sig[self.pix] ** 2; A, c, chi2 = self._lin(self.M, d, w)
        j = np.unravel_index(np.argmin(chi2), chi2.shape); x, y = self.gx[j[1]], self.gy[j[0]]
        # parabolic refinement inside the grid
        def refine(g, k, axis):
            if 0 < k < len(g) - 1:
                c0, c1, c2 = (chi2[j[0], k - 1], chi2[j[0], k], chi2[j[0], k + 1]) if axis == 1 else (chi2[k - 1, j[1]], chi2[k, j[1]], chi2[k + 1, j[1]])
                den = c0 - 2 * c1 + c2
                return g[k] + (0.5 * (c0 - c2) / den) * (g[1] - g[0]) if den > 0 else g[k]
            return g[k]
        return dict(x=float(refine(self.gx, j[1], 1)), y=float(refine(self.gy, j[0], 0)), A=float(A[j]), chi2=float(chi2[j]), edge=bool(j[0] in (0, len(self.gy) - 1) or j[1] in (0, len(self.gx) - 1)))
    def fixed(self, D, sig, x, y):
        d = D[self.pix]; w = 1 / sig[self.pix] ** 2; m = self.prf.locate(x, y, self.shape)[self.pix]
        A, c, chi2 = self._lin(m[None, None, :], d, w); return dict(A=float(A[0, 0]), chi2=float(chi2[0, 0]))

def fake_sets(C, n_ep, nsets, T14):
    t = C["t"]; bad = in_tr(t, "A", 3.0) | in_tr(t, "B", 3.0); cand = t[~bad]; sets = []; tries = 0
    while len(sets) < nsets and tries < nsets * 300:
        tries += 1; cs = []
        for _ in range(n_ep * 30):
            c = rng.choice(cand); win = np.abs(t - c) < 2.6 * T14
            if (np.abs(t - c) < 0.4 * T14).sum() < 2 or bad[win].any() or any(abs(c - x) < 6 * T14 for x in cs): continue
            if ((t - c < -0.75 * T14) & (t - c > -2.5 * T14)).sum() < 2 or ((t - c > 0.75 * T14) & (t - c < 2.5 * T14)).sum() < 2: continue
            cs.append(c)
            if len(cs) == n_ep: break
        if len(cs) == n_ep: sets.append(cs)
    return sets

def analyse(name, sig, n_noise=300, n_inj=300):
    C = load(name); srcs = gaia(); e = EPH[sig]; T14 = e["T14"]; other = "B" if sig == "A" else "A"
    sc0 = scene(C, srcs); tg0 = [s for s in sc0 if s["is_target"]][0]
    prf = TESS_PRF(C["cam"], C["ccd"], C["sector"], C["col0"] + tg0["x"], C["row0"] + tg0["y"], localdatadir=os.path.join(AUDD, "prf"))
    astro, sc = astrometry(C, prf, srcs); tg = [s for s in sc if s["is_target"]][0]
    t = C["t"]; allowed = ~in_tr(t, other, 1.0)
    n0 = int(np.floor((t.min() - e["T0"]) / e["P"])); n1 = int(np.ceil((t.max() - e["T0"]) / e["P"]))
    eps = []
    for n in range(n0, n1 + 1):
        tc = e["T0"] + n * e["P"]
        nl = (allowed & (t - tc < -0.75 * T14) & (t - tc > -2.5 * T14)).sum(); nr = (allowed & (t - tc > 0.75 * T14) & (t - tc < 2.5 * T14)).sum()
        if (allowed & (np.abs(t - tc) < 0.4 * T14)).sum() >= 2 and nl + nr >= 4 and max(nl, nr) >= 2:
            eps.append((n, tc))
    if not eps: return dict(name=name, signal=sig, n_epochs=0)
    tcs = [x[1] for x in eps]; Dreal = diff_images(C, C["cube"], tcs, T14, allowed); Dm = Dreal.mean(axis=0)
    shape = Dm.shape; F = Fitter(prf, shape, tg["x"], tg["y"])
    fsets = fake_sets(C, len(Dreal), n_noise, T14)
    fakes = np.array([diff_images(C, C["cube"], s, T14, np.ones(len(t), bool)).mean(axis=0) for s in fsets])
    sig_map = np.std(fakes, axis=0); sig_map = np.where(sig_map > 0, sig_map, np.median(sig_map))
    free = F.free(Dm, sig_map); fixt = F.fixed(Dm, sig_map, tg["x"], tg["y"]); depth = fixt["A"] / astro["flux"]
    # per-epoch free fits (bootstrap over transits as a second error estimate)
    per_ep = [F.free(Di, sig_map * np.sqrt(len(Dreal))) for Di in Dreal]
    # injections on the target (model-free)
    isets = fake_sets(C, len(Dreal), n_inj, T14); ny, nx = shape; yy, xx = np.mgrid[:ny, :nx]; inside = np.hypot(xx - tg["x"], yy - tg["y"]) <= 3
    inj = []
    for s in isets:
        cube_i = C["cube"].copy(); box = np.zeros(len(t), bool)
        for c in s: box |= np.abs(t - c) < T14 / 2
        for k in np.where(box)[0]: cube_i[k][inside] *= (1 - depth)
        Di = diff_images(C, cube_i, s, T14, np.ones(len(t), bool)).mean(axis=0); fr = F.free(Di, sig_map); ft = F.fixed(Di, sig_map, tg["x"], tg["y"])
        inj.append(dict(dx=fr["x"] - tg["x"], dy=fr["y"] - tg["y"], A=ft["A"], chi2_t=ft["chi2"], Di=Di if len(inj) < 150 else None))
    dxy = np.array([[r["dx"], r["dy"]] for r in inj]); mu = dxy.mean(0); cov = np.cov(dxy.T)
    off = np.array([free["x"] - tg["x"], free["y"] - tg["y"]]); z = (off - mu) / np.sqrt(np.diag(cov))
    maha = float(np.sqrt((off - mu) @ np.linalg.inv(cov) @ (off - mu)))
    # neighbours able to host the signal (full-inclusion flux ratio) and their fits / injections
    neigh = []
    for s in sc:
        if s["is_target"]: continue
        fr_ = 10 ** (-0.4 * (s["T"] - tg["T"])); req = depth * (1 + fr_) / fr_; dist = np.hypot(s["x"] - tg["x"], s["y"] - tg["y"]) * PIX
        if dist > 130: continue
        fx = F.fixed(Dm, sig_map, s["x"], s["y"])
        rec = dict(id=s["id"], T=round(s["T"], 2), sep_arcsec=round(dist, 1), req_eclipse=float(req), dchi2_vs_target=fx["chi2"] - fixt["chi2"])
        if req <= 1.0 and len([n for n in neigh if n.get("inj")]) < 6:
            m_nb = prf.locate(s["x"], s["y"], shape) * fr_ * astro["flux"]; m_tg = prf.locate(tg["x"], tg["y"], shape) * astro["flux"]
            ap = np.hypot(xx - tg["x"], yy - tg["y"]) <= 2.0; ecl = depth * (m_tg[ap].sum() + m_nb[ap].sum()) / m_nb[ap].sum()
            dch = []; offs = []
            for sset in isets[:150]:
                cube_i = C["cube"].copy(); box = np.zeros(len(t), bool)
                for c in sset: box |= np.abs(t - c) < T14 / 2
                cube_i[box] -= ecl * m_nb
                Di = diff_images(C, cube_i, sset, T14, np.ones(len(t), bool)).mean(axis=0)
                ft = F.fixed(Di, sig_map, tg["x"], tg["y"]); fn = F.fixed(Di, sig_map, s["x"], s["y"]); fr2 = F.free(Di, sig_map)
                dch.append(fn["chi2"] - ft["chi2"]); offs.append(np.hypot(fr2["x"] - tg["x"], fr2["y"] - tg["y"]) * PIX)
            dch = np.array(dch)
            rec["inj"] = dict(eclipse_needed_in_aperture=float(ecl), dchi2_neighbour_inj_median=float(np.median(dch)), dchi2_neighbour_inj_min=float(dch.min()),
                              frac_neighbour_inj_dchi2_le_real=float(np.mean(dch <= rec["dchi2_vs_target"])), median_offset_arcsec_if_on_neighbour=float(np.median(offs)))
        neigh.append(rec)
    neigh.sort(key=lambda r: r["sep_arcsec"])
    # sky offsets (arcsec, east positive) from the WCS
    p0 = C["wcs"].all_pix2world([[tg["x"] - astro["shift"][0], tg["y"] - astro["shift"][1]]], 0)[0]
    p1 = C["wcs"].all_pix2world([[free["x"] - astro["shift"][0], free["y"] - astro["shift"][1]]], 0)[0]
    dra = (p1[0] - p0[0]) * np.cos(np.radians(p0[1])) * 3600; dde = (p1[1] - p0[1]) * 3600
    # Jacobian for the covariance in sky coordinates
    J = np.zeros((2, 2))
    for k, dv in enumerate(((0.5, 0), (0, 0.5))):
        pp = C["wcs"].all_pix2world([[tg["x"] - astro["shift"][0] + dv[0], tg["y"] - astro["shift"][1] + dv[1]]], 0)[0]
        J[:, k] = [(pp[0] - p0[0]) * np.cos(np.radians(p0[1])) * 3600 / 0.5, (pp[1] - p0[1]) * 3600 / 0.5]
    sky_cov = J @ cov @ J.T; sky_mu = J @ mu
    dumps = [float(x) for x in C["dumps"]]; dump_dist = [float(np.min(np.abs(np.array(dumps) - tc)) * 24) if dumps else None for tc in tcs]
    res = dict(name=name, signal=sig, sector=C["sector"], n_epochs=len(Dreal), epochs=[int(x[0]) for x in eps], cadence_min=C["cad"] * 1440,
               astrometric_shift_px=astro["shift"], target_stamp_xy=[tg["x"], tg["y"]], depth_ppm_prf=depth * 1e6,
               free_fit=dict(dx_px=float(off[0]), dy_px=float(off[1]), dRA_arcsec=float(dra), dDec_arcsec=float(dde), offset_arcsec=float(np.hypot(*off) * PIX), chi2=free["chi2"], edge=free["edge"]),
               target_fixed_chi2=fixt["chi2"], n_fit_pixels=int(F.pix.sum()), n_noise_sets=len(fakes),
               injections=dict(n=len(inj), mean_dx_px=float(mu[0]), mean_dy_px=float(mu[1]), sd_dx_px=float(np.sqrt(cov[0, 0])), sd_dy_px=float(np.sqrt(cov[1, 1])),
                               z_x=float(z[0]), z_y=float(z[1]), mahalanobis=maha, frac_injections_further=float(np.mean(np.hypot(dxy[:, 0], dxy[:, 1]) >= np.hypot(*off))),
                               median_recovered_depth_ppm=float(np.median([r["A"] for r in inj]) / astro["flux"] * 1e6)),
               sky=dict(real_minus_injmean_arcsec=[float(dra - sky_mu[0]), float(dde - sky_mu[1])], cov_arcsec2=sky_cov.tolist()),
               per_epoch_offsets_arcsec=[[float((p["x"] - tg["x"]) * PIX), float((p["y"] - tg["y"]) * PIX)] for p in per_ep],
               bootstrap_over_epochs_sd_arcsec=[float(np.std([p["x"] for p in per_ep]) / np.sqrt(len(per_ep)) * PIX), float(np.std([p["y"] for p in per_ep]) / np.sqrt(len(per_ep)) * PIX)],
               momentum_dumps_btjd=dumps, transit_to_nearest_dump_h=dump_dist, neighbours=neigh[:12])
    return res, Dm, sig_map, tg

if __name__ == "__main__":
    which = sys.argv[1:] or list(SETS)
    out = {}; panels = []
    for name in which:
        for sig in ("A", "B"):
            r = analyse(name, sig)
            if isinstance(r, tuple):
                r, Dm, sig_map, tg = r; panels.append((name, sig, Dm, sig_map, tg, r))
            out[f"{name}_{sig}"] = r
            if "free_fit" in r:
                ff, ij = r["free_fit"], r["injections"]
                print(f"{name:8s} {sig}: n_ep={r['n_epochs']} depth={r['depth_ppm_prf']:.0f} ppm; offset {ff['offset_arcsec']:.1f}\" (dRA {ff['dRA_arcsec']:+.1f}, dDec {ff['dDec_arcsec']:+.1f}); "
                      f"inj sd ({ij['sd_dx_px']*PIX:.1f}, {ij['sd_dy_px']*PIX:.1f})\" z=({ij['z_x']:+.2f}, {ij['z_y']:+.2f}) maha {ij['mahalanobis']:.2f} frac_further {ij['frac_injections_further']:.2f}; "
                      f"dump dist h: {[round(x,1) if x else None for x in r['transit_to_nearest_dump_h']]}", flush=True)
                for nb in r["neighbours"]:
                    if nb.get("inj"): print(f"     {nb['sep_arcsec']}\" T={nb['T']} needs {nb['inj']['eclipse_needed_in_aperture']*100:.0f}% ecl; real dchi2(N-T)={nb['dchi2_vs_target']:.1f}; nb-inj dchi2 median {nb['inj']['dchi2_neighbour_inj_median']:.1f} min {nb['inj']['dchi2_neighbour_inj_min']:.1f}", flush=True)
            else:
                print(name, sig, "no epochs", flush=True)
            json.dump(out, open(os.path.join(OUT, f"a2_localisation_{'_'.join(which)}{'' if CLEAN98 else '_all98'}.json"), "w"), indent=1, default=float)
    if panels:
        fig, axs = plt.subplots(2, len(panels), figsize=(2.6 * len(panels), 5.2), squeeze=False)
        for k, (name, sig, Dm, sm, tg, r) in enumerate(panels):
            ax = axs[0, k]; ax.imshow(Dm, origin="lower", cmap="viridis"); ax.plot(tg["x"], tg["y"], "r+", ms=10)
            ax.plot(r["free_fit"]["dx_px"] + tg["x"], r["free_fit"]["dy_px"] + tg["y"], "wx", ms=8); ax.set_title(f"{name} {sig}: diff", fontsize=8)
            ax = axs[1, k]; ax.imshow(Dm / sm, origin="lower", cmap="RdBu_r", vmin=-4, vmax=4); ax.plot(tg["x"], tg["y"], "k+", ms=10); ax.set_title("diff / sigma", fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, f"a2_diffimages_{'_'.join(which)}{'' if CLEAN98 else '_all98'}.png"), dpi=130)
    print("done")
