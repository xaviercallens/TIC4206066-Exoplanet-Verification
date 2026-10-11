"""T8 check 4: what aperture dependence would an off-target source produce in the package's apertures?
Uses the official TESS PRF (TESS_PRF, MAST start_s0004 files already in opus55 00_data) through the opus55 audit's
read-only PRFModel wrapper, Gaia DR3 + TIC 8.2 cone files from opus55 00_data, and the package's own aperture centres.
Hypotheses: (T) signal on the target; (N1) signal on TIC 4206077 (Gaia DR3 3220388296975320448, 43.3", nearest star
able to produce A); (N2) signal on TIC 4206057 (84", T = 13.06).
Output: fractional depth ratio d(r)/d(2 px) and predicted r1 - r2 for a 470 ppm r = 2 depth (A) or 900 ppm (B).
Writes only to OUT/_work/T8_pixel_recheck/out/."""
import json, os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
W = os.path.dirname(os.path.abspath(__file__))
U3 = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels"
sys.path.insert(0, U3)
import common as U   # opus55 helper, read-only (PYTHONDONTWRITEBYTECODE=1)
sys.path.remove(U3); del sys.modules["common"]
sys.path.insert(0, os.path.join(W, "pkg_copy", "scripts"))
import common as PK   # package copy
ODIR = os.path.join(W, "out")

g = U.load_gaia(180.0)
tgt = g[g["source_id"] == U.TARGET_GAIA].iloc[0]
N1, N2 = 3220388296975320448, 3220385243255181184
res = {"prf_files": U.PRFDIR, "gaia_rows": int(len(g))}
for name in ("S98_ffi", "S32_ffi", "S32_tpf2min", "S6_ffi"):
    path = os.path.join(PK.ROOT, PK.DATA[name])
    with fits.open(path) as h:
        h0, h1, h2 = h[0].header, h[1].header, h[2].header
        tmid = np.nanmedian(h[1].data["TIME"]) + 2457000.0
    sector, cam, ccd = int(h0["SECTOR"]), int(h0["CAMERA"]), int(h0["CCD"])
    col0, row0 = int(h1["1CRV4P"]), int(h1["2CRV4P"]); wcs = WCS(h2)
    jyear = 2000.0 + (tmid - 2451545.0) / 365.25
    ra, dec = U.propagate(g, jyear)
    x, y = wcs.all_world2pix(ra, dec, 0)
    _, _, _, c, med = PK.load(name)
    ny, nx = med.shape
    it = int(np.where(g["source_id"].values == U.TARGET_GAIA)[0][0])
    prf = U.PRFModel(cam, ccd, sector, col0 + x[it], row0 + y[it])
    Tt = float(tgt["Tmag_used"])
    imgs = {}
    total = np.zeros((ny, nx))
    for k in range(len(g)):
        if not (-7 < x[k] < nx + 6 and -7 < y[k] < ny + 6): continue
        im = prf.image(x[k], y[k], (ny, nx)) * 10 ** (-0.4 * (float(g["Tmag_used"].iloc[k]) - Tt))
        total += im
        sid = int(g["source_id"].iloc[k])
        if sid in (U.TARGET_GAIA, N1, N2): imgs[sid] = im
    yy, xx = np.mgrid[0:ny, 0:nx]; rr = np.hypot(yy - c[0], xx - c[1])
    out = {"sector": sector, "aperture_centre_pix_rowcol": list(c), "target_xy_pix": [float(x[it]), float(y[it])]}
    for lab, sid in (("T_target", U.TARGET_GAIA), ("N1_TIC4206077_43arcsec", N1), ("N2_TIC4206057_84arcsec", N2)):
        if sid not in imgs: out[lab] = "outside stamp"; continue
        k = int(np.where(g["source_id"].values == sid)[0][0])
        dE, dN = U.offsets_arcsec(float(ra[it]), float(dec[it]), ra[k], dec[k])
        frac = {}
        for r in (0.5, 1.0, 2.0, 3.0):
            ap = rr <= r
            frac[str(r)] = float(imgs[sid][ap].sum() / total[ap].sum())   # depth per unit eclipse depth of that star
        ratio = {r: frac[r] / frac["2.0"] for r in frac}
        out[lab] = {"sep_arcsec": float(np.hypot(dE[0], dN[0])), "dE_dN_arcsec": [float(dE[0]), float(dN[0])],
                    "Tmag": float(g["Tmag_used"].iloc[k]),
                    "flux_fraction_of_star_in_ap_over_total": frac, "depth_ratio_vs_r2": ratio,
                    "pred_r1_minus_r2_ppm_if_r2depth_470": float(470 * (ratio["1.0"] - 1)),
                    "pred_r0.5_minus_r2_ppm_if_r2depth_470": float(470 * (ratio["0.5"] - 1)),
                    "pred_r1_minus_r3_ppm_if_r2depth_900": float(900 * (ratio["1.0"] - ratio["3.0"])),
                    "pred_r2_minus_r3_ppm_if_r2depth_900": float(900 * (1 - ratio["3.0"])),
                    "eclipse_depth_needed_for_470ppm_in_r2": float(470e-6 / frac["2.0"])}
    res[name] = out
    print(name, json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk in ("sep_arcsec", "depth_ratio_vs_r2", "pred_r1_minus_r2_ppm_if_r2depth_470", "eclipse_depth_needed_for_470ppm_in_r2")}) for k, v in out.items()}, default=float))
json.dump(res, open(os.path.join(ODIR, "t8_prf_offtarget_signature.json"), "w"), indent=1, default=float)
