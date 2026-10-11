"""Step 5: independent aperture photometry from the TESScut SPOC-FFI cubes (S6, S32, S98) and from the SPOC 2-min
target pixel file (S32), using one identical method in every sector.

Method (all sectors):
  * cadence mask: lightkurve TessQualityFlags.DEFAULT_BITMASK (17087 = attitude tweak, safe mode, coarse point, earth
    point, argabrightening, desat, manual exclude, impulsive outlier, bad calibration). Straylight (2048) is NOT
    removed here but carried as a column (TESS SDPDD: informational flag), so its effect can be tested later.
  * target pixel position: Gaia DR3 position propagated with proper motion to the mid-sector epoch -> WCS of the
    cutout; compared with the flux-weighted centroid of the median image (WCS sanity check).
  * background (FFI cutouts only; TESScut FFI cutouts are NOT background subtracted): per-cadence median of pixels
    whose time-median flux is in the lowest 30 % AND that are > 5 px from the target.  SPOC TPF FLUX is already
    background-subtracted (SDPDD), so no extra background is removed there.
  * apertures: circles of radius 1.0, 1.5, 2.0, 2.5, 3.0 px centred on the target position (pixel centre inside).
Output: out/phot/<name>.npz with time, fluxes per aperture, bkg, centroids, quality, straylight flag.
"""
import csv
import json
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, save_json  # noqa

DEFAULT_BITMASK = 17087
STRAY = 2048
RADII = [1.0, 1.5, 2.0, 2.5, 3.0]
os.makedirs(os.path.join(AUD, "out", "phot"), exist_ok=True)

g = list(csv.DictReader(open(os.path.join(AUD, "data/gaia/gaia_dr3_source.csv"))))[0]
RA0, DE0 = float(g["ra"]), float(g["dec"])
PMRA, PMDE = float(g["pmra"]), float(g["pmdec"])


def radec_at(btjd):
    yr = Time(btjd + 2457000.0, format="jd", scale="tdb").jyear
    dt = yr - 2016.0
    return RA0 + PMRA * dt / 3.6e6 / np.cos(np.radians(DE0)), DE0 + PMDE * dt / 3.6e6, yr


def phot(name, path, kind):
    h = fits.open(path)
    d = h[1].data
    t = np.asarray(d["TIME"], float)
    fl = np.asarray(d["FLUX"], float)
    fe = np.asarray(d["FLUX_ERR"], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    ok = np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1) & ((q & DEFAULT_BITMASK) == 0)
    t, fl, fe, q = t[ok], fl[ok], fe[ok], q[ok]
    w = WCS(h[2].header)
    ny, nx = fl.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    tmid = float(np.median(t))
    ra, de, yr = radec_at(tmid)
    tx, ty = w.all_world2pix([[ra, de]], 0)[0]
    med = np.nanmedian(fl, axis=0)
    info = dict(name=name, file=os.path.relpath(path, AUD), kind=kind, n_cad=int(ok.sum()), n_rejected=int((~ok).sum()),
                epoch_jyear=yr, target_radec_epoch=[ra, de], target_px_wcs=[float(tx), float(ty)])
    if kind == "ffi":
        far = np.hypot(xx - tx, yy - ty) > 5
        faint = med < np.nanpercentile(med[far], 30)
        bpix = far & faint
        bkg = np.nanmedian(fl[:, bpix], axis=1)
        info["n_bkg_pixels"] = int(bpix.sum())
    else:
        bkg = np.zeros(len(t))
        bpix = None
    flb = fl - bkg[:, None, None]
    # flux-weighted centroid of the median background-subtracted image in a 5x5 box around the WCS position
    medb = np.nanmedian(flb, axis=0)
    box = (np.abs(xx - np.round(tx)) <= 2) & (np.abs(yy - np.round(ty)) <= 2)
    im = np.where(box, np.clip(medb, 0, None), 0)
    cx, cy = (im * xx).sum() / im.sum(), (im * yy).sum() / im.sum()
    info["median_image_centroid_px"] = [float(cx), float(cy)]
    info["wcs_minus_centroid_px"] = [float(tx - cx), float(ty - cy)]
    fluxes, errs, npix = {}, {}, {}
    for r in RADII:
        ap = np.hypot(xx - tx, yy - ty) <= r
        fluxes[r] = flb[:, ap].sum(axis=1)
        errs[r] = np.sqrt((fe[:, ap] ** 2).sum(axis=1))
        npix[r] = int(ap.sum())
    # per-cadence flux-weighted centroid (5x5 box) for centroid-motion tests
    bb = np.where(box[None], np.clip(flb, 0, None), 0)
    s = bb.sum(axis=(1, 2))
    ccx = (bb * xx[None]).sum(axis=(1, 2)) / s
    ccy = (bb * yy[None]).sum(axis=(1, 2)) / s
    pc1 = np.asarray(d["POS_CORR1"], float)[ok] if "POS_CORR1" in d.columns.names else np.zeros(len(t))
    pc2 = np.asarray(d["POS_CORR2"], float)[ok] if "POS_CORR2" in d.columns.names else np.zeros(len(t))
    tc = np.asarray(d["TIMECORR"], float)[ok] if "TIMECORR" in d.columns.names else np.zeros(len(t))
    cad = np.asarray(d["CADENCENO"])[ok]
    out = dict(time=t, cadenceno=cad, timecorr=tc, bkg=bkg, cx=ccx, cy=ccy, pos_corr1=pc1, pos_corr2=pc2, quality=q,
               straylight=((q & STRAY) > 0))
    for r in RADII:
        out[f"flux_r{r}"] = fluxes[r]
        out[f"err_r{r}"] = errs[r]
    np.savez_compressed(os.path.join(AUD, "out", "phot", f"{name}.npz"), **out)
    # simple precision metric per aperture: robust std of 1-h binned, 0.5-d-median-normalised flux
    prec = {}
    for r in RADII:
        f = fluxes[r] / np.nanmedian(fluxes[r])
        # remove slow trend with a running median of 0.5 d
        order = np.argsort(t)
        tt, ff = t[order], f[order]
        trend = np.array([np.median(ff[np.abs(tt - x) < 0.25]) for x in tt[:: max(1, len(tt) // 2000)]])
        tr = np.interp(tt, tt[:: max(1, len(tt) // 2000)], trend)
        res = ff / tr - 1
        nb = max(1, int(round(1 / 24 / np.median(np.diff(tt)))))
        binned = np.array([res[i:i + nb].mean() for i in range(0, len(res) - nb, nb)])
        prec[str(r)] = dict(npix=npix[r], rms_1h_ppm=float(1.4826 * np.median(np.abs(binned - np.median(binned))) * 1e6),
                            median_flux=float(np.nanmedian(fluxes[r])))
    info["aperture_precision"] = prec
    info["straylight_frac"] = float(np.mean(out["straylight"]))
    info["cadence_min"] = float(np.median(np.diff(t)) * 1440)
    return info


res = []
import glob  # noqa

for s in (6, 32, 98):
    p = glob.glob(os.path.join(AUD, f"data/tesscut/tess-s{s:04d}-*_astrocut.fits"))[0]
    res.append(phot(f"tesscut_s{s}", p, "ffi"))
tpf = glob.glob(os.path.join(AUD, "data/mast/tess2020324010417-s0032-*_tp.fits"))[0]
res.append(phot("spoc_tpf_s32", tpf, "tpf"))
tsp = glob.glob(os.path.join(AUD, "data/mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits"))[0]
res.append(phot("tessspoc_tp_s6", tsp, "tpf"))
tsp32 = glob.glob(os.path.join(AUD, "data/mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_tp.fits"))[0]
res.append(phot("tessspoc_tp_s32", tsp32, "tpf"))
save_json(res, "out/s05_photometry.json")
for r in res:
    print(r["name"], r["n_cad"], "cad", round(r["cadence_min"], 2), "min; WCS-centroid px", [round(x, 3) for x in r["wcs_minus_centroid_px"]],
          "stray", round(r["straylight_frac"], 3))
    print("   ", {k: (v["npix"], round(v["rms_1h_ppm"])) for k, v in r["aperture_precision"].items()})
