"""Shared inputs for the TIC 4206066 pixel re-check. Read-only on project data; writes only inside this package."""
import hashlib, json, os
import numpy as np
from astropy.io import fits

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(PKG))            # tess-jev
OUT, PLOTS = os.path.join(PKG, "out"), os.path.join(PKG, "plots")

# Ephemeris: publication audit (addenda/tic4206066_PUBLICATION_AUDIT/out/ephem_current.json), BTJD.
EPHEM_FILE = "addenda/tic4206066_PUBLICATION_AUDIT/out/ephem_current.json"
_e = json.load(open(os.path.join(ROOT, EPHEM_FILE)))
EPH = {s: (_e[s]["P"], _e[s]["T0"], _e[s]["T14"]) for s in ("A", "B")}

DATA = {
    "S6_ffi":      "addenda/deep_followup_tic4206066/pixels/tess-s0006-1-1_80.859830_-1.326960_21x21_astrocut.fits",
    "S32_ffi":     "addenda/deep_followup_tic4206066/pixels/tess-s0032-1-2_80.859830_-1.326960_21x21_astrocut.fits",
    "S32_tpf2min": "addenda/deep_followup_tic4206066/lc/tess2020324010417-s0032-0000000004206066-0200-s_tp.fits",
    "S98_ffi":     "addenda/deep_followup_tic4206066/pixels/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits",
}
RADII = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def phase_dist(t, P, T0):
    return ((t - T0 + 0.5 * P) % P) - 0.5 * P


QLP_S98 = "addenda/deep_followup_tic4206066/lc/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"
TESSCUT_SHIFT_S = {"S6_ffi": -11.78, "S32_ffi": -3.28, "S98_ffi": 0.49}   # as in rnaas_tic4206066_pubprep/analysis/a2_prf_localisation.py L40
LK_DEFAULT, STRAY = 17087, 2048                                            # same cadence selection as pubprep a2 (L39, L73-L83)


def load(name, raw=False):
    """Cadence selection as in pubprep a2: SPOC default bitmask removed; in S98 also stray-light bit 2048 and the
    QLP 2^29 Earth-scattered-light intervals. TESScut cubes are NOT background-subtracted (BACKAPP = 0, FLUX_BKG = 0):
    a per-cadence background (median of faint pixels > 5 px from the target, below the 30th percentile of the median
    image, as in pubprep a2) is subtracted and returned. SPOC TPF FLUX is already background-subtracted; its
    FLUX_BKG summed over the stamp is returned as the background regressor.
    raw=True reproduces the first-pass selection (QUALITY == 0, no background subtraction) used by s02."""
    d = fits.open(os.path.join(ROOT, DATA[name]))[1].data
    t = d["TIME"].astype(float); F = d["FLUX"].astype(float); q = d["QUALITY"].astype(np.int64)
    B = d["FLUX_BKG"].astype(float)
    fin = np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
    if raw:
        good = fin & (q == 0)
    else:
        good = fin & ((q & LK_DEFAULT) == 0)
        if name.startswith("S98"):
            good &= (q & STRAY) == 0
            dq = fits.open(os.path.join(ROOT, QLP_S98))[1].data
            tq = np.asarray(dq["TIME"], float); qq = np.asarray(dq["QUALITY"]).astype(np.int64); m = np.isfinite(tq); tq, qq = tq[m], qq[m]
            j = np.clip(np.searchsorted(tq, t), 1, len(tq) - 1); j = np.where(np.abs(tq[j-1] - t) < np.abs(tq[j] - t), j - 1, j)
            good &= ~((np.abs(tq[j] - t) < 60 / 86400) & ((qq[j] & (1 << 29)) != 0))
    t, F, B = t[good], F[good], B[good]
    med = np.nanmedian(F, 0); ny, nx = med.shape
    cy, cx = np.unravel_index(np.nanargmax(med[ny//2-2:ny//2+3, nx//2-2:nx//2+3]), (5, 5))
    c = (int(cy + ny//2 - 2), int(cx + nx//2 - 2))
    if raw:
        return t, F, np.nansum(B.reshape(len(t), -1), 1), c, med
    if name.endswith("_ffi"):
        t = t + TESSCUT_SHIFT_S[name] / 86400
        yy, xx = np.mgrid[0:ny, 0:nx]; far = np.hypot(yy - c[0], xx - c[1]) > 5
        bp = far & (med < np.percentile(med[far], 30))
        bkg = np.median(F[:, bp], axis=1); F = F - bkg[:, None, None]
        med = np.nanmedian(F, 0)
    else:
        bkg = np.nansum(B.reshape(len(t), -1), 1)
    return t, F, bkg, c, med


def radii_grid(shape, c):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.hypot(yy - c[0], xx - c[1])


def centroid(F, c):
    sub = F[:, c[0]-2:c[0]+3, c[1]-2:c[1]+3]; tot = sub.sum((1, 2)); gy, gx = np.mgrid[-2:3, -2:3]
    return (sub * gx).sum((1, 2)) / tot, (sub * gy).sum((1, 2)) / tot
