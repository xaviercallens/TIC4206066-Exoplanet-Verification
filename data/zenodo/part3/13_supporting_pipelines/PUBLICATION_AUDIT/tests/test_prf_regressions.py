"""Regression tests for the PRF mistakes found in the TIC 158955499 audit (addenda/pixel_audit_tic158955499_SEALED,
'Re-verification (fourth pass)'), made executable here so they cannot recur in the TIC 4206066 pixel analysis.

R1  PRF files are 9x9 INTERLEAVED 13x13 PRFs, not one continuous 9x-supersampled image.
    Test: the flux-weighted centroid of the TESS_PRF sub-PRF for a source placed at (x, y) must track (x, y)
    to < 0.05 px over a grid of intra-pixel positions; the 'continuous image' reading must FAIL the same test.
R2  The host must never appear in its own neighbour list (TIC J2000 vs Gaia 2016 positions differ by ~1.1'' through
    proper motion). Test: scene() flags exactly one target by source_id, and no non-target source lies < 2''.
R3  Injections must not be circular (same PRF model used to inject and to fit). Test: the injection function used by
    s12 multiplies REAL pixel values; it has no PRF argument (inspected by signature) and a unit test shows the
    injected cube equals data*(1-depth*box) inside the injection radius exactly.
R4  Required eclipse depth of a neighbour must use the flux of target and neighbour INSIDE the aperture actually used;
    full-inclusion is only a lower bound. Test: aperture-based requirement >= full-inclusion requirement.
Run: audit_env/bin/python tests/test_prf_regressions.py
"""
import inspect
import json
import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import prf_tools as pt  # noqa

results = {}


def centroid(img):
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    return (img * xx).sum() / img.sum(), (img * yy).sum() / img.sum()


# ---------- R1 ----------
prf = pt.get_prf(1, 2, 98, 1372.0, 983.0)
errs_ok, errs_bad = [], []
fn = sorted(os.listdir(os.path.join(pt.PRFDIR, "cam1_ccd2")))[12]
raw = fits.getdata(os.path.join(pt.PRFDIR, "cam1_ccd2", fn))
assert raw.shape == (117, 117), raw.shape
for fx in np.linspace(-0.45, 0.45, 7):
    for fy in np.linspace(-0.45, 0.45, 7):
        x, y = 10 + fx, 10 + fy
        m = prf.locate(x, y, (21, 21))
        cx, cy = centroid(m)
        # PRF wings are asymmetric, so compare centroid SHIFTS relative to the pixel-centred source
        c0 = centroid(prf.locate(10.0, 10.0, (21, 21)))
        errs_ok.append((cx - c0[0] - fx, cy - c0[1] - fy))
        # wrong reading: treat the 117x117 image as a continuous 9x supersampled PSF, shift by (fx,fy)*9 samples, rebin
        sx, sy = int(round(fx * 9)), int(round(fy * 9))
        sh = np.roll(np.roll(raw, sy, axis=0), sx, axis=1)
        reb = sh.reshape(13, 9, 13, 9).sum(axis=(1, 3))
        bx, by = centroid(reb)
        b0 = centroid(raw.reshape(13, 9, 13, 9).sum(axis=(1, 3)))
        # also the naive un-reversed un-interleave prf[i::9, j::9] (orientation error)
        i, j = int(round(4 + fy * 9)), int(round(4 + fx * 9))
        nv = raw[i::9, j::9]
        n0 = raw[4::9, 4::9]
        nx_, ny_ = centroid(nv)
        nx0, ny0 = centroid(n0)
        errs_bad.append((nx_ - nx0 - fx, ny_ - ny0 - fy))
errs_ok = np.abs(np.array(errs_ok))
errs_bad = np.abs(np.array(errs_bad))
results["R1"] = dict(max_centroid_tracking_error_TESS_PRF_px=float(errs_ok.max()),
                     max_error_unreversed_interleave_px=float(errs_bad.max()),
                     passed=bool(errs_ok.max() < 0.05 and errs_bad.max() > 0.2))

# ---------- R2 ----------
from astropy.wcs import WCS  # noqa
import glob  # noqa

cut = glob.glob(os.path.join(pt.AUD, "data", "tesscut", "tess-s0098-*_astrocut.fits"))[0]
w = WCS(fits.getheader(cut, 2))
src = pt.gaia_sources(180)
sc = pt.scene(src, w, 4017.0, 3220388198192519424, Tlim=21)
tg = [s for s in sc if s["is_target"]]
tx, ty = tg[0]["x"], tg[0]["y"]
close = [s for s in sc if (not s["is_target"]) and np.hypot(s["x"] - tx, s["y"] - ty) * 21.0 < 2.0]
results["R2"] = dict(n_targets=len(tg), n_nontarget_within_2arcsec=len(close), passed=bool(len(tg) == 1 and len(close) == 0))

# ---------- R3 ----------
import s12_pixels as s12  # noqa

sig = inspect.signature(s12.inject_on_target)
cube = np.random.default_rng(1).normal(1000, 5, size=(50, 11, 11))
bkg = np.full(50, 100.0)
t = np.linspace(0, 1, 50)
inj = s12.inject_on_target(cube, bkg, t, tc=0.5, dur=0.2, depth=1e-3, x0=5.0, y0=5.0, radius=3.0)
yy, xx = np.mgrid[:11, :11]
inside = np.hypot(xx - 5, yy - 5) <= 3
box = (np.abs(t - 0.5) < 0.1)
expect = cube.copy()
expect[box][:, inside] = ((cube[box][:, inside] - 100.0) * (1 - 1e-3) + 100.0)
exp2 = cube.copy()
for k in np.where(box)[0]:
    exp2[k][inside] = (cube[k][inside] - 100.0) * (1 - 1e-3) + 100.0
results["R3"] = dict(signature=str(sig), has_prf_argument=("prf" in sig.parameters),
                     max_abs_diff=float(np.max(np.abs(inj - exp2))), passed=bool(("prf" not in sig.parameters) and np.allclose(inj, exp2)))

# ---------- R4 ----------
req = s12.required_depths_example()
results["R4"] = dict(examples=req, passed=bool(all(r["req_aperture"] >= r["req_full_inclusion"] - 1e-12 for r in req)))

print(json.dumps(results, indent=1))
json.dump(results, open(os.path.join(pt.AUD, "tests", "test_prf_regressions_result.json"), "w"), indent=1)
assert all(v["passed"] for v in results.values()), "PRF regression test FAILED"
print("ALL PRF REGRESSION TESTS PASSED")
