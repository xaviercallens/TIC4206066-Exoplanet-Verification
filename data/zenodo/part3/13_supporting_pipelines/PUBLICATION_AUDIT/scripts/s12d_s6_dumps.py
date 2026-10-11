"""Step 12d: S6 localisation with the momentum-dump coincidence handled.
In S6 the reaction-wheel desaturations occur every 3.125 d and every A transit falls 0.4-3.5 h from one (P_A = 3.1828 d).
A pointing jump at a dump inside the difference-image window produces a dipole that biases the fitted position.
Here: (i) per-pixel baseline = linear in time + a STEP at the dump time when a dump lies inside the window;
(ii) injections are placed at the SAME offsets from dumps as the real transits (tc_fake = dump_j + (tc_real - dump_near)),
avoiding real A/B transits, so real and injected events see the same pointing history.
Also reports the S6 flux-weighted centroid (TESS-SPOC S6 LC MOM_CENTR1/2) around dumps.
"""
import itertools
import json
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import prf_tools as pt  # noqa
import s12_pixels as s12  # noqa
from audit_common import AUD, save_json  # noqa


def dumps_from_cube():
    import glob
    h = fits.open(glob.glob(os.path.join(AUD, "data/tesscut/tess-s0006-*_astrocut.fits"))[0])
    t = np.asarray(h[1].data["TIME"], float)
    q = np.asarray(h[1].data["QUALITY"]).astype(np.int64)
    d = np.sort(t[((q & 32) > 0) & np.isfinite(t)]) + lcio.TESSCUT_SHIFT_S[6] / 86400
    out = []
    for x in d:
        if not out or x - out[-1] > 0.1:
            out.append(x)
    return np.array(out)


def diff_image_step(C, cube, tcs, T14, okmask, dumps):
    t = C["t"]
    ny, nx = cube.shape[1:]
    D = []
    for tc in tcs:
        i_in = okmask & (np.abs(t - tc) < 0.4 * T14)
        i_o = okmask & (np.abs(t - tc) > 0.75 * T14) & (np.abs(t - tc) < 2.5 * T14)
        if i_in.sum() < 2 or i_o.sum() < 4:
            continue
        cols = [np.ones(i_o.sum()), t[i_o] - tc]
        dd = [x for x in dumps if abs(x - tc) < 2.5 * T14]
        stepcols_in = []
        for x in dd:
            s_o = (t[i_o] > x).astype(float)
            if 0 < s_o.sum() < len(s_o):
                cols.append(s_o)
                stepcols_in.append(float(np.mean(t[i_in] > x)))
        X = np.vstack(cols).T
        coef, *_ = np.linalg.lstsq(X, cube[i_o].reshape(i_o.sum(), -1), rcond=None)
        base = coef[0] + sum(coef[2 + k] * stepcols_in[k] for k in range(len(stepcols_in))) if stepcols_in else coef[0]
        D.append(base.reshape(ny, nx) - cube[i_in].mean(axis=0))
    return np.array(D)


def main(sig="A"):
    C = s12.load_cube("S6_ffi")
    dumps = dumps_from_cube()
    shape = C["cube"].shape[1:]
    sc0 = pt.scene(C["sources"], C["wcs"], C["tmid"], s12.TARGET_ID)
    tg0 = [s for s in sc0 if s["is_target"]][0]
    prf = pt.get_prf(C["cam"], C["ccd"], C["sector"], C["col0"] + tg0["x"], C["row0"] + tg0["y"])
    astro, srcs, _ = s12.astrometric_solution(C, prf)
    tg = [s for s in srcs if s["is_target"]][0]
    e = s12.EPH[sig]
    other = "B" if sig == "A" else "A"
    okc = ~lcio.in_transit(C["t"], other, 1.0, s12.EPH)
    eps = s12.epochs_for(C, sig)
    tcs = np.array([x[1] for x in eps])
    ofs = np.array([tc - dumps[np.argmin(np.abs(dumps - tc))] for tc in tcs])
    Dr = diff_image_step(C, C["cube"], tcs, e["T14"], okc, dumps)
    allok = np.ones(len(C["t"]), bool)
    # Dump-matched fake epochs are impossible here: every S6 dump already has an A transit within 3.5 h, so any fake
    # placed at the same offset from a dump overlaps a real transit window. Instead: random fake epochs (as in s12),
    # processed with the SAME step-baseline code; the step only activates when a dump falls inside a window.
    cand_sets = s12.random_epoch_sets(C, sig, len(tcs), 300, np.random.default_rng(3))
    fakes = np.array([diff_image_step(C, C["cube"], s, e["T14"], allok, dumps).mean(0) for s in cand_sets])
    sig_map = np.std(fakes, axis=0)
    # also the real data WITHOUT the step (as in s12) for comparison
    Dr_nostep = s12.diff_image(C, C["cube"], tcs, e["T14"], okc)
    G = s12.GridFitter(prf, shape, tg["x"], tg["y"])
    fr = G.fit(Dr.mean(0), sig_map)
    fr_nostep = G.fit(Dr_nostep.mean(0), sig_map)
    depth = G.fit_fixed(Dr.mean(0), sig_map, tg["x"], tg["y"])["A"] / astro["target_flux"]
    offs = []
    for s in cand_sets:
        ci = s12.inject_on_target(C["cube"] + C["bkg"][:, None, None], C["bkg"], C["t"], s, e["T14"], depth, tg["x"], tg["y"], 30.0) \
            - C["bkg"][:, None, None]
        Di = diff_image_step(C, ci, s, e["T14"], allok, dumps)
        if len(Di) != len(ofs):
            continue
        f = G.fit(Di.mean(0), sig_map)
        offs.append((f["x"] - tg["x"], f["y"] - tg["y"]))
    o = np.array(offs)
    real = np.array([fr["x"] - tg["x"], fr["y"] - tg["y"]])
    mu, cov = o.mean(0), np.cov(o.T)
    res = dict(signal=sig, dumps=dumps.tolist(), transit_minus_dump_h=(ofs * 24).tolist(), n_epochs=len(tcs),
               n_dump_matched_sets=len(o), real_px=real.tolist(), inj_mean_px=mu.tolist(), inj_sd_px=np.sqrt(np.diag(cov)).tolist(),
               z=[float((real[i] - mu[i]) / np.sqrt(cov[i, i])) for i in (0, 1)],
               chi2_2dof=float((real - mu) @ np.linalg.inv(cov) @ (real - mu)), depth_ppm=float(depth * 1e6),
               real_px_without_step=[fr_nostep["x"] - tg["x"], fr_nostep["y"] - tg["y"]],
               note="injections at random epochs (dump-matched sets impossible: every dump has an A transit within 3.5 h)")
    return res


if __name__ == "__main__":
    out = {s: main(s) for s in ("A",)}
    save_json(out, "out/s12d_s6_dumps.json")
    print(json.dumps(out, indent=1)[:3000])
