"""Step 12b: independent routes for the pixel localisation (the PRF fit of s12 shows a +x / east offset of 0.3-0.5 px for
signal A in S6 and S98 at ~2.5-3 sigma against on-target injections). Three further tests, all injection-calibrated with
the SAME fake epoch sets (model-free on-target injections into real pixels):
  C1  flux-weighted centroid of a fixed r = 2 px aperture: in-transit minus local out-of-transit, averaged over epochs
  C2  empirical-PRF fit: difference image = A * (local out-of-transit image shifted by dx, dy) + c (cubic spline shift)
  C3  per-epoch PRF fits (free position) of the real transits vs per-epoch fits of injected transits (consistency of
      individual events; shows whether one or two events drive the offset)
Also reports the scatter-light proxy (background level) during the real in-transit windows vs random windows.
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import shift as ndshift
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import prf_tools as pt  # noqa
import s12_pixels as s12  # noqa
from audit_common import save_json  # noqa


def centroid_shift(C, cube, tcs, T14, okmask, x0, y0, r=2.0):
    t = C["t"]
    ny, nx = cube.shape[1:]
    yy, xx = np.mgrid[:ny, :nx]
    ap = np.hypot(xx - x0, yy - y0) <= r
    F = cube[:, ap].sum(axis=1)
    cx = (cube[:, ap] * xx[ap]).sum(axis=1) / F
    cy = (cube[:, ap] * yy[ap]).sum(axis=1) / F
    out = []
    for tc in tcs:
        i_in = okmask & (np.abs(t - tc) < 0.4 * T14)
        i_o = okmask & (np.abs(t - tc) > 0.75 * T14) & (np.abs(t - tc) < 2.5 * T14)
        if i_in.sum() < 2 or i_o.sum() < 4:
            continue
        X = np.vstack([np.ones(i_o.sum()), t[i_o] - tc]).T
        bx = np.linalg.lstsq(X, cx[i_o], rcond=None)[0][0]
        by = np.linalg.lstsq(X, cy[i_o], rcond=None)[0][0]
        bf = np.linalg.lstsq(X, F[i_o], rcond=None)[0][0]
        out.append((cx[i_in].mean() - bx, cy[i_in].mean() - by, 1 - F[i_in].mean() / bf))
    return np.array(out)


def emp_prf_fit(D, O, box):
    def res(p):
        m = p[0] * ndshift(O, (p[2], p[1]), order=3, mode="nearest") + p[3]
        return ((D - m)[box]).ravel()
    s = least_squares(res, [D[box].sum() / O[box].sum(), 0.0, 0.0, 0.0], bounds=([-1, -2, -2, -np.inf], [1, 2, 2, np.inf]))
    return s.x


def run(name, sig, nsets=200, seed=11):
    C = s12.load_cube(name)
    shape = C["cube"].shape[1:]
    sc0 = pt.scene(C["sources"], C["wcs"], C["tmid"], s12.TARGET_ID)
    tg0 = [s for s in sc0 if s["is_target"]][0]
    prf = pt.get_prf(C["cam"], C["ccd"], C["sector"], C["col0"] + tg0["x"], C["row0"] + tg0["y"])
    astro, srcs, _ = s12.astrometric_solution(C, prf)
    tg = [s for s in srcs if s["is_target"]][0]
    e = s12.EPH[sig]
    other = "B" if sig == "A" else "A"
    okc = ~lcio.in_transit(C["t"], other, 1.0, s12.EPH)
    allok = np.ones(len(C["t"]), bool)
    eps = s12.epochs_for(C, sig)
    tcs = [x[1] for x in eps]
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    box = np.hypot(xx - tg["x"], yy - tg["y"]) <= 3.5
    # real
    cs_real = centroid_shift(C, C["cube"], tcs, e["T14"], okc, tg["x"], tg["y"])
    Dreal = s12.diff_image(C, C["cube"], tcs, e["T14"], okc)
    O = np.nanmedian(C["cube"][~(lcio.in_transit(C["t"], "A", 1.5, s12.EPH) | lcio.in_transit(C["t"], "B", 1.5, s12.EPH))], axis=0)
    emp_real = emp_prf_fit(Dreal.mean(0), O, box)
    depth = cs_real[:, 2].mean() if len(cs_real) else 5e-4
    # per-epoch PRF fits (real)
    fsets = s12.random_epoch_sets(C, sig, len(tcs), 300, np.random.default_rng(seed))
    fakes = np.array([s12.diff_image(C, C["cube"], s, e["T14"], allok).mean(0) for s in fsets])
    sig_map = np.std(fakes, axis=0)
    fakes1 = np.array([s12.diff_image(C, C["cube"], [s[0]], e["T14"], allok)[0] for s in fsets])
    sig1 = np.std(fakes1, axis=0)
    G = s12.GridFitter(prf, shape, tg["x"], tg["y"])
    per_ep = []
    for (n, tc), De in zip(eps, Dreal):
        f = G.fit(De, sig1)
        per_ep.append(dict(epoch=int(n), dx=f["x"] - tg["x"], dy=f["y"] - tg["y"], A=f["A"]))
    # injections: same fake sets for all three statistics
    isets = s12.random_epoch_sets(C, sig, len(tcs), nsets, np.random.default_rng(seed + 1))
    inj_c, inj_e, inj_1 = [], [], []
    for s in isets:
        ci = s12.inject_on_target(C["cube"] + C["bkg"][:, None, None], C["bkg"], C["t"], s, e["T14"], depth, tg["x"], tg["y"], 30.0) \
            - C["bkg"][:, None, None]
        c_ = centroid_shift(C, ci, s, e["T14"], allok, tg["x"], tg["y"])
        inj_c.append(c_[:, :2].mean(0))
        Di = s12.diff_image(C, ci, s, e["T14"], allok)
        inj_e.append(emp_prf_fit(Di.mean(0), O, box)[1:3])
        f1 = G.fit(Di[0], sig1)
        inj_1.append((f1["x"] - tg["x"], f1["y"] - tg["y"]))
    inj_c, inj_e, inj_1 = np.array(inj_c), np.array(inj_e), np.array(inj_1)

    def z(real, arr):
        return [float((real[i] - arr[:, i].mean()) / arr[:, i].std()) for i in (0, 1)]

    real_c = cs_real[:, :2].mean(0)
    res = dict(name=name, signal=sig, n_epochs=len(tcs), depth_aperture_r2=float(depth),
               C1_centroid=dict(real_mpx=[float(real_c[0] * 1e3), float(real_c[1] * 1e3)],
                                inj_mean_mpx=[float(inj_c[:, 0].mean() * 1e3), float(inj_c[:, 1].mean() * 1e3)],
                                inj_sd_mpx=[float(inj_c[:, 0].std() * 1e3), float(inj_c[:, 1].std() * 1e3)], z=z(real_c, inj_c)),
               C2_empirical_prf=dict(real_px=[float(emp_real[1]), float(emp_real[2])],
                                     inj_mean_px=[float(inj_e[:, 0].mean()), float(inj_e[:, 1].mean())],
                                     inj_sd_px=[float(inj_e[:, 0].std()), float(inj_e[:, 1].std())],
                                     z=z(emp_real[1:3], inj_e)),
               C3_per_epoch=dict(real=per_ep, inj_single_sd_px=[float(inj_1[:, 0].std()), float(inj_1[:, 1].std())],
                                 inj_single_mean_px=[float(inj_1[:, 0].mean()), float(inj_1[:, 1].mean())],
                                 real_mean_px=[float(np.mean([p["dx"] for p in per_ep])), float(np.mean([p["dy"] for p in per_ep]))],
                                 real_median_px=[float(np.median([p["dx"] for p in per_ep])), float(np.median([p["dy"] for p in per_ep]))],
                                 chi2_x_vs_inj=float(np.sum(((np.array([p["dx"] for p in per_ep]) - inj_1[:, 0].mean()) / inj_1[:, 0].std()) ** 2)),
                                 chi2_y_vs_inj=float(np.sum(((np.array([p["dy"] for p in per_ep]) - inj_1[:, 1].mean()) / inj_1[:, 1].std()) ** 2))))
    return res


if __name__ == "__main__":
    out = {}
    for name in sys.argv[1:] or ["S6_ffi", "S32_ffi", "S32_tpf2min", "S98_ffi"]:
        for sig in ("A", "B"):
            r = run(name, sig)
            out[f"{name}_{sig}"] = r
            c1, c2, c3 = r["C1_centroid"], r["C2_empirical_prf"], r["C3_per_epoch"]
            print(f"{name:12s} {sig} n={r['n_epochs']}: C1 centroid z=({c1['z'][0]:+.2f},{c1['z'][1]:+.2f}) real {c1['real_mpx']} inj {c1['inj_mean_mpx']}±{c1['inj_sd_mpx']}; "
                  f"C2 empPRF z=({c2['z'][0]:+.2f},{c2['z'][1]:+.2f}) real {c2['real_px']}; "
                  f"C3 per-epoch chi2 x {c3['chi2_x_vs_inj']:.1f} y {c3['chi2_y_vs_inj']:.1f} (n={r['n_epochs']}) median {c3['real_median_px']}", flush=True)
            save_json(out, f"out/s12b_centroid_{'_'.join(sys.argv[1:]) or 'all'}.json")
