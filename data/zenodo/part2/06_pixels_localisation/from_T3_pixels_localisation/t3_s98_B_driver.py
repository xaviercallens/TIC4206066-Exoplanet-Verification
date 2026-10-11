"""T3 minimal rerun for signal B in S98: audit recipe (s12_pixels/s12c, imported read-only) vs mask, with a
jackknife over the 3 clean events; real free-fit offset only (no injections). Compare: audit stored
real_sky_arcsec (+8.33, +0.10) (s12c clean_S98.B) and U S98_clean29_B (+1.19, -1.02)."""
import os, sys, json
import numpy as np
W = "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation"
sys.path.insert(0, "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_PUBLICATION_AUDIT/scripts")
from astropy.io import fits
import lcio, prf_tools as pt, s12_pixels as s12
from s12c_clean_localisation import qlp_clean_mask
AUD_INJ_MEAN_B = np.array([0.0008333333333227344, -0.011333333333343945])
def u_clean29_mask(C):
    d = fits.getdata(os.path.join(lcio.MAST, lcio.PRODUCTS["QLP_S98_SYS"][0]), 1)
    tq, q = np.array(d["TIME"], float), np.array(d["QUALITY"], np.int64); ok = np.isfinite(tq); tq, q = tq[ok], q[ok]
    t = C["t"] - lcio.TESSCUT_SHIFT_S[98] / 86400.0
    idx = np.clip(np.searchsorted(tq, t), 1, len(tq) - 1); j = np.where(np.abs(tq[idx] - t) < np.abs(tq[idx - 1] - t), idx, idx - 1)
    return ~((np.abs(tq[j] - t) < 150.0 / 86400) & (((q[j] & 2048) > 0) | ((q[j] & (1 << 29)) > 0)))
res = {}
for mk in ("frozen", "clean29"):
    C = s12.load_cube("S98_ffi"); keep = qlp_clean_mask(C["t"]) if mk == "frozen" else u_clean29_mask(C)
    for k in ("t", "cube", "bkg", "pc1", "pc2"): C[k] = C[k][keep]
    shape = C["cube"].shape[1:]
    sc0 = pt.scene(C["sources"], C["wcs"], C["tmid"], s12.TARGET_ID); tg0 = [s for s in sc0 if s["is_target"]][0]
    prf = pt.get_prf(C["cam"], C["ccd"], C["sector"], C["col0"] + tg0["x"], C["row0"] + tg0["y"])
    astro, srcs, _ = s12.astrometric_solution(C, prf); tg = [s for s in srcs if s["is_target"]][0]
    e = s12.EPH["B"]; okc = ~lcio.in_transit(C["t"], "A", 1.0, s12.EPH)
    eps = s12.epochs_for(C, "B"); tcs = [x[1] for x in eps]
    fsets = s12.random_epoch_sets(C, "B", len(tcs), 300, np.random.default_rng(5))
    sig_map = np.std(np.array([s12.diff_image(C, C["cube"], s, e["T14"], np.ones(len(C["t"]), bool)).mean(0) for s in fsets]), 0)
    G = s12.GridFitter(prf, shape, tg["x"], tg["y"])
    w = C["wcs"]; p0 = w.all_pix2world([[tg["x"], tg["y"]]], 0)[0]
    px = w.all_pix2world([[tg["x"] + 1, tg["y"]]], 0)[0]; py = w.all_pix2world([[tg["x"], tg["y"] + 1]], 0)[0]; cosd = np.cos(np.radians(p0[1]))
    J = np.array([[(px[0] - p0[0]) * cosd * 3600, (py[0] - p0[0]) * cosd * 3600], [(px[1] - p0[1]) * 3600, (py[1] - p0[1]) * 3600]])
    def fit(tl):
        D = s12.diff_image(C, C["cube"], tl, e["T14"], okc).mean(0); fr = G.fit(D, sig_map)
        r = np.array([fr["x"] - tg["x"], fr["y"] - tg["y"]]); return r, J @ (r - AUD_INJ_MEAN_B)
    r, s = fit(tcs); jk = [dict(dropped_tc=tcs[k], sky=fit([t for i, t in enumerate(tcs) if i != k])[1].tolist()) for k in range(len(tcs))]
    res[mk] = dict(epochs=[int(x[0]) for x in eps], tcs=tcs, real_px=r.tolist(), sky_minus_audit_injmean=s.tolist(), jackknife=jk)
    print(mk, "epochs", res[mk]["epochs"], "real_px", np.round(r, 3), "sky", np.round(s, 2), "jk", [(round(j["dropped_tc"], 2), np.round(j["sky"], 1).tolist()) for j in jk], flush=True)
json.dump(res, open(os.path.join(W, "t3_s98_B_driver.json"), "w"), indent=1, default=float)
