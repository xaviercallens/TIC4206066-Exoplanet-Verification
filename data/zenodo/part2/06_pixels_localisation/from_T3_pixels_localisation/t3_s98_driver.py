"""T3 minimal rerun: why does the PUBLICATION_AUDIT S98 clean 3.18 d localisation (9 events) give (+1.2, -8.4)"
while U (opus55 03_pixels) gives (+5.6, -0.3)" for the SAME 9 events and identical pixel data?

Only the real-data free-fit offset is recomputed (no injections). Audit code is imported read-only from
tic4206066_PUBLICATION_AUDIT/scripts (PYTHONDONTWRITEBYTECODE=1); U code read-only from opus55.../03_pixels.
Nothing is written outside this work dir.
Variants (audit recipe = s12_pixels.load_cube/epochs_for/diff_image/GridFitter, s12c noise map seed 5):
  AUD_frozen     : audit clean mask (QLP FROZEN_QLP_MASK + nearest QLP cadence within 2 min)  -> must reproduce real_px (0.05,-0.45)
  AUD_clean29    : same recipe, U's clean29 mask (TESScut 2048 + QLP 2048 + QLP 2^29, 150 s match)
  AUD_frozen_jk  : jackknife over the 9 events
  AUD_frozen_unif: uniform pixel weights instead of the empirical noise map
Variants (U recipe = pixloc.Dataset + run_real.analyse):
  U_clean29      : U default -> must reproduce (+5.65,-0.32)
  U_frozenmask   : U recipe with the audit's FROZEN clean mask
"""
import os, sys, json
import numpy as np
W = "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation"
AS = "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_PUBLICATION_AUDIT/scripts"
US = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels"
from astropy.io import fits
res = {}

# ------------------------------------------------------------------ audit recipe
sys.path.insert(0, AS)
import lcio, prf_tools as pt, s12_pixels as s12
from s12c_clean_localisation import qlp_clean_mask

def u_clean29_mask(C):
    d = fits.getdata(os.path.join(lcio.MAST, lcio.PRODUCTS["QLP_S98_SYS"][0]), 1)
    tq, q = np.array(d["TIME"], float), np.array(d["QUALITY"], np.int64); ok = np.isfinite(tq); tq, q = tq[ok], q[ok]
    t = C["t"] - lcio.TESSCUT_SHIFT_S[98] / 86400.0   # U matches on raw TESScut time
    idx = np.clip(np.searchsorted(tq, t), 1, len(tq) - 1); j = np.where(np.abs(tq[idx] - t) < np.abs(tq[idx - 1] - t), idx, idx - 1)
    close = np.abs(tq[j] - t) < 150.0 / 86400
    bad = close & (((q[j] & 2048) > 0) | ((q[j] & (1 << 29)) > 0))
    return ~bad

def audit_like(mask_kind, drop=None, uniform=False, tag=""):
    C = s12.load_cube("S98_ffi")
    keep = qlp_clean_mask(C["t"]) if mask_kind == "frozen" else u_clean29_mask(C)
    for k in ("t", "cube", "bkg", "pc1", "pc2"):
        C[k] = C[k][keep]
    shape = C["cube"].shape[1:]
    sc0 = pt.scene(C["sources"], C["wcs"], C["tmid"], s12.TARGET_ID); tg0 = [s for s in sc0 if s["is_target"]][0]
    prf = pt.get_prf(C["cam"], C["ccd"], C["sector"], C["col0"] + tg0["x"], C["row0"] + tg0["y"])
    astro, srcs, _ = s12.astrometric_solution(C, prf); tg = [s for s in srcs if s["is_target"]][0]
    e = s12.EPH["A"]; okc = ~lcio.in_transit(C["t"], "B", 1.0, s12.EPH)
    eps = s12.epochs_for(C, "A"); tcs = [x[1] for x in eps]
    fsets = s12.random_epoch_sets(C, "A", len(tcs), 300, np.random.default_rng(5))
    fakes = np.array([s12.diff_image(C, C["cube"], s, e["T14"], np.ones(len(C["t"]), bool)).mean(0) for s in fsets])
    sig_map = np.ones(shape) if uniform else np.std(fakes, 0)
    G = s12.GridFitter(prf, shape, tg["x"], tg["y"])
    w = C["wcs"]; p0 = w.all_pix2world([[tg["x"], tg["y"]]], 0)[0]
    px = w.all_pix2world([[tg["x"] + 1, tg["y"]]], 0)[0]; py = w.all_pix2world([[tg["x"], tg["y"] + 1]], 0)[0]; cosd = np.cos(np.radians(p0[1]))
    J = np.array([[(px[0] - p0[0]) * cosd * 3600, (py[0] - p0[0]) * cosd * 3600], [(px[1] - p0[1]) * 3600, (py[1] - p0[1]) * 3600]])
    def fit(tlist):
        D = s12.diff_image(C, C["cube"], tlist, e["T14"], okc).mean(0); fr = G.fit(D, sig_map)
        real = np.array([fr["x"] - tg["x"], fr["y"] - tg["y"]]); return real, J @ real
    real, sky = fit(tcs)
    out = dict(n_cadences=int(len(C["t"])), epochs=[int(x[0]) for x in eps], tcs=tcs, real_px=real.tolist(), sky_raw_arcsec=sky.tolist(),
               audit_inj_mean_px=[0.030666666666655906, -0.03883333333334387],
               sky_minus_audit_injmean=(J @ (real - np.array([0.030666666666655906, -0.03883333333334387]))).tolist())
    if drop == "jk":
        jk = []
        for k in range(len(tcs)):
            r2, s2 = fit([t for i, t in enumerate(tcs) if i != k]); jk.append(dict(dropped_tc=tcs[k], sky_raw=s2.tolist()))
        out["jackknife"] = jk
    print(tag, mask_kind, "uniform" if uniform else "noisemap", "n_cad", len(C["t"]), "epochs", out["epochs"], "real_px", np.round(real, 3), "sky", np.round(sky, 2), flush=True)
    if drop == "jk":
        for j in out["jackknife"]: print("   drop", round(j["dropped_tc"], 3), np.round(j["sky_raw"], 2), flush=True)
    return out

res["AUD_frozen"] = audit_like("frozen", drop="jk", tag="AUD")
res["AUD_clean29"] = audit_like("clean29", drop="jk", tag="AUD")
res["AUD_frozen_unif"] = audit_like("frozen", uniform=True, tag="AUD")
json.dump(res, open(os.path.join(W, "t3_s98_driver.json"), "w"), indent=1, default=float)

# ------------------------------------------------------------------ U recipe
for m in ("lcio", "prf_tools", "s12_pixels", "s12c_clean_localisation", "audit_common"):
    sys.modules.pop(m, None)
sys.path.remove(AS); sys.path.insert(0, US)
import common as ucommon
ucommon.OUT = W    # never write into U
import pixloc, run_real
from pixloc import Dataset

def u_run(maskkind):
    ds = Dataset("S98_ffi", "clean29")
    if maskkind == "frozen":
        # rebuild the good mask as: U default bitmask + TESScut 2048 + audit FROZEN QLP clean mask (2-min match)
        c = ds.cube; good = c.good(ucommon.DEFAULT_BITMASK) & ((c.quality & 2048) == 0)
        d = fits.getdata(ucommon.FILES["QLP98"], 1); tq, q = np.array(d["TIME"], float), np.array(d["QUALITY"], np.int64)
        badq = ((q & 1610619071) != 0) | ~np.isfinite(tq); tqq = np.where(np.isfinite(tq), tq, -1e9)
        t = c.time + 0.49 / 86400
        idx = np.clip(np.searchsorted(tqq, t), 1, len(tqq) - 1); j = np.where(np.abs(tqq[idx] - t) < np.abs(tqq[idx - 1] - t), idx, idx - 1)
        keep = ~badq[j] & (np.abs(tqq[j] - t) < 2.0 / 1440)
        good &= keep & np.isfinite(c.time)
        ds.good = good; ds.t = c.time[good]; ds.F = c.flux[good]; ds.medimg = np.nanmedian(ds.F, axis=0)
    sf = ds.scene_fit()
    evs = [e for e in ds.event_windows("A", min_side_frac=run_real.MIN_SIDE) if e["ok"]]
    out, *_ = run_real.analyse(ds, "A", sf, events=evs, want_nb=False)
    jk = []
    for k in range(len(evs)):
        o2, *_ = run_real.analyse(ds, "A", sf, events=[e for i, e in enumerate(evs) if i != k], want_nb=False)
        jk.append(dict(dropped_tc=evs[k]["tc"], off=[o2["off_E"], o2["off_N"]]))
    r = dict(n_cadences=int(len(ds.t)), event_tc=[e["tc"] for e in evs], off=[out["off_E"], out["off_N"]], snr=out["snr_target"], jackknife=jk)
    print("U", maskkind, "n_cad", len(ds.t), "events", np.round(r["event_tc"], 3), "off", np.round(r["off"], 2), flush=True)
    for j in jk: print("   drop", round(j["dropped_tc"], 3), np.round(j["off"], 2), flush=True)
    return r

res["U_clean29"] = u_run("clean29")
res["U_frozenmask"] = u_run("frozen")
json.dump(res, open(os.path.join(W, "t3_s98_driver.json"), "w"), indent=1, default=float)
print("done")
