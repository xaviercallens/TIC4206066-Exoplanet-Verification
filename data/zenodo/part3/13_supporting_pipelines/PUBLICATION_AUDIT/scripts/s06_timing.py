"""Step 6: time-system audit.
(a) Which time system every product declares (TIMESYS/TIMEREF/TASSIGN/BJDREFI/BJDREFF/TIMEUNIT/TIMEPIXR).
(b) TESScut TIME carries the mission barycentric correction for the CCD centre (astrocut docs). The target-specific
    correction differs by (r_obs . (n_target - n_ccdcentre))/c. We compute the barycentric correction for the target
    and compare with the TIMECORR column (applied by SPOC for the CCD centre) -> expected TESScut offset.
(c) Cadence-by-cadence comparison of TESScut TIME with the product TIME (TESS-SPOC S6/S32 FFI LCs, QLP S6/S32/S98).
    The product TIME should equal TESScut TIME + (corr_target - corr_centre) if both use target-specific corrections.
QLP had documented barycentric errors in S74-79 and S99-104 (tess.mit.edu/qlp); S98 must be checked explicitly.
"""
import glob
import os
import sys

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris
from astropy.io import fits
from astropy.time import Time

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, save_json  # noqa

solar_system_ephemeris.set("builtin")
MAST = os.path.join(AUD, "data", "mast")
target = SkyCoord(80.85982665567 * u.deg, -1.32696401034 * u.deg)
geo = EarthLocation.from_geocentric(0 * u.m, 0 * u.m, 0 * u.m)
out = {"declared": {}, "compare": {}}

prods = {
    "SPOC_2min_S32": "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits",
    "TESS-SPOC_S6": "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits",
    "TESS-SPOC_S32": "hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits",
    "QLP_S6": "hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits",
    "QLP_S32": "hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits",
    "QLP_S98": "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits",
    "eleanor_S6": "hlsp_gsfc-eleanor-lite_tess_ffi_s0006-0000000004206066_tess_v1.0_lc.fits",
    "TGLC_S6": "hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits",
    "T16_S6": "hlsp_t16_tess_ffi_s0006-cam1-ccd1-03220388198192519424_tess_v01_llc.fits",
    "TASOC_S6": "hlsp_tasoc_tess_ffi_tic00004206066-s0006-cam1-ccd1-c1800_tess_v05_cbv-lc.fits",
    "TARS_S6": "hlsp_tars_tess_ffi_s0006-0000000004206066_tess_v01_lc.fits",
    "TARS_S32": "hlsp_tars_tess_ffi_s0032-0000000004206066_tess_v01_lc.fits",
}
for k, fn in prods.items():
    h = fits.open(os.path.join(MAST, fn))
    hd = {**h[0].header, **h[1].header}
    out["declared"][k] = {kk: (str(hd.get(kk)) if hd.get(kk) is not None else None) for kk in
                          ["TIMESYS", "TIMEREF", "TASSIGN", "BJDREFI", "BJDREFF", "TIMEUNIT", "TIMEPIXR", "TELAPSE",
                           "TSTART", "TSTOP", "TIMEDEL", "CALIB", "PROCVER", "DATA_REL", "TICVER", "CROWDSAP",
                           "FLFRCSAP", "BESTAP", "DETMET", "CAMERA", "CCD"]}
    out["declared"][k]["columns"] = h[1].data.columns.names
    tu = h[1].columns["TIME"].unit if "TIME" in h[1].data.columns.names else None
    out["declared"][k]["TIME_unit"] = str(tu)

# TESScut vs products
for s in (6, 32, 98):
    p = glob.glob(os.path.join(AUD, f"data/tesscut/tess-s{s:04d}-*_astrocut.fits"))[0]
    h = fits.open(p)
    t = np.asarray(h[1].data["TIME"], float)
    tc = np.asarray(h[1].data["TIMECORR"], float)
    ok = np.isfinite(t)
    t, tc = t[ok], tc[ok]
    # barycentric correction for the target with the geocentre as observer (TESS-Earth offset < ~2 light-seconds)
    tt = Time(t + 2457000.0 - tc, format="jd", scale="tdb")  # approx. spacecraft time (TDB) before correction
    ltt = tt.light_travel_time(target, kind="barycentric", location=geo).to_value(u.day)
    dcorr = (ltt - tc) * 86400.0  # target correction minus applied (CCD-centre) correction [s]
    info = dict(n=int(len(t)), timecorr_applied_s=[float(tc.min() * 86400), float(tc.max() * 86400)],
                target_minus_applied_corr_s=dict(min=float(dcorr.min()), max=float(dcorr.max()), median=float(np.median(dcorr))))
    comps = {}
    for k, fn in prods.items():
        if not k.endswith(f"_S{s}"):
            continue
        d = fits.open(os.path.join(MAST, fn))[1].data
        if "TIME" not in d.columns.names:
            continue
        tp = np.asarray(d["TIME"], float)
        tp = tp[np.isfinite(tp)]
        if k == "SPOC_2min_S32":
            continue
        idx = np.clip(np.searchsorted(t, tp), 1, len(t) - 1)
        near = np.where(np.abs(t[idx] - tp) < np.abs(t[idx - 1] - tp), idx, idx - 1)
        dt_s = (tp - t[near]) * 86400.0
        good = np.abs(dt_s) < 120
        dexp = dcorr[near][good]
        comps[k] = dict(n_matched=int(good.sum()), product_minus_tesscut_s=dict(median=float(np.median(dt_s[good])), std=float(np.std(dt_s[good])),
                                                                               min=float(dt_s[good].min()), max=float(dt_s[good].max())),
                        expected_if_target_specific_s=float(np.median(dexp)),
                        residual_s=dict(median=float(np.median(dt_s[good] - dexp)), p01=float(np.percentile(dt_s[good] - dexp, 1)),
                                        p99=float(np.percentile(dt_s[good] - dexp, 99))))
        # time dependence of the residual (split sector into 4 quarters) - catches drifting barycentric errors
        q = np.array_split(np.arange(good.sum()), 4)
        comps[k]["residual_by_quarter_s"] = [float(np.median((dt_s[good] - dexp)[i])) for i in q]
    out["compare"][f"S{s}"] = dict(tesscut=info, products=comps)

# SPOC 2-min vs TESScut S32: 2-min timestamps are target-specific; compare interpolated
d = fits.open(os.path.join(MAST, prods["SPOC_2min_S32"]))[1].data
t2 = np.asarray(d["TIME"], float)
t2 = t2[np.isfinite(t2)]
h = fits.open(glob.glob(os.path.join(AUD, "data/tesscut/tess-s0032-*_astrocut.fits"))[0])
t10 = np.asarray(h[1].data["TIME"], float)
# FFI mid-times should coincide with the mid-time of the 5 two-min cadences they span (SPOC: FFIs are co-adds of 2-s
# frames; 2-min and FFI exposures are synchronous). Nearest-neighbour offsets modulo 2 min tell the phase.
idx = np.clip(np.searchsorted(t10, t2), 1, len(t10) - 1)
dd = np.minimum(np.abs(t10[idx] - t2), np.abs(t10[idx - 1] - t2)) * 86400
out["compare"]["S32_2min_vs_tesscut"] = dict(min_abs_offset_s=float(dd.min()), note="2-min grid vs 10-min FFI grid; offsets are k*120 s + const")
save_json(out, "out/s06_timing.json")
for s in ("S6", "S32", "S98"):
    c = out["compare"][s]
    print(s, "TESScut target-minus-applied corr (s):", {k: round(v, 2) for k, v in c["tesscut"]["target_minus_applied_corr_s"].items()})
    for k, v in c["products"].items():
        print("   ", k, "n", v["n_matched"], "prod-tesscut median", round(v["product_minus_tesscut_s"]["median"], 2), "s; expected",
              round(v["expected_if_target_specific_s"], 2), "; residual median", round(v["residual_s"]["median"], 2),
              "[", round(v["residual_s"]["p01"], 2), ",", round(v["residual_s"]["p99"], 2), "] by quarter", [round(x, 2) for x in v["residual_by_quarter_s"]])
for k, v in out["declared"].items():
    print(k, {kk: v[kk] for kk in ["TIMESYS", "TIMEREF", "BJDREFI", "BJDREFF", "TIMEUNIT", "CALIB", "CROWDSAP", "FLFRCSAP"]})
