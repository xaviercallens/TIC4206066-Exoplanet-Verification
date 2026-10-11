"""Step 15: stellar variability / rotation / activity.
TESS: Lomb-Scargle of the NON-cotrended fluxes (SPOC/TESS-SPOC SAP_FLUX, QLP SAP_FLUX, own TESScut r=2.5 px) per sector,
0.3-14 d, with the spectral window and a bootstrap false-alarm level (1000 shuffles of the flux values).
ASAS-SN Sky Patrol v2 (g and V; file retrieved 2026-09-25T03:20:56Z by the deep follow-up with pyasassn, re-used here
with sha256 recorded; quality 'G' only; per-camera zero points removed): LS 1-100 d, window, bootstrap FAP.
Harmonic bookkeeping vs P_A, P_B and the 1-d alias of each peak. Catalogue rotation periods (Gao+2025, Schmitt+2026,
Wang+2025, TSS_ARC, TARS) are listed from the priority re-check and tested here against the data.
X-ray: eRASS1 1eRASS J052326.6-011938 (DET_LIKE_0 = 8.04, 4.5 counts, 4.8'' from the target).
"""
import csv
import json
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
from audit_common import AUD, OLD_DEEP, save_json, sha256_file  # noqa

rng = np.random.default_rng(15)
EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))


def ls_analysis(t, f, pmin, pmax, nboot=1000, label=""):
    m = np.isfinite(t) & np.isfinite(f)
    t, f = t[m], f[m]
    f = f / np.nanmedian(f) - 1
    freq = np.linspace(1 / pmax, 1 / pmin, 20000)
    ls = LombScargle(t, f)
    pw = ls.power(freq)
    win = LombScargle(t, np.ones_like(t), fit_mean=False, center_data=False).power(freq)
    order = np.argsort(pw)[::-1]
    peaks = []
    for j in order:
        P = 1 / freq[j]
        if all(abs(P / p["P"] - 1) > 0.05 for p in peaks):
            peaks.append(dict(P=float(P), power=float(pw[j]), window_power_at_P=float(win[j])))
        if len(peaks) >= 5:
            break
    boot = []
    for _ in range(nboot):
        boot.append(LombScargle(t, rng.permutation(f)).power(freq).max())
    fap_lvl = dict(p01=float(np.percentile(boot, 99)), p001=float(np.percentile(boot, 99.9)))
    for p in peaks:
        p["bootstrap_fap"] = float(np.mean(np.array(boot) >= p["power"]))
        p["rel"] = [f"{s}:{lab}" for s, P0 in (("A", EPH["A"]["P"]), ("B", EPH["B"]["P"])) for lab, fac in (("P", 1), ("2P", 2), ("P/2", .5), ("3P", 3), ("4P", 4), ("5P", 5))
                    if abs(p["P"] / (P0 * fac) - 1) < 0.03]
    amp = float(np.sqrt(peaks[0]["power"] * 2) * np.std(f)) if peaks else None
    return dict(label=label, n=int(len(t)), baseline_d=float(t.max() - t.min()), peaks=peaks, bootstrap_levels=fap_lvl,
                approx_semi_amplitude_ppm=amp * 1e6 if amp else None, rms_ppm=float(np.std(f) * 1e6))


out = {"tess": [], "asassn": {}}
for key in ("TESSSPOC_S6_SAP", "SPOC2min_S32_SAP", "QLP_S6_SAP", "QLP_S32_SAP", "QLP_S98_SAP"):
    lc = lcio.load_product(key, "frozen" if "S98" in key else "default")
    t, f = lc["t"], lc["f"]
    if lc["cad"] < 0.004:        # bin to 30 min for speed
        e = np.arange(t.min(), t.max() + 1 / 48, 1 / 48)
        i = np.digitize(t, e)
        t = np.array([t[i == k].mean() for k in np.unique(i)])
        f = np.array([f[i == k].mean() for k in np.unique(i)])
    r = ls_analysis(t, f, 0.3, 14.0, label=key)
    out["tess"].append(r)
    print(key, [(round(p["P"], 2), round(p["power"], 3), p["bootstrap_fap"], p["rel"]) for p in r["peaks"][:3]], "rms", round(r["rms_ppm"]))
for nm, r_ in (("tesscut_s6", 2.5), ("tesscut_s32", 2.5), ("tesscut_s98", 2.5)):
    lc = lcio.load_own(nm, r_, drop_stray=True)
    t, f = lc["t"], lc["f"]
    e = np.arange(t.min(), t.max() + 1 / 48, 1 / 48)
    i = np.digitize(t, e)
    t = np.array([t[i == k].mean() for k in np.unique(i)])
    f = np.array([f[i == k].mean() for k in np.unique(i)])
    r = ls_analysis(t, f, 0.3, 14.0, label=nm + "_own_raw")
    out["tess"].append(r)
    print(nm, [(round(p["P"], 2), round(p["power"], 3), p["bootstrap_fap"], p["rel"]) for p in r["peaks"][:3]], "rms", round(r["rms_ppm"]))

# ASAS-SN
fn = os.path.join(OLD_DEEP, "queries", "asassn_skypatrol_v2_lc.csv")
out["asassn"]["file"] = os.path.relpath(fn, AUD)
out["asassn"]["sha256"] = sha256_file(fn)
rows = list(csv.DictReader(open(fn)))
for band in ("g", "V"):
    R = [r for r in rows if r["phot_filter"] == band and r["quality"] == "G" and float(r["mag_err"]) < 0.05]
    if len(R) < 30:
        continue
    t = np.array([float(r["jd"]) for r in R]) - 2457000
    mag = np.array([float(r["mag"]) for r in R])
    cam = np.array([r["camera"] for r in R])
    for c in np.unique(cam):                      # remove per-camera zero-point offsets
        mag[cam == c] -= np.median(mag[cam == c])
    clip = np.abs(mag - np.median(mag)) < 5 * 1.4826 * np.median(np.abs(mag - np.median(mag)))
    f = 10 ** (-0.4 * mag[clip])
    res = ls_analysis(t[clip], f, 1.0, 100.0, nboot=300, label=f"ASAS-SN {band}")
    res["n_cameras"] = int(len(np.unique(cam)))
    res["jd_range"] = [float(t.min()), float(t.max())]
    res["mad_mmag"] = float(1.4826 * np.median(np.abs(mag[clip] - np.median(mag[clip]))) * 1e3)
    out["asassn"][band] = res
    print("ASAS-SN", band, res["n"], [(round(p["P"], 2), round(p["power"], 4), p["bootstrap_fap"]) for p in res["peaks"][:4]], "MAD mmag", round(res["mad_mmag"], 1))

# catalogue periods to test (from the priority re-check, priority_matches.json)
out["catalogue_periods_reported"] = {"Gao+2025 ROT": 8.1838, "Schmitt+2026": 8.1106, "Wang+2025": 7.9764, "TSS_ARC mirror (TESS-SVC)": 7.9834,
                                     "TARS (S32, doubled 15.85)": 7.93}
# X-ray
out["xray"] = dict(source="1eRASS J052326.6-011938", sep_arcsec=4.797, det_like_0=8.0396881, net_counts=4.49, flux_cgs=3.36859955e-14,
                   flux_lower_err_cgs=1.43962075e-14, distance_pc=35.58)
Lx = 4 * np.pi * (35.58 * 3.0857e18) ** 2 * 3.36859955e-14
for L in (0.090, 0.107):
    out["xray"][f"log_Lx_Lbol_Lbol{L}"] = float(np.log10(Lx / (L * 3.828e33)))
out["xray"]["Lx_erg_s"] = float(Lx)
save_json(out, "out/s15_activity.json")
print(json.dumps(out["xray"], indent=1))
