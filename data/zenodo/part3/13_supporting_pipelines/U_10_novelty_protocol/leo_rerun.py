"""Audit 11: independent re-derivation of the Sector 98 detection statistics and LEO-Vetter flux-level
test outcomes for TIC 4206066.  Own code (no project/addenda code imported).  Uses only:
  * the public QLP S98 light curve downloaded fresh from MAST (sha256 recorded),
  * wotan biweight (Hippke et al. 2019), astropy BoxLeastSquares (Kovacs et al. 2002),
  * the public leo-vetter 1.2.0 package (Kunimoto et al. 2025) installed into the audit env,
  * TIC 8.2 stellar parameters fetched fresh from MAST.
The quality-bit choice follows the manuscript's description (S98 scattered-light intervals removed) and the
documented QLP recommendations; it is a method choice, not copied code.
Output: leo_rerun.json
"""
import os, json, hashlib, datetime, warnings
import numpy as np, requests
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from wotan import flatten
warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw"); os.makedirs(RAW, exist_ok=True)
LOG = os.path.join(HERE, "queries.log")
URL = "https://archive.stsci.edu/hlsps/qlp/s0098/0000/0000/0420/6066/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"
now = lambda: datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
out = {}

fn = os.path.join(RAW, "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")
if not os.path.exists(fn):
    r = requests.get(URL, timeout=300, headers={"User-Agent": "novelty-audit/1.0"})
    open(fn, "wb").write(r.content)
    with open(LOG, "a") as f:
        f.write("\t".join([now(), "MAST QLP S98 light curve TIC 4206066", URL, str(r.status_code), str(len(r.content)),
                           hashlib.sha256(r.content).hexdigest(), "downloaded for independent MES/LEO re-run"]) + "\n")
out["lc_url"] = URL
out["lc_sha256"] = hashlib.sha256(open(fn, "rb").read()).hexdigest()
out["lc_sha256_matches_pipeline_record"] = out["lc_sha256"] == "a034136697b4f75765a91503107f55a455b365d19cd270ed633845efcdcb188e"

h = fits.open(fn)
d = h[1].data
t = np.array(d["TIME"], float); q = np.array(d["QUALITY"], int)
f = np.array(d["SYS_RM_FLUX"], float); fe = np.array(d["DET_FLUX_ERR"], float)
# QLP-update page recommends excluding bits 3, 6, 12, 31 (values 4, 32, 2048, 2^30); SPOC/TICA bits 1,2,8,128 (attitude tweak,
# safe mode, Earth point, manual exclude) and bit 30 (2^29, the scattered-light half-orbit flag masked per the manuscript text).
bad = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 2048 | 4096 | (1 << 29) | (1 << 30)
good = np.isfinite(t) & np.isfinite(f) & ((q & bad) == 0)
if np.isfinite(fe[good]).mean() < 0.9:
    fe = np.full_like(f, np.nan)
t, f, fe = t[good], f[good], fe[good]
f_n = f / np.nanmedian(f)
trend = flatten(t, f_n, method="biweight", window_length=1.0, break_tolerance=0.1, edge_cutoff=0.25, return_trend=True)[1]
fl = f_n / trend
ok = np.isfinite(fl)
t, fl, fe = t[ok], fl[ok], (fe[ok] / np.nanmedian(f)) / trend[ok] if np.isfinite(fe).any() else None
mad = 1.4826 * np.median(np.abs(fl - np.median(fl)))
keep = fl < 1 + 4 * mad            # upper-only clipping (flares/cosmics), never touches dips
t, fl = t[keep], fl[keep]
if fe is not None:
    fe = fe[keep]
    if not np.all(np.isfinite(fe)):
        fe = np.full_like(fl, 1.4826 * np.median(np.abs(np.diff(fl))) / np.sqrt(2))
else:
    fe = np.full_like(fl, 1.4826 * np.median(np.abs(np.diff(fl))) / np.sqrt(2))
out["n_cadences_used"] = int(len(t)); out["robust_std_ppm"] = float(mad * 1e6)

# ---- own BLS search (two iterations) ------------------------------------------------------------
bls_res = []
mask = np.ones(len(t), bool)
for it in range(2):
    bls = BoxLeastSquares(t[mask], fl[mask], fe[mask])
    periods = np.exp(np.linspace(np.log(0.5), np.log(25.0), 60000))
    durs = np.array([0.75, 1.0, 1.5, 1.75, 2.0, 2.5, 3.0]) / 24
    pg = bls.power(periods, durs, objective="snr")
    i = int(np.argmax(pg.power))
    P, T0, D = float(pg.period[i]), float(pg.transit_time[i]), float(pg.duration[i])
    # refine on a fine grid around the peak
    pf = np.linspace(P * 0.998, P * 1.002, 4001)
    dfine = np.arange(0.5, 4.01, 0.25) / 24
    pg2 = bls.power(pf, dfine, objective="snr")
    j = int(np.argmax(pg2.power))
    P, T0, D = float(pg2.period[j]), float(pg2.transit_time[j]), float(pg2.duration[j])
    bls_res.append(dict(iteration=it + 1, period=P, t0_btjd=T0, duration_hr=D * 24, depth_ppm=float(pg2.depth[j] * 1e6),
                        bls_snr=float(pg2.power[j])))
    ph = ((t - T0 + 0.5 * P) % P) - 0.5 * P
    mask &= np.abs(ph) > D           # mask +-1 duration around found transits
out["own_bls"] = bls_res

# ---- LEO-Vetter 1.2.0 flux-level metrics and tests -------------------------------------------
from astroquery.mast import Catalogs
from leo_vetter.main import TCELightCurve
from leo_vetter import thresholds as TH
from leo_vetter.stellar import quadratic_ldc
tic = Catalogs.query_criteria(catalog="Tic", ID=4206066)[0]
star = {k: float(tic[k]) if tic[k] is not None and np.isfinite(float(tic[k])) else np.nan for k in
        ("rad", "e_rad", "mass", "e_mass", "Teff", "e_Teff", "logg", "rho", "e_rho", "Tmag")}
if not np.isfinite(star["e_mass"]):
    star["e_mass"] = 0.1 * star["mass"]
    star["e_mass_note"] = "TIC e_mass missing -> 10% adopted"
u1, u2 = quadratic_ldc(star["Teff"], star["logg"])
star["u1"], star["u2"] = float(u1), float(u2)
out["star"] = star
FA = [TH.weak, TH.invalid_transits, TH.bad_shape, TH.non_unique, TH.chases, TH.dmm, TH.single_event, TH.bad_fit,
      TH.sinusoidal, TH.unphysical_duration, TH.asymmetric, TH.chi, TH.data_gapped]
FP = [TH.odd_even, TH.vshaped, TH.large, TH.secondary]

cases = {
    "A_pipeline_ephemeris": (3.182945191903069, 3998.7348383728113, 1.8 / 24),
    "B_pipeline_ephemeris": (11.132975671662091, 4000.8681717061445, 1.6 / 24),
}
for r in bls_res:
    tag = "A" if abs(r["period"] / 3.1829 - 1) < 0.01 else ("B" if abs(r["period"] / 11.133 - 1) < 0.01 else "other")
    cases[f"{tag}_own_bls_iter{r['iteration']}"] = (r["period"], r["t0_btjd"], r["duration_hr"] / 24)
# the pipeline searched each signal after masking the other signal's transits; do the same for the LEO run
res = {}
# pipeline order: B was signal 1 (searched on the unmasked light curve), A was signal 2 (searched after masking B)
for k in list(cases):
    cases[k + "_nomask"] = cases[k]
for name, (P, T0, D) in cases.items():
    tt, ff, ee = t, fl, fe
    other = None
    if name.endswith("_nomask"):
        pass
    elif name.startswith("A"):
        other = (11.132975671662091, 4000.8681717061445, 1.6 / 24)
    elif name.startswith("B"):
        other = (3.182945191903069, 3998.7348383728113, 1.8 / 24)
    if other is not None:
        pho = ((tt - other[1] + 0.5 * other[0]) % other[0]) - 0.5 * other[0]
        m = np.abs(pho) > other[2]
        tt, ff, ee = tt[m], ff[m], ee[m]
    try:
        tlc = TCELightCurve(4206066, tt, ff, ff, ee, P, T0, D, planetno=1)
        tlc.compute_flux_metrics(star, verbose=False)
        M = tlc.metrics
        fa = {fn.__name__: bool(fn(M, TH._default_thresholds)[0]) for fn in FA}
        fp = {fn.__name__: bool(fn(M, TH._default_thresholds)[0]) for fn in FP}
        res[name] = dict(P=P, T0=T0, dur_hr=D * 24, MES=float(M["MES"]), new_MES=float(M.get("new_MES", np.nan)),
                         N_transit=int(M["N_transit"]), new_N_transit=int(M.get("new_N_transit", -1)), depth_ppm=float(M["dep"] * 1e6),
                         FA_failed=[k for k, v in fa.items() if v], FP_failed=[k for k, v in fp.items() if v],
                         n_FA_tests=len(fa), n_FP_tests=len(fp),
                         key_metrics={k: float(M[k]) for k in ("SHP", "DMM", "CHI", "max_SES", "sig_dep", "sig_sec", "Fred", "FA1", "FA2",
                                                              "sig_pri", "sig_ter", "sig_pos", "mean_chases", "med_chases", "n_in", "sig_w", "sig_r", "Rp", "transit_b", "transit_RpRs")
                                      if k in M})
    except Exception as e:
        import traceback
        res[name] = dict(error=traceback.format_exc()[-1500:])
out["leo"] = res
json.dump(out, open(os.path.join(HERE, "leo_rerun.json"), "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str)[:6000])
