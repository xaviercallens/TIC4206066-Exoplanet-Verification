"""Audit 8: run TRICERATOPS 1.1.0 on my own folds, recording per-draw STP/PTP weights.
Run with trice_env/bin/python, cwd = this directory (TRICERATOPS writes 4206066_gaia_background.csv to cwd).
usage: f02_trice_run.py VARIANT SIGNAL FOLD NRUNS [--N 1000000] [--M 0.606 --R 0.623 --Teff 4030 --plx 28.075]
        [--ap spoc|3x3|5x5|1x1] [--exptime_min 6] [--sigma emp|formal] [--molusc path] [--seed0 100] [--sectors 32,98]
"""
import argparse, json, os, pickle, sys, time, warnings, hashlib
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__))
SCR = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache"
os.makedirs(SCR, exist_ok=True)
import numpy as np
import lightkurve
lightkurve.conf.cache_dir = SCR
import triceratops.triceratops as tt
import triceratops.marginal_likelihoods as ml
from importlib.metadata import version

ap = argparse.ArgumentParser()
ap.add_argument("variant"); ap.add_argument("signal"); ap.add_argument("fold"); ap.add_argument("nruns", type=int)
ap.add_argument("--N", type=int, default=1_000_000)
ap.add_argument("--M", type=float, default=0.606); ap.add_argument("--R", type=float, default=0.623)
ap.add_argument("--Teff", type=float, default=4030.0); ap.add_argument("--plx", type=float, default=28.075)
ap.add_argument("--ap", default="spoc"); ap.add_argument("--exptime_min", type=float, default=6.0)
ap.add_argument("--sigma", default="emp"); ap.add_argument("--molusc", default=None)
ap.add_argument("--seed0", type=int, default=100); ap.add_argument("--sectors", default="32,98")
ap.add_argument("--tdepth_ppm", type=float, default=None)
ap.add_argument("--no_record", action="store_true")
a = ap.parse_args()
sectors = [int(s) for s in a.sectors.split(",")]
OUT = os.path.join(H, "runs", a.variant); os.makedirs(OUT, exist_ok=True)

# ---------- target (cached per sector set) ----------
pk = os.path.join(H, "runs", f"target_{'_'.join(map(str, sectors))}.pkl")
if os.path.exists(pk):
    target = pickle.load(open(pk, "rb"))
else:
    t0 = time.time()
    target = tt.target(ID=4206066, sectors=np.array(sectors), search_radius=10, lightkurve_cache_dir=SCR)
    pickle.dump(target, open(pk, "wb"))
    target.stars.to_csv(os.path.join(H, "runs", f"stars_{'_'.join(map(str, sectors))}_TIC.csv"), index=False)
    print("target init %.0fs, %d stars, bkg file %s" % (time.time() - t0, len(target.stars), target.trilegal_fname), flush=True)
for k, v in dict(mass=a.M, rad=a.R, Teff=a.Teff, plx=a.plx).items():
    target.update_star(4206066, k, v)

# ---------- apertures ----------
def aperture(k, sector):
    tc = np.round(target.pix_coords[k][0])     # target pixel (col,row) in this sector
    if a.ap == "spoc" and sector == 32:
        # SPOC S32 optimal aperture from the 2-min TPF (APERTURE bit 2); CRVAL1P=1368 (col), CRVAL2P=939 (row)
        from astropy.io import fits
        h = fits.open(os.path.join(H, "..", "00_data", "mast", "tess2020324010417-s0032-0000000004206066-0200-s_tp.fits"))
        img = h[2].data; c0, r0 = h[2].header["CRVAL1P"], h[2].header["CRVAL2P"]
        rr, cc = np.nonzero((img & 2) > 0)
        return np.c_[c0 + cc, r0 + rr].astype(float)
    n = {"spoc": 3, "3x3": 3, "5x5": 5, "1x1": 1}[a.ap]   # S98: QLP BESTAP header = '3x3 square'
    h_ = n // 2
    return np.array([[tc[0] + i, tc[1] + j] for i in range(-h_, h_ + 1) for j in range(-h_, h_ + 1)], float)
aps = [aperture(k, s) for k, s in enumerate(sectors)]

lc = np.loadtxt(os.path.join(H, "folds", a.fold + ".csv"), delimiter=",")
fm = json.load(open(os.path.join(H, "folds", "folds_meta.json")))["folds"][a.fold]
sigma = fm["sigma_empirical_oot"] if a.sigma == "emp" else fm["sigma_formal_mean"]
P = {"A": 3.1827856, "B": 11.13274}[a.signal]
tdepth = (a.tdepth_ppm or fm["depth_core_ppm"]) * 1e-6
target.calc_depths(tdepth=tdepth, all_ap_pixels=aps)
stars = target.stars.copy()

# ---------- recording hooks ----------
REC = {"scen": None}
orig = {n: getattr(ml, n) for n in ("sample_q_companion", "sample_ecc", "_log_mean_exp")}
def q_hook(x, M_s):
    q = orig["sample_q_companion"](x, M_s)
    if REC["scen"] in ("STP", "PTP") and "q" not in REC[REC["scen"]]: REC[REC["scen"]]["q"] = np.array(q, copy=True)
    return q
def e_hook(x, planet, P_orb):
    e = orig["sample_ecc"](x, planet, P_orb)
    if REC["scen"] in ("STP", "PTP") and "e" not in REC[REC["scen"]]: REC[REC["scen"]]["e"] = np.array(e, copy=True)
    return e
def lme_hook(logw, *, N_total):
    if REC["scen"] in ("STP", "PTP") and "logw" not in REC[REC["scen"]]: REC[REC["scen"]]["logw"] = np.array(logw, copy=True)
    return orig["_log_mean_exp"](logw, N_total=N_total)
ml.sample_q_companion = q_hook; ml.sample_ecc = e_hook; ml._log_mean_exp = lme_hook
def wrap(name, fn):
    def w(*args, **kw):
        REC["scen"] = name; REC[name] = {}
        try: return fn(*args, **kw)
        finally: REC["scen"] = None
    return w
# pandas-3 copy-on-write shim: TRICERATOPS 1.1.0 does `qs_comp = df[...]["mass ratio"].values; qs_comp[qs_comp < 0.1/M_s] = ...`,
# which raises "assignment destination is read-only" under pandas >= 3. Serve writable numpy copies instead (same numbers).
if a.molusc is not None:
    import types, pandas as _pd
    class _Frame:
        def __init__(self, d): self.d = d
        def __getitem__(self, k):
            if isinstance(k, str): return types.SimpleNamespace(values=np.array(self.d[k], copy=True))
            return _Frame({c: v[k] for c, v in self.d.items()})
    def _read_csv(f):
        df = _pd.read_csv(f); return _Frame({c: df[c].to_numpy(copy=True) for c in df.columns})
    ml.read_csv = _read_csv
if not a.no_record:
    tt.lnZ_STP = wrap("STP", ml.lnZ_STP); tt.lnZ_PTP = wrap("PTP", ml.lnZ_PTP)

def molusc_q(M_s):
    import pandas as pd
    d = pd.read_csv(a.molusc); A_ = d["semi-major axis(AU)"].values; E_ = d["eccentricity"].values
    q = d[A_ * (1 - E_) > 10]["mass ratio"].values.copy(); q[q < 0.1 / M_s] = 0.1 / M_s
    return np.pad(q, (0, a.N - len(q)))

def summarize(scen, M_s):
    r = REC.get(scen, {})
    if "logw" not in r: return None
    q = r.get("q") if a.molusc is None else molusc_q(M_s)
    lw = r["logw"]; fin = np.isfinite(lw)
    if fin.sum() == 0: return dict(n_finite=0)
    mx = lw[fin].max(); w = np.zeros_like(lw); w[fin] = np.exp(lw[fin] - mx); W = w.sum()
    m = q * M_s; e = r["e"]
    out = dict(n_finite=int(fin.sum()), lnZ_check=float(mx + np.log(W / a.N)))
    edges_m = [0.0, 0.1, 0.2, 0.3, 0.38, 0.45, 0.54, 0.607]
    out["wfrac_mass_bins"] = {f"{edges_m[i]:.2f}-{edges_m[i+1]:.2f}": float(w[(m >= edges_m[i]) & (m < edges_m[i+1])].sum() / W) for i in range(len(edges_m) - 1)}
    edges_e = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    out["wfrac_ecc_bins"] = {f"{edges_e[i]:.2f}-{edges_e[i+1]:.2f}": float(w[(e >= edges_e[i]) & (e < edges_e[i+1])].sum() / W) for i in range(len(edges_e) - 1)}
    out["wmean_mass"] = float((w * m).sum() / W); out["wmean_ecc"] = float((w * e).sum() / W)
    out["wfrac_mass_gt_0.38"] = float(w[m > 0.38].sum() / W)
    out["wfrac_mass_gt_0.54"] = float(w[m > 0.54].sum() / W)
    out["wfrac_ecc_lt_0.1"] = float(w[e < 0.1].sum() / W)
    out["wfrac_mass_le_0.38_ecc_lt_0.1"] = float(w[(m <= 0.38) & (e < 0.1)].sum() / W)
    # evidence if companions > 0.38 Msun are excluded (manuscript's own overluminosity limit)
    keep = (m <= 0.38)
    out["lnZ_mass_le_0.38"] = float(mx + np.log(max(w[keep].sum(), 1e-300) / a.N))
    # evidence if additionally restricted to near-circular orbits e<0.1
    keep2 = keep & (e < 0.1)
    out["lnZ_mass_le_0.38_e_lt_0.1"] = float(mx + np.log(max(w[keep2].sum(), 1e-300) / a.N))
    out["lnZ_e_lt_0.1"] = float(mx + np.log(max(w[e < 0.1].sum(), 1e-300) / a.N))
    for thr in (0.30, 0.45, 0.54):
        out[f"lnZ_mass_le_{thr:.2f}"] = float(mx + np.log(max(w[m <= thr].sum(), 1e-300) / a.N))
    return out

def fpp_from_lnZ(lnZ):
    lz = np.array(lnZ, float); p = np.exp(lz - np.nanmax(lz[np.isfinite(lz)])); p[~np.isfinite(lz)] = 0; p /= p.sum()
    return float(1 - (p[0] + p[3] + p[9])), p

res_path = os.path.join(OUT, "runs.json")
runs = json.load(open(res_path))["runs"] if os.path.exists(res_path) else []
done = {r["seed"] for r in runs}
for i in range(a.nruns):
    seed = a.seed0 + i
    if seed in done: continue
    np.random.seed(seed); t1 = time.time()
    target.calc_probs(time=lc[:, 0], flux_0=lc[:, 1], flux_err_0=float(sigma), P_orb=P, N=a.N, parallel=True,
                      verbose=0, exptime=a.exptime_min / 1440, molusc_file=a.molusc, dilution_corrected=True)
    pr = target.probs
    row = dict(seed=seed, FPP=float(target.FPP), NFPP=float(target.NFPP), degenerate=bool(target.FPP_degenerate),
               probs={f"{s}_{int(idn)}" if k >= 15 else s: float(p) for k, (s, idn, p) in enumerate(zip(pr.scenario, pr.ID, pr.prob)) if p > 0},
               lnZ=[float(x) for x in target.lnZ], seconds=time.time() - t1)
    row["best_TP"] = dict(R_p=float(pr.R_p[0]), b=float(pr.b[0]), inc=float(pr.inc[0]), ecc=float(pr.ecc[0]))
    row["best_STP"] = dict(M_s=float(pr.M_s[6]), R_s=float(pr.R_s[6]), R_p=float(pr.R_p[6]), ecc=float(pr.ecc[6]), b=float(pr.b[6]))
    if not a.no_record:
        sS = summarize("STP", a.M); sP = summarize("PTP", a.M)
        row["STP_draws"] = sS; row["PTP_draws"] = sP
        lz = list(target.lnZ)
        if sS and "lnZ_mass_le_0.38" in sS:
            lz2 = lz.copy(); lz2[6] = sS["lnZ_mass_le_0.38"]
            if sP and "lnZ_mass_le_0.38" in sP: lz2[3] = sP["lnZ_mass_le_0.38"]
            row["FPP_comp_mass_le_0.38"] = fpp_from_lnZ(lz2)[0]
            lz3 = lz.copy(); lz3[6] = sS["lnZ_mass_le_0.38_e_lt_0.1"]
            if sP and "lnZ_mass_le_0.38" in sP: lz3[3] = sP["lnZ_mass_le_0.38"]
            row["FPP_comp_mass_le_0.38_STP_e_lt_0.1"] = fpp_from_lnZ(lz3)[0]
            for thr in ("0.30", "0.45", "0.54"):
                if f"lnZ_mass_le_{thr}" in sS:
                    lz5 = lz.copy(); lz5[6] = sS[f"lnZ_mass_le_{thr}"]
                    if sP and f"lnZ_mass_le_{thr}" in sP: lz5[3] = sP[f"lnZ_mass_le_{thr}"]
                    row[f"FPP_comp_mass_le_{thr}"] = fpp_from_lnZ(lz5)[0]
            lz4 = lz.copy(); lz4[6] = sS["lnZ_e_lt_0.1"]
            row["FPP_STP_e_lt_0.1"] = fpp_from_lnZ(lz4)[0]
    runs.append(row)
    f = np.array([r["FPP"] for r in runs])
    summ = dict(variant=a.variant, signal=a.signal, fold=a.fold, fold_meta=fm, args=vars(a), sigma_used=sigma, tdepth=tdepth,
                apertures=[x.tolist() for x in aps], triceratops=version("triceratops"), bkg_file=target.trilegal_fname,
                stars_used=stars[["ID", "Tmag", "sep (arcsec)", "mass", "rad", "Teff", "plx", "fluxratio", "tdepth"]].head(15).to_dict("list"),
                n_runs=len(runs), FPP_mean=float(f.mean()), FPP_sd=float(f.std(ddof=1)) if len(f) > 1 else None,
                FPP_min=float(f.min()), FPP_max=float(f.max()), NFPP_mean=float(np.mean([r["NFPP"] for r in runs])), runs=runs)
    json.dump(summ, open(res_path, "w"), indent=1)
    print(a.variant, seed, "FPP=%.4f NFPP=%.2e P_STP=%.4f  (%.0fs)" % (row["FPP"], row["NFPP"], row["probs"].get("STP", 0), row["seconds"]),
          ("FPP(M<=0.38)=%.4f" % row.get("FPP_comp_mass_le_0.38", np.nan)), flush=True)
print("DONE", a.variant, "FPP %.4f +- %.4f" % (f.mean(), f.std(ddof=1) if len(f) > 1 else 0))
