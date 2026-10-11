"""Step 13b: TRICERATOPS 1.1.0 runs from this audit's clean inputs (run with the separate trice_env).
usage: s13_fpp_run.py SIGNAL VARIANT FOLD NRUNS [N]
  SIGNAL  A | B
  VARIANT repro   : previous configuration (Teff 4080 K, plx 28.075, exptime 200 s, all TIC neighbours)
          base    : Mann+2015/2019 M, R; Teff 4030 K; Gaia DR3 parallax corrected for the Lindegren+2021 zero point;
                    exptime = 6 min (= bin width of the fold)
          loc     : base + remove every neighbour that could host the signal (all lie >= 43'' away and are excluded by
                    the pixel localisation, s12/s12c, and by the aperture-flux requirement > 100 %)
          molusc5 : loc + MOLUSC survivors (age 5 Gyr, this audit's run)   [TRICERATOPS keeps a(1-e) > 10 AU only]
          molusc1 : loc + MOLUSC survivors (age 1 Gyr)
          tic     : loc but TIC 8.2 stellar parameters (M 0.63, R 0.648, Teff 4010, plx 28.0487)
  FOLD    S32+S98 | S32only | S98only | S6only
Every run: np.random.seed(1000+i); calc_probs(N draws, parallel=True); FPP/NFPP and the full scenario table saved.
FPP = 1 - (P_TP + P_PTP + P_DTP), NFPP = sum over nearby-star scenarios (triceratops.py lines 1697-1701).
"""
import datetime
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")
import numpy as np  # noqa
import triceratops  # noqa
import triceratops.triceratops as tr  # noqa

AUD = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
IN = os.path.join(AUD, "out", "fpp", "inputs")
sig, var, fold, nrun = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
N = int(sys.argv[5]) if len(sys.argv) > 5 else 1_000_000
tag = f"{sig}_{var}_{fold}"
outd = os.path.join(AUD, "out", "fpp", "runs", tag)
os.makedirs(outd, exist_ok=True)
EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
aps_json = json.load(open(os.path.join(IN, "apertures.json")))
DEPTH = {"A": 450e-6, "B": 800e-6}[sig]   # conservative (lower) depths -> more neighbours count as possible hosts
SECTORS = [6, 32, 98]
STAR = {"repro": dict(mass=0.6064, rad=0.6230, Teff=4080.0, plx=28.075),
        "tic": dict(mass=0.630, rad=0.648133, Teff=4010.0, plx=28.0487)}.get(var, dict(mass=0.6062, rad=0.6229, Teff=4030.0, plx=28.1083))
if fold == "S6only":
    EXPT = 30 / 1440
elif var == "repro":
    EXPT = {"S32+S98": 200 / 86400, "S32only": 120 / 86400, "S98only": 200 / 86400}[fold]
else:
    EXPT = 6 / 1440
molusc = {"molusc5": "/private/tmp/claude-501/-Users-paulsmith-space--3/21ae3b23-e2c8-49b9-927a-27d6e2d5b6a3/scratchpad/molusc_age5_kept.csv",
          "molusc1": "/private/tmp/claude-501/-Users-paulsmith-space--3/21ae3b23-e2c8-49b9-927a-27d6e2d5b6a3/scratchpad/molusc_age1_kept.csv"}.get(var)

t0 = time.time()
cache = "/private/tmp/claude-501/-Users-paulsmith-space--3/21ae3b23-e2c8-49b9-927a-27d6e2d5b6a3/scratchpad/trice_lk_cache"  # large TESScut cache kept out of the audit folder
os.makedirs(cache, exist_ok=True)
target = tr.target(ID=4206066, sectors=np.array(SECTORS), search_radius=10, lightkurve_cache_dir=cache)
log = dict(tag=tag, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), triceratops_version=getattr(triceratops, "__version__", "1.1.0"),
           N=N, exptime_d=EXPT, depth_for_calc_depths=DEPTH, star_params=STAR, n_tic_stars=int(len(target.stars)),
           tic_target_row_before=target.stars.iloc[0][["ID", "Tmag", "mass", "rad", "Teff", "plx"]].to_dict(),
           background_file=str(target.trilegal_fname), molusc_file=molusc)
for k, v in STAR.items():
    target.update_star(4206066, k, v)
aps = [np.array(aps_json[str(s)]["pixels"]) for s in SECTORS]
target.calc_depths(tdepth=DEPTH, all_ap_pixels=aps)
cap = target.stars[target.stars.tdepth > 0][["ID", "Tmag", "sep (arcsec)", "PA (E of N)", "fluxratio", "tdepth"]]
log["possible_hosts_before_removal"] = cap.to_dict("records")
if var in ("loc", "molusc5", "molusc1", "tic"):
    drop = [int(i) for i in cap["ID"].values[1:]]
    if drop:
        target.remove_star(np.array(drop))
    log["removed_neighbours"] = drop
lc = np.loadtxt(os.path.join(IN, f"{sig}_{fold}_binned.csv"), delimiter=",")
tt, ff, ee = lc[:, 0], lc[:, 1], lc[:, 2]
runs = []
for i in range(nrun):
    np.random.seed(1000 + i)
    t1 = time.time()
    target.calc_probs(time=tt, flux_0=ff, flux_err_0=float(np.mean(ee)), P_orb=EPH[sig]["P"], N=N, parallel=True, verbose=0,
                      exptime=EXPT, molusc_file=molusc, dilution_corrected=True)
    pr = target.probs.copy()
    pr.to_csv(os.path.join(outd, f"probs_run{i:02d}.csv"), index=False)
    runs.append(dict(run=i, seed=1000 + i, FPP=float(target.FPP), NFPP=float(target.NFPP), degenerate=bool(getattr(target, "FPP_degenerate", False)),
                     seconds=time.time() - t1, P_TP=float(pr.prob[0]), P_PTP=float(pr.prob[3]), P_STP=float(pr.prob[6]), P_DTP=float(pr.prob[9]),
                     P_EB_all=float(pr.prob[[1, 2, 4, 5, 7, 8, 10, 11, 13, 14]].sum())))
    print(tag, i, "FPP=%.4g NFPP=%.3g (%.0fs)" % (target.FPP, target.NFPP, time.time() - t1), flush=True)
    json.dump(dict(log=log, runs=runs), open(os.path.join(outd, "summary.json"), "w"), indent=1, default=str)
f = np.array([r["FPP"] for r in runs])
n = np.array([r["NFPP"] for r in runs])
log.update(finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), total_s=time.time() - t0,
           FPP_mean=float(f.mean()), FPP_std=float(f.std(ddof=1)), FPP_median=float(np.median(f)), FPP_min=float(f.min()),
           FPP_max=float(f.max()), FPP_sem=float(f.std(ddof=1) / np.sqrt(len(f))), NFPP_mean=float(n.mean()), NFPP_max=float(n.max()),
           n_degenerate=int(sum(r["degenerate"] for r in runs)), n_saturated=int((f > 0.99).sum()))
json.dump(dict(log=log, runs=runs), open(os.path.join(outd, "summary.json"), "w"), indent=1, default=str)
print("DONE", tag, "FPP %.4g ± %.2g (sd), sem %.2g; range %.3g-%.3g; NFPP mean %.2g" % (f.mean(), f.std(ddof=1), log["FPP_sem"], f.min(), f.max(), n.mean()))
