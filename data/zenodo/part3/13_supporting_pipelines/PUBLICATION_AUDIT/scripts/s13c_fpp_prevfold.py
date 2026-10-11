"""Step 13c: does the earlier B fold reproduce FPP_B = 0.146, and is the 0.146 -> 0.18 change due to the fold?
usage: s13c_fpp_prevfold.py VARIANT NRUNS [N]      (run with the separate trice_env, after s13_fpp_run.py)
  prevfold_repro : exact earlier configuration (deep_followup_tic4206066/fpp/run_fpp.py, variant A, fold S32+S98):
                   earlier fold file and ephemeris, calc_depths depth 783 ppm, Teff 4080 K, plx 28.075, exptime 200 s,
                   all TIC neighbours kept.
  prevfold_base  : this audit's 'base' configuration (s13_fpp_run.py) with only the fold swapped for the earlier one.
Together with B_repro (audit fold + earlier configuration) and B_base (audit fold + audit configuration) this separates
the effect of the fold from the effect of the stellar parameters / exposure time.
The earlier input files are only read.
"""
import datetime
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")
import numpy as np  # noqa
import triceratops.triceratops as tr  # noqa

AUD = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OLD = os.path.abspath(os.path.join(AUD, "..", "deep_followup_tic4206066", "fpp"))
var, nrun = sys.argv[1], int(sys.argv[2])
N = int(sys.argv[3]) if len(sys.argv) > 3 else 1_000_000
tag = f"B_{var}_S32+S98"
outd = os.path.join(AUD, "out", "fpp", "runs", tag)
os.makedirs(outd, exist_ok=True)
old_meta = json.load(open(os.path.join(OLD, "inputs", "inputs_meta.json")))["signals"]["p11.13"]
old_fold = os.path.join(OLD, "inputs", "p11.13_S32+S98_binned.csv")
aps_json = json.load(open(os.path.join(OLD, "inputs", "apertures.json")))
EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
SECTORS = [6, 32, 98]
if var == "prevfold_repro":
    STAR, EXPT, DEPTH, P = dict(mass=0.6064, rad=0.6230, Teff=4080.0, plx=28.075), 200 / 86400, old_meta["depth"], old_meta["P"]
elif var == "prevfold_base":
    STAR, EXPT, DEPTH, P = dict(mass=0.6062, rad=0.6229, Teff=4030.0, plx=28.1083), 6 / 1440, 800e-6, EPH["B"]["P"]
else:
    sys.exit("unknown variant")

t0 = time.time()
cache = "/private/tmp/claude-501/-Users-paulsmith-space--3/21ae3b23-e2c8-49b9-927a-27d6e2d5b6a3/scratchpad/trice_lk_cache"
target = tr.target(ID=4206066, sectors=np.array(SECTORS), search_radius=10, lightkurve_cache_dir=cache)
log = dict(tag=tag, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), N=N, exptime_d=EXPT,
           depth_for_calc_depths=DEPTH, P_orb=P, star_params=STAR, fold_file=old_fold, n_tic_stars=int(len(target.stars)))
for k, v in STAR.items():
    target.update_star(4206066, k, v)
target.calc_depths(tdepth=DEPTH, all_ap_pixels=[np.array(aps_json[str(s)]["pixels"]) for s in SECTORS])
log["possible_hosts"] = target.stars[target.stars.tdepth > 0][["ID", "Tmag", "sep (arcsec)", "tdepth"]].to_dict("records")
lc = np.loadtxt(old_fold, delimiter=",")
tt, ff, ee = lc[:, 0], lc[:, 1], lc[:, 2]
runs = []
for i in range(nrun):
    np.random.seed(1000 + i)
    t1 = time.time()
    target.calc_probs(time=tt, flux_0=ff, flux_err_0=float(np.mean(ee)), P_orb=P, N=N, parallel=True, verbose=0,
                      exptime=EXPT, molusc_file=None, dilution_corrected=True)
    pr = target.probs.copy()
    pr.to_csv(os.path.join(outd, f"probs_run{i:02d}.csv"), index=False)
    runs.append(dict(run=i, seed=1000 + i, FPP=float(target.FPP), NFPP=float(target.NFPP), seconds=time.time() - t1,
                     degenerate=bool(getattr(target, "FPP_degenerate", False)),
                     P_TP=float(pr.prob[0]), P_PTP=float(pr.prob[3]), P_STP=float(pr.prob[6]), P_DTP=float(pr.prob[9]),
                     P_EB_all=float(pr.prob[[1, 2, 4, 5, 7, 8, 10, 11, 13, 14]].sum())))
    print(tag, i, "FPP=%.4g NFPP=%.3g (%.0fs)" % (target.FPP, target.NFPP, time.time() - t1), flush=True)
    json.dump(dict(log=log, runs=runs), open(os.path.join(outd, "summary.json"), "w"), indent=1, default=str)
f = np.array([r["FPP"] for r in runs])
nf = np.array([r["NFPP"] for r in runs])
old_runs = json.load(open(os.path.join(OLD, "runs", "p11.13_A_S32+S98", "summary.json")))["runs"]
old_f = np.array([r["FPP"] for r in old_runs])
log.update(finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), total_s=time.time() - t0,
           FPP_mean=float(f.mean()), FPP_std=float(f.std(ddof=1)), FPP_min=float(f.min()), FPP_max=float(f.max()),
           NFPP_mean=float(nf.mean()), NFPP_max=float(nf.max()), n_degenerate=int(sum(r["degenerate"] for r in runs)),
           n_saturated=int((f > 0.99).sum()),
           earlier_FPP_mean=float(old_f.mean()),
           max_abs_diff_vs_earlier_same_seed=float(np.max(np.abs(f - old_f[:len(f)]))) if var == "prevfold_repro" else None)
json.dump(dict(log=log, runs=runs), open(os.path.join(outd, "summary.json"), "w"), indent=1, default=str)
print("DONE", tag, "FPP %.4g ± %.2g (sd); earlier %.4g" % (f.mean(), f.std(ddof=1), old_f.mean()),
      "max |diff| same seed %.2g" % log["max_abs_diff_vs_earlier_same_seed"] if var == "prevfold_repro" else "")
