"""Verifier: direct TRICERATOPS 1.1.0 runs reproducing MOLUSC mode WITHOUT the molusc_file code path (no pandas shim):
sample_q_companion -> MOLUSC survivor mass ratios (clamped at 0.1 Msun, zero-padded to N, as marginal_likelihoods.py
l.949-951 does) and lnprior_bound_TP/EB -> 0 (as l.1038), applied only inside the bound-companion scenarios.
variants: coded = survivors with a(1-e) > 10 AU (as TRICERATOPS does); gaiaonly = all MOLUSC survivors (no 10 AU cut).
usage: v07_molusc_direct.py VARIANT SIGNAL FOLD NRUNS   (trice_env, cwd=verify/)"""
import json, os, pickle, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
import triceratops.triceratops as tt
import triceratops.marginal_likelihoods as ml
var, sig, fold, nruns = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
N, M = 1_000_000, 0.606
F = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_run/fresh5_kept.csv"
d = pd.read_csv(F)
sel = (d["semi-major axis(AU)"].to_numpy() * (1 - d["eccentricity"].to_numpy()) > 10) if var == "coded" else np.ones(len(d), bool)
q = d["mass ratio"].to_numpy()[sel].copy(); q[q < 0.1 / M] = 0.1 / M
qpad = np.pad(q, (0, N - len(q)))
T = pickle.load(open(os.path.join(H, "target_v.pkl"), "rb"))
for k, v in dict(mass=M, rad=0.623, Teff=4012.0, plx=28.075).items():
    T.update_star(4206066, k, v)
c98 = np.round(T.pix_coords[1][0])
aps = [np.asarray(T._v_spoc_aps[0], float), np.array([[c98[0] + i, c98[1] + j] for i in (-1, 0, 1) for j in (-1, 0, 1)], float)]
fm = json.load(open(os.path.join(H, "folds", "folds_meta.json")))["folds"][fold]
lc = np.loadtxt(os.path.join(H, "folds", fold + ".csv"), delimiter=",")
T.calc_depths(tdepth=fm["depth_ppm"] * 1e-6, all_ap_pixels=aps)
st = {"on": False}
q0, lb_tp, lb_eb = ml.sample_q_companion, ml.lnprior_bound_TP, ml.lnprior_bound_EB
ml.sample_q_companion = lambda x, M_s: qpad.copy() if st["on"] else q0(x, M_s)
ml.lnprior_bound_TP = lambda M_s, plx, seps: np.zeros(len(seps)) if st["on"] else lb_tp(M_s, plx, seps)
ml.lnprior_bound_EB = lambda M_s, plx, seps: np.zeros(len(seps)) if st["on"] else lb_eb(M_s, plx, seps)
def wrap(fn):
    def w(*a, **k):
        st["on"] = True
        try: return fn(*a, **k)
        finally: st["on"] = False
    return w
for nm in ("lnZ_PTP", "lnZ_PEB", "lnZ_STP", "lnZ_SEB"):
    setattr(tt, nm, wrap(getattr(ml, nm)))
out_path = os.path.join(H, "runs_v", f"molusc_{var}_{sig}_{fold}.json"); rows = []
for i in range(nruns):
    np.random.seed(700 + i); t1 = time.time()
    T.calc_probs(time=lc[:, 0], flux_0=lc[:, 1], flux_err_0=float(fm["sigma_oot"]), P_orb={"A": 3.1827856, "B": 11.13274}[sig],
                 N=N, parallel=True, verbose=0, exptime=fm["binmin"] / 1440.0, dilution_corrected=True)
    rows.append(dict(seed=700 + i, FPP=float(T.FPP), NFPP=float(T.NFPP), STP=float(T.probs.prob[6]), PTP=float(T.probs.prob[3]),
                     TP=float(T.probs.prob[0]), n_q=int(len(q)), sec=time.time() - t1))
    print(var, sig, rows[-1], flush=True)
    json.dump(dict(variant=var, signal=sig, fold=fold, molusc_file=F, runs=rows, FPP_mean=float(np.mean([r["FPP"] for r in rows]))),
              open(out_path, "w"), indent=1)
