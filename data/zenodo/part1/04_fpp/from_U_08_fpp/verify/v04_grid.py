"""Verifier: mean likelihood of the STP and PTP scenarios at FIXED companion mass, L(m) = E[L | m],
by calling TRICERATOPS 1.1.0's own lnZ_STP / lnZ_PTP with sample_q_companion replaced by a constant and the
companion prior lnprior_bound_TP replaced by 0 (so lnZ returned = ln E[L|m]).  Other draws (Rp, inc, ecc, w) are
TRICERATOPS' own priors.  usage: v04_grid.py SIGNAL FOLD [--N 4e5]
"""
import argparse, json, os, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
import triceratops.marginal_likelihoods as ml

p = argparse.ArgumentParser(); p.add_argument("signal"); p.add_argument("fold"); p.add_argument("--N", type=float, default=4e5)
p.add_argument("--Teff", type=float, default=4012.0); p.add_argument("--tag", default=""); p.add_argument("--seedoff", type=int, default=1000)
a = p.parse_args(); N = int(a.N)
M, R, plx = 0.606, 0.623, 28.075
fm = json.load(open(os.path.join(H, "folds", "folds_meta.json")))["folds"][a.fold]
lc = np.loadtxt(os.path.join(H, "folds", a.fold + ".csv"), delimiter=",")
P = {"A": 3.1827856, "B": 11.13274}[a.signal]
cur = {"q": None}
ml.sample_q_companion = lambda x, M_s: np.full(len(x), cur["q"])
ml.lnprior_bound_TP = lambda M_s, plx, seps: np.zeros(len(seps))
masses = list(np.round(np.arange(0.10, 0.581, 0.02), 3)) + [0.59, 0.60, 0.606]
out_path = os.path.join(H, "runs_v", f"grid_{a.signal}_{a.fold}{a.tag}.json")
res = json.load(open(out_path)) if os.path.exists(out_path) else {"rows": []}
done = {r["m"] for r in res["rows"]}
for i, m in enumerate(masses):
    if m in done:
        continue
    cur["q"] = m / M
    row = {"m": m}
    for scen, fn in (("STP", ml.lnZ_STP), ("PTP", ml.lnZ_PTP)):
        np.random.seed(a.seedoff + i); t1 = time.time()
        r = fn(lc[:, 0], lc[:, 1], float(fm["sigma_oot"]), P, M, R, a.Teff, 0.0, plx, None, "TESS", N, True, "TESS",
               False, fm["binmin"] / 1440.0, 20, None)
        row["lnL_" + scen] = float(r["lnZ"]); row["sec_" + scen] = time.time() - t1
        if scen == "STP":
            row["best_ecc"] = float(r["ecc"][0]); row["n_top_ecc_gt_0.1"] = int(np.sum(r["ecc"] > 0.1))
    res["rows"].append(row)
    res.update(signal=a.signal, fold=a.fold, N=N, M=M, R=R, Teff=a.Teff, plx=plx)
    json.dump(res, open(out_path, "w"), indent=1)
    print(a.signal, "m=%.3f lnL_STP=%.2f lnL_PTP=%.2f (%.0f+%.0fs) best_ecc=%.2f" % (
        m, row["lnL_STP"], row["lnL_PTP"], row["sec_STP"], row["sec_PTP"], row["best_ecc"]), flush=True)
print("DONE grid", a.signal)
