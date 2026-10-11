"""Verifier: TRICERATOPS 1.1.0 runs on own folds and own target (trice_env, cwd=verify/).
modes:
  default         : standard calc_probs
  trunc<Mmax>     : companion-mass prior truncated at Mmax (Msun) for the bound-companion scenarios
                    (PTP, PEB, STP, SEB): q drawn from TRICERATOPS' own sample_q_companion by rejection, then
                    ln P(M<=Mmax) added to those lnZ (so the prior mass is removed, not renormalised).
usage: v03_runs.py MODE SIGNAL FOLD NRUNS [--N 1e6] [--seed0 500] [--Teff 4012]
"""
import argparse, json, os, pickle, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
import triceratops.triceratops as tt
import triceratops.marginal_likelihoods as ml
from triceratops.priors import sample_q_companion as q_orig

p = argparse.ArgumentParser()
p.add_argument("mode"); p.add_argument("signal"); p.add_argument("fold"); p.add_argument("nruns", type=int)
p.add_argument("--N", type=float, default=1e6); p.add_argument("--seed0", type=int, default=500)
p.add_argument("--Teff", type=float, default=4012.0); p.add_argument("--M", type=float, default=0.606)
p.add_argument("--R", type=float, default=0.623); p.add_argument("--plx", type=float, default=28.075)
a = p.parse_args(); N = int(a.N)

T = pickle.load(open(os.path.join(H, "target_v.pkl"), "rb"))
for k, v in dict(mass=a.M, rad=a.R, Teff=a.Teff, plx=a.plx).items():
    T.update_star(4206066, k, v)
# apertures: S32 = SPOC optimal aperture returned by TRICERATOPS' own get_spoc_apertures(); S98 = QLP BESTAP 3x3
c98 = np.round(T.pix_coords[1][0])
ap98 = np.array([[c98[0] + i, c98[1] + j] for i in (-1, 0, 1) for j in (-1, 0, 1)], float)
aps = [np.asarray(T._v_spoc_aps[0], float), ap98]

fm = json.load(open(os.path.join(H, "folds", "folds_meta.json")))["folds"][a.fold]
lc = np.loadtxt(os.path.join(H, "folds", a.fold + ".csv"), delimiter=",")
binmin = fm["binmin"]
T.calc_depths(tdepth=fm["depth_ppm"] * 1e-6, all_ap_pixels=aps)
P = {"A": 3.1827856, "B": 11.13274}[a.signal]

state = {"trunc": False}
Mmax = float(a.mode[5:]) if a.mode.startswith("trunc") else None
if Mmax is not None:
    # probability mass of the default q prior below Mmax (large sample, fixed seed, own RNG)
    rs = np.random.RandomState(12345)
    qq = q_orig(rs.rand(4_000_000), a.M)
    frac_le = float(np.mean(qq * a.M <= Mmax))

    def q_patched(x, M_s):
        if not state["trunc"]:
            return q_orig(x, M_s)
        n = len(x); out = np.empty(0)
        while out.size < n:
            qd = q_orig(np.random.rand(2 * n), M_s)
            out = np.concatenate([out, qd[qd * M_s <= Mmax]])
        return out[:n]
    ml.sample_q_companion = q_patched

    def wrap(fn):
        def w(*args, **kw):
            state["trunc"] = True
            try:
                return fn(*args, **kw)
            finally:
                state["trunc"] = False
        return w
    for nm in ("lnZ_PTP", "lnZ_PEB", "lnZ_STP", "lnZ_SEB"):
        setattr(tt, nm, wrap(getattr(ml, nm)))

out_path = os.path.join(H, "runs_v", f"{a.mode}_{a.signal}_{a.fold}.json")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
res = json.load(open(out_path)) if os.path.exists(out_path) else {"runs": []}
for i in range(a.nruns):
    seed = a.seed0 + i
    if any(r["seed"] == seed for r in res["runs"]):
        continue
    np.random.seed(seed); t1 = time.time()
    T.calc_probs(time=lc[:, 0], flux_0=lc[:, 1], flux_err_0=float(fm["sigma_oot"]), P_orb=P, N=N, parallel=True,
                 verbose=0, exptime=binmin / 1440.0, dilution_corrected=True)
    lnZ = np.array(T.lnZ, float); scen = list(T.probs.scenario)
    if Mmax is not None:
        for j in (3, 4, 5, 6, 7, 8):
            assert scen[j] in ("PTP", "PEB", "PEBx2P", "STP", "SEB", "SEBx2P"), scen[j]
            lnZ[j] = lnZ[j] + np.log(frac_le)
    fin = np.isfinite(lnZ); pr = np.zeros_like(lnZ); pr[fin] = np.exp(lnZ[fin] - lnZ[fin].max()); pr /= pr.sum()
    assert scen[0] == "TP" and scen[3] == "PTP" and scen[9] == "DTP"
    FPP = float(1 - pr[0] - pr[3] - pr[9]); NFPP = float(pr[15:].sum())
    row = dict(seed=seed, FPP=FPP, NFPP=NFPP, FPP_trice=float(T.FPP),
               P={s: float(pr[k]) for k, s in enumerate(scen[:15]) if pr[k] > 1e-8}, lnZ=lnZ.tolist(), sec=time.time() - t1)
    res["runs"].append(row)
    f = np.array([r["FPP"] for r in res["runs"]])
    res.update(mode=a.mode, signal=a.signal, fold=a.fold, N=N, M=a.M, R=a.R, Teff=a.Teff, plx=a.plx, sigma=fm["sigma_oot"],
               depth_ppm=fm["depth_ppm"], exptime_min=binmin, aps=[x.tolist() for x in aps],
               frac_le=None if Mmax is None else frac_le, FPP_mean=float(f.mean()),
               FPP_sd=float(f.std(ddof=1)) if f.size > 1 else None, FPP_min=float(f.min()), FPP_max=float(f.max()))
    json.dump(res, open(out_path, "w"), indent=1)
    print(a.mode, a.signal, seed, "FPP=%.4f (trice %.4f) NFPP=%.1e STP=%.4f PTP=%.4f TP=%.4f %.0fs" % (
        FPP, T.FPP, NFPP, row["P"].get("STP", 0), row["P"].get("PTP", 0), row["P"].get("TP", 0), row["sec"]), flush=True)
print("DONE", a.mode, a.signal, a.fold, "FPP %.4f +- %s" % (res["FPP_mean"], res["FPP_sd"]))
