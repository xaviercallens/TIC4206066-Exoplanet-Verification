"""R5: direct TRICERATOPS 1.1.0 runs on the ADOPTED fold (the manuscript headline fold: tic4206066_PUBLICATION_AUDIT
out/fpp/inputs/{A,B}_S32+S98_binned.csv, 6-min bins, sigma = mean bin error, exptime 6 min, star M 0.6062 R 0.6229
Teff 4030 plx 28.1083, dilution_corrected=True, N = 1e6, seeds 1000+i as in s13_fpp_run.py).
Target object: the verifier's TRICERATOPS target pickle (08_fpp/verify/target_v.pkl; read only), star updated.

modes
  default                 standard calc_probs
  ol<M>_all               overluminosity limit M applied at every separation <= 2.2" (= verifier trunc<M>)
  ol<M>_bl0.2             overluminosity limit applied only to Gaia-blended companions (projected sep < 0.2", Way+2026)
  gaia_ol<M>_bl0.2        as above + MOLUSC (Wood+2021) RUWE and Gaia-contrast survival at all separations
  gaia_only               MOLUSC RUWE + Gaia-contrast survival only
  default_teff4012        default, but Teff 4012 K (ARIADNE value in the manuscript) instead of 4030 K
Implementation: in PTP/PEB/STP/SEB the companion mass is drawn from TRICERATOPS' own sample_q_companion and accepted
with probability s(m) = sum_a dF(a) S(m,a) / f_comp(<2.2") (S from r03: TRICERATOPS default separation prior x MOLUSC
survival x overluminosity rule); ln E_p[s] is then added to lnZ of those scenarios (prior mass removed, not renormalised).
usage: r05_direct_runs.py MODE SIGNAL NRUNS      (trice_env, cwd = this folder)
"""
import json, os, pickle, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
import triceratops.triceratops as tt
import triceratops.marginal_likelihoods as ml
from triceratops.priors import sample_q_companion as q_orig, lnprior_bound_TP

mode, sig, nrun = sys.argv[1], sys.argv[2], int(sys.argv[3])
AUD = os.path.abspath(os.path.join(H, "..", "..", "..", "tic4206066_PUBLICATION_AUDIT"))
V = os.path.abspath(os.path.join(H, "..", "..", "08_fpp", "verify"))
STAR = dict(mass=0.6062, rad=0.6229, Teff=4012.0 if mode == "default_teff4012" else 4030.0, plx=28.1083)
M = STAR["mass"]
EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
lc = np.loadtxt(os.path.join(AUD, "out", "fpp", "inputs", f"{sig}_S32+S98_binned.csv"), delimiter=",")
T = pickle.load(open(os.path.join(V, "target_v.pkl"), "rb"))
for k, v in STAR.items():
    T.update_star(4206066, k, v)
c98 = np.round(T.pix_coords[1][0])
aps = [np.asarray(T._v_spoc_aps[0], float), np.array([[c98[0] + i, c98[1] + j] for i in (-1, 0, 1) for j in (-1, 0, 1)], float)]
T.calc_depths(tdepth={"A": 450e-6, "B": 800e-6}[sig], all_ap_pixels=aps)

# ---------------- survival s(m) ----------------
s_of_m = None
if not mode.startswith("default"):
    plx = STAR["plx"]; d_pc = 1000 / plx
    df = pd.read_csv("/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/r5_cache/r5_all.csv",
                     usecols=["mass ratio", "semi-major axis(AU)", "Projected Separation(AU)", "Full Rejected"])
    mm = df["mass ratio"].to_numpy() * 0.606; aa = df["semi-major axis(AU)"].to_numpy()
    rho = df["Projected Separation(AU)"].to_numpy() / d_pc; rej_m = df["Full Rejected"].to_numpy().astype(bool)
    a_edges = np.geomspace(3.0, 2.2 * d_pc, 31)
    Fa = np.exp(lnprior_bound_TP(M, plx, a_edges / d_pc)); Fa[~np.isfinite(Fa)] = 0.0; dF = np.diff(Fa); f22 = Fa[-1]
    m_edges = np.round(np.arange(0.10, 0.6061, 0.02), 3); m_edges[-1] = 0.6061; m_c = 0.5 * (m_edges[1:] + m_edges[:-1])
    if mode.startswith("ol") and mode.endswith("_all"):
        Mlim = float(mode[2:].split("_")[0]); s_of_m = lambda m: (m <= Mlim).astype(float)
    else:
        use_gaia = mode.startswith("gaia")
        Mlim = float(mode.split("ol")[1].split("_")[0]) if "ol" in mode else 9.0
        rho_b = float(mode.split("_bl")[1]) if "_bl" in mode else 0.0
        rej = (rej_m if use_gaia else np.zeros(len(df), bool)) | ((rho < rho_b) & (mm > Mlim))
        sel = (mm >= 0.1) & (aa >= a_edges[0]) & (aa < a_edges[-1])
        im = np.clip(np.digitize(mm[sel], m_edges) - 1, 0, len(m_c) - 1); ia = np.clip(np.digitize(aa[sel], a_edges) - 1, 0, len(dF) - 1)
        n = np.zeros((len(m_c), len(dF))); s = np.zeros_like(n)
        np.add.at(n, (im, ia), 1); np.add.at(s, (im, ia), (~rej[sel]).astype(float))
        S = np.where(n > 0, s / np.maximum(n, 1), 1.0); sm = (S @ dF) / f22
        s_of_m = lambda m, sm=sm: np.interp(m, m_c, sm)
    rs = np.random.RandomState(12345); m_big = q_orig(rs.rand(4_000_000), M) * M
    Es = float(np.mean(s_of_m(m_big)))
state = {"on": False}
def q_patched(x, M_s):
    if not state["on"]:
        return q_orig(x, M_s)
    n = len(x); out = np.empty(0)
    while out.size < n:
        qd = q_orig(np.random.rand(2 * n), M_s)
        acc = np.random.rand(qd.size) < s_of_m(qd * M_s)
        out = np.concatenate([out, qd[acc]])
    return out[:n]
if s_of_m is not None:
    ml.sample_q_companion = q_patched
    def wrap(fn):
        def w(*args, **kw):
            state["on"] = True
            try: return fn(*args, **kw)
            finally: state["on"] = False
        return w
    for nm in ("lnZ_PTP", "lnZ_PEB", "lnZ_STP", "lnZ_SEB"):
        setattr(tt, nm, wrap(getattr(ml, nm)))

os.makedirs(os.path.join(H, "runs"), exist_ok=True)
out_path = os.path.join(H, "runs", f"{mode}_{sig}_adoptedfold.json")
res = json.load(open(out_path)) if os.path.exists(out_path) else {"runs": []}
for i in range(nrun):
    seed = 1000 + i
    if any(r["seed"] == seed for r in res["runs"]): continue
    np.random.seed(seed); t1 = time.time()
    T.calc_probs(time=lc[:, 0], flux_0=lc[:, 1], flux_err_0=float(np.mean(lc[:, 2])), P_orb=EPH[sig]["P"], N=1_000_000,
                 parallel=True, verbose=0, exptime=6 / 1440, dilution_corrected=True)
    lnZ = np.array(T.lnZ, float); scen = list(T.probs.scenario)
    assert scen[0] == "TP" and scen[3] == "PTP" and scen[6] == "STP" and scen[9] == "DTP"
    if s_of_m is not None:
        for j in (3, 4, 5, 6, 7, 8):
            lnZ[j] += np.log(Es)
    fin = np.isfinite(lnZ); pr = np.zeros_like(lnZ); pr[fin] = np.exp(lnZ[fin] - lnZ[fin].max()); pr /= pr.sum()
    row = dict(seed=seed, FPP=float(1 - pr[0] - pr[3] - pr[9]), NFPP=float(pr[15:].sum()), FPP_trice_raw=float(T.FPP),
               P={s: float(pr[k]) for k, s in enumerate(scen[:15]) if pr[k] > 1e-8}, lnZ=lnZ.tolist(), sec=time.time() - t1)
    res["runs"].append(row)
    f = np.array([r["FPP"] for r in res["runs"]])
    res.update(mode=mode, signal=sig, fold=f"PUBLICATION_AUDIT {sig}_S32+S98_binned.csv (6-min)", star=STAR, exptime_min=6, N=1_000_000,
               E_s=None if s_of_m is None else Es, FPP_mean=float(f.mean()), FPP_sd=float(f.std(ddof=1)) if f.size > 1 else None,
               FPP_min=float(f.min()), FPP_max=float(f.max()), NFPP_max=float(max(r["NFPP"] for r in res["runs"])))
    json.dump(res, open(out_path, "w"), indent=1)
    print(mode, sig, seed, "FPP=%.4f NFPP=%.1e STP=%.4f PTP=%.4f %.0fs" % (row["FPP"], row["NFPP"], row["P"].get("STP", 0), row["P"].get("PTP", 0), row["sec"]), flush=True)
print("DONE", mode, sig, "FPP %.4f +- %s  E_s=%s" % (res["FPP_mean"], res["FPP_sd"], res["E_s"]))
