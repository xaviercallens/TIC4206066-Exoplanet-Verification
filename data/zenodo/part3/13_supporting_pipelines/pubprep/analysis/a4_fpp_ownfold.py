"""TRICERATOPS 1.1.0 with an independently built fold (third fold, to measure the fold sensitivity of the FPP).
Run with the TRICERATOPS environment (numpy 1.23):  trice_env/bin/python a4_fpp_ownfold.py A|B NRUNS
Fold: my detrended light curves (a1_timing_depths: SPOC S32 2-min PDCSAP + QLP S98 frozen mask; biweight 0.75 d with
both signals masked), the other signal removed (+-1 T14), folded on the linear ephemeris of a1 (injection-error fit),
+-0.25 d, 6-min inverse-variance bins. Star: Mann M, R; Teff 4030 K; ZP-corrected parallax (as the audit baseline).
Apertures: the audit's apertures.json (SPOC optimal apertures S6/S32, QLP 3x3 S98), which equal the deep follow-up's.
Neighbours are kept (baseline variant). Seeds 2000+i.
"""
import json, os, sys, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "out", "fpp"); os.makedirs(OUT, exist_ok=True)
args = [a for a in sys.argv[1:] if not a.startswith("--")]; sig, nrun = args[0], int(args[1]); N = int(args[2]) if len(args) > 2 else 1_000_000
if "--build" in sys.argv or not os.path.exists(os.path.join(OUT, f"fold_{sig}.csv")):
    import a1_timing_depths as a1
    R = json.load(open(os.path.join(H, "out", "a1_results.json")))
    key = "A__pipe" if sig == "A" else "B_mann__pipe"
    eph = R["signals"][sig][key]["ephemeris"]["e_dt_inj_min"]; P, T0 = eph["P"], eph["T0"]; T14 = R["signals"][sig][key]["template_T14_h"] / 24
    other = "B" if sig == "A" else "A"
    T, F, E = [], [], []
    for name in ("S32", "S98"):
        L = a1.detrend(a1.load_lc(name)); keep = ~a1.in_tr(L["t"], other, 1.0)
        ph = (((L["t"] - T0) / P + 0.5) % 1 - 0.5) * P; w = keep & (np.abs(ph) < 0.25)
        T.append(ph[w]); F.append(L["fd"][w]); E.append(np.full(w.sum(), L["sig"]))
    t, f, e = np.concatenate(T), np.concatenate(F), np.concatenate(E); edges = np.arange(-0.25, 0.25 + 1e-9, 6 / 1440); idx = np.digitize(t, edges)
    rows = []
    for i in range(1, len(edges)):
        m = idx == i
        if m.sum() < 3: continue
        ww = 1 / e[m] ** 2; rows.append((t[m].mean(), np.sum(ww * f[m]) / ww.sum(), 1 / np.sqrt(ww.sum())))
    np.savetxt(os.path.join(OUT, f"fold_{sig}.csv"), np.array(rows), delimiter=",", header="t_from_mid_d,flux,flux_err")
    json.dump(dict(P=P, T0=T0, T14_d=T14, n_raw=int(len(t)), n_bins=len(rows)), open(os.path.join(OUT, f"fold_{sig}_meta.json"), "w"), indent=1)
    print("fold built", sig, len(rows), "bins"); sys.exit(0)
import triceratops, triceratops.triceratops as tr
AUDD = os.path.abspath(os.path.join(H, "..", "..", "tic4206066_PUBLICATION_AUDIT"))
aps_json = json.load(open(os.path.join(AUDD, "out", "fpp", "inputs", "apertures.json")))
meta = json.load(open(os.path.join(OUT, f"fold_{sig}_meta.json")))
cache = "/private/tmp/claude-501/-Users-paulsmith-space--3/21ae3b23-e2c8-49b9-927a-27d6e2d5b6a3/scratchpad/trice_lk_cache"
target = tr.target(ID=4206066, sectors=np.array([6, 32, 98]), search_radius=10, lightkurve_cache_dir=cache)
for k, v in dict(mass=0.6062, rad=0.6229, Teff=4030.0, plx=28.1083).items(): target.update_star(4206066, k, v)
target.calc_depths(tdepth={"A": 450e-6, "B": 800e-6}[sig], all_ap_pixels=[np.array(aps_json[str(s)]["pixels"]) for s in (6, 32, 98)])
lc = np.loadtxt(os.path.join(OUT, f"fold_{sig}.csv"), delimiter=","); runs = []
for i in range(nrun):
    np.random.seed(2000 + i); t1 = time.time()
    target.calc_probs(time=lc[:, 0], flux_0=lc[:, 1], flux_err_0=float(np.mean(lc[:, 2])), P_orb=meta["P"], N=N, parallel=True, verbose=0, exptime=6 / 1440, dilution_corrected=True)
    pr = target.probs; runs.append(dict(seed=2000 + i, FPP=float(target.FPP), NFPP=float(target.NFPP), P_STP=float(pr.prob[6]), seconds=time.time() - t1))
    print(sig, i, "FPP=%.4g NFPP=%.3g" % (target.FPP, target.NFPP), flush=True)
    f = np.array([r["FPP"] for r in runs])
    json.dump(dict(signal=sig, fold=meta, N=N, star=dict(mass=0.6062, rad=0.6229, Teff=4030.0, plx=28.1083), runs=runs, FPP_mean=float(f.mean()), FPP_sd=float(f.std(ddof=1)) if len(f) > 1 else None,
                   FPP_min=float(f.min()), FPP_max=float(f.max()), triceratops=__import__('importlib.metadata').metadata.version('triceratops')), open(os.path.join(OUT, f"fpp_ownfold_{sig}.json"), "w"), indent=1)
print("DONE", sig, f.mean(), f.std(ddof=1))
