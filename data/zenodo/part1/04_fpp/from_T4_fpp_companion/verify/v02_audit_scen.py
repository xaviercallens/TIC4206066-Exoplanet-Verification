"""Verifier: scenario totals from the PUBLICATION_AUDIT probs_run*.csv files (read-only) and deep_followup runs:
EB total, STP/FPP, BTP+BEB, NFPP max; to test C085 and the 'STP >= 99.66%' / 'BTP+BEB <= 1.2e-4' claims beyond U08."""
import glob, os, json
import numpy as np, pandas as pd
A = "/Users/paulsmith/space <3/tess-jev/addenda"
out = {}
for d in sorted(glob.glob(A + "/tic4206066_PUBLICATION_AUDIT/out/fpp/runs/A_*")):
    eb, stp, btpbeb, nf, fpp = [], [], [], [], []
    for f in sorted(glob.glob(d + "/probs_run*.csv")):
        p = pd.read_csv(f); s = dict(zip(p.scenario, p.prob))
        F = 1 - s["TP"] - s["PTP"] - s["DTP"]
        fpp.append(F); eb.append(sum(v for k, v in s.items() if "EB" in k)); stp.append(s["STP"] / F)
        btpbeb.append(s.get("BTP", 0) + s.get("BEB", 0) + s.get("BEBx2P", 0)); nf.append(sum(v for k, v in s.items() if k.startswith("N")))
    out[os.path.basename(d)] = dict(n=len(fpp), FPP_mean=float(np.mean(fpp)), EB_mean=float(np.mean(eb)), EB_max=float(np.max(eb)),
                                    STP_over_FPP_mean=float(np.mean(stp)), STP_over_FPP_min=float(np.min(stp)), BTP_BEB_mean=float(np.mean(btpbeb)), BTP_BEB_max=float(np.max(btpbeb)), NFPP_max=float(np.max(nf)))
json.dump(out, open("v02_audit_scen.json", "w"), indent=1)
for k, v in out.items(): print(k, {a: (round(b, 7) if isinstance(b, float) else b) for a, b in v.items()})
