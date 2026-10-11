"""C085: total probability of all EB scenarios (keys containing 'EB') and scenario shares, from stored U/08_fpp run files.
Also records mean P(STP), P(DTP), P(BTP), P(BEB) for signal A default runs (what TRICERATOPS covers)."""
import json, glob, os
import numpy as np
U = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066"
out = {}
for f in sorted(glob.glob(U + "/08_fpp/runs/*/runs.json")):
    d = json.load(open(f)); runs = d["runs"]
    eb = [sum(v for k, v in r["probs"].items() if "EB" in k) for r in runs]
    g = lambda key: float(np.mean([r["probs"].get(key, 0.0) for r in runs]))
    nb = [sum(v for k, v in r["probs"].items() if k.startswith("N")) for r in runs]
    out[d["variant"]] = dict(n=len(runs), P_EB_all_mean=float(np.mean(eb)), P_EB_all_max=float(np.max(eb)),
                             STP=g("STP"), PTP=g("PTP"), DTP=g("DTP"), BTP=g("BTP"), BEB=g("BEB"), TP=g("TP"),
                             nearby_all_mean=float(np.mean(nb)), STP_over_FPP=float(np.mean([r["probs"]["STP"] / r["FPP"] for r in runs])))
json.dump(out, open(os.path.dirname(os.path.abspath(__file__)) + "/t4_06_eb_totals.json", "w"), indent=1)
for k, v in out.items():
    print(f"{k:22s} n={v['n']:2d} EBall={v['P_EB_all_mean']:.2e} (max {v['P_EB_all_max']:.1e}) STP={v['STP']:.4f} STP/FPP={v['STP_over_FPP']:.4f} DTP={v['DTP']:.4f} BTP={v['BTP']:.1e} BEB={v['BEB']:.1e} Nall={v['nearby_all_mean']:.1e}")
