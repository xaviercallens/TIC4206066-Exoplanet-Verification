"""T1_A_timing: which events move T0/P between the pubprep 20-set and the U 19-set (U verifier times, v2_real_events.csv)."""
import json, os, numpy as np, pandas as pd
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r1_ephemerides.py")).read().split("res = {}")[0])
v2 = pd.read_csv(f"{U}/02_timing_A/verify/v2_real_events.csv").set_index(["sector", "epoch"])
v4 = json.load(open(f"{U}/02_timing_A/verify/v4_results.json"))["injections"]
s19 = [(6,-795),(6,-794),(6,-793),(6,-791),(6,-790),(6,-789),(32,-572),(32,-571),(32,-570),(32,-569),(32,-568),(98,0),(98,1),(98,4),(98,5),(98,9),(98,10),(98,13),(98,14)]
out = {}
for R in (60, 144):
    sd = {s: v4[f"S{s}_R{R}"]["sd"] for s in (6, 32, 98)}
    for name, S in (("19", s19), ("19+E-1", s19+[(98,-1)]), ("19+E-573", s19+[(32,-573)]), ("19-E-795", [k for k in s19 if k!=(6,-795)]),
                    ("pkg20", [k for k in s19 if k!=(6,-795)]+[(98,-1),(32,-573)]), ("21_all", s19+[(98,-1),(32,-573)])):
        d = v2.loc[S].reset_index()
        r = wls(d.epoch, d.tc + d[f"dt{R}"]/1440, d.sector.map(sd)/1440)
        out[f"R{R}:{name}"] = {k: r[k] for k in ("n","P","eP","T0_BJD","eT0","chi2","dof","corr_T0_P")}
        print(f"R{R} {name:9s} n={r['n']} P={r['P']:.7f}+-{r['eP']:.1e} T0={r['T0_BJD']:.5f}+-{r['eT0']:.5f} chi2={r['chi2']:.2f}/{r['dof']} corr={r['corr_T0_P']:.2f}")
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r2_event_leverage.json"), "w"), indent=1)
