"""T1_A_timing: S98-only linear ephemerides from stored per-transit times (C045), with and without the one-sided E-1."""
import json, os, numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(here, "r1_ephemerides.py")).read().split("res = {}")[0])
Pfin, ePfin = 3.182785, 0.000006
out = {}
pt = pd.read_csv(f"{PUB}/a1_per_transit.csv"); pa = pt[(pt.signal=="A")&(pt.photometry=="pipe")&(pt.sector==98)]
pa = pa[(pa.coverage>=0.7)&(pa.snr_white>=2)&(~pa.edge)]
for name, d in (("pubprep_S98_9", pa), ("pubprep_S98_8_noEm1", pa[pa.epoch!=-1])):
    r = wls(d.epoch, d.tc_pred + d.dt_min/1440, np.full(len(d), 10.066410592769511/1440)); out[name] = r
v2 = pd.read_csv(f"{U}/02_timing_A/verify/v2_real_events.csv"); v4 = json.load(open(f"{U}/02_timing_A/verify/v4_results.json"))["injections"]
s = v2[v2.sector==98]
for R in (60, 144):
    for name, d in ((f"Uver_S98_9_R{R}", s), (f"Uver_S98_8_noEm1_R{R}", s[s.epoch!=-1])):
        out[name] = wls(d.epoch, d.tc + d[f"dt{R}"]/1440, np.full(len(d), v4[f"S98_R{R}"]["sd"]/1440))
for k, r in out.items():
    r["z_vs_final"] = (r["P"] - Pfin) / np.hypot(r["eP"], ePfin)
    print(f"{k:24s} n={r['n']} P={r['P']:.6f}+-{r['eP']:.6f} z={r['z_vs_final']:+.2f} chi2={r['chi2']:.1f}/{r['dof']}")
json.dump(out, open(os.path.join(here, "r3_s98_only.json"), "w"), indent=1)
