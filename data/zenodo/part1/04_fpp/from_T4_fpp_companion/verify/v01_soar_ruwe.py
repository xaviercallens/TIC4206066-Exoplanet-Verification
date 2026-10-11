"""Verifier: how many SOAR speckle binaries (R4 calibration table, read-only) with one companion at 0.2-2" and dI<=2
have RUWE <= target RUWE 1.2143 (i.e. would NOT be rejected by a RUWE-above-observed rule)?"""
import pandas as pd, json, numpy as np
U = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066"
d = pd.read_csv(U + "/referee/resolve_R4/out/r4_02_binaries_flags.csv")
out = {}
for lo, hi, dm in ((0.2, 2.0, 2.0), (0.2, 2.0, 1.76), (0.2, 0.7, 2.0), (0.7, 2.0, 2.0)):
    s = d[(d.ncomp == 1) & (d.sep_min >= lo) & (d.sep_min <= hi) & (d.dm_of_closest <= dm)]
    out[f"sep{lo}-{hi}_dI<={dm}"] = dict(n=int(len(s)), n_ruwe_le_target=int((s.ruwe <= 1.2143).sum()),
        n_ruwe_le_1p4=int((s.ruwe <= 1.4).sum()), ruwe_min=float(s.ruwe.min()) if len(s) else None,
        seps_with_ruwe_le_target=sorted(s.sep_min[s.ruwe <= 1.2143].round(2).tolist()),
        n_Cstar_flag=int(s.f_Cstar_gt_target_plus_2sig.sum()))
json.dump(out, open("v01_soar_ruwe.json", "w"), indent=1); print(json.dumps(out, indent=1))
