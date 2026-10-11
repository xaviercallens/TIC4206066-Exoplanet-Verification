"""Step 13d: markdown table of every TRICERATOPS run set (out/fpp/runs/*/summary.json) for PUBLICATION_TRUTH.md section 7.
Values are copied by code, not by hand. Writes out/fpp/fpp_table.md and prints it."""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD  # noqa

LABEL = {
    "repro": "earlier configuration (Teff 4080 K, uncorrected ϖ), audit fold",
    "base": "audit baseline: Mann M, R; Teff 4030 K; ZP-corrected ϖ",
    "loc": "baseline − neighbours excluded by the pixel localisation",
    "molusc5": "loc + MOLUSC companion limits (5 Gyr)",
    "molusc1": "loc + MOLUSC companion limits (1 Gyr)",
    "tic": "loc with TIC 8.2 stellar parameters",
    "prevfold_repro": "earlier configuration AND earlier fold (exact earlier run)",
    "prevfold_base": "audit baseline with the earlier fold",
}
ORDER = ["repro", "base", "loc", "molusc5", "molusc1", "tic", "prevfold_repro", "prevfold_base"]
rows = []
for d in sorted(glob.glob(os.path.join(AUD, "out", "fpp", "runs", "*"))):
    p = os.path.join(d, "summary.json")
    if not os.path.exists(p):
        continue
    x = json.load(open(p))
    L, R = x["log"], x["runs"]
    if "FPP_mean" not in L:
        continue
    sig, rest = os.path.basename(d).split("_", 1)
    var, fold = rest.rsplit("_", 1)
    stp = np.mean([r["P_STP"] for r in R])
    rows.append(dict(sig=sig, var=var, fold=fold, n=len(R), N=L["N"], mean=L["FPP_mean"], sd=L["FPP_std"], lo=L["FPP_min"], hi=L["FPP_max"],
                     nfpp=L.get("NFPP_mean", float("nan")), stp=stp, sat=L.get("n_saturated", 0), exp=L["exptime_d"] * 1440))
key = lambda r: (r["sig"], 0 if r["fold"] == "S32+S98" else 1, ORDER.index(r["var"]) if r["var"] in ORDER else 99, r["fold"])
rows.sort(key=key)
lines = ["| signal | variant | fold | exptime | FPP mean ± sd | range | NFPP | STP share of FPP | runs × draws |", "|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    lines.append(f"| {r['sig']} | {LABEL.get(r['var'], r['var'])} | {r['fold']} | {r['exp']:.3g} min | **{r['mean']:.4f} ± {r['sd']:.4f}** | {r['lo']:.4f}–{r['hi']:.4f} | "
                 f"{r['nfpp']:.1e} | {r['stp'] / r['mean']:.2f} | {r['n']} × {r['N']:.0e} |")
md = "\n".join(lines)
open(os.path.join(AUD, "out", "fpp", "fpp_table.md"), "w").write(md + "\n")
print(md)
