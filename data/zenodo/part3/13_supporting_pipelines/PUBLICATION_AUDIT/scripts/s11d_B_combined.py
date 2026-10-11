"""Step 11d: indicative overall false-alarm probability of signal B from independent tests on independent data
(Fisher's method). Inputs are read from the audit outputs:
  S98 alone  : s07b BLS null (circular-shift), per S98 product/mask
  S6         : s11c distribution-free rank test (pipeline photometry)
  S32        : s11b period+epoch-randomised null (pipeline photometry)
Caveats: the S98 value depends on the product/mask (itself a choice); Fisher's method assumes each p is uniform under
the null; the S6 and S32 tests are conditional on the S98 ephemeris but use independent data."""
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, save_json  # noqa

J = lambda p: json.load(open(os.path.join(AUD, p)))
s98 = {r["name"]: r["fap"] for r in J("out/s07b_s98_B_significance.json")["results"]}
p6 = [r for r in J("out/s11c_single_box_census.json")["results"] if r["name"] == "pipe S6"][0]["rank_test"]["p_rank"]
p32 = [r for r in J("out/s11b_nullB_periods.json")["results"] if r["name"] == "pipe S32 only"][0]["fap_empirical"]
out = dict(inputs=dict(S98_alone=s98, S6_rank=p6, S32_period_epoch_null=p32), combos={})
for name, p98 in s98.items():
    for label, ps in (("S98+S6", [p98, p6]), ("S98+S6+S32", [p98, p6, p32])):
        chi2 = -2 * np.sum(np.log(ps))
        p = float(stats.chi2.sf(chi2, 2 * len(ps)))
        out["combos"][f"{name} | {label}"] = dict(p=p, sigma_one_sided=float(stats.norm.isf(p)))
        print(f"{name:24s} {label:11s} p = {p:.1e} ({stats.norm.isf(p):.1f} sigma)")
save_json(out, "out/s11d_B_combined.json")
