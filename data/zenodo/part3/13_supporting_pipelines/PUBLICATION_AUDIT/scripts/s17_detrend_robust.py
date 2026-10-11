"""Step 17: detrending robustness. For both signals and the pipeline set, repeat the per-event depth measurement
(s09 method) and the empirical-error box S/N (s10 method) with: wotan biweight windows 0.5 / 0.75 / 1.0 / 1.5 d,
wotan 'cosine' (robust sum of cosines, window 1.0 d), and NO global detrending (local linear baseline in each event
window only = s09 fit alone). Reports per-sector depths and total S/N for each choice.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import s09_events as ev  # noqa
from audit_common import AUD, save_json  # noqa

EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
raw = [lcio.load_product("TESSSPOC_S6_PDC"), lcio.load_product("SPOC2min_S32_PDC"), lcio.load_product("QLP_S98_SYS", "frozen")]
variants = [("biweight", 0.5), ("biweight", 0.75), ("biweight", 1.0), ("biweight", 1.5), ("cosine", 1.0), ("none", None)]
out = {}
for meth, win in variants:
    tag = f"{meth}_{win}"
    lcs = [lcio.detrend(lc, window=win or 1.0, method=meth, eph=EPH) for lc in raw]
    out[tag] = {}
    for sig in ("A", "B"):
        rows = ev.analyse(sig, lcs, tag)
        summ = {lc["key"]: ev.summarise(rows, lc["key"]) for lc in lcs}
        R = [r for r in rows if "depth_ppm" in r]
        d = np.array([r["depth_ppm"] for r in R])
        e = np.array([r["depth_err_ppm"] for r in R])
        w = 1 / e ** 2
        tot = dict(mean=float(np.sum(w * d) / w.sum()), err=float(1 / np.sqrt(w.sum())), snr=float(np.sum(w * d) / np.sqrt(w.sum())),
                   n=int(len(R)), chi2=float(np.sum(((d - np.sum(w * d) / w.sum()) / e) ** 2)))
        out[tag][sig] = dict(per_sector={k: (None if v is None else dict(n=v["n_events"], depth=v["mean_depth_ppm"], err=v["err_ppm"]))
                                         for k, v in summ.items()}, all=tot)
        ps = "  ".join(f"S{lc['sector']}:{summ[lc['key']]['mean_depth_ppm']:.0f}±{summ[lc['key']]['err_ppm']:.0f}" for lc in lcs if summ[lc['key']])
        print(f"{tag:14s} {sig}: all {tot['mean']:.0f}±{tot['err']:.0f} (S/N {tot['snr']:.1f}, chi2 {tot['chi2']:.1f}/{tot['n']-1})  {ps}", flush=True)
save_json(out, "out/s17_detrend_robust.json")
