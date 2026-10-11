"""Audit 9: percentile of the target's RUWE, AEN, ipd_gof_harmonic_amplitude and RV statistics among comparison stars."""
import json, sys
from pathlib import Path
import numpy as np
from astropy.table import Table
HERE = Path(__file__).resolve().parent
t = Table.read(HERE / "data" / sys.argv[1]); g = Table.read(HERE / "data/gaia_dr3_source.ecsv")[0]
t = t[t["source_id"] != g["source_id"]]
out = {"n": len(t)}
for c in ["ruwe", "astrometric_excess_noise", "ipd_gof_harmonic_amplitude", "rv_renormalised_gof", "rv_amplitude_robust"]:
    x = np.array(t[c], float); x = x[np.isfinite(x)]
    out[c] = dict(target=float(g[c]), n=int(x.size), median=float(np.median(x)), pct_below_target=float(np.mean(x < g[c]) * 100),
                  p95=float(np.percentile(x, 95)))
x = np.array(t["rv_chisq_pvalue"], float); x = x[np.isfinite(x)]
out["rv_chisq_pvalue"] = dict(target=float(g["rv_chisq_pvalue"]), n=int(x.size), frac_le_target=float(np.mean(x <= g["rv_chisq_pvalue"])),
                              frac_lt_0.01=float(np.mean(x < 0.01)))
print(json.dumps(out, indent=1)); json.dump(out, open(HERE / "out/ruwe_rv_percentiles.json", "w"), indent=1)
