"""T6_star rerun 2: percentile of the target's Gaia DR3 RUWE (1.2142985) among comparable stars.

Re-evaluates the two stored samples (audit 828 stars; U 383 stars) and one fresh Gaia DR3 TAP query that also
returns source_id and binarity flags, so that the target can be excluded and an 'apparently single' subset
(non_single_star = 0, ipd_frac_multi_peak = 0, no other flags) can be formed.
Output: out/t6_02_ruwe_pct.json, data/ruwe_sample_fresh.csv
"""
import json, datetime, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/Users/paulsmith/space <3/tess-jev")
H = Path(__file__).resolve().parent; (H / "out").mkdir(exist_ok=True); (H / "data").mkdir(exist_ok=True)
TID, TR = 3220388198192519424, 1.2142985


def stats(r, label):
    r = np.asarray(r, float); r = r[np.isfinite(r)]
    return dict(sample=label, n=int(r.size), median=float(np.median(r)),
                pct_strictly_below=float(np.mean(r < TR) * 100), pct_le=float(np.mean(r <= TR) * 100),
                frac_gt_1p4=float(np.mean(r > 1.4)),
                pct_among_ruwe_lt_1p4=float(np.mean(r[r < 1.4] < TR) * 100), n_lt_1p4=int((r < 1.4).sum()))


res = {"target_ruwe": TR, "samples": []}
a = pd.read_csv(ROOT / "addenda/tic4206066_PUBLICATION_AUDIT/data/gaia/ruwe_comparison_sample.csv")
res["audit_contains_target_value"] = bool(np.any(np.isclose(a.ruwe, TR, atol=1e-7)))
res["samples"].append(stats(a.ruwe, "audit s14: G10-11, BP-RP 1.60-1.90, plx>15, 5p only (stored, target not removed)"))
u = pd.read_csv(ROOT / "addenda/opus55_final_audit_tic4206066/09_companion/verify/v4_ruwe_sample.csv")
res["U_contains_target_value"] = bool(np.any(np.isclose(u.ruwe, TR, atol=1e-7)))
res["samples"].append(stats(u.ruwe, "U v4: G10.3-10.7, BP-RP 1.6-1.9, plx>15 (stored, target not removed)"))

q = ("SELECT source_id, phot_g_mean_mag, bp_rp, parallax, ruwe, astrometric_params_solved, non_single_star, "
     "ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, duplicated_source, rv_renormalised_gof, rv_chisq_pvalue, rv_nb_transits "
     "FROM gaiadr3.gaia_source WHERE phot_g_mean_mag BETWEEN 10.0 AND 11.0 AND bp_rp BETWEEN 1.60 AND 1.90 "
     "AND parallax > 15 AND ruwe IS NOT NULL")
fresh = None
try:
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = -1
    t = Gaia.launch_job(q).get_results().to_pandas()
    t.columns = [c.lower() for c in t.columns]
    p = H / "data/ruwe_sample_fresh.csv"; t.to_csv(p, index=False)
    res["fresh_query"] = dict(url="https://gea.esac.esa.int/tap-server/tap", adql=q, nrows=len(t),
                              sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                              utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"))
    fresh = t
except Exception as e:  # noqa
    res["fresh_query"] = f"failed: {e!r}"

if fresh is not None:
    f = fresh[fresh.source_id != TID]
    res["fresh_target_found"] = bool((fresh.source_id == TID).any())
    res["samples"].append(stats(f.ruwe, "fresh: G10-11, BP-RP 1.60-1.90, plx>15, all solutions, target removed"))
    res["samples"].append(stats(f[f.astrometric_params_solved == 31].ruwe, "fresh: as above, 5p only"))
    single = f[(f.non_single_star == 0) & (f.ipd_frac_multi_peak == 0) & (~f.duplicated_source.astype(bool))]
    res["samples"].append(stats(single.ruwe, "fresh: apparently single (non_single_star=0, ipd_frac_multi_peak=0, not duplicated)"))
    nb = f[(np.abs(f.phot_g_mean_mag - 10.513) < 0.2) & (np.abs(f.bp_rp - 1.739) < 0.1)]
    res["samples"].append(stats(nb.ruwe, "fresh: narrow box G 10.31-10.71, BP-RP 1.64-1.84"))
    # Gaia RV statistics among the same stars (for context only)
    rv = f[np.isfinite(f.rv_renormalised_gof) & (f.rv_nb_transits >= 10)]
    res["rv_context"] = dict(n=int(len(rv)), pct_gof_below_2p9986=float(np.mean(rv.rv_renormalised_gof < 2.998642) * 100),
                             frac_p_lt_0p01=float(np.mean(rv.rv_chisq_pvalue < 0.01)),
                             frac_katz_variable=float(np.mean((rv.rv_chisq_pvalue < 0.01) & (rv.rv_renormalised_gof > 4))),
                             pct_p_below_target=float(np.mean(rv.rv_chisq_pvalue <= 4.0627176e-05) * 100))
json.dump(res, open(H / "out/t6_02_ruwe_pct.json", "w"), indent=1)
print(json.dumps(res, indent=1))
