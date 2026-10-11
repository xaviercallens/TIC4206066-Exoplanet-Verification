"""T6_star rerun 4: large Monte Carlo of rho* = M/R^3 (Mann+2015 R(M_Ks) with 2.89% scatter; Mann+2019 Table 6 n=5
M(M_Ks) with 2.4% total error from Table 7), to decide the printed rounding of the upper error (+0.24 vs +0.25).
Coefficients as used in opus55 06_stellar/s01_relations.py and audit s04_stellar.py (checked there against arXiv sources).
M_Ks ~ N(4.938170, 0.029064) (opus55 06_stellar/out/s01_relations.json MK). Output: out/t6_04_density_mc.json"""
import json
from pathlib import Path
import numpy as np
H = Path(__file__).resolve().parent
out = []
for seed in (1, 2, 3):
    rng = np.random.default_rng(seed); N = 20_000_000
    MK = rng.normal(4.938170, 0.029064, N)
    R = (1.9515 - 0.3520 * MK + 0.01680 * MK**2) * (1 + rng.normal(0, 0.0289, N))
    x = MK - 7.5
    logM = -0.642 - 0.208 * x - 8.43e-4 * x**2 + 7.87e-3 * x**3 + 1.42e-4 * x**4 - 2.13e-4 * x**5
    M = 10**logM * (1 + rng.normal(0, 0.024, N))
    rho = M / R**3
    lo, med, hi = np.percentile(rho, [15.8655, 50, 84.1345])
    lg = np.log10(6.674e-8 * M * 1.98847e33 / (R * 6.957e10) ** 2)
    out.append(dict(seed=seed, rho_median=float(med), plus=float(hi - med), minus=float(med - lo),
                    R_med=float(np.median(R)), R_std=float(R.std()), M_med=float(np.median(M)), M_std=float(M.std()),
                    logg_med=float(np.median(lg)), logg_std=float(lg.std())))
json.dump(out, open(H / "out/t6_04_density_mc.json", "w"), indent=1)
for o in out: print({k: round(v, 5) for k, v in o.items()})
