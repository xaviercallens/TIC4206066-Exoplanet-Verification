"""Planet radii (A from own chains; B from the k values quoted in the task/manuscript sources) and a summary table of
all joint-fit variants -> out/summary_fits.csv, out/radii.json."""
import glob, json
import numpy as np
from common import *

rng = np.random.default_rng(3)
rows = []
for fn in sorted(glob.glob(os.path.join(OUT, "fit_*.json"))):
    d = json.load(open(fn))
    q = lambda key: d[key]
    rows.append(dict(tag=d["tag"], n_steps=d["n_steps"], n_over_tau=round(d["n_over_tau"], 1), npts=d["npts"],
                     k_med=q("k")[2], k_lo=q("k")[2] - q("k")[1], k_hi=q("k")[3] - q("k")[2],
                     b_med=q("b")[2], b_lo=q("b")[2] - q("b")[1], b_hi=q("b")[3] - q("b")[2],
                     T14_med=q("T14_h")[2], T14_err=(q("T14_h")[3] - q("T14_h")[1]) / 2,
                     Rp_med=q("Rp_Re")[2], Rp_err=(q("Rp_Re")[3] - q("Rp_Re")[1]) / 2,
                     rho_med=q("rho")[2], rho_q97_7=d["rho_upper"]["q97_7"], rho_q95=d["rho_upper"]["q95"],
                     depth_mid_ppm=q("depth_mid_model_ppm")[2], u1=q("u1")[2], u2=q("u2")[2]))
import csv
with open(os.path.join(OUT, "summary_fits.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow({kk: (round(v, 5) if isinstance(v, float) else v) for kk, v in r.items()})
for r in rows:
    print(f"{r['tag']:38s} k={r['k_med']:.5f} -{r['k_lo']:.5f}/+{r['k_hi']:.5f} b={r['b_med']:.2f} "
          f"T14={r['T14_med']:.3f}+-{r['T14_err']:.3f} Rp={r['Rp_med']:.3f}+-{r['Rp_err']:.3f} rho97.7={r['rho_q97_7']:.2f} "
          f"dmid={r['depth_mid_ppm']:.0f} n/tau={r['n_over_tau']}")


def rp(kmed, ksig=0.0, n=400000):
    kk = rng.normal(kmed, ksig, n) if ksig > 0 else np.full(n, kmed)
    r = kk * rng.normal(R_STAR, E_R_STAR, n) * RSUN_OVER_REARTH
    return dict(median=float(np.median(r)), sd=float(np.std(r)), nominal=float(kmed * R_STAR * RSUN_OVER_REARTH))


radii = dict(RSUN_OVER_REARTH=RSUN_OVER_REARTH,
             A_manuscript_k=rp(0.0212, 0.0007),
             B_density_prior_k0p03305=rp(0.03305, 0.00147),   # 0.00147 = half 68% width of the earlier chain (diagnostic)
             B_density_prior_k0p03305_Ronly=rp(0.03305),
             B_free_density_k0p0296=rp(0.0296))
print(json.dumps(radii, indent=1))
json.dump(radii, open(os.path.join(OUT, "radii.json"), "w"), indent=1)
