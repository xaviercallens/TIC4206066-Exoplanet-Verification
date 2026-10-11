"""Verifier: secondary-eclipse scan (windows need >= 0.15 d baseline on both sides) over phase 0.06-0.94 (box width T14 = 1.89 h), flux relative to local quadratic
baselines (same segment), both signals excluded; significance vs the robust scatter of the scan itself."""
import numpy as np, json, sys
sys.argv = ["x", "0.0212", "0.35", "2.51", "0.316", "0.313", "SYS_RM_FLUX", "scan"]
src = open("v04_tests.py").read()
pre = src.split("R = dict(label=LABEL")[0]
exec(pre)
exec("def boxflux" + src.split("def boxflux")[1].split("v05, e05 = boxflux")[0])
TWO_SIDED = True
phs = np.arange(0.06, 0.9401, 0.004)
vals = np.array([boxflux(p, T14)[0] for p in phs])
ok = np.isfinite(vals)
sd = 1.4826 * np.median(np.abs(vals[ok] - np.median(vals[ok])))
j = np.nanargmin(vals)
out = dict(width_h=T14 * 24, n_phases=int(ok.sum()), robust_sd_ppm=float(sd), min_flux_ppm=float(vals[j]), phase_min=float(phs[j]),
           snr_min=float(vals[j] / sd), max_flux_ppm=float(np.nanmax(vals)), flux_at_0p5=float(vals[np.argmin(np.abs(phs - 0.5))]))
rng2 = np.random.default_rng(5)
v05, e05 = boxflux(0.5, T14)
rnd = np.array([boxflux(p, T14)[0] for p in rng2.uniform(0.1, 0.9, 300) if abs(p - 0.5) > 0.04]); rnd = rnd[np.isfinite(rnd)]
out["phase0p5_two_sided"] = dict(flux_ppm=float(v05), sig_formal_ppm=float(e05), sig_emp_std_ppm=float(np.std(rnd)), n_random=len(rnd))
print(out)
np.savetxt(os.path.join(OUT, "secscan_T14.txt"), np.c_[phs, vals], header="phase flux_ppm (negative = dip)")
json.dump(out, open(os.path.join(OUT, "secscan_T14.json"), "w"), indent=1)
