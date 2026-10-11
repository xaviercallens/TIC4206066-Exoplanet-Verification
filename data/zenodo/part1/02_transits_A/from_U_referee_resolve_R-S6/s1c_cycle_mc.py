"""R-S6 (a): probability of choosing the wrong S32<->S98 cycle count from in-sector timing alone, with per-event timing
errors drawn from the (heavy-tailed) injection error pools of s1 (selected injections), for the real epoch sets."""
import json
import numpy as np
from rs6_common import OUT, os
R = json.load(open(os.path.join(OUT, "s1_timing_R60.json")))
rng = np.random.default_rng(3)
cc = R["cycle_count"]
pools = {s: np.loadtxt(os.path.join(OUT, f"s1_inj_S{s}_R60.txt")) for s in (32, 98)}
pools = {s: p[p[:, 1] == 1, 0] / 1440.0 for s, p in pools.items()}
sig = {s: R["sigma_used_min"][str(s)] / 1440.0 for s in (32, 98)}
ev = {s: np.array([e["epoch"] for e in R["events"] if e["sector"] == s and e["sel"]], float) for s in (32, 98)}
P = cc["P_S98only"] * 0 + 3.1827828
spacing = cc["alias_spacing_d"]
N = 200000
wrong = 0
dev = []
for _ in range(N):
    est = []
    for s in (32, 98):
        E = ev[s] - ev[s].mean()
        y = rng.choice(pools[s], len(E))
        X = np.vstack([np.ones_like(E), E]).T
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        est.append((c[1], sig[s] / np.sqrt(np.sum(E ** 2))))
    w = np.array([1 / e[1] ** 2 for e in est]); pe = np.array([e[0] for e in est])
    d = np.sum(w * pe) / np.sum(w)
    dev.append(d)
dev = np.array(dev)
out = dict(n_mc=N, alias_spacing_d=spacing, frac_wrong_cycle=float(np.mean(np.abs(dev) > spacing / 2)),
           sd_P_insector=float(np.std(dev)), q99_9=float(np.percentile(np.abs(dev), 99.9)), max_abs_dev=float(np.max(np.abs(dev))))
print(out)
json.dump(out, open(os.path.join(OUT, "s1c_cycle_mc.json"), "w"), indent=1)
