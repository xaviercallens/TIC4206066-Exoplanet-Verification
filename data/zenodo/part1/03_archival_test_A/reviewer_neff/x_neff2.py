"""Reviewer experiment 2: N_eff of the per-sector / joint OPT null over the full fake-period range U(2.5, 4.5),
measured directly as n_samples / integrated autocorrelation time on a log-spaced period grid; main detrending
(edge_cutoff 0) and the 0.375 d edge-cut variant."""
import json, time, sys
import numpy as np
src = open('s2_archival_fap.py').read().split('# ------------------------------------------------------------------ observed')[0]
EDGE = float(sys.argv[1]); TAG = sys.argv[2]
sys.argv = ['x', '0.75', '10', 'xneff2', str(EDGE)]
ns = {}; exec(src, ns)
D, stats_for, SIG, ellipse_grid = ns['D'], ns['stats_for'], ns['SIG'], ns['ellipse_grid']
sT, sP, rho = SIG['adopted']; EG = ellipse_grid(sT, sP, rho)
def iact(x):
    x = np.asarray(x, float); x = x - x.mean(); n = len(x)
    ac = np.correlate(x, x, 'full')[n-1:] / (np.arange(n, 0, -1) * x.var())
    s = 0.0
    for k in range(1, n // 2):
        if ac[k] <= 0: break
        s += ac[k]
    return 1 + 2 * s
t0 = time.time()
n = 5000
P0s = np.exp(np.linspace(np.log(2.5), np.log(4.5), n))
rng = np.random.default_rng(7)
rows = []
for P0 in P0s:
    anc = {s: D[s]['tmid'] for s in (6, 32)}          # fixed anchors: per-sector opt is anchor-independent
    st = stats_for(P0, anc, sT, sP, rho, True, EG)
    rows.append((st['opt_S6'], st['opt_S32'], st['opt_joint']))
rows = np.array(rows)
out = dict(edge_cutoff_d=EDGE, n=n, P_range=[2.5, 4.5], runtime_s=time.time() - t0)
for j, lab in enumerate(('opt_S6', 'opt_S32', 'opt_joint_fixed_delta')):
    tau = iact(rows[:, j])
    # also block-max check: split into 20 blocks; report max per block spread
    out[lab] = dict(iact_samples=float(tau), N_eff=float(n / tau), bound_k0=float(1 / (n / tau + 1)),
                    mean=float(rows[:, j].mean()), sd=float(rows[:, j].std()), max=float(rows[:, j].max()),
                    p99=float(np.percentile(rows[:, j], 99)))
np.save(f'x_neff2_{TAG}_rows.npy', rows); np.save(f'x_neff2_{TAG}_P0.npy', P0s)
json.dump(out, open(f'x_neff2_{TAG}.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
