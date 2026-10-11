"""Verifier: S32 2-min POS_CORR in-transit displacement for 3.18 d transits vs random epochs (own code)."""
import numpy as np, json
from astropy.io import fits
d = fits.getdata("/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/tess2020324010417-s0032-0000000004206066-0200-s_tp.fits", 1)
BITS = 1|2|4|8|16|32|128|512|2048|4096|16384
t = np.array(d["TIME"], float); q = d["QUALITY"]; p1 = np.array(d["POS_CORR1"], float); p2 = np.array(d["POS_CORR2"], float)
ok = np.isfinite(t) & ((q & BITS) == 0) & np.isfinite(p1) & np.isfinite(p2)
t, p1, p2 = t[ok], p1[ok], p2[ok]
P, T0, T14 = 3.1827856, 3998.7375, 1.87/24; W = 4/24
PB, T0B, T14B = 11.13274, 4000.865, 1.65/24
def inB(tt): ph = ((tt-T0B)/PB+.5)%1-.5; return np.abs(ph*PB) < T14B/2+0.5/24
def inA(tt): ph = ((tt-T0)/P+.5)%1-.5; return np.abs(ph*P) < T14/2+0.5/24
def shift(tc, arr):
    s = (np.abs(t-tc) <= T14/2+W) & ~inB(t)
    tt = t[s]; y = arr[s]; box = (np.abs(tt-tc) < T14/2).astype(float)
    if box.sum() < 20 or (tt < tc-T14/2).sum() < 30 or (tt > tc+T14/2).sum() < 30: return None
    X = np.vstack([np.ones_like(tt), tt-tc, box]).T
    return np.linalg.lstsq(X, y, rcond=None)[0][-1]
n = np.arange(np.ceil((t.min()-T0)/P), np.floor((t.max()-T0)/P)+1); tcs = T0+n*P
real = [(tc, shift(tc, p1), shift(tc, p2)) for tc in tcs]
real = [r for r in real if r[1] is not None]
m1 = np.mean([r[1] for r in real]); m2 = np.mean([r[2] for r in real]); N = len(real)
rng = np.random.default_rng(1); rnd = []
for k in range(2000):
    vals = []
    while len(vals) < N:
        tc = rng.uniform(t.min(), t.max())
        if inA(np.array([tc]))[0] or np.any(np.abs(tcs-tc) < 0.3): continue
        a, b = shift(tc, p1), shift(tc, p2)
        if a is None: continue
        vals.append((a, b))
    rnd.append(np.mean(vals, 0))
rnd = np.array(rnd)
res = dict(n=N, per_event=real, mean_pc1_mpx=m1*1e3, mean_pc2_mpx=m2*1e3, rnd_sd_mpx=(rnd.std(0)*1e3).tolist(),
           rnd_mean_mpx=(rnd.mean(0)*1e3).tolist(), frac_rnd_abs_ge_pc2=float(np.mean(np.abs(rnd[:,1]-rnd[:,1].mean()) >= abs(m2-rnd[:,1].mean()))))
print(json.dumps(res, indent=1, default=float)); json.dump(res, open("v2_poscorr.json","w"), indent=1, default=float)
