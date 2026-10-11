"""V05: what is the S/N-10 box at BTJD ~2200.0 in the manuscript's S32 null?  Quality bits of SPOC 2-min S32 near the
sector end and around the other strong null boxes (2187.17, 2193.28)."""
import numpy as np, json, vcommon as v
d, t, q = v.load("s32")
f = np.asarray(d["PDCSAP_FLUX"], float)
fin = np.isfinite(t) & np.isfinite(f)
LK = 17087   # manuscript mask (a1_timing_depths.LK_DEFAULT)
out = {}
for tc in (2187.167, 2193.284, 2200.000):
    w = fin & (np.abs(t - tc) < 0.2)
    bits = {}
    for b in range(17):
        n = int(np.sum((q[w] >> b) & 1))
        if n: bits[2 ** b] = n
    passLK = w & ((q & LK) == 0); passV = w & ((q & v.SPOC_BAD) == 0)
    out[str(tc)] = dict(n=int(w.sum()), bits=bits, n_pass_manuscript_mask=int(passLK.sum()), n_pass_verifier_mask=int(passV.sum()),
                        last_time_finite=float(t[fin].max()))
    print(tc, out[str(tc)])
w = fin & (t > 2199.0)
print("t>2199: n", w.sum(), "unique q", np.unique(q[w])[:20])
json.dump(out, open(v.HERE + "/out/v05_s32_edge.json", "w"), indent=1)
