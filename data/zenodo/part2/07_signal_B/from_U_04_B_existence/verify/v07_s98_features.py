"""V07: strongest single 1.6-h / 3-h box dips in S98 (q==0, own biweight, A removed) for DET_FLUX vs SYS_RM_FLUX:
identifies the feature that inflates the SYS_RM null."""
import numpy as np, json, vcommon as v
from v02_link import Grid
out = {}
for col in ("DET_FLUX", "SYS_RM_FLUX"):
    d, t, q = v.load("qlp98")
    f = np.asarray(d[col], float); ok = np.isfinite(t) & np.isfinite(f) & (q == 0); t, f = t[ok], f[ok]
    gi = np.where(np.diff(t) > 0.5)[0]; ed = np.r_[0, gi + 1, len(t)]
    for a, b in zip(ed[:-1], ed[1:]): f[a:b] /= np.median(f[a:b])
    A = v.inA(t); fd = v.biweight(t, f, A, 0.75); k = np.isfinite(fd) & ~A; t, fd = t[k], fd[k]
    s = v.rstd(fd); k = fd < 1 + 5 * s; t, fd = t[k], fd[k]
    g = Grid(t, fd); z = g.dep / g.err; okz = np.isfinite(z); x, zz = g.g[okz], z[okz]
    boxes = []; bl = np.zeros(len(x), bool)
    for i in np.argsort(-zz):
        if bl[i]: continue
        boxes.append((round(float(x[i]), 3), round(float(zz[i]), 2), round(float(g.dep[okz][i] * 1e6))))
        bl |= np.abs(x - x[i]) < v.DUR
        if len(boxes) >= 8: break
    # also most negative (bumps)
    out[col] = dict(rstd_ppm=s * 1e6, beta=g.beta, top_dips=boxes, n=len(t))
    print(col, json.dumps(out[col]))
json.dump(out, open(v.HERE + "/out/v07_s98_features.json", "w"), indent=1)
