"""Verifier diagnostic: depth of S98 epoch -1 (segment 3995.161-3995.606, steep scattered-light ramp, transit at the
segment end) and of all S98 epochs under different detrendings (own implementations):
  (a) wotan biweight 0.75 d, transit points removed BEFORE detrending and trend np.interp'ed under the mask (+-1.5 T14),
      break_tolerance 0.3 d   [mimics the earlier pipeline's mask handling]
  (b) wotan biweight 0.75 d with wotan's own mask (+-1 T14), break_tolerance 0.5 d  [auditor 'base']
  (c) local quadratic, same segment only (v02)  (d) local linear over +-0.3 d ignoring gaps  [~auditor 'raw']
Depth = least-squares scale of the manuscript template (k 0.0212, b 0.35, rho 2.51) + constant, |dt| < 0.25 d."""
import numpy as np, json
from wotan import flatten
from vcommon import *
t, f, cad, ss = load(98)
K, B, RHO = 0.0212, 0.35, 2.51
a = ars(RHO, PA); T14 = t14(PA, K, B, a)
dA = t - (T0A + np.round((t - T0A) / PA) * PA); dB = t - (T0B + np.round((t - T0B) / PB) * PB)
eA = np.round((t - T0A) / PA).astype(int)
mA = np.abs(dA) < 1.5 * T14A; mB = np.abs(dB) < 1.5 * T14B
msk = mA | mB
tr_a = flatten(t[~msk], f[~msk], method="biweight", window_length=0.75, break_tolerance=0.3, return_trend=True)[1]
ok = np.isfinite(tr_a)
tr_a = np.interp(t, t[~msk][ok], tr_a[ok])
mk2 = (np.abs(dA) < T14A) | (np.abs(dB) < T14B)
tr_b = flatten(t, f, method="biweight", window_length=0.75, break_tolerance=0.5, return_trend=True, mask=mk2)[1]
res = {}
for e in (-1, 0, 1, 4, 5, 8, 9, 10, 13, 14):
    tc = T0A + e * PA
    w = (np.abs(t - tc) < 0.25) & (np.abs(dB) > T14B)
    if (w & (np.abs(t - tc) < T14 / 2)).sum() < 3:
        continue
    tmp = 1 - Model(t[w], cad, ss)(tc, PA, K, B, a, 0.316, 0.313)
    dm = tmp.max()
    row = {}
    for lab, fd in (("a_interp_under_mask", f[w] / tr_a[w]), ("b_wotan_mask", f[w] / tr_b[w])):
        X = np.c_[tmp, np.ones_like(tmp)]
        c = np.linalg.lstsq(X, 1 - fd, rcond=None)[0]; row[lab] = c[0] * dm * 1e6
    # (c) local quadratic, same segment, OOT |dt|>0.75 T14A, +-0.45 d
    w2 = (np.abs(t - tc) < 0.45) & segment_of(t, tc) & (np.abs(dB) > T14B)
    oot = w2 & (np.abs(t - tc) > 0.75 * T14A)
    x = t - tc
    c2 = np.polyfit(x[oot], f[oot], 2); fd2 = f / np.polyval(c2, x)
    w3 = w2 & (np.abs(x) < 0.25)
    tmp3 = 1 - Model(t[w3], cad, ss)(tc, PA, K, B, a, 0.316, 0.313)
    row["c_localquad"] = np.linalg.lstsq(np.c_[tmp3, np.ones_like(tmp3)], 1 - fd2[w3], rcond=None)[0][0] * dm * 1e6
    # (d) raw linear + template jointly over +-0.3 d, ignoring gaps
    w4 = (np.abs(x) < 0.3) & (np.abs(dB) > T14B)
    tmp4 = 1 - Model(t[w4], cad, ss)(tc, PA, K, B, a, 0.316, 0.313)
    row["d_rawlinear"] = np.linalg.lstsq(np.c_[tmp4, np.ones_like(tmp4), x[w4]], 1 - f[w4], rcond=None)[0][0] * dm * 1e6
    res[e] = {k: round(float(v)) for k, v in row.items()}
    print(e, res[e])
json.dump(dict(template_mid_depth_ppm=float(dm * 1e6), epochs=res), open(os.path.join(OUT, "s98_epoch_detrend_diag.json"), "w"), indent=1)
