"""What is the 'gap' after the S98 3.18 d event at BTJD 3995.555? Raw TESScut TIME/QUALITY and QLP QUALITY (read-only)."""
import numpy as np, json
from astropy.io import fits
R = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/"
h = fits.getdata(R + "tesscut/21x21/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits", 1, memmap=True)
t = np.array(h["TIME"], float); q = np.array(h["QUALITY"], np.int64)
m = (t > 3995.40) & (t < 3995.95)
d = fits.getdata(R + "mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", 1)
tq = np.array(d["TIME"], float); qq = np.array(d["QUALITY"], np.int64); mq = (tq > 3995.40) & (tq < 3995.95)
out = dict(tesscut_n=int(m.sum()), tesscut_time_finite=int(np.isfinite(t[m]).sum()),
           tesscut_quality_values={str(k): int(v) for k, v in zip(*np.unique(q[m], return_counts=True))},
           tesscut_segments=[], qlp_n=int(mq.sum()), qlp_quality_values={str(k): int(v) for k, v in zip(*np.unique(qq[mq], return_counts=True))})
# segments of constant quality in TESScut
tt, qqq = t[m], q[m]
start = 0
for i in range(1, len(tt) + 1):
    if i == len(tt) or qqq[i] != qqq[start]:
        out["tesscut_segments"].append([float(tt[start]), float(tt[i - 1]), int(qqq[start]), i - start]); start = i
# QLP segments
tq2, q2 = tq[mq], qq[mq]; seg = []; s = 0
for i in range(1, len(tq2) + 1):
    if i == len(tq2) or q2[i] != q2[s]:
        seg.append([float(tq2[s]), float(tq2[i - 1]), int(q2[s]), i - s]); s = i
out["qlp_segments"] = seg
dt = np.diff(tq2); out["qlp_time_gaps_gt_0.01d"] = [[float(tq2[i]), float(tq2[i + 1])] for i in np.where(dt > 0.01)[0]]
print(json.dumps(out, indent=1))
json.dump(out, open("/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation/t3_gap3995.json", "w"), indent=1)
