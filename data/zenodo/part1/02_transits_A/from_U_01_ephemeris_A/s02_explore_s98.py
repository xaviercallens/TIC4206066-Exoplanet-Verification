"""Explore QLP S98 quality bits: contiguous intervals of bits 11 (2048), 29 (2^29), 30 (2^30) vs background."""
import numpy as np, json
from astropy.io import fits
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
h = fits.open('../00_data/mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits')
d = h[1].data
t = np.array(d['TIME'], float); q = np.array(d['QUALITY'], np.int64)
def intervals(mask):
    out = []; i = 0; n = len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1] and (t[j + 1] - t[j]) < 0.01: j += 1
            out.append([round(t[i], 4), round(t[j], 4), int(j - i + 1)]); i = j + 1
        else: i += 1
    return out
res = {}
for name, bit in [('b11_2048_straylight', 11), ('b29_undocumented', 29), ('b30_qlp_bad', 30), ('b5_dump', 5), ('b2_coarse', 2)]:
    iv = intervals(((q >> bit) & 1).astype(bool))
    res[name] = iv
    print(name, len(iv), 'intervals; long ones:', [x for x in iv if x[2] > 20])
json.dump(res, open('out_s02_s98_quality_intervals.json', 'w'), indent=1)
fig, ax = plt.subplots(5, 1, figsize=(16, 13), sharex=True)
for a, col in zip(ax[:4], ['SAP_FLUX', 'SYS_RM_FLUX', 'DET_FLUX', 'SAP_BKG']):
    y = np.array(d[col], float)
    good = q == 0
    a.plot(t[~good], y[~good], ',', color='C3'); a.plot(t[good], y[good], ',', color='k'); a.set_ylabel(col)
    if col != 'SAP_BKG':
        med = np.nanmedian(y[good]); a.set_ylim(med - 0.006 * abs(med), med + 0.004 * abs(med))
for bit, yv, c in [(11, 1, 'C1'), (29, 2, 'C2'), (30, 3, 'C3'), (5, 4, 'C4')]:
    m = ((q >> bit) & 1).astype(bool); ax[4].plot(t[m], np.full(m.sum(), yv), '|', color=c)
ax[4].set_yticks([1, 2, 3, 4]); ax[4].set_yticklabels(['2048', '2^29', '2^30', 'dump'])
P, T0 = 3.1827856, 3998.7375
for a in ax:
    for n in range(-4, 16): a.axvline(T0 + n * P, color='C0', lw=0.5, alpha=0.5)
plt.tight_layout(); plt.savefig('fig_s02_s98_quality.png', dpi=90)
