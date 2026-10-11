"""S98 depth vs detrending variant and per-epoch (diagnostic for the k discrepancy)."""
import numpy as np, json
src = open('s05_depth_check.py').read().split('for s, cols in')[0]
exec(src)
from astropy.io import fits
res = {}
for v in ['base', 'det98', 'win05', 'b29', 'mask05']:
    for s in (6, 32, 98):
        z = np.load(f'out/lc_{v}_S{s}.npz'); t, f = z['t'], z['f']
        dA, _ = phase_dt(t, P, T0); dB, _ = phase_dt(t, P_B, T0_B)
        m = (np.abs(dA) < 0.3) & (np.abs(dB) > T14_B_H / 24)
        r = fit_local(t[m], f[m], float(z['cad']), s, order=0)
        res[f'{v}_S{s}'] = [float(x) for x in r]
        print(v, s, [round(x, 1) for x in r])
d = fits.open(os.path.join(DATA, FILES[98][0]))[1].data
q = np.asarray(d['QUALITY']).astype(np.int64); t0 = np.asarray(d['TIME'], float); f0 = np.asarray(d['SYS_RM_FLUX'], float)
g = (q == 0) & np.isfinite(f0); tr, fr = t0[g], f0[g] / np.nanmedian(f0[g])
z = np.load('out/lc_base_S98.npz'); t, f = z['t'], z['f']
ep_rows = []
for e in range(-1, 15):
    tc = T0 + e * P
    m = (np.abs(t - tc) < 0.3) & (np.abs(phase_dt(t, P_B, T0_B)[0]) > T14_B_H / 24)
    mr = (np.abs(tr - tc) < 0.3) & (np.abs(phase_dt(tr, P_B, T0_B)[0]) > T14_B_H / 24)
    if (np.abs(t[m] - tc) < T14 / 2).sum() < 3:
        continue
    a = fit_local(t[m], f[m], 200., 98, order=0); r = fit_local(tr[mr], fr[mr], 200., 98, order=1)
    ep_rows.append(dict(epoch=e, wotan=a[0], raw_poly1=r[0], sig=a[1], n=int(m.sum())))
    print(e, 'wotan', round(a[0]), 'rawpoly1', round(r[0]), '+-', round(a[1]), 'n', m.sum())
json.dump(dict(variants=res, s98_epochs=ep_rows), open('out/detrend_variants.json', 'w'), indent=1)
