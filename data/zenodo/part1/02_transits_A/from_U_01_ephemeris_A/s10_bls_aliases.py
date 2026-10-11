"""Route 2a: BLS (astropy BoxLeastSquares, likelihood objective) on the combined 3-sector light curve and on
sector subsets, fine period grid around 3.1828 d covering the S32<->S98 (+-1, +-2 cycles) and S6<->S32 aliases.
Light curves detrended WITHOUT any A mask (wotan biweight 0.75 d) so that all trial periods are treated
equally; B in-transit points removed."""
import numpy as np, json
from astropy.timeseries import BoxLeastSquares
from a_common import *
from a_global import b_in

def get(key, s98_mask='clean'):
    d = load(key, s98_mask=s98_mask)
    flat, trend = detrend(d['t'], d['f'], np.zeros(len(d['t']), bool), window=0.75, break_tol=0.1)
    ok = np.isfinite(flat) & ~b_in(d['t'])
    t = d['t'][ok]; y = flat[ok]
    sig = 1.4826 * np.median(np.abs(np.diff(y))) / np.sqrt(2)
    # clip upward outliers only (flares/cosmics)
    good = y < 1 + 5 * sig
    return t[good], y[good], np.full(good.sum(), sig)

D = {s: get(s) for s in (6, 32, 98)}
np.savez('lc/nomask_detrended.npz', **{f't{s}': D[s][0] for s in D}, **{f'y{s}': D[s][1] for s in D}, **{f'e{s}': D[s][2] for s in D})

def run(secs, pmin, pmax, dp, durs=(1.5 / 24, 1.8 / 24, 2.1 / 24)):
    t = np.concatenate([D[s][0] for s in secs]); y = np.concatenate([D[s][1] for s in secs]); e = np.concatenate([D[s][2] for s in secs])
    bls = BoxLeastSquares(t, y, dy=e)
    P = np.arange(pmin, pmax, dp)
    r = bls.power(P, list(durs), objective='likelihood', oversample=20)
    return P, np.array(r.power), r

res = {}
# combined, fine grid
P, pw, r = run((6, 32, 98), 3.1628, 3.2028, 2e-6)
dchi = 2 * pw
j = int(np.argmax(pw))
res['combined'] = dict(P_best=float(P[j]), dchi2_best=float(dchi[j]), T0_best=float(r.transit_time[j]),
                       duration_h=float(r.duration[j] * 24), depth_ppm=float(r.depth[j] * 1e6), depth_snr=float(r.depth_snr[j]))
# local maxima within the grid, top 15
from scipy.signal import argrelmax
idx = argrelmax(dchi, order=200)[0]
idx = idx[np.argsort(dchi[idx])[::-1]][:15]
res['combined']['top_local_maxima'] = [dict(P=float(P[i]), dchi2=float(dchi[i]), ddchi2_vs_best=float(dchi[j] - dchi[i])) for i in idx]
np.savez_compressed('out_s10_bls_combined_periodogram.npz', P=P, dchi2=dchi)
# expected alias positions
Pb = P[j]
T98 = MS_T0A; T32 = MS_T0A - 571 * MS_PA; T6 = MS_T0A - 791 * MS_PA
alias = {}
for k in (-2, -1, 1, 2):
    P98_32 = (T98 - T32) / (571 + k)
    P32_6 = (T32 - T6) / (220 + k)
    P98_6 = (T98 - T6) / (791 + k)
    def near(Pa):
        w = np.abs(P - Pa) < 2e-4
        return float(dchi[w].max()) if w.any() else None
    alias[k] = dict(P_S98_S32=P98_32, dchi2_comb_at=near(P98_32), P_S32_S6=P32_6, dchi2_comb_at_S32S6=near(P32_6),
                    P_S98_S6=P98_6, dchi2_comb_at_S98S6=near(P98_6))
res['alias_positions'] = alias
res['alias_spacing'] = dict(S98_S32=float(MS_PA / 571), S32_S6=float(MS_PA / 220), S98_S6=float(MS_PA / 791))
# subsets
for name, secs, lo, hi, dp in [('S32+S98', (32, 98), 3.1628, 3.2028, 2e-6), ('S6+S98', (6, 98), 3.1628, 3.2028, 2e-6),
                               ('S6+S32', (6, 32), 3.1628, 3.2028, 5e-6),
                               ('S98', (98,), 3.10, 3.27, 2e-5), ('S32', (32,), 3.05, 3.32, 5e-5), ('S6', (6,), 3.0, 3.4, 5e-5)]:
    P2, pw2, r2 = run(secs, lo, hi, dp)
    d2 = 2 * pw2; j2 = int(np.argmax(d2))
    idx2 = argrelmax(d2, order=50)[0]; idx2 = idx2[np.argsort(d2[idx2])[::-1]][:8]
    res[name] = dict(P_best=float(P2[j2]), dchi2_best=float(d2[j2]), T0_best=float(r2.transit_time[j2]),
                     depth_ppm=float(r2.depth[j2] * 1e6),
                     dchi2_at_ms_P=float(d2[np.argmin(np.abs(P2 - MS_PA))]),
                     top_local_maxima=[dict(P=float(P2[i]), dchi2=float(d2[i])) for i in idx2])
    np.savez_compressed(f"out_s10_bls_{name.replace('+', '_')}_periodogram.npz", P=P2, dchi2=d2)
    print(name, res[name]['P_best'], res[name]['dchi2_best'], [(round(m['P'], 6), round(m['dchi2'], 1)) for m in res[name]['top_local_maxima'][:5]])
print('combined', json.dumps(res['combined'], indent=1))
print(json.dumps(res['alias_positions'], indent=1))
json.dump(res, open('out_s10_bls_aliases.json', 'w'), indent=1)

import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(3, 1, figsize=(12, 10))
ax[0].plot(P, dchi, lw=0.5); ax[0].set_title('combined S6+S32+S98 BLS (dchi2 = 2 x likelihood power)')
for k, a in alias.items(): ax[0].axvline(a['P_S98_S32'], color='r', ls=':', lw=0.8)
ax[0].axvline(MS_PA, color='g', lw=0.8)
w = np.abs(P - MS_PA) < 0.001
ax[1].plot(P[w], dchi[w], lw=0.7); ax[1].axvline(MS_PA, color='g', lw=0.8); ax[1].set_title('zoom +-0.001 d')
for name in ['S32+S98', 'S6+S98', 'S6+S32']:
    z = np.load(f"out_s10_bls_{name.replace('+', '_')}_periodogram.npz"); ax[2].plot(z['P'], z['dchi2'], lw=0.5, label=name)
ax[2].legend(); ax[2].axvline(MS_PA, color='g', lw=0.8)
plt.tight_layout(); plt.savefig('fig_s10_bls_aliases.png', dpi=90)
