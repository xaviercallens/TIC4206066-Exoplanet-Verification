"""Diagnostic figure for the archival FAP of signal A (reads s1/s2 outputs; recomputes only the folded curves)."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fap_common import *

s1 = json.load(open('s1_s98_prediction.json')); A = s1['adopted']; NREF = A['NREF']
R = json.load(open('a_fap_main.json')); z = np.load('_diag_main.npz')
P, T = A['P'], A['Tref']
NMID = {6: -792, 32: -569}
col = {6: '#1f5fbf', 32: '#c0392b', 98: '#1e8449'}

fig, ax = plt.subplots(2, 2, figsize=(11, 7.2))
taus = z['taus']; tw = (taus + P / 2) % P - P / 2; o = np.argsort(tw)
for i, s in enumerate((6, 32)):
    a = ax[0, i]
    zz = z['z6'] if s == 6 else z['z32']
    k = NMID[s] - NREF
    sig_s = np.sqrt(A['sT'] ** 2 + k ** 2 * A['sP'] ** 2 + 2 * k * A['corr'] * A['sT'] * A['sP'])
    sl = s1['likelihood']; sig_l = np.sqrt(sl['sT'] ** 2 + k ** 2 * sl['sP'] ** 2 + 2 * k * sl['corr'] * sl['sT'] * sl['sP'])
    a.axvspan(-min(sig_l, P / 2) * 24, min(sig_l, P / 2) * 24, color='#f6e3b4', label=f'likelihood 1$\\sigma$ = {sig_l*24:.0f} h')
    a.plot(tw[o] * 24, zz[o], color=col[s], lw=1)
    tp = T + k * P
    tf = FINAL_T0 + NMID[s] * FINAL_P; td = DISC_T0 + NMID[s] * DISC_P
    for x, c, lab in ((0, 'k', 'S98-only prediction centre'), ((td - tp), '#8e44ad', 'discovery-record centre'),
                      ((tf - tp), '#d35400', 'multi-sector ephemeris')):
        xw = ((x + P / 2) % P - P / 2) * 24
        a.axvline(xw, color=c, ls='--', lw=1, label=lab)
    n1max = R['N1_fixed_period_phase_null'][f'S{s}_null_max']
    a.axhline(n1max, color='#117a65', lw=0.8, ls=':', label='max S/N away from true phase')
    a.set_xlim(-P / 2 * 24, P / 2 * 24)
    a.set_xlabel('epoch offset from S98-only prediction (h), period fixed at %.5f d' % P)
    a.set_ylabel('box S/N (white units)')
    a.set_title(f'Sector {s}: adopted 1$\\sigma$ = {sig_s*24:.0f} h ($\\gg$ T14 = 1.87 h)', fontsize=10)
    a.set_ylim(top=a.get_ylim()[1] * 1.45); a.legend(fontsize=6.5, loc='upper right', ncol=2, framealpha=0.9)

a = ax[1, 0]
N2 = R['N2']['N2']['adopted']
bins = np.linspace(0, 15, 76)
a.hist(z['null_opt_joint'], bins=bins, color='#5dade2', label='null: optimised joint (fake periods, N=%d)' % len(z['null_opt_joint']))
a.hist(z['null_fixed_joint'][np.isfinite(z['null_fixed_joint'])], bins=np.linspace(-6, 15, 106), histtype='step',
       color='#2e4053', label='null: fixed-epoch joint')
a.axvline(N2['opt_joint']['obs'], color='#c0392b', lw=2, label='observed optimised joint S/N %.1f' % N2['opt_joint']['obs'])
a.axvline(N2['fixed_joint']['obs'], color='#8e44ad', lw=2, ls='--', label='observed fixed-epoch joint S/N %.1f' % N2['fixed_joint']['obs'])
a.set_yscale('log'); a.set_xlabel('joint S6+S32 box S/N'); a.set_ylabel('trials'); a.legend(fontsize=7)

a = ax[1, 1]
rec = R['recovered']['adopted']
off = 0
for s in (98, 6, 32):
    t, f, _ = load(s)
    fl = detrend(t, f, b_mask(t, T14B)); ok = np.isfinite(fl) & ~b_mask(t); t, fl = t[ok], fl[ok]
    if s == 98:
        Pp, Tc, lab = P, T, 'S98 (at S98-only ephemeris)'
    else:
        Pp, Tc, lab = rec['P_rec'], rec[f'tc_S{s}_mid'], f'S{s} (at recovered epoch)'
    ph = phase_dist(t, Pp, Tc) * 24
    sel = np.abs(ph) < 6
    b = np.arange(-6, 6.01, 0.5); idx = np.digitize(ph[sel], b)
    yb = np.array([np.mean(fl[sel][idx == j]) if np.sum(idx == j) > 2 else np.nan for j in range(1, len(b))])
    a.plot(0.5 * (b[1:] + b[:-1]), (yb - 1) * 1e6 - off, 'o-', ms=3, color=col[s], label=lab)
    off += 900
a.axvspan(-0.935, 0.935, color='#fdebd0')
a.set_xlabel('hours from mid-transit'); a.set_ylabel('flux - 1 (ppm), offset'); a.legend(fontsize=7)
fig.suptitle('TIC 4206066 signal A: archival recovery vs the S98-only prediction', fontsize=11)
fig.tight_layout()
fig.savefig('a_fap_diagnostic.png', dpi=90)
print('saved')
