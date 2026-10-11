"""Diagnostic plot: every A epoch with data (baseline variant): raw flux + trend (top) and detrended flux (bottom)."""
import numpy as np, csv
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from a_common import *
from a_timing import template, tpl_eval
x = np.load('template_params.npy'); T0g, Pg = x[0], x[1]
rows = [r for r in csv.DictReader(open('out_s06_per_transit_all_variants.csv')) if r['variant'] == 'base' and r['opt'] == 'fixed' and r['tc'] not in ('', 'nan')]
tags = {6: 'S6_w0.75', 32: 'S32_w0.75', 98: 'S98_w0.75_clean_SYS_RM_FLUX'}
Z = {s: np.load(f'lc/{tags[s]}.npz') for s in tags}
raw98 = load(98, s98_mask='keepSL')
n = len(rows); ncol = 6; nrow = int(np.ceil(n / ncol))
fig, axs = plt.subplots(nrow * 2, ncol, figsize=(24, 5.2 * nrow))
for i, r in enumerate(rows):
    s = int(r['sector']); z = Z[s]; tp = float(r['tc_pred_global']); tc = float(r['tc'])
    ax1 = axs[2 * (i // ncol), i % ncol]; ax2 = axs[2 * (i // ncol) + 1, i % ncol]
    w = np.abs(z['t'] - tp) < 0.6
    t = z['t'][w]; b = z['inB'][w]
    ax1.plot((t - tp) * 24, z['f'][w], '.', ms=2, color='0.5'); ax1.plot((t - tp) * 24, z['trend'][w], 'r-')
    if s == 98:
        ww = np.abs(raw98['t'] - tp) < 0.6
        ax1.plot((raw98['t'][ww] - tp) * 24, raw98['f'][ww], '.', ms=1, color='C1', alpha=0.4, zorder=0)
    ax1.set_title(f"S{s} E{r['epoch']} cov {r['coverage']} snr {float(r['snr']):.1f}", fontsize=9)
    y = z['flat'][w]
    ax2.plot((t[~b] - tp) * 24, y[~b], '.', ms=2, color='0.6'); ax2.plot((t[b] - tp) * 24, y[b], 'x', ms=2, color='C2')
    # bins of 20 min
    be = np.arange(-0.6, 0.6001, 20 / 1440.)
    idx = np.digitize(t[~b] - tp, be)
    bc = [(0.5 * (be[k - 1] + be[k]) * 24, np.mean(y[~b][idx == k])) for k in range(1, len(be)) if (idx == k).sum() >= 2]
    if bc: bc = np.array(bc); ax2.plot(bc[:, 0], bc[:, 1], 'ko', ms=3)
    tpl = template(z['key'].item() if not z['key'].item().isdigit() else int(z['key'].item()), x)
    tt = np.linspace(tp - 0.6, tp + 0.6, 2000)
    ax2.plot((tt - tp) * 24, tpl_eval(tpl, tt - tp), 'b-', lw=1, label='pred')
    ax2.plot((tt - tp) * 24, tpl_eval(tpl, tt - tc), 'r-', lw=1, label='meas')
    ax2.set_ylim(1 - 0.0025, 1 + 0.0015); ax2.axvline(0, color='b', lw=0.5)
plt.tight_layout(); plt.savefig('fig_s06_events_baseline.png', dpi=70)
