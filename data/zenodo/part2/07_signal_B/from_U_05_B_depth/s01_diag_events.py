"""Diagnostic plot: every B event in the three baseline products (detrended, biweight 0.75 d), with A transit positions,
quality-flagged cadences, and the fixed-shape template at the manuscript depths. Output plots/s01_events.png"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from b_run import *
tpl = Template(dict(k=0.033, b=0.87, rho=2.51, u1=0.32, u2=0.31))
P = {6: 'S6_tessspoc_pdc', 32: 'S32_spoc2m_pdc', 98: 'S98_qlp_sysrm'}
fig, axs = plt.subplots(3, 3, figsize=(15, 11)); axs = axs.ravel(); k = 0
for sec in (6, 32, 98):
    lc_all = get_lc(P[sec], spoc_mask=0 if sec != 98 else MASK_S98_BASE)
    lcd = detrend(lc_all)
    evs = events_in(lcd, EPH_B, T14=tpl.t14_d, half_window=0.3, min_cov=0.25, drop_A_hw=1e-9)
    for ev in evs:
        ax = axs[k]; k += 1
        w = np.abs(lcd['t'] - ev['tc']) < 0.3
        tt = (lcd['t'][w] - ev['tc']) * 24; ff = (lcd['fd'][w] - 1) * 1e6; qq = lcd['q'][w]
        ax.plot(tt, ff, '.', ms=2, color='0.6' if sec != 6 else 'k')
        bad = (qq & 17087) != 0
        ax.plot(tt[bad], ff[bad], 'rx', ms=8, label='SPOC default-mask flagged')
        if sec != 6:
            nb = int(round(lcd['cad'] and (20 / 1440) / lcd['cad']))
            if nb > 1:
                o = np.argsort(tt); n = len(o) // nb * nb
                tb = tt[o][:n].reshape(-1, nb).mean(1); fb = ff[o][:n].reshape(-1, nb).mean(1)
                ax.plot(tb, fb, 'o', color='C0', ms=4, label='20-min bins')
        tm = np.linspace(-7.2, 7.2, 600) / 24
        ax.plot(tm * 24, -1000 * tpl.unit(tm, lcd['cad']), 'C3-', lw=1, label='template 1000 ppm')
        dA = phase_dt(lcd['t'][w], EPH_A)
        inA = np.abs(dA) < EPH_A['T14'] / 2
        if inA.any():
            ax.axvspan(tt[inA].min(), tt[inA].max(), color='orange', alpha=0.2, label='A in transit')
            ax.axvspan(tt[np.abs(dA) < EPH_A['T14']].min(), tt[np.abs(dA) < EPH_A['T14']].max(), color='orange', alpha=0.08, label='A mask +-T14')
        ax.axvspan(-tpl.t14_d * 12, tpl.t14_d * 12, color='C3', alpha=0.06)
        ax.set_title(f"S{sec} epoch {ev['epoch']} tc={ev['tc']:.3f} cov={ev['cov']:.2f}", fontsize=9)
        ax.set_ylim(-3500, 2500); ax.set_xlabel('hours from predicted mid-transit'); ax.set_ylabel('ppm')
        ax.legend(fontsize=6, loc='lower left')
for j in range(k, len(axs)):
    axs[j].axis('off')
plt.tight_layout(); plt.savefig(os.path.join(PLOTS, 's01_events.png'), dpi=110)
print('saved')
