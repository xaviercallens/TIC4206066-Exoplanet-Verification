"""O-C diagram at the manuscript ephemeris: own per-transit times (baseline variant), injection errors
(16-84% half width), sector-level shifts (s11), own global fits (LSQ, emcee)."""
import numpy as np, json, csv
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from a_common import *
rows = [r for r in csv.DictReader(open('out_s06_per_transit_all_variants.csv')) if r['variant'] == 'base' and r.get('tc', '') not in ('', 'nan')]
inj = json.load(open('out_s07_injection_summary_base.json'))
sig = {s: inj[f'base_S{s}']['cov07']['p16_p84_half'] for s in (6, 32, 98)}
cc = json.load(open('out_s11_cycle_count.json'))['sector_shift_at_ms_ephemeris']
lsq = json.load(open('out_s05_globalfit_lsq.json'))['raw_lin']['params']
try:
    em = json.load(open('out_s09_emcee_long.json'))
except Exception:
    em = json.load(open('out_s09_emcee_base.json'))
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True, gridspec_kw=dict(width_ratios=[1, 1, 1.6]))
for ax, s in zip(axs, (6, 32, 98)):
    for r in [r for r in rows if int(r['sector']) == s]:
        n = int(r['epoch']); tc = float(r['tc']); oc = (tc - (MS_T0A + n * MS_PA)) * 1440
        good = float(r['coverage']) >= 0.7 and r['overlap_B'] == 'False' and float(r['snr']) >= 2 and r['edge'] == 'False'
        ax.errorbar(n, oc, yerr=sig[s], fmt='o' if good else 's', mfc='C0' if good else 'none', color='C0' if good else 'C3', ms=5, capsize=2)
    ns = np.array([int(r['epoch']) for r in rows if int(r['sector']) == s])
    x = np.linspace(ns.min() - 0.5, ns.max() + 0.5, 50)
    for P, T0, c, lab in [(lsq['P'][0], lsq['T0'][0], 'C2', 'own LSQ'), (em['P']['median'], em['T0']['median'], 'C1', 'own emcee')]:
        ax.plot(x, (T0 + x * P - (MS_T0A + x * MS_PA)) * 1440, color=c, lw=1, label=lab)
    sh = cc[str(s)]
    ax.axhspan(sh['shift_min'] - sh['lo_min'], sh['shift_min'] + sh['hi_min'], color='C4', alpha=0.3, label='sector shift fit')
    ax.axhline(0, color='k', lw=0.5)
    ax.set_title(f'S{s}'); ax.set_xlabel('epoch (manuscript numbering)')
axs[0].set_ylabel('O-C vs manuscript ephemeris [min]'); axs[0].set_ylim(-90, 150); axs[2].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig_s14_oc.png', dpi=90)
