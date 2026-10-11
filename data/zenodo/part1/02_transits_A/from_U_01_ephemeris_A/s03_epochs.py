"""Step 3: enumerate every predicted Signal-A epoch in S6, S32, S98 (manuscript ephemeris) with coverage,
in-transit cadences, proximity to gaps/momentum dumps/scattered light, and overlap with Signal B."""
import numpy as np, json, csv
from a_common import *

rows = []
for key in [6, 32, 98]:
    variants = ['clean', 'keepSL'] if key == 98 else ['clean']
    dd = {v: load(key, s98_mask=v) for v in variants}
    d = dd['clean']
    tall = d['t_all'][np.isfinite(d['t_all'])]
    cad = d['exptime']
    nexp = T14A / cad
    for n in epochs_in(tall.min() - T14A, tall.max() + T14A, MS_PA, MS_T0A):
        tc = MS_T0A + n * MS_PA
        r = dict(sector=key, epoch=int(n), tc_pred=round(tc, 5))
        inwin_all = np.abs(tall - tc) < T14A / 2
        r['n_in_all_cadences'] = int(inwin_all.sum())
        for v in variants:
            dv = dd[v]
            ing = np.abs(dv['t'] - tc) < T14A / 2
            r[f'n_in_{v}'] = int(ing.sum())
            r[f'cov_{v}'] = round(float(ing.sum() / nexp), 3)
            # wider window coverage (±1.5 T14) for baseline
            r[f'n_win3_{v}'] = int((np.abs(dv['t'] - tc) < 1.5 * T14A).sum())
        if key == 98:
            ta = d['t_all']; fin = np.isfinite(ta)
            w = fin & (np.abs(ta - tc) < T14A / 2)
            q = d['q_all'][w]
            r['frac_in_bit2048'] = round(float(((q & STRAY) != 0).mean()) if w.sum() else np.nan, 3)
            r['frac_in_bit2^29'] = round(float(((q & QLP_BIT29) != 0).mean()) if w.sum() else np.nan, 3)
            r['frac_in_bit2^30'] = round(float(((q & QLP_BIT30) != 0).mean()) if w.sum() else np.nan, 3)
        # dumps
        if len(d['dump_t']):
            dt = d['dump_t'] - tc
            j = np.argmin(np.abs(dt)); r['dump_nearest_h'] = round(float(dt[j] * 24), 2)
        # gaps in the good (clean) data: nearest edge of a gap > 3 cadences (or > 0.05 d)
        tg = d['t']; gap = np.diff(tg) > max(3 * cad, 0.02)
        edges = np.concatenate([tg[:-1][gap], tg[1:][gap], [tg.min(), tg.max()]])
        de = edges - tc; j = np.argmin(np.abs(de)); r['gap_edge_nearest_h'] = round(float(de[j] * 24), 2)
        # B overlap
        nb = np.round((tc - MS_T0B) / MS_PB); tb = MS_T0B + nb * MS_PB
        r['dt_to_B_h'] = round(float((tb - tc) * 24), 2)
        r['overlap_B'] = bool(abs(tb - tc) < 0.5 * (T14A + T14B) + 0.25 / 24)
        rows.append(r)

keys = []
for r in rows:
    for k in r:
        if k not in keys: keys.append(k)
with open('out_s03_epochs.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in rows]
for r in rows: print(r)
# summary counts
summ = {}
for key in [6, 32, 98]:
    rr = [r for r in rows if r['sector'] == key]
    summ[key] = dict(n_pred=len(rr),
                     n_any_data=sum(r['n_in_all_cadences'] > 0 for r in rr),
                     n_cov_clean_gt0=sum(r['n_in_clean'] > 0 for r in rr),
                     n_cov_clean_ge05=sum(r['cov_clean'] >= 0.5 for r in rr),
                     n_cov_clean_ge07=sum(r['cov_clean'] >= 0.7 for r in rr),
                     n_overlap_B=sum(r['overlap_B'] for r in rr))
print(json.dumps(summ, indent=1))
json.dump(dict(rows=rows, summary=summ), open('out_s03_epochs.json', 'w'), indent=1)
