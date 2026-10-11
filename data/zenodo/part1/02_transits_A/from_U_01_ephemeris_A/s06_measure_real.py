"""Route 1 (part 1): per-transit mid-times of every predicted A epoch that has data, baseline + sensitivity
variants. Predictions (search-window centres) from the own global LSQ ephemeris (s05, raw_lin)."""
import numpy as np, json, csv
from a_common import *
from a_timing import template, measure
from a_global import b_in

x = np.load('template_params.npy')
T0g, Pg = x[0], x[1]
epochs = json.load(open('out_s03_epochs.json'))['rows']
B = {6: 'S6_w0.75', 32: 'S32_w0.75', 98: 'S98_w0.75_clean_SYS_RM_FLUX'}
def V(**kw):
    d = dict(tags=dict(B), flux='f', order=1, depth_free=False, scale=1.0); d.update(kw); return d
VARIANTS = {
 'base': V(),
 'quad': V(order=2),
 'free': V(depth_free=True),
 'dur0.9': V(scale=0.9),
 'dur1.1': V(scale=1.1),
 'detr_w0.75': V(flux='flat', order=0),
 'detr_w0.5': V(flux='flat', order=0, tags={6: 'S6_w0.5', 32: 'S32_w0.5', 98: 'S98_w0.5_clean_SYS_RM_FLUX'}),
 'detr_w1.0': V(flux='flat', order=0, tags={6: 'S6_w1.0', 32: 'S32_w1.0', 98: 'S98_w1.0_clean_SYS_RM_FLUX'}),
 'S98keepSL': V(tags={6: B[6], 32: B[32], 98: 'S98_w0.75_keepSL_SYS_RM_FLUX'}),
 'S98no30': V(tags={6: B[6], 32: B[32], 98: 'S98_w0.75_no30_SYS_RM_FLUX'}),
 'S98SAP': V(tags={6: B[6], 32: B[32], 98: 'S98_w0.75_clean_SAP_FLUX'}),
 'S98DET': V(tags={6: B[6], 32: B[32], 98: 'S98_w0.75_clean_DET_FLUX'}),
 'S32ffi': V(tags={6: B[6], 32: 'S32ffi_w0.75', 98: B[98]}),
}
if __name__ == '__main__':
    rows = []
    for vname, v in VARIANTS.items():
        for sec in [6, 32, 98]:
            z = np.load(f"lc/{v['tags'][sec]}.npz"); key = z['key'].item(); key = int(key) if key.isdigit() else key
            t = z['t']; y = z[v['flux']]; sig = float(z['sig_pp']); cad = float(z['exptime'])
            m = ~b_in(t)
            tpl = template(key, x, v['scale'])
            for e in [r for r in epochs if r['sector'] == sec]:
                n = e['epoch']; tp = T0g + n * Pg
                cov = float((np.abs(t[m] - tp) < T14A / 2).sum() * cad / T14A)
                row = dict(variant=vname, sector=sec, epoch=n, tc_pred_global=tp, coverage=round(cov, 3),
                           overlap_B=e['overlap_B'], sig_pp_ppm=round(sig * 1e6, 1))
                r = measure(t[m], y[m], sig, tpl, tp, depth_free=v['depth_free'], order=v['order'], exptime=cad) if cov > 0 else None
                if r is not None:
                    row.update(r); row['oc_min'] = (r['tc'] - tp) * 1440; row['err_formal_min'] = r['err_formal'] * 1440
                rows.append(row)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open('out_s06_per_transit_all_variants.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in rows]
    for r in rows:
        if r['variant'] in ('base',) and 'tc' in r:
            print(r['variant'], r['sector'], r['epoch'], 'cov', r['coverage'], 'oc %.1f +- %.1f' % (r['oc_min'], r['err_formal_min']),
                  'snr %.2f' % r['snr'], 'covtc %.2f' % r['cov_at_tc'], 'pre/post h %.1f/%.1f' % (r['base_pre_h'], r['base_post_h']),
                  'edge', r['edge'], 'B', r['overlap_B'])
    # compact table of O-C per variant
    import collections
    tab = collections.OrderedDict()
    for r in rows:
        if 'tc' in r:
            tab.setdefault((r['sector'], r['epoch']), {})[r['variant']] = round(r['oc_min'], 1)
    print('sector epoch ' + ' '.join(VARIANTS.keys()))
    for k, v in tab.items():
        print(k, [v.get(n, None) for n in VARIANTS])
