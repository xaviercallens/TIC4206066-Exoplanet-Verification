"""Route 1 (part 2): injection-recovery calibration of per-transit timing errors (own implementation).
For each sector, the global-fit template is injected (multiplicatively) into the real raw normalised light
curve at random times away from real A and B transits (> 0.45 d), and the mid-time is re-measured with exactly
the same routine as the real transits (window centred on the injected time). >= 400 injections per sector."""
import numpy as np, json, csv, sys
from a_common import *
from a_timing import template, measure
from a_global import b_in
from s06_measure_real import VARIANTS

x = np.load('template_params.npy'); T0g, Pg = x[0], x[1]
NINJ = int(sys.argv[1]) if len(sys.argv) > 1 else 400
variants = sys.argv[2].split(',') if len(sys.argv) > 2 else ['base']
rng = np.random.default_rng(20260925)
out = []
summ = {}
for vname in variants:
    v = VARIANTS[vname]
    for sec in [6, 32, 98]:
        z = np.load(f"lc/{v['tags'][sec]}.npz"); key = z['key'].item(); key = int(key) if key.isdigit() else key
        t = z['t']; y = z[v['flux']]; sig = float(z['sig_pp']); cad = float(z['exptime'])
        m = ~b_in(t); t = t[m]; y = y[m]
        tpl = template(key, x, v['scale'])
        from a_timing import tpl_eval
        rec = []; ntry = 0
        while len(rec) < NINJ and ntry < 50 * NINJ:
            ntry += 1
            tc = rng.uniform(t.min(), t.max())
            # away from real A and B transits
            dA = np.abs((tc - T0g + 0.5 * Pg) % Pg - 0.5 * Pg); dB = np.abs((tc - MS_T0B + 0.5 * MS_PB) % MS_PB - 0.5 * MS_PB)
            if dA < 0.45 or dB < 0.45:
                continue
            cov = (np.abs(t - tc) < T14A / 2).sum() * cad / T14A
            if cov < 0.5:
                continue
            yi = y * tpl_eval(tpl, t - tc) if v['flux'] == 'f' else y + (tpl_eval(tpl, t - tc) - 1.0)
            r = measure(t, yi, sig, tpl, tc, depth_free=v['depth_free'], order=v['order'], exptime=cad)
            if r is None:
                continue
            r.update(dict(variant=vname, sector=sec, tc_inj=tc, cov_inj=cov, dt_min=(r['tc'] - tc) * 1440,
                          err_formal_min=r['err_formal'] * 1440))
            rec.append(r)
        dt = np.array([r['dt_min'] for r in rec]); snr = np.array([r['snr'] for r in rec])
        cov = np.array([r['cov_inj'] for r in rec])
        def stats(sel):
            d = dt[sel]
            if len(d) < 10: return None
            return dict(n=int(len(d)), mad_sigma=float(1.4826 * np.median(np.abs(d - np.median(d)))), median=float(np.median(d)),
                        p16_p84_half=float(0.5 * (np.percentile(d, 84) - np.percentile(d, 16))),
                        std_all=float(np.std(d)), frac_gt30min=float(np.mean(np.abs(d) > 30)),
                        frac_gt60min=float(np.mean(np.abs(d) > 60)),
                        std_clip5=float(np.std(d[np.abs(d - np.median(d)) < 5 * 1.4826 * np.median(np.abs(d - np.median(d)))])))
        s = dict(all_cov05=stats(cov >= 0.5), cov07=stats(cov >= 0.7), cov07_snr2=stats((cov >= 0.7) & (snr >= 2)),
                 cov09=stats(cov >= 0.9), cov05_07=stats((cov >= 0.5) & (cov < 0.9)),
                 median_formal_err_min=float(np.median([r['err_formal_min'] for r in rec])),
                 n_tries=ntry)
        summ[f'{vname}_S{sec}'] = s
        out += rec
        print(vname, sec, json.dumps(s))
tag = '_'.join(variants)
json.dump(summ, open(f'out_s07_injection_summary_{tag}.json', 'w'), indent=1)
keys = list(out[0].keys())
with open(f'out_s07_injections_{tag}.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in out]
