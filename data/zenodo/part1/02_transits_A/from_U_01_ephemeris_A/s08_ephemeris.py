"""Route 1 (part 3): weighted linear ephemeris from per-transit mid-times with injection-calibrated errors.
Several event-selection rules and error metrics; per-sector mean O-C at the manuscript ephemeris."""
import numpy as np, json, csv
from scipy.optimize import least_squares
from a_common import *

rows = list(csv.DictReader(open('out_s06_per_transit_all_variants.csv')))
inj = json.load(open('out_s07_injection_summary_base.json'))
def f(r, k):
    try: return float(r[k])
    except: return np.nan

SIG = {m: {s: inj[f'base_S{s}']['cov07'][m] for s in (6, 32, 98)} for m in ('p16_p84_half', 'std_all', 'mad_sigma', 'std_clip5')}

def linfit(n, t, e, nref=0.0):
    """weighted LSQ t = T0 + P (n - nref)."""
    A = np.vstack([np.ones_like(n), n - nref]).T
    W = 1 / e ** 2
    C = np.linalg.inv(A.T @ (A * W[:, None]))
    b = C @ (A.T @ (W * t))
    r = t - A @ b
    chi2 = float((r ** 2 * W).sum()); dof = len(t) - 2
    return dict(T0=float(b[0]), P=float(b[1]), eT0=float(np.sqrt(C[0, 0])), eP=float(np.sqrt(C[1, 1])),
                chi2=chi2, dof=dof, rms_oc_min=float(np.sqrt(np.mean(r ** 2)) * 1440), oc_min=(r * 1440).tolist(),
                cov=C.tolist())

def huberfit(n, t, e):
    def res(p): return (t - p[0] - p[1] * n) / e
    r = least_squares(res, [MS_T0A, MS_PA], loss='huber', f_scale=1.5, x_scale=[1e-3, 1e-6])
    J = r.jac; C = np.linalg.inv(J.T @ J)
    return dict(T0=float(r.x[0]), P=float(r.x[1]), eT0=float(np.sqrt(C[0, 0])), eP=float(np.sqrt(C[1, 1])))

def select(variant, rule):
    out = []
    for r in rows:
        if r['variant'] != variant or r.get('tc', '') in ('', 'nan'):
            continue
        sec = int(r['sector']); ep = int(r['epoch']); cov = f(r, 'coverage'); snr = f(r, 'snr')
        ovB = r['overlap_B'] == 'True'; edge = r['edge'] == 'True'
        ok = {
          'all23': cov >= 0.5,
          'T1_cov07_noB_snr2': cov >= 0.7 and not ovB and snr >= 2 and not edge,
          'ms20_like': cov >= 0.7 and not ovB and snr >= 2 and not edge and not (sec == 6 and ep == -795),
          'T2_cov07_noB_nosnrcut': cov >= 0.7 and not ovB and not edge,
          'T3_T1_minus_S32m573': cov >= 0.7 and not ovB and snr >= 2 and not edge and not (sec == 32 and ep == -573),
          'T4_T1_basepost1h': cov >= 0.7 and not ovB and snr >= 2 and not edge and f(r, 'base_post_h') >= 1 and f(r, 'base_pre_h') >= 1,
        }[rule]
        if ok:
            out.append((sec, ep, f(r, 'tc'), f(r, 'err_formal') * 1440))
    return out

res = {'injection_sigma_min': SIG}
for variant in ['base', 'quad', 'free', 'dur0.9', 'dur1.1', 'detr_w0.75', 'detr_w0.5', 'detr_w1.0', 'S98keepSL', 'S98no30', 'S98SAP', 'S98DET', 'S32ffi']:
    for rule in ['all23', 'T1_cov07_noB_snr2', 'ms20_like', 'T2_cov07_noB_nosnrcut', 'T3_T1_minus_S32m573', 'T4_T1_basepost1h']:
        ev = select(variant, rule)
        if len(ev) < 5: continue
        n = np.array([e[1] for e in ev], float); t = np.array([e[2] for e in ev])
        secs = np.array([e[0] for e in ev])
        for metric in ['p16_p84_half', 'std_all', 'formal']:
            if metric == 'formal':
                e = np.array([max(x[3], 0.5) for x in ev]) / 1440
            else:
                e = np.array([SIG[metric][s] for s in secs]) / 1440
            fit = linfit(n, t, e)
            # centre epoch
            C = np.array(fit['cov']); nc = -C[0, 1] / C[1, 1]
            fitc = dict(n_center=float(nc), T0_center=float(fit['T0'] + nc * fit['P']),
                        eT0_center=float(np.sqrt(C[0, 0] + 2 * nc * C[0, 1] + nc ** 2 * C[1, 1])))
            hub = huberfit(n, t, e)
            key = f'{variant}|{rule}|{metric}'
            per = {s: int((secs == s).sum()) for s in (6, 32, 98)}
            # per-sector weighted mean O-C at manuscript ephemeris
            oc_ms = (t - (MS_T0A + MS_PA * n)) * 1440
            sec_oc = {}
            for s in (6, 32, 98):
                w = secs == s
                if w.sum():
                    ww = 1 / (e[w] * 1440) ** 2
                    sec_oc[s] = [float((oc_ms[w] * ww).sum() / ww.sum()), float(1 / np.sqrt(ww.sum())),
                                 float(np.median(oc_ms[w]))]
            res[key] = dict(N=len(ev), N_per_sector=per, epochs=[(int(a), int(b)) for a, b, _, _ in ev],
                            P=fit['P'], eP=fit['eP'], T0=fit['T0'], eT0=fit['eT0'], chi2=fit['chi2'], dof=fit['dof'],
                            rms_oc_min=fit['rms_oc_min'], oc_min=fit['oc_min'], center=fitc, huber=hub,
                            sector_mean_oc_at_ms_ephem_min=sec_oc,
                            dP_vs_ms_sigma=float((fit['P'] - MS_PA) / np.hypot(fit['eP'], MS_ePA)),
                            dT0_vs_ms_min=float((fit['T0'] - MS_T0A) * 1440))
for k, v in res.items():
    if k == 'injection_sigma_min': print(k, v); continue
    if k.split('|')[0] in ('base',) or k.endswith('p16_p84_half'):
        print(f"{k:55s} N={v['N']:2d} {v['N_per_sector']} P={v['P']:.7f}+-{v['eP']:.7f} T0={v['T0']:.4f}+-{v['eT0']:.4f} "
              f"chi2={v['chi2']:.1f}/{v['dof']} rms={v['rms_oc_min']:.1f} Huber P={v['huber']['P']:.7f} T0={v['huber']['T0']:.4f} "
              f"secOC={ {s: [round(a, 1) for a in x[:2]] for s, x in v['sector_mean_oc_at_ms_ephem_min'].items()} }")
json.dump(res, open('out_s08_ephemeris_fits.json', 'w'), indent=1)
