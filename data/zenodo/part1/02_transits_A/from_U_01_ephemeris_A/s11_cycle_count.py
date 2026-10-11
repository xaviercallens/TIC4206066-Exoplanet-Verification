"""Step 6: cycle count S6 -> S32 -> S98.
(a) sector-only periods (template fit, shape fixed at the global LSQ fit, per-window linear baselines profiled;
    errors: dchi2 = 1 with beta-inflated sigma, and bootstrap over transits);
(b) propagation of each single-sector / two-sector ephemeris to the other sectors: predicted-time sigma vs P;
(c) global-fit dchi2 of alternative cycle counts (S98<->S32 offsets k = -3..3), P,T0 re-optimised locally;
(d) per-sector mean transit-time offset at the manuscript P, T0 (sector-level shift fit)."""
import numpy as np, json
from scipy.optimize import minimize_scalar, least_squares
from a_common import *
from a_global import build_windows, window_chi2, null_chi2, beta_factor, TAGS_BASE

x = np.load('template_params.npy'); _, _, K, Bb, ARS, Q1, Q2 = x
BETA = {s: max(1.0, beta_factor(TAGS_BASE[s])[0]) for s in (6, 32, 98)}
def W_for(secs, P=MS_PA, T0=MS_T0A, inflate=True):
    W = build_windows(tags={s: TAGS_BASE[s] for s in secs}, P=P, T0=T0, order=1)
    if inflate:
        for w in W: w['sig'] *= BETA[w['sector']]
    return W

def fit_PT(W, P0, T00, nref):
    """fit T0 (at epoch nref) and P with shape fixed; return values, errors (from dchi2=1 curvature), chi2."""
    def res(p):
        return window_chi2(W, p[0] - nref * p[1], p[1], K, Bb, ARS, Q1, Q2, return_resid=True)[1]
    r = least_squares(res, [T00 + nref * P0, P0], x_scale=[1e-3, 1e-5], diff_step=[1e-8, 1e-9])
    J = r.jac; C = np.linalg.inv(J.T @ J)
    return r.x, np.sqrt(np.diag(C)), float((r.fun ** 2).sum()), C

out = dict(beta=BETA)
# (a) single sectors and pairs
ref_epoch = {6: -792, 32: -570, 98: 5}      # near each sector's middle (manuscript epoch numbering)
single = {}
for secs in [(6,), (32,), (98,), (32, 98), (6, 32), (6, 98), (6, 32, 98)]:
    W = W_for(secs)
    nref = int(np.round(np.mean([w['epoch'] for w in W])))
    p, e, chi2, C = fit_PT(W, MS_PA, MS_T0A, nref)
    # bootstrap over transits
    rng = np.random.default_rng(1)
    bs = []
    for i in range(200):
        idx = rng.integers(0, len(W), len(W))
        Wb = [W[j] for j in idx]
        if len(set(w['epoch'] for w in Wb)) < 3: continue
        try:
            pb, _, _, _ = fit_PT(Wb, p[1], p[0] - nref * p[1], nref)
            bs.append(pb)
        except Exception:
            pass
    bs = np.array(bs)
    single['+'.join(f'S{s}' for s in secs)] = dict(n_transits=len(W), epochs=[w['epoch'] for w in W], nref=nref,
        Tref=float(p[0]), eTref=float(e[0]), P=float(p[1]), eP=float(e[1]), chi2=chi2,
        dchi2_vs_null=float(null_chi2(W) - chi2),
        P_boot_std=float(np.std(bs[:, 1])) if len(bs) > 10 else None, Tref_boot_std=float(np.std(bs[:, 0])) if len(bs) > 10 else None,
        P_boot_16_84_half=float(0.5 * np.diff(np.percentile(bs[:, 1], [15.87, 84.13]))[0]) if len(bs) > 10 else None,
        cov=C.tolist())
    s = single['+'.join(f'S{s}' for s in secs)]
    print(secs, 'N', len(W), 'P=%.6f +- %.6f (boot %.6f)' % (s['P'], s['eP'], s['P_boot_std'] or np.nan), 'Tref=%.4f +- %.4f' % (s['Tref'], s['eTref']), 'nref', nref)
out['single_sector_fits'] = single

# (b) propagation
prop = {}
tgt = {6: -792, 32: -570, 98: 5}
for src in ['S98', 'S32', 'S6', 'S32+S98', 'S6+S32', 'S6+S98']:
    s = single[src]; C = np.array(s['cov'])
    for dst in (6, 32, 98):
        if f'S{dst}' in src.split('+'): continue
        dn = tgt[dst] - s['nref']
        tpred = s['Tref'] + dn * s['P']
        sig = np.sqrt(C[0, 0] + 2 * dn * C[0, 1] + dn ** 2 * C[1, 1])
        sig_boot = np.sqrt(s['Tref_boot_std'] ** 2 + (dn * s['P_boot_std']) ** 2) if s['P_boot_std'] else None
        tms = MS_T0A + tgt[dst] * MS_PA
        prop[f'{src}->S{dst}'] = dict(n_cycles=int(dn), t_pred=float(tpred), sigma_pred_d=float(sig),
                                      sigma_pred_boot_d=float(sig_boot) if sig_boot else None,
                                      P_over_sigma=float(MS_PA / sig), P_over_sigma_boot=float(MS_PA / sig_boot) if sig_boot else None,
                                      offset_vs_ms_ephem_d=float(tpred - tms), offset_in_sigma=float((tpred - tms) / sig),
                                      alias_dP_for_pm1_cycle=float(MS_PA / abs(dn)))
        print(src, '->', dst, prop[f'{src}->S{dst}'])
out['propagation'] = prop

# (c) alternative cycle counts in the global fit
alt = {}
W_all_ms = W_for((6, 32, 98))
T98 = MS_T0A; T32 = MS_T0A - 571 * MS_PA
for k in range(-3, 4):
    Pk = (T98 - T32) / (571 + k)
    Wk = W_for((6, 32, 98), P=Pk, T0=T98)
    # local optimisation of P and T0 (at epoch 0 of S98 numbering, relative to this alias)
    p, e, chi2, C = fit_PT(Wk, Pk, T98, 0)
    alt[k] = dict(P_alias_start=float(Pk), P_fit=float(p[1]), T0_fit=float(p[0]), chi2=chi2, chi2_null=float(null_chi2(Wk)),
                  n_windows=len(Wk), n_per_sector={s: sum(w['sector'] == s for w in Wk) for s in (6, 32, 98)})
    alt[k]['dchi2_vs_null'] = alt[k]['chi2_null'] - chi2
best = max(alt.values(), key=lambda a: a['dchi2_vs_null'])['dchi2_vs_null']
for k in alt:
    alt[k]['ddchi2_vs_best'] = best - alt[k]['dchi2_vs_null']
    print('k', k, 'P', round(alt[k]['P_fit'], 6), 'dchi2', round(alt[k]['dchi2_vs_null'], 1), 'ddchi2', round(alt[k]['ddchi2_vs_best'], 1), alt[k]['n_per_sector'])
out['alternative_cycle_counts_S98_S32'] = alt

# (d) per-sector time shift at the manuscript ephemeris (P, T0 fixed; shape fixed), beta-inflated sigma
shifts = {}
for sec in (6, 32, 98):
    W = W_for((sec,))
    grid = np.arange(-60, 60.01, 0.25) / 1440.
    chi = np.array([window_chi2(W, MS_T0A + g, MS_PA, K, Bb, ARS, Q1, Q2) for g in grid])
    j = int(np.argmin(chi)); cmin = chi[j]
    inside = grid[chi - cmin <= 1.0]
    shifts[sec] = dict(shift_min=float(grid[j] * 1440), lo_min=float((grid[j] - inside.min()) * 1440),
                       hi_min=float((inside.max() - grid[j]) * 1440), n_transits=len(W),
                       shift_over_sigma=float(grid[j] / (0.5 * (inside.max() - inside.min()))))
    # also with the manuscript T0 uncertainty and P uncertainty propagated
    n_mid = np.mean([w['epoch'] for w in W])
    s_ms = np.sqrt(MS_eT0A ** 2 + (n_mid * MS_ePA) ** 2)
    shifts[sec]['ms_ephem_sigma_at_sector_min'] = float(s_ms * 1440)
    tot = np.hypot(0.5 * (inside.max() - inside.min()), s_ms)
    shifts[sec]['shift_over_total_sigma'] = float(grid[j] / tot)
    print('sector', sec, shifts[sec])
out['sector_shift_at_ms_ephemeris'] = shifts
json.dump(out, open('out_s11_cycle_count.json', 'w'), indent=1, default=float)
