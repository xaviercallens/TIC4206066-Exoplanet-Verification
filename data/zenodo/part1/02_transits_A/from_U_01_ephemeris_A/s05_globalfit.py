"""Route 2b (part 1): global batman fit (least squares) of the 23 A windows (coverage >= 0.5, clean masks),
P and T0 free, per-window linear baselines profiled out (raw normalised flux). Variants: quadratic baseline,
detrended flux + constant baseline, S98 with scattered-light intervals kept."""
import numpy as np, json
from scipy.optimize import least_squares
from a_common import *
from a_global import build_windows, window_chi2, null_chi2, TAGS_BASE

NAMES = ['T0', 'P', 'k', 'b', 'aRs', 'q1', 'q2']

def resid(p, W):
    T0, P, k, b, aRs, q1, q2 = p
    if not (0 < q1 < 1 and 0 < q2 < 1 and 0 <= b < 1 + k and aRs > 1.5):
        return np.full(sum(len(w['t']) for w in W), 1e3)
    return window_chi2(W, *p, return_resid=True)[1]

def fit(W, p0):
    lo = [p0[0] - 0.05, p0[1] - 0.002, 0.005, 0.0, 3.0, 0.001, 0.001]
    hi = [p0[0] + 0.05, p0[1] + 0.002, 0.08, 0.95, 40.0, 0.999, 0.999]
    r = least_squares(resid, p0, args=(W,), bounds=(lo, hi), x_scale=[1e-3, 1e-6, 1e-3, 0.1, 1, 0.1, 0.1],
                      diff_step=[1e-7, 1e-9, 1e-4, 1e-3, 1e-4, 1e-3, 1e-3])
    J = r.jac; cov = np.linalg.pinv(J.T @ J)
    return r.x, cov, float((r.fun ** 2).sum()), len(r.fun)

def summarize(W, x, cov, chi2, n):
    e = np.sqrt(np.diag(cov))
    out = dict(n_windows=len(W), n_per_sector={s: sum(w['sector'] == s for w in W) for s in (6, 32, 98)},
               epochs={s: [w['epoch'] for w in W if w['sector'] == s] for s in (6, 32, 98)},
               n_points=n, chi2=chi2, chi2_null=float(null_chi2(W)),
               params={nm: [float(v), float(s)] for nm, v, s in zip(NAMES, x, e)},
               corr_T0_P=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])))
    out['dchi2_transit'] = out['chi2_null'] - chi2
    nc = -cov[0, 1] / cov[1, 1]
    out['n_center'] = float(nc)
    out['T0_at_center'] = [float(x[0] + nc * x[1]), float(np.sqrt(cov[0, 0] + 2 * nc * cov[0, 1] + nc ** 2 * cov[1, 1]))]
    k, b, aRs, P = x[2], x[3], x[4], x[1]
    out['T14_h'] = float(P / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / (aRs * np.sin(np.arccos(b / aRs)))) * 24)
    u1, u2 = q_to_u(x[5], x[6]); out['u1u2'] = [float(u1), float(u2)]
    return out

res = {}
p0 = [MS_T0A, MS_PA, 0.021, 0.35, 12.4, 0.36, 0.3]
for name, kw in [('raw_lin', dict(order=1)), ('raw_quad', dict(order=2)), ('detr_const', dict(order=0, flux='flat')),
                 ('raw_lin_S98keepSL', dict(order=1, tags={6: 'S6_w0.75', 32: 'S32_w0.75', 98: 'S98_w0.75_keepSL_SYS_RM_FLUX'})),
                 ('raw_lin_S32ffi', dict(order=1, tags={6: 'S6_w0.75', 32: 'S32ffi_w0.75', 98: 'S98_w0.75_clean_SYS_RM_FLUX'}))]:
    W = build_windows(**kw)
    x, cov, chi2, n = fit(W, p0)
    x2, cov2, chi22, _ = fit(W, [MS_T0A + 0.003, MS_PA - 0.00001, 0.02, 0.2, 13.0, 0.4, 0.35])
    s = summarize(W, x, cov, chi2, n)
    s['restart'] = dict(chi2=chi22, T0=float(x2[0]), P=float(x2[1]))
    res[name] = s
    if name == 'raw_lin':
        np.save('template_params.npy', x)
    print(name, 'N', s['n_windows'], s['n_per_sector'], 'P', s['params']['P'], 'T0', s['params']['T0'], 'dchi2', round(s['dchi2_transit'], 1),
          'T14h', round(s['T14_h'], 3), 'k', s['params']['k'][0], 'b', s['params']['b'][0], 'restart dchi2', round(chi22 - chi2, 3))
json.dump(res, open('out_s05_globalfit_lsq.json', 'w'), indent=1)
