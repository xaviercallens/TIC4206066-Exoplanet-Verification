"""Route 2b (part 1): global batman fit of all three sectors, P and T0 free (least squares).
Data: own detrended light curves (wotan biweight 0.75 d, A+B masked), points within +-0.35 d of A epochs,
B in-transit points removed. Per-point sigma = per-sector point-to-point scatter.
Parameters: T0 (epoch 0, BTJD), P, k=Rp/R*, b, a/R*, q1, q2 (Kipping 2013), circular orbit."""
import numpy as np, json, sys
from scipy.optimize import least_squares
from a_common import *

def load_set(tags, P=MS_PA, T0=MS_T0A, half=0.35):
    S = []
    for tag in tags:
        z = np.load(f'lc/{tag}.npz'); key = z['key'].item(); key = int(key) if key.isdigit() else key
        t = z['t']; fl = z['flat']; inB = z['inB']
        sel = transit_mask(t, P, T0, half) & ~inB
        S.append(dict(tag=tag, key=key, t=t[sel], y=fl[sel], sig=float(z['sig_pp'])))
    return S

def model_all(p, S):
    T0, P, k, b, aRs, q1, q2 = p
    u1, u2 = q_to_u(q1, q2)
    return [batman_lc(s['t'], s['key'], T0, P, k, aRs, b, u1, u2) for s in S]

def resid(p, S):
    T0, P, k, b, aRs, q1, q2 = p
    if not (0 < q1 < 1 and 0 < q2 < 1 and 0 <= b < 1 + k and aRs > 1.5):
        return np.concatenate([np.full(len(s['t']), 1e3) for s in S])
    m = model_all(p, S)
    return np.concatenate([(s['y'] - mm) / s['sig'] for s, mm in zip(S, m)])

def fit(S, p0):
    lo = [p0[0] - 0.05, p0[1] - 0.002, 0.005, 0.0, 3.0, 0.001, 0.001]
    hi = [p0[0] + 0.05, p0[1] + 0.002, 0.08, 0.95, 40.0, 0.999, 0.999]
    r = least_squares(resid, p0, args=(S,), bounds=(lo, hi), x_scale=[1e-3, 1e-6, 1e-3, 0.1, 1, 0.1, 0.1],
                      diff_step=[1e-7, 1e-9, 1e-4, 1e-3, 1e-4, 1e-3, 1e-3])
    J = r.jac; cov = np.linalg.pinv(J.T @ J)
    chi2 = float((r.fun ** 2).sum()); n = len(r.fun)
    return r.x, np.sqrt(np.diag(cov)), chi2, n, cov

if __name__ == '__main__':
    base = ['S6_w0.75', 'S32_w0.75', 'S98_w0.75_clean_SYS_RM_FLUX']
    S = load_set(base)
    p0 = [MS_T0A, MS_PA, 0.021, 0.35, 12.4, 0.36, 0.3]
    x, e, chi2, n, cov = fit(S, p0)
    # refit from a second start to check convergence
    x2, e2, chi22, _, _ = fit(S, [MS_T0A + 0.003, MS_PA - 0.00001, 0.02, 0.2, 13.0, 0.4, 0.35])
    names = ['T0', 'P', 'k', 'b', 'aRs', 'q1', 'q2']
    res = dict(data=base, n_points=n, chi2=chi2, params={nm: [float(v), float(s)] for nm, v, s in zip(names, x, e)},
               restart_params={nm: float(v) for nm, v in zip(names, x2)}, restart_chi2=chi22,
               corr_T0_P=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])))
    u1, u2 = q_to_u(x[5], x[6]); res['u1u2'] = [float(u1), float(u2)]
    # depth and T14 of the template
    k, b, aRs, P = x[2], x[3], x[4], x[1]
    inc = np.arccos(b / aRs)
    T14 = P / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / (aRs * np.sin(inc)))
    res['T14_h'] = float(T14 * 24)
    tt = np.linspace(-0.1, 0.1, 2001)
    import batman
    lc = batman_lc(tt + x[0], 32, x[0], P, k, aRs, b, u1, u2)
    res['depth_ppm_center'] = float((1 - lc.min()) * 1e6)
    # weighted-centre epoch: T0 at epoch n_c minimizing correlation
    nc = -cov[0, 1] / cov[1, 1]
    res['n_center'] = float(nc)
    res['T0_at_center'] = [float(x[0] + nc * x[1]), float(np.sqrt(cov[0, 0] + 2 * nc * cov[0, 1] + nc ** 2 * cov[1, 1]))]
    print(json.dumps(res, indent=1))
    json.dump(res, open('out_s05_globalfit_lsq.json', 'w'), indent=1)
    np.save('template_params.npy', x)
