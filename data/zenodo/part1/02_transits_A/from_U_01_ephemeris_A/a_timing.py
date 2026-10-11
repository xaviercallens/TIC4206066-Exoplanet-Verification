"""Per-transit mid-time measurement with a fixed-shape template (own implementation).
Template = batman model from the global fit, integrated over each sector's exposure time, tabulated on a fine
grid of (t - tc) and interpolated.
Model per transit (linear in the nuisance parameters for a given tc):
   fixed depth : y = poly(t) * template(t - tc)
   depth free  : y = poly(t) + A * (template(t - tc) - 1)
poly = polynomial of order `order` in (t - tc_pred) fitted jointly (order 0 on detrended flux; order 1-2 on raw
normalised flux). The data window is |t - tc_pred| < WIN_HALF restricted to the continuous data segment
(no gap > GAP_SPLIT) that contains the prediction."""
import numpy as np
from a_common import *

GRID_HALF = 0.10            # search window for tc: +-0.10 d (+-144 min) around the prediction
GRID_STEP = 0.5 / 1440.     # 0.5 min
WIN_HALF = 0.30             # data window +-0.30 d (7.2 h)
GAP_SPLIT = 0.1             # d; S98 downlink gaps are ~0.21 d
_TPL = {}


def template(key, x, scale=1.0):
    kk = (key, tuple(np.round(x, 10)), scale)
    if kk in _TPL:
        return _TPL[kk]
    T0, P, k, b, aRs, q1, q2 = x
    u1, u2 = q_to_u(q1, q2)
    dt = np.arange(-0.6, 0.6, 2.0 / 86400.)
    if scale == 1.0:
        lc = batman_lc(dt + T0, key, T0, P, k, aRs, b, u1, u2)
    else:
        import batman
        fine = np.arange(-0.8, 0.8, 2.0 / 86400.)
        p = batman.TransitParams(); p.t0 = T0; p.per = P; p.rp = k; p.a = aRs
        p.inc = np.degrees(np.arccos(b / aRs)); p.ecc = 0; p.w = 90; p.limb_dark = 'quadratic'; p.u = [u1, u2]
        un = batman.TransitModel(p, fine / scale + T0).light_curve(p)
        nbox = max(1, int(round(EXPTIME[key] / (2.0 / 86400.))))
        conv = np.convolve(un, np.ones(nbox) / nbox, mode='same')
        lc = np.interp(dt, fine, conv)
    _TPL[kk] = (dt, lc)
    return _TPL[kk]


def tpl_eval(tpl, dt):
    g, lc = tpl
    return np.interp(dt, g, lc, left=1.0, right=1.0)


def segment_window(t, tc_pred, win_half=WIN_HALF, gap=GAP_SPLIT):
    """Boolean selection: points within win_half of tc_pred and in the same continuous segment as tc_pred.
    If tc_pred falls inside a gap, the window keeps both neighbouring segments (the event is then flagged
    by its coverage anyway)."""
    sel = np.abs(t - tc_pred) < win_half
    if sel.sum() < 2:
        return sel
    ts = t[sel]
    idx = np.where(np.diff(ts) > gap)[0]
    if len(idx) == 0:
        return sel
    edges = np.concatenate([[ts[0]], ts[idx + 1]])     # segment starts
    ends = np.concatenate([ts[idx], [ts[-1]]])          # segment ends
    keep = np.zeros(len(ts), bool)
    inside = [(a <= tc_pred <= b) for a, b in zip(edges, ends)]
    if any(inside):
        j = inside.index(True)
        keep = (ts >= edges[j]) & (ts <= ends[j])
    else:   # prediction in a gap: keep nearest segment on each side
        keep[:] = True
    out = np.zeros(len(t), bool); out[np.where(sel)[0][keep]] = True
    return out


def measure(t, y, sig, tpl, tc_pred, depth_free=False, order=0, grid_half=GRID_HALF, step=GRID_STEP,
            exptime=None, t14=T14A):
    sel = segment_window(t, tc_pred)
    tt = t[sel]; yy = y[sel]
    if len(tt) < 8:
        return None
    grid = tc_pred + np.arange(-grid_half, grid_half + step / 2, step)
    s2 = sig ** 2
    xs = (tt - tc_pred) / WIN_HALF
    Pm = np.vstack([xs ** j for j in range(order + 1)])            # (order+1, npts)
    M = tpl_eval(tpl, tt[None, :] - grid[:, None])                 # (ngrid, npts)
    chi = np.empty(len(grid))
    for i in range(len(grid)):
        m = M[i]
        if depth_free:
            X = np.vstack([Pm, m - 1.0]).T
        else:
            X = (Pm * m[None, :]).T
        beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
        r = yy - X @ beta
        chi[i] = (r @ r) / s2
    j = int(np.argmin(chi))
    tc = grid[j]
    if 0 < j < len(grid) - 1:
        c_0, c_1, c_2 = chi[j - 1], chi[j], chi[j + 1]
        den = c_0 - 2 * c_1 + c_2
        if den > 0:
            tc = grid[j] + 0.5 * step * (c_0 - c_2) / den
    cmin = chi.min()
    lo = j
    while lo > 0 and chi[lo] - cmin < 1: lo -= 1
    hi = j
    while hi < len(grid) - 1 and chi[hi] - cmin < 1: hi += 1
    err = 0.5 * (grid[hi] - grid[lo])
    edge = (j <= 2) or (j >= len(grid) - 3)
    # depth-free fit at best tc -> S/N, depth
    m = tpl_eval(tpl, tt - tc)
    X = np.vstack([Pm, m - 1.0]).T
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    C = np.linalg.pinv(X.T @ X) * s2
    snr = beta[-1] / np.sqrt(C[-1, -1]) if C[-1, -1] > 0 else np.nan
    # null (no transit) chi2 with the same baseline polynomial
    Xn = Pm.T
    bn, *_ = np.linalg.lstsq(Xn, yy, rcond=None); rn = yy - Xn @ bn
    chi_null = (rn @ rn) / s2
    # coverage at the measured time and baseline on each side
    cad = exptime if exptime is not None else np.median(np.diff(tt))
    inn = np.abs(tt - tc) < t14 / 2
    cov = inn.sum() * cad / t14
    pre = tt[tt < tc - t14 / 2]; post = tt[tt > tc + t14 / 2]
    base_pre_h = (tc - t14 / 2 - pre.min()) * 24 if len(pre) else 0.0
    base_post_h = (post.max() - tc - t14 / 2) * 24 if len(post) else 0.0
    return dict(tc=float(tc), err_formal=float(err), edge=bool(edge), snr=float(snr), depth_scale=float(beta[-1]),
                dchi2_vs_null=float(chi_null - cmin), n_win=int(len(tt)), chi2min=float(cmin),
                cov_at_tc=float(cov), base_pre_h=float(base_pre_h), base_post_h=float(base_post_h))
