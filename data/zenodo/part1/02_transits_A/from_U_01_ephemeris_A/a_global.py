"""Joint transit model of all selected A windows with per-window polynomial baselines profiled out
(linear least squares for each trial set of transit parameters). Raw normalised flux (PDCSAP / QLP SYS_RM_FLUX)."""
import numpy as np
from a_common import *
from a_timing import segment_window

TAGS_BASE = {6: 'S6_w0.75', 32: 'S32_w0.75', 98: 'S98_w0.75_clean_SYS_RM_FLUX'}


def b_in(t, margin_min=10.0):
    return transit_mask(t, MS_PB, MS_T0B, 0.5 * T14B + margin_min / 1440.)


def build_windows(tags=TAGS_BASE, P=MS_PA, T0=MS_T0A, epochs=None, min_cov=0.5, order=1, flux='f',
                  win_half=0.30):
    """epochs: dict sector -> list of epochs to use (None: all with coverage >= min_cov)."""
    W = []
    for sec, tag in tags.items():
        z = np.load(f'lc/{tag}.npz'); key = z['key'].item(); key = int(key) if key.isdigit() else key
        t_all = z['t']; y_all = z[flux]; sig = float(z['sig_pp']); cad = float(z['exptime'])
        keepB = ~b_in(t_all)
        t_all = t_all[keepB]; y_all = y_all[keepB]
        ns = epochs_in(t_all.min() - T14A, t_all.max() + T14A, P, T0)
        for n in ns:
            tp = T0 + n * P
            if epochs is not None and n not in epochs.get(sec, []):
                continue
            cov = (np.abs(t_all - tp) < T14A / 2).sum() * cad / T14A
            if epochs is None and cov < min_cov:
                continue
            sel = segment_window(t_all, tp, win_half)
            t = t_all[sel]; y = y_all[sel]
            xs = (t - tp) / win_half
            Pm = np.vstack([xs ** j for j in range(order + 1)])
            W.append(dict(sector=sec, key=key, epoch=int(n), tp=tp, t=t, y=y, sig=sig, Pm=Pm, cov=cov))
    return W


def window_chi2(W, T0, P, k, b, aRs, q1, q2, return_resid=False, tt_shift=None):
    u1, u2 = q_to_u(q1, q2)
    chi = 0.0; res = []
    for w in W:
        m = batman_lc(w['t'], w['key'], T0, P, k, aRs, b, u1, u2)
        X = (w['Pm'] * m[None, :]).T
        beta, *_ = np.linalg.lstsq(X, w['y'], rcond=None)
        r = (w['y'] - X @ beta) / w['sig']
        chi += r @ r
        if return_resid:
            res.append(r)
    return (chi, np.concatenate(res)) if return_resid else chi


def null_chi2(W):
    chi = 0.0
    for w in W:
        X = w['Pm'].T
        beta, *_ = np.linalg.lstsq(X, w['y'], rcond=None)
        r = (w['y'] - X @ beta) / w['sig']; chi += r @ r
    return chi


def beta_factor(tag):
    """Red-noise beta (time-averaging method; Pont et al. 2006, Winn et al. 2008) of the out-of-transit
    detrended flux, median over bins of 30-120 min."""
    z = np.load(f'lc/{tag}.npz'); t = z['t']; y = z['flat']; cad = float(z['exptime'])
    oot = ~(transit_mask(t, MS_PA, MS_T0A, T14A) | transit_mask(t, MS_PB, MS_T0B, T14B))
    t = t[oot]; y = y[oot] - 1; s1 = np.std(y)
    betas = []
    for binmin in [30, 45, 60, 90, 120]:
        nb = binmin / 1440. / cad
        if nb < 2: continue
        edges = np.arange(t.min(), t.max() + binmin / 1440., binmin / 1440.)
        idx = np.digitize(t, edges)
        means = []; counts = []
        for k in np.unique(idx):
            w = idx == k
            if w.sum() >= 0.8 * nb: means.append(y[w].mean()); counts.append(w.sum())
        if len(means) < 10: continue
        M = np.mean(counts); expected = s1 / np.sqrt(M) * np.sqrt(len(means) / (len(means) - 1))
        betas.append(np.std(means) / expected)
    return (float(np.median(betas)) if betas else 1.0), betas
