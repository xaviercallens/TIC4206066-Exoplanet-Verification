"""Step 1: the S98-only prediction for signal A, redone on Sector 98 QLP data alone.

Rule fixed before running (written here before the first execution):
  * data: QLP SYS_RM_FLUX, QUALITY == 0; wotan biweight 0.75 d; masked in the trend: B (reconciled ephemeris,
    +-(T14B/2 + 15 min)) and A at the DISCOVERY-record ephemeris (the only A information that existed before the
    archival look), half-width T14 = 1.87 h. B in-transit cadences are dropped from the A fit.
  * model: symmetric trapezoid (no limb darkening, no external shape prior), shape (depth, T14, T12/T14) fitted
    on S98 alone; then the ephemeris (P, T_ref) is profiled with that shape fixed.
  * uncertainty: the larger of (a) the beta-inflated likelihood width (beta = time-averaging red-noise factor)
    and (b) the event bootstrap (resampling the S98 transits); correlation taken from (a).
  * centre: the best (P, T_ref) of the joint S98 likelihood.
The discovery record (P = 3.182945, no stated error) is propagated as a cross-check with the same sigma.
Nothing from S6, S32 or the multi-sector ephemeris enters this step."""
import json, sys
import numpy as np
from scipy.optimize import least_squares
from fap_common import *

rng = np.random.default_rng(98)
t, f, dumps98 = load(98)
mask = (np.abs(phase_dist(t, DISC_P, DISC_T0)) < T14A) | b_mask(t, T14B)
flat = detrend(t, f, mask)
ok = np.isfinite(flat) & ~b_mask(t)
t, y = t[ok], flat[ok]
oot = np.abs(phase_dist(t, DISC_P, DISC_T0)) > T14A
sig = 1.4826 * np.median(np.abs(np.diff(y[oot]))) / np.sqrt(2)


def beta_time_avg(tt, r, bins_min=(30, 45, 60, 90, 120)):
    s1 = np.std(r); out = []
    for bm in bins_min:
        edges = np.arange(tt.min(), tt.max() + bm / 1440, bm / 1440)
        idx = np.digitize(tt, edges)
        sums = np.bincount(idx, r); cnt = np.bincount(idx)
        good = cnt >= max(3, 0.5 * np.median(cnt[cnt > 0]))
        m = sums[good] / cnt[good]; M = good.sum(); N = np.median(cnt[good])
        out.append(np.std(m) / (s1 / np.sqrt(N) * np.sqrt(M / (M - 1))))
    return float(max(1.0, np.median(out))), [float(x) for x in out]


beta, beta_list = beta_time_avg(t[oot], y[oot] - 1)


def trap(dt, depth, T14, rho):
    a = np.abs(dt); T12 = rho * T14; half = 0.5 * T14
    m = np.zeros_like(dt)
    full = a <= half - T12
    ramp = (a > half - T12) & (a < half)
    m[full] = depth
    m[ramp] = depth * (half - a[ramp]) / T12
    return 1 - m


# events: every discovery-ephemeris epoch with data within +-0.3 d
n_all = np.round((t - DISC_T0) / DISC_P).astype(int)
near = np.abs(t - (DISC_T0 + n_all * DISC_P)) < 0.30
epochs = []
for n in np.unique(n_all[near]):
    sel = near & (n_all == n)
    cov = (np.abs(t[sel] - (DISC_T0 + n * DISC_P)) < T14A / 2).sum() * EXPTIME[98] / T14A
    if cov >= 0.5:
        epochs.append(int(n))
NREF = int(np.round(np.mean(epochs)))            # reference epoch near the centre of the S98 events
TREF0 = DISC_T0 + NREF * DISC_P
W = near & np.isin(n_all, epochs)
tw, yw, nw = t[W], y[W], n_all[W]


def resid(p):
    Tref, P, depth, T14, rho, c = p
    tc = Tref + (nw - NREF) * P
    return (yw - c * trap(tw - tc, depth, T14, rho)) / sig


best = None
for P0 in np.arange(3.1800, 3.1860, 0.0004):
    for dT in (-0.01, 0.0, 0.01):
        p0 = [TREF0 + dT, P0, 5e-4, 1.8 / 24, 0.15, 1.0]
        r = least_squares(resid, p0, bounds=([TREF0 - 0.1, 3.176, 0, 0.8 / 24, 0.02, 0.99],
                                             [TREF0 + 0.1, 3.190, 3e-3, 4 / 24, 0.5, 1.01]), x_scale='jac')
        if best is None or r.cost < best.cost:
            best = r
Tref_b, P_b, depth_b, T14_b, rho_b, c_b = best.x
chi2_min = 2 * best.cost
chi2_null = float(np.sum(((yw - np.median(yw)) / sig) ** 2))

# per-event chi2 curves as a function of that event's mid-time (shape and depth fixed at the global best)
grid_dt = np.arange(-0.25, 0.25 + 1e-9, 0.5 / 1440)
curves, ev_tc, ev_etc = {}, {}, {}
for n in epochs:
    s = nw == n; tt, yy = tw[s], yw[s]
    tc0 = Tref_b + (n - NREF) * P_b
    c2 = np.array([np.sum(((yy - c_b * trap(tt - (tc0 + d), depth_b, T14_b, rho_b)) / sig) ** 2) for d in grid_dt])
    curves[n] = (tc0 + grid_dt, c2)
    loc = np.abs(grid_dt) <= 0.1                           # per-event timing searched within +-2.4 h
    j = int(np.argmin(np.where(loc, c2, np.inf))); ev_tc[n] = float(tc0 + grid_dt[j])
    inside = grid_dt[loc & (c2 <= c2[j] + beta ** 2)]      # beta-inflated delta-chi2 = 1 (outer envelope)
    ev_etc[n] = float(max(0.5 * (inside.max() - inside.min()), 1 / 1440))


def surface(ev_list, Pgrid, Tgrid):
    S = np.zeros((len(Pgrid), len(Tgrid)))
    for n in ev_list:
        tg, c2 = curves[n]
        tcs = Tgrid[None, :] + (n - NREF) * Pgrid[:, None]
        S += np.interp(tcs, tg, c2, left=c2[0] + 1e3, right=c2[-1] + 1e3)
    return S


Pgrid = np.arange(3.1780, 3.1880, 4e-6)
Tgrid = Tref_b + np.arange(-0.03, 0.03 + 1e-9, 0.25 / 1440)
S = surface(epochs, Pgrid, Tgrid)
i, j = np.unravel_index(np.argmin(S), S.shape)
L = np.exp(-(S - S.min()) / (2 * beta ** 2)); L /= L.sum()
PP, TT = np.meshgrid(Pgrid, Tgrid, indexing='ij')
mP, mT = float((L * PP).sum()), float((L * TT).sum())
cPP = float((L * (PP - mP) ** 2).sum()); cTT = float((L * (TT - mT) ** 2).sum()); cPT = float((L * (PP - mP) * (TT - mT)).sum())
lik = dict(P_best=float(Pgrid[i]), Tref_best=float(Tgrid[j]), P_mean=mP, Tref_mean=mT, sP=np.sqrt(cPP), sT=np.sqrt(cTT),
           corr=cPT / np.sqrt(cPP * cTT))
# marginal P quantiles
mp = L.sum(1); cdf = np.cumsum(mp)
lik['P_q16_50_84'] = [float(np.interp(q, cdf, Pgrid)) for q in (0.1587, 0.5, 0.8413)]

# event bootstrap
boot = []
for b in range(1000):
    ev = list(rng.choice(epochs, size=len(epochs), replace=True))
    if len(set(ev)) < 3:
        continue
    Sb = surface(ev, Pgrid, Tgrid)
    ib, jb = np.unravel_index(np.argmin(Sb), Sb.shape)
    boot.append((Pgrid[ib], Tgrid[jb]))
boot = np.array(boot)
np.save('s1_boot_P_Tref.npy', boot.astype(np.float64))
bs = dict(n=len(boot), sP=float(np.std(boot[:, 0])), sT=float(np.std(boot[:, 1])),
          P_q16_84=[float(x) for x in np.percentile(boot[:, 0], [15.87, 84.13])],
          corr=float(np.corrcoef(boot.T)[0, 1]))

# per-transit linear ephemeris (WLS, errors rescaled by sqrt(chi2/dof) if > 1)
nn = np.array(epochs, float); tcs = np.array([ev_tc[n] for n in epochs]); ee = np.array([ev_etc[n] for n in epochs])
A = np.vstack([np.ones_like(nn), nn - NREF]).T
Wt = 1 / ee ** 2
C = np.linalg.inv(A.T @ (A * Wt[:, None])); bb = C @ (A.T @ (Wt * tcs))
res = tcs - A @ bb; chi2t = float(np.sum(res ** 2 * Wt)); dof = len(nn) - 2
C *= max(1.0, chi2t / dof)
pt = dict(P=float(bb[1]), Tref=float(bb[0]), sP=float(np.sqrt(C[1, 1])), sT=float(np.sqrt(C[0, 0])),
          corr=float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])), chi2=chi2t, dof=dof,
          oc_min={int(n): round(float(r * 1440), 2) for n, r in zip(nn, res)},
          etc_min={int(n): round(float(e * 1440), 2) for n, e in zip(nn, ee)})

# adopted prediction (rule in the docstring)
P_pred, T_pred = lik['P_best'], lik['Tref_best']
sP = max(lik['sP'], bs['sP']); sT = max(lik['sT'], bs['sT']); rho_PT = lik['corr']
cov = np.array([[sT ** 2, rho_PT * sT * sP], [rho_PT * sT * sP, sP ** 2]])


def propagate(n, P, T, cv):
    k = n - NREF
    return float(T + k * P), float(np.sqrt(cv[0, 0] + k ** 2 * cv[1, 1] + 2 * k * cv[0, 1]))


# epochs of S6 / S32 (relative to the discovery numbering, E0 = BTJD 3998.73)
out = dict(rule=__doc__, n_points=int(len(tw)), sig_pp_ppm=float(sig * 1e6), beta=beta, beta_bins=beta_list,
           epochs_S98=epochs, NREF=NREF,
           global_fit=dict(Tref=float(Tref_b), P=float(P_b), depth_ppm=float(depth_b * 1e6), T14_h=float(T14_b * 24),
                           T12_over_T14=float(rho_b), offset=float(c_b), dchi2_vs_flat=float(chi2_null - chi2_min)),
           likelihood=lik, bootstrap=bs, per_transit=pt,
           adopted=dict(P=P_pred, Tref=T_pred, sP=sP, sT=sT, corr=rho_PT, NREF=NREF,
                        T0_epoch0=float(T_pred - NREF * P_pred)),
           discovery=dict(P=DISC_P, T0=DISC_T0),
           S98_only_elsewhere=dict(pertransit_T1=[3.18164, 3.18206], LEO=[3.18220, 0.00031],
                                   audit_template_s11=[3.1827602, 0.000235], audit_BLS_s12=3.1827078))
prop = {}
for name, P, T, cv in [('adopted', P_pred, T_pred, cov),
                       ('discovery_centre_adopted_sigma', DISC_P, DISC_T0 + NREF * DISC_P, cov)]:
    d = {}
    for sec, (lo, hi) in {'S6': (1468.28, 1490.04), 'S32': (2174.22, 2200.24)}.items():
        n0 = int(np.ceil((lo - T) / P)) + NREF; n1 = int(np.floor((hi - T) / P)) + NREF
        rows = []
        for n in range(n0, n1 + 1):
            tp, sp = propagate(n, P, T, cv)
            tf = FINAL_T0 + n * FINAL_P
            rows.append(dict(epoch=n, t_pred=round(tp, 5), sigma_h=round(sp * 24, 2),
                             sigma_over_T14=round(sp / T14A, 2), final_minus_pred_h=round((tf - tp) * 24, 2),
                             final_minus_pred_sigma=round((tf - tp) / sp, 2)))
        d[sec] = rows
    prop[name] = d
out['propagation'] = prop
json.dump(out, open('s1_s98_prediction.json', 'w'), indent=1)
print(f"S98 points {len(tw)}  sig_pp {sig*1e6:.0f} ppm  beta {beta:.2f} {beta_list}")
print('epochs', epochs, 'NREF', NREF)
print('global fit', out['global_fit'])
print('likelihood', lik)
print('bootstrap', bs)
print('per-transit', {k: v for k, v in pt.items() if k not in ('oc_min',)})
print('ADOPTED', out['adopted'])
for name in prop:
    for sec in ('S6', 'S32'):
        r = prop[name][sec]; mid = r[len(r) // 2]
        print(name, sec, 'epochs', r[0]['epoch'], '..', r[-1]['epoch'], 'mid', mid)
