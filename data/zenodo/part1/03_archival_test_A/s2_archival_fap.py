"""Step 2-5: false-alarm probability of the archival (S6, S32) recovery of signal A at the S98-only prediction.

Statistic (identical for data and every null trial): box of full width T14 = 1.87 h at every predicted mid-time in the
sector; S/N = mean(1 - flat in boxes) * sqrt(N_in) / sigma_pp (sigma_pp = robust point-to-point white scatter).
The detrending is the same for all trials: wotan biweight 0.75 d with ONLY signal B masked (A is not masked anywhere,
because masking it would favour the true phase). B in-transit cadences (+-(T14B/2 + 15 min)) are dropped.
Joint statistic S6+S32: inverse-variance-weighted box depth over both sectors, S/N = sum(r/s^2)/sqrt(sum(1/s^2)).

Predictions (from s1_s98_prediction.json, S98 data alone): P_pred, T_ref (epoch NREF) and covariance; propagated with
the real cycle lever arms k_s = n_s - NREF (S6 mid epoch -792, S32 mid epoch -569).

Tests
 MAIN (fixed): P = P_pred and epoch = predicted centre, nothing fitted.
 OPT  (optimised): max S/N over (dT, dP) inside the 3-sigma prediction ellipse (per sector: |tau_s| <= 3 sigma_s,
      |dP| <= 3 sigma_P; joint: Mahalanobis <= 3 in (dT, dP)). tau is taken modulo the trial period.
 WPRED (prediction-weighted, no maximisation): log sum_grid w(dT,dP) exp(max(S/N,0)^2/2), w = Gaussian prediction.
Nulls
 N1 fixed period, same (unmasked) data: S/N at every phase offset whose boxes stay >= T14 from the true transits
    (only ~39 independent phases per sector; used for MAIN only).
 N2 fake periods P' ~ U(2.5, 4.5) d on the same data with the true A transits removed (+-T14 around the multi-sector
    ephemeris, used ONLY to remove the signal from the null), random anchors, and EXACTLY the same prediction
    machinery (same sigma_s, lever arms, dP grid, window and weights). N2r: same with P' ~ U(3.10, 3.26) d (brackets
    the 3.125 d S6 momentum-dump cadence)."""
import json, sys, time
import numpy as np
from scipy.stats import norm, chi2 as chi2d, genpareto
import fap_common
from fap_common import *

WINDOW = float(sys.argv[1]) if len(sys.argv) > 1 else 0.75
NPER = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
TAG = sys.argv[3] if len(sys.argv) > 3 else 'main'
fap_common.EDGE = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0   # wotan edge_cutoff variant
rng = np.random.default_rng(2026)
s1 = json.load(open('s1_s98_prediction.json'))
AD = s1['adopted']; NREF = AD['NREF']
P_pred, T_ref = AD['P'], AD['Tref']
SIG = {'adopted': (AD['sT'], AD['sP'], AD['corr']),
       'likelihood': (s1['likelihood']['sT'], s1['likelihood']['sP'], s1['likelihood']['corr'])}
NMID = {6: -792, 32: -569}
HALF = T14A / 2
STEP = 5 / 1440                      # phase grid step (5 min)
NDP = 13                             # dP grid points across +-3 sigma_P

# ------------------------------------------------------------------ data
D = {}
for s in (6, 32):
    t, f, dumps = load(s)
    flat = detrend(t, f, b_mask(t, T14B), window=WINDOW)
    ok = np.isfinite(flat) & ~b_mask(t)
    t, flat = t[ok], flat[ok]
    r = np.median(flat) - flat                                   # positive = dip
    sig = 1.4826 * np.median(np.abs(np.diff(r))) / np.sqrt(2)
    intrue = np.abs(phase_dist(t, FINAL_P, FINAL_T0)) < T14A     # true A transits (+-1.87 h)
    D[s] = dict(t=t, r=r, sig=sig, dumps=dumps, tmid=0.5 * (t.min() + t.max()),
                cr=np.concatenate([[0], np.cumsum(r)]),
                tm=t[~intrue], crm=np.concatenate([[0], np.cumsum(r[~intrue])]))


def boxsums(ts, cr, centers):
    lo = np.searchsorted(ts, centers - HALF); hi = np.searchsorted(ts, centers + HALF)
    return (cr[hi] - cr[lo]), (hi - lo)


def train(s, P, anchor, taus, masked=False):
    """D-sum and N for box trains of period P with phase offsets taus (array) relative to anchor."""
    d = D[s]; ts, cr = (d['tm'], d['crm']) if masked else (d['t'], d['cr'])
    kmin = int(np.floor((ts.min() - anchor - taus.max()) / P)) - 1
    kmax = int(np.ceil((ts.max() - anchor - taus.min()) / P)) + 1
    k = np.arange(kmin, kmax + 1)
    c = anchor + taus[..., None] + k * P
    S, N = boxsums(ts, cr, c)
    return S.sum(-1), N.sum(-1)


def snr(s, S, N):
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(N > 0, S / np.maximum(N, 1) * np.sqrt(N) / D[s]['sig'], 0.0)


def lever(s):
    return NMID[s] - NREF


def table(s, P0, anchor, sP, masked):
    """S/N and (S,N) on grid (dP, tau) with tau in [0, P) step STEP; returns taus, dPs, arrays."""
    dPs = np.linspace(-3 * sP, 3 * sP, NDP)
    taus = np.arange(0, P0, STEP)
    Sg = np.empty((NDP, len(taus))); Ng = np.empty((NDP, len(taus)))
    for i, dp in enumerate(dPs):
        # the anchor is the predicted mid-sector transit; dP also shifts it by lever*dP (handled via tau)
        Sg[i], Ng[i] = train(s, P0 + dp, anchor, taus, masked)
    return taus, dPs, Sg, Ng


def ellipse_grid(sT, sP, rho):
    # grid fine enough that the S6 box train moves by <= 1 phase step between neighbours
    nT = max(15, int(np.ceil(6 * sT / STEP)) + 1)
    nP = int(np.ceil(6 * sP * abs(lever(6)) / STEP)) + 1
    dT = np.linspace(-3 * sT, 3 * sT, nT); dP = np.linspace(-3 * sP, 3 * sP, nP)
    TT, PP = np.meshgrid(dT, dP, indexing='ij')
    C = np.array([[sT ** 2, rho * sT * sP], [rho * sT * sP, sP ** 2]]); Ci = np.linalg.inv(C)
    m2 = Ci[0, 0] * TT ** 2 + 2 * Ci[0, 1] * TT * PP + Ci[1, 1] * PP ** 2
    inside = m2 <= 9.0
    w = np.exp(-0.5 * m2); w[~inside] = 0; w /= w.sum()
    return TT[inside], PP[inside], w[inside]


def lookup(taus, dPs, arr, tau_unwrapped, dP, P0):
    it = np.rint((tau_unwrapped % P0) / STEP).astype(int) % len(taus)
    ip = np.clip(np.rint((dP - dPs[0]) / (dPs[1] - dPs[0])).astype(int), 0, len(dPs) - 1)
    return arr[ip, it]


def stats_for(P0, anchors, sT, sP, rho, masked, EG):
    """All statistics for one (real or fake) prediction: anchors = {s: predicted mid-sector time}."""
    TT, PPd, w = EG
    out = {}
    tabs = {s: table(s, P0, anchors[s], sP, masked) for s in (6, 32)}
    tot_S = 0; tot_V = 0; SS = {}; VV = {}
    for s in (6, 32):
        taus, dPs, Sg, Ng = tabs[s]
        k = lever(s)
        tau_u = TT + k * PPd
        S_ = lookup(taus, dPs, Sg, tau_u, PPd, P0); N_ = lookup(taus, dPs, Ng, tau_u, PPd, P0)
        z = snr(s, S_, N_)
        SS[s] = S_ / D[s]['sig'] ** 2; VV[s] = N_ / D[s]['sig'] ** 2
        # fixed: tau = 0, dP = 0
        i0 = NDP // 2
        out[f'fixed_S{s}'] = float(snr(s, Sg[i0, 0], Ng[i0, 0]))
        # per-sector optimised: |tau_s| <= 3 sigma_s (marginal) and |dP| <= 3 sigma_P
        sig_s = np.sqrt(sT ** 2 + k ** 2 * sP ** 2 + 2 * k * rho * sT * sP)
        if 3 * sig_s >= P0 / 2:
            win = np.ones(len(taus), bool)
        else:
            tw = (taus + P0 / 2) % P0 - P0 / 2
            win = np.abs(tw) <= 3 * sig_s
        Z = snr(s, Sg, Ng)
        out[f'opt_S{s}'] = float(Z[:, win].max())
        # per-sector prediction-weighted
        out[f'wpred_S{s}'] = float(np.log(np.sum(w * np.exp(0.5 * np.maximum(z, 0) ** 2))))
        out[f'_zgrid_S{s}'] = z
    zj = (SS[6] + SS[32]) / np.sqrt(np.maximum(VV[6] + VV[32], 1e-30))
    out['opt_joint'] = float(zj.max())
    out['wpred_joint'] = float(np.log(np.sum(w * np.exp(0.5 * np.maximum(zj, 0) ** 2))))
    i0 = NDP // 2
    S6f = tabs[6][2][i0, 0] / D[6]['sig'] ** 2; S32f = tabs[32][2][i0, 0] / D[32]['sig'] ** 2
    V6f = tabs[6][3][i0, 0] / D[6]['sig'] ** 2; V32f = tabs[32][3][i0, 0] / D[32]['sig'] ** 2
    out['fixed_joint'] = float((S6f + S32f) / np.sqrt(V6f + V32f))
    jb = int(np.argmax(zj)); out['_best'] = (float(TT[jb]), float(PPd[jb])); out['_zj'] = zj
    out['_tabs'] = tabs
    return out


def anchors_real(P0, Tr):
    return {s: Tr + lever(s) * P0 for s in (6, 32)}


# ------------------------------------------------------------------ observed
t0 = time.time()
res = dict(window_d=WINDOW, edge_cutoff_d=fap_common.EDGE, n_fake_periods=NPER, statistic=__doc__.split('Tests')[0].strip(),
           sigma_pp_ppm={f'S{s}': float(D[s]['sig'] * 1e6) for s in (6, 32)})
obs = {}; EGs = {}
for vname, (sT, sP, rho) in SIG.items():
    EGs[vname] = ellipse_grid(sT, sP, rho)
    obs[vname] = stats_for(P_pred, anchors_real(P_pred, T_ref), sT, sP, rho, False, EGs[vname])
# discovery-centred fixed test (its own P and T0, nothing fitted)
disc = {}
for s in (6, 32):
    tc = DISC_T0 + NMID[s] * DISC_P
    S_, N_ = train(s, DISC_P, tc, np.array([0.0]))
    disc[f'S{s}'] = float(snr(s, S_, N_)[0])
S6_, N6_ = train(6, DISC_P, DISC_T0 + NMID[6] * DISC_P, np.array([0.0]))
S32_, N32_ = train(32, DISC_P, DISC_T0 + NMID[32] * DISC_P, np.array([0.0]))
disc['joint'] = float(((S6_ / D[6]['sig'] ** 2 + S32_ / D[32]['sig'] ** 2) /
                       np.sqrt(N6_ / D[6]['sig'] ** 2 + N32_ / D[32]['sig'] ** 2))[0])
# true-ephemeris S/N (for reference only)
true_snr = {}
for s in (6, 32):
    S_, N_ = train(s, FINAL_P, FINAL_T0 + NMID[s] * FINAL_P, np.array([0.0]))
    true_snr[f'S{s}'] = float(snr(s, S_, N_)[0]); true_snr[f'depth_S{s}_ppm'] = float(S_[0] / N_[0] * 1e6)
    true_snr[f'Nin_S{s}'] = int(N_[0])

# ------------------------------------------------------------------ N1: fixed-period phase null (MAIN only)
n1 = {}
for s in (6, 32):
    anc = T_ref + lever(s) * P_pred
    taus = np.arange(0, P_pred, STEP)
    S_, N_ = train(s, P_pred, anc, taus)
    z = snr(s, S_, N_)
    ttrue = FINAL_T0 + np.round((anc - FINAL_T0) / FINAL_P) * FINAL_P
    off = ((anc + taus - ttrue) + P_pred / 2) % P_pred - P_pred / 2
    far = np.abs(off) >= T14A
    for lab, val in (('adopted_centre', obs['adopted'][f'fixed_S{s}']), ('discovery_centre', disc[f'S{s}'])):
        k_ = int(np.sum(z[far] >= val))
        n1[f'S{s}_{lab}'] = dict(obs=val, k=k_, N_grid=int(far.sum()), N_indep=int(far.sum() * STEP / T14A),
                                 fap=(k_ + 1) / (far.sum() + 1))
    n1[f'S{s}_null_sd'] = float(np.std(z[far])); n1[f'S{s}_null_max'] = float(z[far].max())
    n1[f'S{s}_true_phase_snr'] = float(z[~far].max())

# ------------------------------------------------------------------ N2: fake-period nulls
def run_null(Plo, Phi, nper, vnames):
    keys = ['fixed_S6', 'fixed_S32', 'fixed_joint', 'opt_S6', 'opt_S32', 'opt_joint',
            'wpred_S6', 'wpred_S32', 'wpred_joint']
    out = {v: {k: [] for k in keys} for v in vnames}
    for i in range(nper):
        P0 = rng.uniform(Plo, Phi)
        anc = {s: D[s]['tmid'] + rng.uniform(0, P0) for s in (6, 32)}
        for v in vnames:
            sT, sP, rho = SIG[v]
            st = stats_for(P0, anc, sT, sP, rho, True, EGs[v])
            for k in keys:
                out[v][k].append(st[k])
    return {v: {k: np.array(x) for k, x in out[v].items()} for v in vnames}


def fap(null, x):
    null = null[np.isfinite(null)]
    k = int(np.sum(null >= x)); N = len(null)
    return dict(obs=float(x), k=k, N=N, fap=(k + 1) / (N + 1), fap_count=k / N,
                sigma=sigma_from_p((k + 1) / (N + 1)))


def gauss_tail(null, x):
    null = null[np.isfinite(null)]
    mu, sd = float(np.mean(null)), float(np.std(null))
    return dict(mu=mu, sd=sd, z=(x - mu) / sd, p=float(norm.sf((x - mu) / sd)))


def gpd_tail(null, x, q=0.95):
    """peaks-over-threshold extrapolation: generalized Pareto fitted to exceedances over the q-quantile"""
    null = null[np.isfinite(null)]
    u = np.quantile(null, q); exc = null[null > u] - u
    if len(exc) < 10:
        return dict(threshold=float(u), n_exc=int(len(exc)), p=float('nan'), sigma=float('nan'))
    c, loc, scale = genpareto.fit(exc, floc=0)
    p = float((1 - q) * genpareto.sf(x - u, c, 0, scale)) if x > u else float(np.mean(null >= x))
    endpoint = float(u + scale / abs(c)) if c < 0 else float('inf')
    return dict(threshold=float(u), shape=float(c), scale=float(scale), n_exc=int(len(exc)), endpoint=endpoint,
                p=p, sigma=(sigma_from_p(p) if p > 0 else None),
                note=('observed beyond the fitted finite endpoint (shape < 0): extrapolated p = 0, not meaningful'
                      if p == 0 else ''))


def n_eff_periods(Plo, Phi):
    """independent box-train periods per sector in [Plo, Phi]: train shifts by T14 over the baseline"""
    return {s: float((D[s]['t'].max() - D[s]['t'].min()) / T14A * np.log(Phi / Plo)) for s in (6, 32)}


nulls = {'N2': run_null(2.5, 4.5, NPER, list(SIG)), 'N2r': run_null(3.10, 3.26, max(NPER // 5, 500), list(SIG))}
res['runtime_s'] = time.time() - t0
R = {}
PR = {'N2': (2.5, 4.5), 'N2r': (3.10, 3.26)}
for nn, NU in nulls.items():
    R[nn] = {}
    ne = n_eff_periods(*PR[nn])
    for v in SIG:
        R[nn][v] = {}
        sT, sP, rho = SIG[v]
        # independent relative S6-S32 phase positions inside the joint window
        rel = min(P_pred, 6 * abs(lever(32) - lever(6)) * sP) / T14A
        neff = {'S6': ne[6], 'S32': ne[32], 'joint': min(ne[6], ne[32]) * max(1.0, rel)}
        for k, arr in NU[v].items():
            x = obs[v][k]
            e = fap(arr, x)
            e['n_nan'] = int(np.sum(~np.isfinite(arr)))
            if k.startswith('fixed'):
                e['gauss_extrap'] = gauss_tail(arr, x)
            else:
                e['N_eff_indep'] = float(min(e['N'], neff[k.split('_')[1]]))
                e['fap_conservative'] = (e['k'] + 1) / (e['N_eff_indep'] + 1)
                e['sigma_conservative'] = sigma_from_p(e['fap_conservative'])
                e['gpd_extrap'] = gpd_tail(arr, x)
            a = arr[np.isfinite(arr)]
            e['null_p50_p99_max'] = [float(np.percentile(a, 50)), float(np.percentile(a, 99)), float(a.max())]
            R[nn][v][k] = e
        # Fisher combination of the two per-sector p-values (independent data)
        for kind in ('fixed', 'opt', 'wpred'):
            vals = {'counting': ('fap', 'fap'), 'conservative': ('fap_conservative', 'fap_conservative')}
            for lab in ('counting', 'conservative', 'gpd'):
                if kind == 'fixed' and lab != 'counting':
                    continue
                if lab == 'gpd':
                    p6 = R[nn][v][f'{kind}_S6']['gpd_extrap']['p']; p32 = R[nn][v][f'{kind}_S32']['gpd_extrap']['p']
                    if not (p6 > 0 and p32 > 0):
                        continue
                else:
                    p6 = R[nn][v][f'{kind}_S6'][vals[lab][0]]; p32 = R[nn][v][f'{kind}_S32'][vals[lab][0]]
                X = -2 * (np.log(max(p6, 1e-300)) + np.log(max(p32, 1e-300))); pf = float(chi2d.sf(X, 4))
                R[nn][v][f'{kind}_fisher_{lab}'] = dict(p6=p6, p32=p32, X=float(X), p=pf, sigma=sigma_from_p(pf))
res['observed'] = {v: {k: x for k, x in obs[v].items() if not k.startswith('_')} for v in obs}
res['discovery_centre_fixed_snr'] = disc
res['true_ephemeris_snr_reference_only'] = true_snr
res['N1_fixed_period_phase_null'] = n1
res['N2'] = R

# ------------------------------------------------------------------ recovered ephemeris, depth ratio, cycle count
rec = {}
for v in SIG:
    dT, dP = obs[v]['_best']; sT, sP, rho = SIG[v]
    Pr = P_pred + dP; Tr = T_ref + dT
    rec[v] = dict(dT_min=dT * 1440, dP=dP, dT_sigma=dT / sT, dP_sigma=dP / sP, P_rec=Pr,
                  tc_S6_mid=Tr + lever(6) * Pr, tc_S32_mid=Tr + lever(32) * Pr)
    for s in (6, 32):
        tp = T_ref + lever(s) * P_pred; k = lever(s)
        sig_s = np.sqrt(sT ** 2 + k ** 2 * sP ** 2 + 2 * k * rho * sT * sP)
        # nearest-alias recovered time and its offset from the prediction
        dtr = rec[v][f'tc_S{s}_mid'] - tp
        rec[v][f'S{s}_offset_h'] = dtr * 24; rec[v][f'S{s}_offset_sigma'] = dtr / sig_s
        tf = FINAL_T0 + NMID[s] * FINAL_P
        rec[v][f'S{s}_rec_minus_final_min'] = (rec[v][f'tc_S{s}_mid'] - tf) * 1440
    rec[v]['P_rec_minus_final_sigmaPpred'] = (Pr - FINAL_P) / sP
res['recovered'] = rec
res['prediction_vs_final'] = dict(P_pred_minus_final=P_pred - FINAL_P,
                                  in_sigma_adopted=(P_pred - FINAL_P) / SIG['adopted'][1],
                                  in_sigma_likelihood=(P_pred - FINAL_P) / SIG['likelihood'][1])

# depths (box 1.87 h at the recovered archival phase; S98 at its own S98-only ephemeris; same detrending)
t98, f98, _ = load(98)
fl98 = detrend(t98, f98, b_mask(t98, T14B), window=WINDOW)
ok = np.isfinite(fl98) & ~b_mask(t98); t98, r98 = t98[ok], np.median(fl98[ok]) - fl98[ok]
in98 = np.abs(phase_dist(t98, P_pred, T_ref)) < HALF
dep = dict(S98_ppm=float(r98[in98].mean() * 1e6),
           S98_err_ppm=float(1.4826 * np.median(np.abs(np.diff(r98))) / np.sqrt(2) / np.sqrt(in98.sum()) * 1e6))
v = 'adopted'; Pr = rec[v]['P_rec']
for s in (6, 32):
    tc = rec[v][f'tc_S{s}_mid']; inn = np.abs(phase_dist(D[s]['t'], Pr, tc)) < HALF
    dep[f'S{s}_ppm'] = float(D[s]['r'][inn].mean() * 1e6)
    dep[f'S{s}_err_white_ppm'] = float(D[s]['sig'] / np.sqrt(inn.sum()) * 1e6)
    dep[f'ratio_S{s}_to_S98'] = dep[f'S{s}_ppm'] / dep['S98_ppm']
dep['note'] = ('box depths, biweight %.2f d, A unmasked in the trend in all sectors (same for all trials); white '
               'errors understate red noise. Reconciled fitted depths: S6 521+-53, S32 476+-75, S98 496+-48 ppm '
               '(ratios 1.05, 0.96).' % WINDOW)
dep['ratio_reconciled'] = dict(S6=521 / 496, S32=476 / 496)


def depth_variant(t_, f_, P_, T_, mask_A, drop_epochs=(), cad=None):
    m = b_mask(t_, T14B) | (np.abs(phase_dist(t_, P_, T_)) < T14A if mask_A else False)
    fl = detrend(t_, f_, m, window=WINDOW)
    ok_ = np.isfinite(fl) & ~b_mask(t_)
    tt, rr = t_[ok_], np.median(fl[ok_]) - fl[ok_]

    def boxdepth(T__):
        ph = phase_dist(tt, P_, T__); n = np.round((tt - T__) / P_).astype(int)
        inn = np.abs(ph) < HALF
        keep = np.zeros_like(inn)
        for e in np.unique(n[inn]):                 # events with >= 50 % box coverage only
            sel = inn & (n == e)
            if sel.sum() * cad / T14A >= 0.5 and e not in drop_epochs:
                keep |= sel
        return rr[keep].mean() if keep.any() else np.nan
    d0 = boxdepth(T_)
    # empirical (red-noise) error: SD of the same box-train depth at phases >= 2 T14 from the transit
    offs = np.arange(2 * T14A, P_ - 2 * T14A, T14A)
    dd = np.array([boxdepth(T_ + o) for o in offs])
    return float(d0 * 1e6), float(np.nanstd(dd) * 1e6), len(offs)


dv = {}
for mA in (False, True):
    lab = 'A_masked' if mA else 'A_unmasked'
    t98_, f98_, _ = load(98)
    d98, e98, n98 = depth_variant(t98_, f98_, P_pred, T_ref, mA, cad=EXPTIME[98])
    # S98 epoch -1 (NREF-relative -7) sits at a data-gap edge (reconciled facts exclude it from the depth)
    d98b, e98b, _ = depth_variant(t98_, f98_, P_pred, T_ref, mA, drop_epochs=(-1 - NREF,), cad=EXPTIME[98])
    dv[lab] = dict(S98=[d98, e98], S98_without_Em1=[d98b, e98b])
    for s in (6, 32):
        ts_, fs_, _ = load(s)
        d_, e_, n_ = depth_variant(ts_, fs_, Pr, rec[v][f'tc_S{s}_mid'], mA, cad=EXPTIME[s])
        dv[lab][f'ratio_S{s}_to_S98_without_Em1'] = [d_ / d98b, d_ / d98b * np.hypot(e_ / d_, e98b / d98b)]
        dv[lab][f'S{s}'] = [d_, e_]
        rat = d_ / d98; erat = rat * np.hypot(e_ / d_, e98 / d98)
        dv[lab][f'ratio_S{s}_to_S98'] = [rat, erat]
    dv[lab]['n_null_phases'] = n98
dep['box_depth_variants_ppm_[value,empirical_err]'] = dv
res['depths'] = dep

# cycle count from S98 + one archival sector
boot = np.load('s1_boot_P_Tref.npy')
cc = {}
for s in (6, 32):
    tc = rec['adopted'][f'tc_S{s}_mid']
    n_best = int(np.round((T_ref - tc) / P_pred))
    rows = []
    for dn in (-2, -1, 0, 1, 2):
        N = n_best + dn; Pn = (T_ref - tc) / N
        rows.append(dict(cycles=N, P=Pn, z_adopted=(Pn - P_pred) / SIG['adopted'][1],
                         z_likelihood=(Pn - P_pred) / SIG['likelihood'][1]))
    Ps = np.array([r['P'] for r in rows])
    nearest = np.argmin(np.abs(boot[:, 0][:, None] - Ps[None, :]), axis=1)
    for i, r in enumerate(rows):
        r['boot_frac_nearest'] = float(np.mean(nearest == i))
    cc[f'S98+S{s}'] = dict(alias_spacing_dP=float(P_pred / abs(n_best)), rows=rows,
                           unique_selected_cycles=int(n_best), final_ephemeris_cycles=int(-NMID[s] + NREF))
res['cycle_count'] = cc

# ------------------------------------------------------------------ dump check (S6)
d6 = D[6]; dumps = np.sort(d6['dumps'])
cent = np.arange(d6['t'].min() + HALF, d6['t'].max() - HALF, STEP)
S_, N_ = boxsums(d6['tm'], d6['crm'], cent)
okc = N_ >= 3
z1 = np.where(okc, S_ / np.maximum(N_, 1) * np.sqrt(N_) / d6['sig'], np.nan)
ddump = np.min(np.abs(cent[:, None] - dumps[None, :]), axis=1) * 24
near = okc & (ddump >= 0.4) & (ddump <= 3.6); farm = okc & (ddump > 6)
res['S6_dump_check'] = dict(
    n_dumps=int(len(dumps)), dump_spacing_d=float(np.median(np.diff(dumps))),
    single_box_snr_near_dumps_0p4_3p6h=dict(mean=float(np.nanmean(z1[near])), sd=float(np.nanstd(z1[near])),
                                            p95=float(np.nanpercentile(z1[near], 95)), n=int(near.sum())),
    single_box_snr_far_gt6h=dict(mean=float(np.nanmean(z1[farm])), sd=float(np.nanstd(z1[farm])),
                                 p95=float(np.nanpercentile(z1[farm], 95)), n=int(farm.sum())),
    note='single T14 boxes on S6 data with true A transits removed; S/N in white-noise units')

json.dump(res, open(f'a_fap_{TAG}.json', 'w'), indent=1, default=float)
np.savez_compressed(f'_diag_{TAG}.npz',
                    z6=snr(6, *train(6, P_pred, T_ref + lever(6) * P_pred, np.arange(0, P_pred, STEP))),
                    z32=snr(32, *train(32, P_pred, T_ref + lever(32) * P_pred, np.arange(0, P_pred, STEP))),
                    taus=np.arange(0, P_pred, STEP),
                    null_opt_joint=nulls['N2']['adopted']['opt_joint'],
                    null_opt_joint_lik=nulls['N2']['likelihood']['opt_joint'],
                    null_fixed_joint=nulls['N2']['adopted']['fixed_joint'],
                    null_wpred_joint=nulls['N2']['adopted']['wpred_joint'])
print(json.dumps({k: res[k] for k in ('observed', 'discovery_centre_fixed_snr', 'true_ephemeris_snr_reference_only',
                                      'N1_fixed_period_phase_null', 'recovered', 'depths', 'cycle_count',
                                      'S6_dump_check', 'prediction_vs_final', 'runtime_s')}, indent=1, default=float))
for nn in R:
    for v in R[nn]:
        for k, e in R[nn][v].items():
            print(nn, v, k, {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in e.items()
                             if kk in ('obs', 'k', 'N', 'fap', 'sigma', 'p')},
                  e.get('fap_conservative', ''), e.get('gpd_extrap', e.get('gauss_extrap', '')))
