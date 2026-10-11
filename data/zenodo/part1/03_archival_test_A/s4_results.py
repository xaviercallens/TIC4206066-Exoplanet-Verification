"""Assemble a_fap_results.json from s1 (S98-only prediction) and the s2 runs (main + detrending variants)."""
import json
import numpy as np

s1 = json.load(open('s1_s98_prediction.json'))
runs = {t: json.load(open(f'a_fap_{t}.json')) for t in ('main', 'edge0375', 'w050', 'w100')}
M = runs['main']
T14 = 1.87


def g(nn, v, k, run='main'):
    return runs[run]['N2'][nn][v][k]


def brief(e):
    out = dict(obs=round(e['obs'], 3), k=e['k'], N=e['N'], fap_counting=e['fap'], sigma_counting=round(e['sigma'], 2),
               null_median_p99_max=[round(x, 2) for x in e['null_p50_p99_max']])
    if 'N_eff_indep' in e:
        out.update(N_eff_independent=round(e['N_eff_indep']), fap_conservative=e['fap_conservative'],
                   sigma_conservative=round(e['sigma_conservative'], 2))
        gp = e['gpd_extrap']
        out['gpd_tail_extrapolation'] = dict(p=gp.get('p'), shape=gp.get('shape'), endpoint=gp.get('endpoint'),
                                             note=gp.get('note', ''))
    if 'gauss_extrap' in e:
        out['gauss_extrapolation_p'] = e['gauss_extrap']['p']
    return out


prop = {}
for name in ('adopted', 'discovery_centre_adopted_sigma'):
    prop[name] = {}
    for sec in ('S6', 'S32'):
        rows = s1['propagation'][name][sec]
        mid = [r for r in rows if r['epoch'] in (-792, -569)][0]
        prop[name][sec] = dict(epochs=[rows[0]['epoch'], rows[-1]['epoch']], mid_epoch=mid['epoch'],
                               t_pred_mid_BTJD=mid['t_pred'], sigma_h=mid['sigma_h'],
                               sigma_over_T14=mid['sigma_over_T14'],
                               multisector_minus_pred_h=mid['final_minus_pred_h'],
                               multisector_minus_pred_sigma=mid['final_minus_pred_sigma'],
                               all_epochs=[dict(epoch=r['epoch'], t_pred=r['t_pred']) for r in rows])
lik = s1['likelihood']
sig_lik = {}
for sec, k in (('S6', -798), ('S32', -575)):
    s = np.sqrt(lik['sT'] ** 2 + k ** 2 * lik['sP'] ** 2 + 2 * k * lik['corr'] * lik['sT'] * lik['sP'])
    sig_lik[sec] = dict(sigma_h=round(float(s * 24), 2), sigma_over_T14=round(float(s * 24 / T14), 2))

res = dict(
    target='TIC 4206066 (StKM 1-561), signal A (3.18 d)',
    question='How often would noise produce a dip this significant where the S98-only prediction points, in the '
             'archival S6 (TESS-SPOC 30 min PDCSAP) and S32 (SPOC 2 min PDCSAP) data?',
    chronology=dict(
        discovery_record='candidates/tic4206066_p3.182945/candidate.json, S98 QLP search at 2026-09-25 00:09 CEST: '
                         'P = 3.182945191903069 d, T0 = BTJD 3998.73484, no period uncertainty stated (grid step ~2e-4 d)',
        first_archival_look='04:45-04:50 CEST (A03 harmonic audit QLP downloads; pool multisector pass 1). Every '
                            'S98-only fit other than the discovery record was made after this, by analysts who had '
                            'seen archival data; the S98-only fit here (s1) is likewise post hoc in time but uses S98 '
                            'data only and a rule written before it was run.'),
    prediction=dict(
        method='QLP SYS_RM_FLUX QUALITY==0, biweight 0.75 d (A masked at the discovery ephemeris, B at the '
               'reconciled one), trapezoid fitted on S98 alone, (P, T_ref) profile likelihood with beta-inflated '
               'errors and a 9-event bootstrap; adopted sigma = larger of the two',
        P=s1['adopted']['P'], sigma_P_adopted_bootstrap=s1['adopted']['sP'], sigma_P_likelihood=lik['sP'],
        T_ref_BTJD=s1['adopted']['Tref'], T_ref_epoch=s1['adopted']['NREF'], sigma_T_ref_d=s1['adopted']['sT'],
        corr_P_T=s1['adopted']['corr'], beta_red_noise=s1['beta'],
        S98_events=s1['epochs_S98'], S98_depth_trapezoid_ppm=s1['global_fit']['depth_ppm'],
        S98_T14_h=s1['global_fit']['T14_h'],
        per_transit_linear_P=[s1['per_transit']['P'], s1['per_transit']['sP']],
        per_transit_note='E-1 (data-gap edge) has O-C +107 min; per-transit WLS not adopted',
        other_S98_only_estimates=s1['S98_only_elsewhere'],
        P_pred_minus_multisector=M['prediction_vs_final']['P_pred_minus_final'],
        in_sigma_adopted=M['prediction_vs_final']['in_sigma_adopted'],
        in_sigma_likelihood=M['prediction_vs_final']['in_sigma_likelihood'],
        propagation=prop, propagation_likelihood_sigma=sig_lik,
        verdict=('YES: the propagated 1-sigma timing uncertainty (24.8 h in S6, 17.9 h in S32 with the adopted '
                 'bootstrap sigma; 10.9 h and 7.8 h with the likelihood-only sigma) is 4-13 times the 1.87 h transit '
                 'duration. The +-3 sigma window covers the whole 76 h orbit (adopted) or 62-86 % of it (likelihood). '
                 'The prediction centre misses the archival transits by 12.9 h (S6) and 9.3 h (S32) = 0.5 sigma; '
                 'the discovery-record centre misses by 3.1 h and 2.2 h (0.1 sigma); both misses exceed T14/2.')),
    statistic=M['statistic'],
    nulls=dict(
        N1='fixed period P_pred, same unmasked data, all 5-min phase offsets whose boxes stay >= T14 from the true '
           'transits (872 grid phases, ~38 independent); MAIN (fixed-epoch) test only',
        N2='fake periods U(2.5, 4.5) d, true A transits removed (+-1.87 h), random anchors, identical machinery '
           '(same sigma, lever arms, dP grid, window, weights); 20 000 trials (~164 independent periods in S6, ~195 '
           'in S32; joint ~3 700 with random relative anchors)',
        N2r='as N2 with P\' in U(3.10, 3.26) d, which brackets the 3.125 d S6 momentum-dump cadence; 4 000 trials',
        conservative_rule='fap_conservative = (k+1)/(N_eff+1) with N_eff the number of statistically independent '
                          'null realisations (not the number of trials)'),
    MAIN_fixed_epoch=dict(
        adopted_centre=dict(
            S6=dict(N2=brief(g('N2', 'adopted', 'fixed_S6')), N1=M['N1_fixed_period_phase_null']['S6_adopted_centre']),
            S32=dict(N2=brief(g('N2', 'adopted', 'fixed_S32')), N1=M['N1_fixed_period_phase_null']['S32_adopted_centre']),
            joint=dict(N2=brief(g('N2', 'adopted', 'fixed_joint'))),
            fisher=g('N2', 'adopted', 'fixed_fisher_counting')),
        discovery_centre=dict(snr=M['discovery_centre_fixed_snr'],
                              S6_N1=M['N1_fixed_period_phase_null']['S6_discovery_centre'],
                              S32_N1=M['N1_fixed_period_phase_null']['S32_discovery_centre']),
        verdict='No dip at the predicted times: S/N -0.51 (S6), +0.50 (S32), joint +0.01; FAP 0.61 / 0.25 / 0.34. '
                'The pre-specified fixed-epoch test is uninformative because the prediction is far less precise '
                'than T14; it does NOT support a statement like "FAP X at the S98-predicted times".'),
    OPTIMISED=dict(
        adopted_sigma={k: brief(g('N2', 'adopted', k)) for k in ('opt_S6', 'opt_S32', 'opt_joint', 'wpred_S6',
                                                                   'wpred_S32', 'wpred_joint')},
        adopted_sigma_fisher=dict(conservative=g('N2', 'adopted', 'opt_fisher_conservative'),
                                  counting_if_all_trials_independent=g('N2', 'adopted', 'opt_fisher_counting')),
        likelihood_sigma={k: brief(g('N2', 'likelihood', k)) for k in ('opt_S6', 'opt_S32', 'opt_joint')},
        likelihood_sigma_fisher_conservative=g('N2', 'likelihood', 'opt_fisher_conservative'),
        N2r_dump_bracketing={k: brief(g('N2r', 'adopted', k)) for k in ('opt_S6', 'opt_S32', 'opt_joint')},
        fixed_period_phase_null_zscore=dict(
            S6=dict(true_phase_snr=M['N1_fixed_period_phase_null']['S6_true_phase_snr'],
                    null_sd=M['N1_fixed_period_phase_null']['S6_null_sd'],
                    null_max=M['N1_fixed_period_phase_null']['S6_null_max']),
            S32=dict(true_phase_snr=M['N1_fixed_period_phase_null']['S32_true_phase_snr'],
                     null_sd=M['N1_fixed_period_phase_null']['S32_null_sd'],
                     null_max=M['N1_fixed_period_phase_null']['S32_null_max'])),
        verdict=('Optimising the epoch (and dP) inside the 3-sigma prediction window recovers the transits: S/N 9.8 '
                 '(S6), 9.0 (S32), 13.2 joint, within 7 and 2 min of the multi-sector ephemeris. No null trial '
                 'reaches these values (0 of 20 000 each). Because only ~160-195 null realisations per sector are '
                 'independent, the empirical per-sector FAP is only bounded: <= 6.1e-3 (S6), <= 5.1e-3 (S32), ~2.5-2.6 '
                 'sigma; joint < 3e-4 (~3.5 sigma, N_eff ~3 700); Fisher of the per-sector bounds < 3.5e-4 (3.4 '
                 'sigma). S6 clears the null maximum by a factor 2.1 (9.8 vs 4.6); S32 only just (9.0 vs 8.3, the S32 '
                 'null tail comes from segment-edge detrending artefacts; with a 0.375 d edge cut the S32 null '
                 'maximum drops to 5.2 while the observed stays 8.85). Tail extrapolations (GPD) put the observed '
                 'values beyond the fitted finite endpoint of the null, i.e. formally p ~ 0; not quoted.')),
    detrending_variants={t: dict(window_d=runs[t]['window_d'], edge_cutoff_d=runs[t]['edge_cutoff_d'],
                                 obs_opt=[round(runs[t]['observed']['adopted'][k], 2) for k in
                                          ('opt_S6', 'opt_S32', 'opt_joint')],
                                 obs_fixed=[round(runs[t]['observed']['adopted'][k], 2) for k in
                                            ('fixed_S6', 'fixed_S32', 'fixed_joint')],
                                 null_max_opt=[round(g('N2', 'adopted', k, t)['null_p50_p99_max'][2], 2) for k in
                                               ('opt_S6', 'opt_S32', 'opt_joint')],
                                 k_opt=[g('N2', 'adopted', k, t)['k'] for k in ('opt_S6', 'opt_S32', 'opt_joint')],
                                 N=g('N2', 'adopted', 'opt_S6', t)['N'],
                                 fisher_conservative_p=g('N2', 'adopted', 'opt_fisher_conservative', t)['p'])
                         for t in runs},
    recovered=dict(M['recovered'],
                   note='archival-only best (tau_S6, tau_S32) mapped to (dT, dP) at the S98 reference epoch; dT is '
                        'poorly constrained by archival data alone (a 7 min error in tau_S6 moves dT by ~17 min), '
                        'so the dT sigma is not a meaningful tension'),
    depths=M['depths'],
    cycle_count=dict(M['cycle_count'],
                     verdict=('S98 + S6: 798 cycles; next alias (799) at 2.6 sigma (adopted) / 5.8 sigma '
                              '(likelihood); 16.8 % of the S98 event bootstrap lies nearer the 799 alias. S98 + S32: '
                              '575 cycles; 576 at 3.7 / 8.5 sigma; bootstrap 10.6 %. Unique under the likelihood error, '
                              'not robustly unique under the heavy-tailed S98 bootstrap with one archival sector alone; '
                              'S6+S32 together (223-cycle lever) fix it (reconciled A.cycle_count).')),
    S6_dump_check=M['S6_dump_check'],
    caveats=[
        'The prediction is weak: S98-only P estimates range 3.18164-3.18276 d across methods (this fit 3.18211); '
        'propagated over ~800 cycles the uncertainty is 4-13 x T14, so "at the predicted time" really means "at the '
        'predicted period, anywhere in most of the phase".',
        'All S98-only fits except the discovery record were made after the first archival look (04:45 CEST); this one '
        'too (rule fixed in the script before running, S98 data only).',
        'Per-sector FAPs are upper bounds set by the ~160-195 independent null periods available in one sector; '
        'counting 0/20 000 as independent (5e-5, 3.9 sigma) would overstate the evidence.',
        'S/N is in white-noise units; red noise is handled only through the empirical nulls (N1 null SD 1.06 in S6, '
        '1.41 in S32).',
        'S32 null tail is set by segment-edge detrending artefacts (single box S/N up to 11 at the sector end); the '
        'observed S32 value barely exceeds the null maximum without an edge cut.',
        'S6 momentum dumps (every 3.125 d) lie 0.4-3.5 h from every S6 A transit. Boxes 0.4-3.6 h from dumps (true '
        'transits removed) show mean S/N -0.71 (brightening side), SD 0.93, vs -0.05 / 1.06 far from dumps, so dumps '
        'do not manufacture dips; only ~8 independent dump-adjacent boxes, so this check is weak. The N2r null '
        '(P\' 3.10-3.26 d, including dump-aligned trains) has lower maxima than N2.',
        'Detrending: biweight 0.5/0.75/1.0 d and a 0.375 d edge cut change the optimised S/N by -1.1/+0.7 (joint '
        '12.1-13.9) and never produce a null exceedance.',
        'A is never masked in the archival detrending (same for all trials); this shrinks box depths by ~7-9 % '
        '(S6 422 vs 454 ppm masked).',
        'Depth ratio uses box depths; the S98 epoch -1 event (data-gap edge) drags the S98 box depth from 389 to 299 '
        'ppm; ratios quoted without it.'],
    printable=dict(
        sentence='(false-alarm probability <0.7% in each of S6 and S32, <0.03% jointly, at the S98-only period)',
        words=15,
        do_not_print='"FAP X at the S98-predicted times": at the predicted times there is no dip; the recovery needs '
                     'an epoch search over most of the orbit (look-elsewhere included in the numbers above).'),
    files=['s1_s98_prediction.py', 's1_s98_prediction.json', 's1_s98_prediction.log', 's1_boot_P_Tref.npy',
           'fap_common.py', 's2_archival_fap.py', 'a_fap_main.json', 'a_fap_edge0375.json', 'a_fap_w050.json',
           'a_fap_w100.json', 's2_*.log', '_diag_main.npz', 's3_plot.py', 'a_fap_diagnostic.png', 's4_results.py',
           'a_fap_results.json'],
)
json.dump(res, open('a_fap_results.json', 'w'), indent=1, default=float)
print(json.dumps(res['OPTIMISED']['adopted_sigma_fisher'], indent=1))
print(json.dumps(res['detrending_variants'], indent=1))
