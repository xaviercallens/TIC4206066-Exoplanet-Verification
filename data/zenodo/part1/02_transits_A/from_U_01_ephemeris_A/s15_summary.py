"""Collect headline numbers (manuscript vs recomputed) into out_summary.json and out_numeric_rows.csv."""
import json, csv, os, numpy as np
from a_common import *
L = lambda f: json.load(open(f))
lsq = L('out_s05_globalfit_lsq.json'); eph = L('out_s08_ephemeris_fits.json')
em = L('out_s09_emcee_long.json') if os.path.exists('out_s09_emcee_long.json') else L('out_s09_emcee_base.json')
bls = L('out_s10_bls_aliases.json'); cyc = L('out_s11_cycle_count.json'); fr = L('out_s12_s98_frozen.json'); s110 = L('out_s13_s110.json')
ep = L('out_s03_epochs.json'); bary = L('out_s01b_bary_validation.json'); inj = L('out_s07_injection_summary_base.json')
T1 = eph['base|T1_cov07_noB_snr2|p16_p84_half']; T3 = eph['base|T3_T1_minus_S32m573|p16_p84_half']
M20 = eph['base|ms20_like|std_all']; M20c = eph['base|ms20_like|p16_p84_half']
g = lsq['raw_lin']['params']
S = dict(
 P_route1_T1=[T1['P'], T1['eP']], P_route1_T3=[T3['P'], T3['eP']], P_route1_huber=[T1['huber']['P'], T1['huber']['eP']],
 P_route1_ms20like_std=[M20['P'], M20['eP']], P_lsq=g['P'], P_emcee=[em['P']['median'], em['P']['std']], P_bls=bls['combined']['P_best'],
 T0_route1_T1=[T1['T0'], T1['eT0']], T0_route1_T3=[T3['T0'], T3['eT0']], T0_route1_huber=[T1['huber']['T0'], T1['huber']['eT0']],
 T0_lsq=g['T0'], T0_emcee=[em['T0']['median'], em['T0']['lo'], em['T0']['hi']],
 chi2_T1=[T1['chi2'], T1['dof']], chi2_T3=[T3['chi2'], T3['dof']], chi2_ms20like_std=[M20['chi2'], M20['dof']], chi2_ms20like_core=[M20c['chi2'], M20c['dof']],
 rms_T3=T3['rms_oc_min'], rms_T1=T1['rms_oc_min'],
 sector_shift=cyc['sector_shift_at_ms_ephemeris'], alt_cycles=cyc['alternative_cycle_counts_S98_S32'],
 single=cyc['single_sector_fits'], prop=cyc['propagation'], frozen=fr, s110_counts=s110['counts'],
 epochs_summary=ep['summary'], injection=inj, bary=bary, emcee_file='out_s09_emcee_long.json' if os.path.exists('out_s09_emcee_long.json') else 'out_s09_emcee_base.json',
 emcee_autocorr=em['autocorr'], emcee_T0_center=em['T0_center'], emcee_n_center=em['n_center'])
json.dump(S, open('out_summary.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in S.items() if k.startswith(('P_', 'T0_', 'chi2', 'rms', 'emcee'))}, indent=1, default=float))
