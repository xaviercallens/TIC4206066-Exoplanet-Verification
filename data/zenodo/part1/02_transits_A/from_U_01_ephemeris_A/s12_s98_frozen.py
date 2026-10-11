"""Step 7: the historical 3.182945 d (frozen S98-only BLS) vs the multi-sector ~3.18279 d.
Re-run an S98-only BLS following the frozen search configuration as read from config/prejev_search.json
(QLP SYS_RM_FLUX, bad-bit mask 1610619071, wotan biweight 1.0 d, break 0.1 d, edge cut 0.25 d, upper 4-sigma
clip, 3-cadence bins, likelihood BLS, 3-10 d band durations, df = D_min/(3*baseline*P_max)); own code.
Then S98-only period uncertainty from (i) BLS dchi2 = 1 width, (ii) template fit (s11), (iii) per-transit
linear ephemeris of the S98 events (s06/s07)."""
import numpy as np, json, csv
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from wotan import flatten
from a_common import *

h = fits.open(os.path.join(DATA, FILES[98])); d = h[1].data
t = np.array(d['TIME'], float); f = np.array(d['SYS_RM_FLUX'], float); q = np.array(d['QUALITY'], np.int64)
good = np.isfinite(t) & np.isfinite(f) & ((q & 1610619071) == 0)
t = t[good]; f = f[good] / np.nanmedian(f[good])
flat = flatten(t, f, method='biweight', window_length=1.0, break_tolerance=0.1, edge_cutoff=0.25)
ok = np.isfinite(flat); t = t[ok]; y = flat[ok]
s = 1.4826 * np.median(np.abs(y - np.median(y)))
keep = y < np.median(y) + 4 * s; t = t[keep]; y = y[keep]
# bin by 3 cadences (consecutive)
nb = len(t) // 3
tb = t[:nb * 3].reshape(nb, 3); yb = y[:nb * 3].reshape(nb, 3)
cont = (tb[:, 2] - tb[:, 0]) < 3 * 200 / 86400.
tb = tb[cont].mean(1); yb = yb[cont].mean(1)
eb = np.full(len(tb), 1.4826 * np.median(np.abs(np.diff(yb))) / np.sqrt(2))
baseline = t.max() - t.min()
durs = np.array([0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5]) / 24
df = (0.75 / 24) / (3 * baseline * 10.0)
freqs = np.arange(1 / 3.30, 1 / 3.07, df)
P = np.sort(1 / freqs)
bls = BoxLeastSquares(tb, yb, dy=eb)
r = bls.power(P, durs, objective='likelihood', oversample=5)
j = int(np.argmax(r.power))
out = dict(n_unbinned=int(len(t)), n_binned=int(len(tb)), baseline_d=float(baseline), df=float(df),
           dP_at_3p18=float(3.18 ** 2 * df), coarse_P=float(P[j]), coarse_dur_h=float(r.duration[j] * 24),
           coarse_T0=float(r.transit_time[j]), coarse_depth_ppm=float(r.depth[j] * 1e6), coarse_dchi2=float(2 * r.power[j]))
# refine: fine period grid (1e-6 d) +- 40 coarse steps, durations 0.25 h steps
Pf = np.arange(P[j] - 40 * out['dP_at_3p18'], P[j] + 40 * out['dP_at_3p18'], 1e-6)
dursf = np.arange(1.0, 3.01, 0.25) / 24
rf = bls.power(Pf, dursf, objective='likelihood', oversample=20)
jf = int(np.argmax(rf.power)); dch = 2 * np.array(rf.power)
inside = Pf[dch >= dch[jf] - 1]
out.update(dict(fine_P=float(Pf[jf]), fine_dur_h=float(rf.duration[jf] * 24), fine_T0=float(rf.transit_time[jf]),
                fine_dchi2=float(dch[jf]), fine_dchi2eq1_halfwidth=float(0.5 * (inside.max() - inside.min())),
                fine_dchi2eq1_range=[float(inside.min()), float(inside.max())]))
# unbinned fine BLS as a check
bls_u = BoxLeastSquares(t, y, dy=np.full(len(t), 1.4826 * np.median(np.abs(np.diff(y))) / np.sqrt(2)))
ru = bls_u.power(Pf, dursf, objective='likelihood', oversample=20)
ju = int(np.argmax(ru.power)); dcu = 2 * np.array(ru.power); insu = Pf[dcu >= dcu[ju] - 1]
out.update(dict(unbinned_fine_P=float(Pf[ju]), unbinned_dchi2eq1_halfwidth=float(0.5 * (insu.max() - insu.min()))))

# (ii) template fit S98-only (s11) and (iii) per-transit timing S98-only
cc = json.load(open('out_s11_cycle_count.json'))['single_sector_fits']['S98']
out['template_fit_S98'] = dict(P=cc['P'], eP_dchi2_beta=cc['eP'], eP_boot=cc['P_boot_std'])
rows = [r for r in csv.DictReader(open('out_s06_per_transit_all_variants.csv')) if r['variant'] == 'base' and r['sector'] == '98'
        and r.get('tc', '') not in ('', 'nan') and float(r['coverage']) >= 0.7 and float(r['snr']) >= 2 and r['edge'] == 'False']
inj = json.load(open('out_s07_injection_summary_base.json'))['base_S98']['cov07']
n = np.array([int(r['epoch']) for r in rows], float); tc = np.array([float(r['tc']) for r in rows])
for metric in ['p16_p84_half', 'std_all']:
    e = np.full(len(n), inj[metric] / 1440)
    A = np.vstack([np.ones_like(n), n]).T; C = np.linalg.inv(A.T @ A / e[0] ** 2); b = C @ (A.T @ tc) / e[0] ** 2
    out[f'pertransit_S98_{metric}'] = dict(N=len(n), P=float(b[1]), eP=float(np.sqrt(C[1, 1])))
FROZEN = 3.182945191903069
for k in ['fine_dchi2eq1_halfwidth']:
    pass
ref = dict(ms=MS_PA, own_global=json.load(open('out_s05_globalfit_lsq.json'))['raw_lin']['params']['P'][0])
sig = dict(bls_dchi2eq1=out['fine_dchi2eq1_halfwidth'], template_dchi2_beta=cc['eP'], template_boot=cc['P_boot_std'],
           pertransit_core=out['pertransit_S98_p16_p84_half']['eP'], pertransit_std=out['pertransit_S98_std_all']['eP'])
out['frozen_P'] = FROZEN
out['frozen_minus_ms_d'] = FROZEN - MS_PA
out['difference_in_sigma'] = {k: float((FROZEN - MS_PA) / v) for k, v in sig.items()}
out['own_S98_BLS_minus_frozen_d'] = out['fine_P'] - FROZEN
print(json.dumps(out, indent=1))
json.dump(out, open('out_s12_s98_frozen.json', 'w'), indent=1)
