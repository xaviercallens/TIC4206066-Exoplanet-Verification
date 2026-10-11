"""Step 1b: validate each product's barycentric correction for the TARGET direction.
For SPOC/TESS-SPOC LCs: TIMECORR (target-specific) vs astropy Roemer delay toward the target (Earth position;
TESS-Earth offset < ~1.3 light-s). For QLP (no TIMECORR column): reconstruct spacecraft time from the matching
TESScut cadence (TIME - TIMECORR) and compare QLP TIME with spacecraft time + target Roemer delay."""
import numpy as np, json
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric
import astropy.units as u
from a_common import DATA
import os
tgt = SkyCoord(80.85983 * u.deg, -1.32696 * u.deg)
n_t = tgt.cartesian.xyz.value
c = 299792458.0
def roemer(t_sc_btjd):
    tt = Time(t_sc_btjd + 2457000.0, format='jd', scale='tdb')
    p = get_body_barycentric('earth', tt)
    r = np.vstack([p.x.to(u.m).value, p.y.to(u.m).value, p.z.to(u.m).value])
    return (n_t[:, None] * r).sum(0) / c / 86400.0
def rd(f):
    with fits.open(os.path.join(DATA, f)) as h:
        d = h[1].data
        return {k: np.array(d[k], float) for k in d.columns.names if k in ('TIME', 'TIMECORR')}
out = {}
for name, f in [('S6_tess-spoc', 'mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits'),
                ('S32_spoc_2min', 'mast/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits')]:
    d = rd(f); w = np.isfinite(d['TIME']) & (d['TIMECORR'] != 0)
    tsc = d['TIME'][w] - d['TIMECORR'][w]
    res = (d['TIMECORR'][w] - roemer(tsc)) * 86400
    out[name] = dict(n=int(w.sum()), timecorr_minus_roemer_target_s=dict(median=float(np.median(res)), min=float(res.min()), max=float(res.max())))
for sec, fq, fc, tol in [(98, 'mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits', 'tesscut/21x21/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits', 100),
                         (32, 'mast/hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits', 'tesscut/21x21/tess-s0032-1-2_80.859830_-1.326960_21x21_astrocut.fits', 300),
                         (6, 'mast/hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits', 'tesscut/21x21/tess-s0006-1-1_80.859830_-1.326960_21x21_astrocut.fits', 900)]:
    q = rd(fq); tc_ = rd(fc)
    tq = q['TIME'][np.isfinite(q['TIME'])]
    tt = tc_['TIME']; tcorr = tc_['TIMECORR']; ok = np.isfinite(tt); tt = tt[ok]; tcorr = tcorr[ok]
    idx = np.clip(np.searchsorted(tt, tq), 1, len(tt) - 1)
    j = np.where(np.abs(tq - tt[idx - 1]) < np.abs(tq - tt[idx]), idx - 1, idx)
    m = np.abs(tq - tt[j]) * 86400 < tol
    tsc = tt[j][m] - tcorr[j][m]
    pred = tsc + roemer(tsc)
    res = (tq[m] - pred) * 86400
    out[f'S{sec}_qlp'] = dict(n=int(m.sum()), qlp_minus_spacecraft_plus_target_roemer_s=dict(median=float(np.median(res)), min=float(res.min()), max=float(res.max()), std=float(np.std(res))))
    # tesscut: TIME vs target-corrected
    res2 = (tt - (tt - tcorr + roemer(tt - tcorr))) * 86400
    out[f'S{sec}_tesscut'] = dict(tesscut_minus_target_corrected_s=dict(median=float(np.median(res2)), min=float(res2.min()), max=float(res2.max())))
print(json.dumps(out, indent=1))
json.dump(out, open('out_s01b_bary_validation.json', 'w'), indent=1)
