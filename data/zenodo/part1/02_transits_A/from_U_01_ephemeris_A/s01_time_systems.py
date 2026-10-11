"""Audit 1 / step 1: time systems of every product used for the Signal A ephemeris.
Reads FITS headers (own code), compares timestamps between products, and quantifies the
TESScut (FFI, CCD-reference-pixel barycentric correction) vs target-direction correction."""
import json, numpy as np, os
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric
import astropy.units as u
from astropy.wcs import WCS

D = '../00_data/'
files = {
 'S6_tess-spoc_lc': 'mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits',
 'S32_spoc_2min_lc': 'mast/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits',
 'S32_tess-spoc_lc': 'mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits',
 'S98_qlp': 'mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits',
 'S32_qlp': 'mast/hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits',
 'S6_qlp': 'mast/hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits',
 'S6_tesscut': 'tesscut/21x21/tess-s0006-1-1_80.859830_-1.326960_21x21_astrocut.fits',
 'S32_tesscut': 'tesscut/21x21/tess-s0032-1-2_80.859830_-1.326960_21x21_astrocut.fits',
 'S98_tesscut': 'tesscut/21x21/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits',
}
KEYS = ['TIMESYS','TIMEREF','TASSIGN','BJDREFI','BJDREFF','BJDREFR','TIMEUNIT','TIMEPIXR','TIERRELA',
        'TSTART','TSTOP','TIMEDEL','INT_TIME','NUM_FRM','FRAMETIM','DEADC','BTC_PIX1','BTC_PIX2','BARYCORR',
        'RA_OBJ','DEC_OBJ','CAMERA','CCD','SECTOR','PROCVER','DATA_REL','CALIB','DETMET']
out = {}
for k, f in files.items():
    with fits.open(D + f, memmap=True) as h:
        rec = {}
        for ext in (0, 1):
            for key in KEYS:
                if key in h[ext].header and key not in rec:
                    v = h[ext].header[key]
                    rec[key] = v if not isinstance(v, (np.generic,)) else v.item()
                    rec[key + '_ext'] = ext
        cols = h[1].columns.names
        rec['columns'] = cols
        t = h[1].data['TIME']
        rec['n'] = int(len(t)); rec['time_first'] = float(np.nanmin(t)); rec['time_last'] = float(np.nanmax(t))
        if 'TIMECORR' in cols:
            tc = h[1].data['TIMECORR']
            rec['timecorr_min_d'] = float(np.nanmin(tc)); rec['timecorr_max_d'] = float(np.nanmax(tc))
        # comments of time keywords
        rec['TIME_TUNIT'] = h[1].header.get('TUNIT1')
        out[k] = rec

# --- empirical timestamp comparisons between products of the same sector ---
def load(k, col='TIME'):
    with fits.open(D + files[k]) as h:
        d = h[1].data
        return np.array(d['TIME'], float), (np.array(d['TIMECORR'], float) if 'TIMECORR' in d.columns.names else None), np.array(d['CADENCENO'])

cmp = {}
def nearest_diff(ta, tb, tol):
    ta = ta[np.isfinite(ta)]; tb = tb[np.isfinite(tb)]
    idx = np.searchsorted(tb, ta); idx = np.clip(idx, 1, len(tb) - 1)
    d1 = ta - tb[idx - 1]; d2 = ta - tb[idx]
    d = np.where(np.abs(d1) < np.abs(d2), d1, d2)
    m = np.abs(d) < tol
    return d[m] * 86400.0
for a, b, tol in [('S98_qlp','S98_tesscut',100/86400), ('S32_qlp','S32_tesscut',300/86400), ('S6_qlp','S6_tesscut',900/86400),
                  ('S6_tess-spoc_lc','S6_tesscut',900/86400), ('S32_tess-spoc_lc','S32_tesscut',300/86400)]:
    ta, _, _ = load(a); tb, _, _ = load(b)
    d = nearest_diff(ta, tb, tol)
    cmp[f'{a}_minus_{b}_s'] = dict(n=int(len(d)), median=float(np.median(d)), min=float(d.min()), max=float(d.max()))

# --- barycentric correction: target direction vs TESScut (FFI reference pixel) ---
# TESScut TIME/TIMECORR come from the FFIs, whose barycentric correction is computed for the
# reference pixel BTC_PIX1/BTC_PIX2 of the CCD (SDPDD FFI header). Compute the difference of the
# Roemer delay n.r/c between the target direction and the CCD reference direction.
target = SkyCoord(80.85983*u.deg, -1.32696*u.deg)
bary = {}
for k in ['S6_tesscut','S32_tesscut','S98_tesscut']:
    with fits.open(D + files[k]) as h:
        t = np.array(h[1].data['TIME'], float); tc = np.array(h[1].data['TIMECORR'], float)
        hdr0 = h[0].header; hdr1 = h[1].header; hdr2 = h[2].header if len(h) > 2 else None
    good = np.isfinite(t)
    tt = t[good][::max(1, good.sum() // 200)]
    tcc = tc[good][::max(1, good.sum() // 200)]
    # reference direction for FFI barycentric correction: CCD reference pixel -> use camera/CCD geometry via
    # tess-point to convert (col,row) of the reference pixel to RA/Dec.
    bary[k] = dict(sector=int(hdr0['SECTOR']), camera=int(hdr0['CAMERA']), ccd=int(hdr0['CCD']))
    # Roemer delay for target direction, using Earth-Moon barycentre ~ TESS position (TESS orbit < 0.003 AU -> <1.5 s)
    times = Time(tt - tcc + 2457000.0, format='jd', scale='tdb')   # arrival time at spacecraft (approx)
    pos = get_body_barycentric('earth', times)
    r = np.vstack([pos.x.to(u.m).value, pos.y.to(u.m).value, pos.z.to(u.m).value])
    c = 299792458.0
    n_t = target.cartesian.xyz.value[:, None]
    roemer_target_d = (n_t * r).sum(0) / c / 86400.0
    bary[k]['timecorr_col_minus_roemer_target_s'] = dict(median=float(np.median((tcc - roemer_target_d) * 86400)),
        min=float(np.min((tcc - roemer_target_d) * 86400)), max=float(np.max((tcc - roemer_target_d) * 86400)))
    bary[k]['note'] = ('TIMECORR (FFI, CCD reference pixel) minus Earth-based Roemer delay toward the target; the '
                       'residual contains (i) direction difference target vs CCD reference pixel and (ii) TESS-Earth offset '
                       '(<~1.3 s) and relativistic terms (<~2 ms).')
out['timestamp_comparisons'] = cmp
out['tesscut_barycentric'] = bary
json.dump(out, open('out_s01_time_systems.json', 'w'), indent=1, default=str)
for k, v in out.items():
    print(k, {kk: vv for kk, vv in v.items() if kk not in ('columns',)} if isinstance(v, dict) else v)
