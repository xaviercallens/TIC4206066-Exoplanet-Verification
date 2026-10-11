"""Common code for Audit 1 (Signal A ephemeris) of TIC 4206066.
Written from scratch for this audit; no code or derived tables from earlier analyses are used.
Times are BTJD = BJD_TDB - 2457000 (all products: TIMESYS=TDB, TIMEREF=SOLARSYSTEM, BJDREFI=2457000,
TIMEPIXR=0.5 [SPOC], i.e. mid-exposure)."""
import os
import numpy as np
from astropy.io import fits

os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '00_data')

# ---- manuscript values (the things being audited) ----
MS_PA, MS_ePA = 3.1827856, 0.0000058
MS_T0A, MS_eT0A = 2460998.7375 - 2457000.0, 0.0023      # BTJD 3998.7375
MS_PB, MS_ePB = 11.13274, 0.00003
MS_T0B, MS_eT0B = 2461000.865 - 2457000.0, 0.003        # BTJD 4000.865
T14A = 1.87 / 24.0       # manuscript T14 (used only for masks/coverage windows)
T14B = 1.65 / 24.0       # manuscript 1.6-1.7 h

FILES = {
    6: 'mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits',
    32: 'mast/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits',
    '32ffi': 'mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits',
    98: 'mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits',
}
EXPTIME = {6: 1800.0 / 86400, 32: 120.0 / 86400, '32ffi': 600.0 / 86400, 98: 200.0 / 86400}
SUPERSAMPLE = {6: 15, 32: 1, '32ffi': 5, 98: 3}

QLP_BIT29 = 1 << 29      # undocumented; marks first ~6 d of every S98 orbit (scattered-light half)
QLP_BIT30 = 1 << 30      # QLP quality flag (tess.mit.edu/qlp-update: bit 31 = 2^30)
STRAY = 2048             # bit 12 predicted stray light (SDPDD / QLP-update)
DUMP = 32


def sector_of(key):
    return 32 if key == '32ffi' else key


def load(key, s98_mask='clean', s98_flux='SYS_RM_FLUX'):
    """Return dict with time, flux (normalised), err, and full-cadence info for coverage bookkeeping.
    s98_mask: 'clean'  -> QUALITY == 0 (drops 2048, 2^29, 2^30, dumps, coarse point)
              'keepSL' -> drop only 2^30, dumps (32), coarse point (4): keeps scattered-light intervals
              'no30'   -> drop 2048, 2^29, dumps, coarse; keep 2^30 cadences
    """
    with fits.open(os.path.join(DATA, FILES[key])) as h:
        d = h[1].data
        t = np.array(d['TIME'], float)
        q = np.array(d['QUALITY'], np.int64)
        if key == 98:
            f = np.array(d[s98_flux], float)
            e = np.array(d['DET_FLUX_ERR'], float)
        else:
            f = np.array(d['PDCSAP_FLUX'], float)
            e = np.array(d['PDCSAP_FLUX_ERR'], float)
    fin = np.isfinite(t) & np.isfinite(f)
    if key == 98:
        if s98_mask == 'clean':
            good = fin & (q == 0)
        elif s98_mask == 'keepSL':
            good = fin & ((q & (QLP_BIT30 | DUMP | 4)) == 0)
        elif s98_mask == 'no30':
            good = fin & ((q & (STRAY | QLP_BIT29 | DUMP | 4)) == 0)
        else:
            raise ValueError(s98_mask)
    else:
        good = fin & (q == 0)
    med = np.nanmedian(f[good])
    fn = f / med
    en = e / med
    # QLP DET_FLUX_ERR may be NaN; fall back to point-to-point scatter
    if (not np.all(np.isfinite(en[good]))) or np.nanmedian(en[good]) <= 0:
        s = 1.4826 * np.median(np.abs(np.diff(fn[good]))) / np.sqrt(2)
        en = np.full_like(fn, s)
    # dump times: SPOC/TESS-SPOC bit 6 (value 32) and, for S98, QLP bit 6 (same cadences as TESScut)
    dump_t = t[np.isfinite(t) & ((q & DUMP) != 0)]
    return dict(key=key, sector=sector_of(key), t_all=t, q_all=q, fin_all=fin, good=good,
                t=t[good], f=fn[good], e=en[good], exptime=EXPTIME[key], dump_t=dump_t,
                sl_all=((q & (STRAY | QLP_BIT29)) != 0) if key == 98 else np.zeros_like(fin))


def transit_mask(t, P, T0, half_width):
    ph = (t - T0 + 0.5 * P) % P - 0.5 * P
    return np.abs(ph) < half_width


def detrend(t, f, mask, window=0.75, break_tol=0.1, edge_cutoff=0.0):
    """wotan time-windowed biweight with in-transit points (mask True) excluded from the trend."""
    from wotan import flatten
    flat, trend = flatten(t, f, method='biweight', window_length=window, break_tolerance=break_tol,
                          edge_cutoff=edge_cutoff, return_trend=True, mask=mask)
    return flat, trend


def b_mask(t, half=None):
    half = T14B if half is None else half
    return transit_mask(t, MS_PB, MS_T0B, half)


def a_mask(t, P=MS_PA, T0=MS_T0A, half=None):
    half = T14A if half is None else half
    return transit_mask(t, P, T0, half)


def prepare(key, window=0.75, s98_mask='clean', s98_flux='SYS_RM_FLUX', P=MS_PA, T0=MS_T0A,
            mask_half_A=T14A, mask_half_B=T14B, break_tol=0.1):
    """Load, detrend with A and B transits masked (half width = 1 x T14 -> 2 x T14 total), return data."""
    d = load(key, s98_mask=s98_mask, s98_flux=s98_flux)
    m = a_mask(d['t'], P, T0, mask_half_A) | b_mask(d['t'], mask_half_B)
    flat, trend = detrend(d['t'], d['f'], m, window=window, break_tol=break_tol)
    ok = np.isfinite(flat)
    d['flat'] = flat
    d['trend'] = trend
    d['ok'] = ok
    d['inB'] = b_mask(d['t'], 0.5 * T14B + 0.25 / 24)   # B in-transit (+15 min margin): excluded from A fits
    # robust per-point noise of the flattened out-of-transit data
    oot = ok & ~m
    d['sig_pp'] = 1.4826 * np.median(np.abs(np.diff(flat[oot]))) / np.sqrt(2)
    d['rms_oot'] = 1.4826 * np.median(np.abs(flat[oot] - np.median(flat[oot])))
    return d


# ---------------- transit model ----------------
_models = {}


def batman_lc(t, key, tc, P, k, aRs, b, u1, u2):
    import batman
    p = batman.TransitParams()
    p.t0 = tc; p.per = P; p.rp = k; p.a = aRs
    cosi = b / aRs
    p.inc = np.degrees(np.arccos(np.clip(cosi, 0, 1)))
    p.ecc = 0.0; p.w = 90.0; p.limb_dark = 'quadratic'; p.u = [u1, u2]
    ss = SUPERSAMPLE[key]
    m = batman.TransitModel(p, t, supersample_factor=ss, exp_time=EXPTIME[key]) if ss > 1 else batman.TransitModel(p, t)
    return m.light_curve(p)


def q_to_u(q1, q2):
    s = np.sqrt(q1)
    return 2 * s * q2, s * (1 - 2 * q2)


def epochs_in(t_lo, t_hi, P, T0):
    n0 = int(np.ceil((t_lo - T0) / P)); n1 = int(np.floor((t_hi - T0) / P))
    return np.arange(n0, n1 + 1)
