"""Step 8: Sector 110. (i) is TIC 4206066 on silicon (tess-point 0.9.5.post1, internal pointing table and the
HEASARC S110 pointing as an override), (ii) A and B transit times inside the HEASARC S110 dates
(2026 Oct 31 - Nov 26) with propagated 1-sigma timing uncertainties (manuscript ephemerides), in BJD_TDB and in
UTC at the observer (geocentre; TESS-Earth light time < 1.3 s)."""
import numpy as np, json
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import tess_stars2px as tp
from a_common import *

ra, dec = 80.85983, -1.32696
out = {}
res = tp.tess_stars2px_function_entry(4206066, ra, dec)
out['tesspoint_internal'] = [dict(sector=int(s), camera=int(c), ccd=int(d), col=float(x), row=float(y))
                             for s, c, d, x, y in zip(res[3], res[4], res[5], res[6], res[7])]
sc = res[8]
i110 = list(sc.sectors).index(110)
out['tesspoint_S110_pointing'] = dict(ra=float(sc.ras[i110]), dec=float(sc.decs[i110]), roll=float(sc.rolls[i110]))
out['heasarc_S110_pointing'] = dict(ra=106.07, dec=-23.79, roll=238.70, camera1=[75.88, -2.69, 321.51],
                                    url='https://heasarc.gsfc.nasa.gov/docs/tess/sector110_summary.html')
# override with HEASARC pointing
with open('docs/s110_override.txt', 'w') as fh:
    fh.write('110 106.07 -23.79 238.70\n')
try:
    r2 = tp.tess_stars2px_function_entry(4206066, ra, dec, sectorOverrideFile='docs/s110_override.txt')
    out['tesspoint_heasarc_override'] = [dict(sector=int(s), camera=int(c), ccd=int(d), col=float(x), row=float(y))
                                         for s, c, d, x, y in zip(r2[3], r2[4], r2[5], r2[6], r2[7])]
except Exception as ex:
    out['tesspoint_heasarc_override'] = str(ex)
print(json.dumps({k: v for k, v in out.items()}, indent=1))

target = SkyCoord(ra * u.deg, dec * u.deg)
geo = EarthLocation.from_geocentric(0 * u.m, 0 * u.m, 0 * u.m)
def bjd_to_utc(bjd_tdb):
    """observer (geocentre) UTC for a barycentric event time: iterate t_obs + ltt(t_obs) = BJD."""
    t = Time(bjd_tdb, format='jd', scale='tdb', location=geo)
    for _ in range(3):
        ltt = t.light_travel_time(target, kind='barycentric')
        t = Time(bjd_tdb, format='jd', scale='tdb', location=geo) - ltt
    return t.utc

def transits(P, eP, T0, eT0, lo_utc, hi_utc):
    lo = Time(lo_utc, scale='utc').tdb.jd - 0.01; hi = Time(hi_utc, scale='utc').tdb.jd + 0.01   # widened; exact below
    n0 = int(np.floor((lo - T0) / P)) - 1; n1 = int(np.ceil((hi - T0) / P)) + 1
    rows = []
    for n in range(n0, n1 + 1):
        tb = T0 + n * P
        s = np.sqrt(eT0 ** 2 + (n * eP) ** 2)
        tu = bjd_to_utc(tb)
        rows.append(dict(n=n, bjd_tdb=tb, btjd=tb - 2457000, sigma_min=s * 1440, utc=tu.iso))
    return rows

wins = {'Oct31_00_to_Nov26_00': ('2026-10-31 00:00:00', '2026-11-26 00:00:00'),
        'Oct31_00_to_Nov26_24': ('2026-10-31 00:00:00', '2026-11-27 00:00:00'),
        'Oct31_12_to_Nov26_12': ('2026-10-31 12:00:00', '2026-11-26 12:00:00')}
for name, P, eP, T0, eT0 in [('A', MS_PA, MS_ePA, MS_T0A + 2457000, MS_eT0A), ('B', MS_PB, MS_ePB, MS_T0B + 2457000, MS_eT0B)]:
    rows = transits(P, eP, T0, eT0, '2026-10-30 00:00:00', '2026-11-28 00:00:00')
    out[f'{name}_transits_Oct30_Nov28'] = rows
    for wn, (a, b) in wins.items():
        cnt = sum(1 for r in rows if a <= r['utc'] <= b)
        out.setdefault('counts', {}).setdefault(wn, {})[name] = cnt
    for r in rows:
        print(name, r['n'], r['utc'], 'BTJD %.4f' % r['btjd'], '+- %.1f min' % r['sigma_min'])
print(json.dumps(out['counts'], indent=1))
# also the own ephemeris (global emcee/LSQ) for A
lsq = json.load(open('out_s05_globalfit_lsq.json'))['raw_lin']['params']
rowsA2 = transits(lsq['P'][0], lsq['P'][1], lsq['T0'][0] + 2457000, lsq['T0'][1], '2026-10-30 00:00:00', '2026-11-28 00:00:00')
out['A_transits_own_global_fit'] = rowsA2
out['A_own_minus_ms_S110_min'] = [(r2['bjd_tdb'] - r1['bjd_tdb']) * 1440 for r1, r2 in zip(out['A_transits_Oct30_Nov28'], rowsA2)]
print('own-ms offsets (min):', np.round(out['A_own_minus_ms_S110_min'], 1))
json.dump(out, open('out_s13_s110.json', 'w'), indent=1, default=float)
