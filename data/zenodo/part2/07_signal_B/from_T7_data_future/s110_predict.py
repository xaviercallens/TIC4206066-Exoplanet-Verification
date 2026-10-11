"""T7 minimal S110 rerun: transit mid-times of A and B inside TESS Sector 110 for every candidate ephemeris.
Inputs: ephemerides quoted from U / pubprep / audit (paths in EPH below); official orbit-segment table
web/TESS_orbit_times.csv (tess.mit.edu, fetched 2026-09-25T15:30Z, file 'Created 2026-09-25 by roland').
Two conversions of BJD_TDB to a calendar time are reported:
  utc_naive = BJD_TDB re-expressed on the UTC scale with NO light-travel correction (what audit/pubprep printed)
  utc_obs   = UTC at which the event is seen at the geocentre (TESS within ~1.3 light-s): t + ltt_bary(t) = BJD_TDB.
Output: s110_predict.json (+ stdout log)."""
import json, csv, numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

target = SkyCoord(80.85983 * u.deg, -1.32696 * u.deg)   # Gaia DR3 J2016 position (U/00_data/README.md)
geo = EarthLocation.from_geocentric(0 * u.m, 0 * u.m, 0 * u.m)

def obs_utc(bjd):
    t = Time(bjd, format='jd', scale='tdb', location=geo)
    for _ in range(4):
        ltt = t.light_travel_time(target, kind='barycentric')
        t = Time(bjd, format='jd', scale='tdb', location=geo) - ltt
    return t.utc

EPH = {
 # name: (P, eP, T0_BJD, eT0, source)
 'A_U':       (3.182785, 0.000006, 2460998.736, 0.003, 'U/proposed_revision/manuscript.tex Sect.2'),
 'A_pubprep': (3.1827856, 0.0000058, 2460998.7375, 0.0023, 'U/01_ephemeris_A/a_common.py MS_PA/MS_T0A = pubprep manuscript'),
 'A_audit':   (3.182785293657263, None, 2457000 + 3998.7367930333035, None, 'audit/out/ephem_current.json A (sigma from audit/out/s16_s110.json)'),
 'B_U':       (11.13274, 0.00003, 2461000.865, 0.003, 'U/proposed_revision/manuscript.tex Sect.3'),
 'B_audit':   (11.132745380039134, None, 2457000 + 4000.8658739613074, None, 'audit/out/ephem_current.json B (sigma from audit/out/s16_s110.json)'),
 'B_pubprep': (11.132738735768573, 3.068946752541719e-05, 2457000 + 4000.8650000483804, 0.0030752564161342496,
               'pubprep/analysis/out/a1_results.json signals.B.B_mann__pipe.ephemeris.e_dt_inj_min (scaled errors)'),
}
T14 = {'A': 1.87 / 24, 'B': 1.66 / 24}   # U manuscript T14_A; audit s16 T14_B = 1.658 h

segs = []
for r in csv.DictReader(l for l in open('web/TESS_orbit_times.csv') if not l.startswith('#')):
    if r['Sector'].strip() == '110':
        segs.append((Time(r['Start of Orbit'], scale='utc'), Time(r['End of Orbit'], scale='utc'), int(r['Orbit'])))
print('S110 segments listed:', [(a.iso, b.iso, o) for a, b, o in segs])
listed_end = segs[-1][1]
sector_start = segs[0][0]
# unpublished last segment: assume the 5.0-h gap seen at every other S109/S110 downlink
last_seg_start_assumed = listed_end + 5 * u.hour

out = {'segments_listed': [dict(start=a.iso, end=b.iso, orbit=o) for a, b, o in segs],
       'sector_start_utc': sector_start.iso, 'listed_coverage_end_utc': listed_end.iso,
       'last_segment_start_assumed_utc': last_seg_start_assumed.iso, 'ephemerides': {}}
for name, (P, eP, T0, eT0, src) in EPH.items():
    sig = name[0]
    rows = []
    n0 = int(np.floor((Time('2026-10-30', scale='utc').tdb.jd - T0) / P))
    n1 = int(np.ceil((Time('2026-11-28', scale='utc').tdb.jd - T0) / P))
    for n in range(n0, n1 + 1):
        bjd = T0 + n * P
        naive = Time(bjd, format='jd', scale='tdb').utc
        ob = obs_utc(bjd)
        s_min = None if eP is None else float(np.sqrt(eT0 ** 2 + (n * eP) ** 2) * 1440)
        half = T14[sig] / 2
        ing, egr = ob - half * u.day, ob + half * u.day
        where = 'before_sector' if egr < sector_start else None
        for a, b, o in segs:
            if ing >= a and egr <= b: where = f'in_listed_segment_orbit{o}'
            elif (ing < a < egr) or (ing < b < egr): where = f'partial_edge_orbit{o}'
        if where is None:
            gaps = [(segs[i][1], segs[i + 1][0]) for i in range(len(segs) - 1)]
            for a, b in gaps:
                if a <= ob <= b: where = f'in_gap_{a.iso[:16]}_{b.iso[:16]}'
        if where is None:
            if listed_end < ob < last_seg_start_assumed: where = 'in_assumed_gap_after_listed_end'
            elif ob >= last_seg_start_assumed: where = 'after_listed_coverage (final segment unpublished)'
        rows.append(dict(n=n, bjd_tdb=round(bjd, 6), btjd=round(bjd - 2457000, 6), utc_naive=naive.iso[:19],
                         utc_obs=ob.iso[:19], obs_minus_naive_min=round((ob - naive).to(u.min).value, 2),
                         sigma_min=None if s_min is None else round(s_min, 2),
                         ingress_obs=ing.iso[:16], egress_obs=egr.iso[:16], status=where))
    cnt = {}
    for lab, lo, hi in [('Oct31_00_Nov26_00_obs', '2026-10-31 00:00', '2026-11-26 00:00'),
                        ('Oct31_00_Nov27_00_obs', '2026-10-31 00:00', '2026-11-27 00:00'),
                        ('Oct31_2355_Nov26_24_obs', '2026-10-31 23:55', '2026-11-27 00:00')]:
        cnt[lab] = sum(lo <= r['utc_obs'] <= hi for r in rows)
    cnt['in_listed_segments_full'] = sum(r['status'].startswith('in_listed') for r in rows)
    cnt['in_listed_gaps'] = sum(r['status'].startswith('in_gap') for r in rows)
    out['ephemerides'][name] = dict(P=P, eP=eP, T0_bjd=T0, eT0=eT0, source=src, events=rows, counts=cnt)
    print('\n==', name, src)
    for r in rows: print(r['n'], r['utc_obs'], '(naive', r['utc_naive'] + ')', r['sigma_min'], r['status'])
    print(cnt)
json.dump(out, open('s110_predict.json', 'w'), indent=1)
