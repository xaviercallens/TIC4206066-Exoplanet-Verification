"""Size of the barycentric light-time term at the S110 B epochs (explains C093 minute-level differences)."""
import json
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
target = SkyCoord(80.85983 * u.deg, -1.32696 * u.deg); geo = EarthLocation.from_geocentric(0 * u.m, 0 * u.m, 0 * u.m)
out = {}
for n in (31, 32, 33):
    bjd = 2461000.865 + n * 11.13274
    naive = Time(bjd, format="jd", scale="tdb").utc.iso        # TDB->UTC only, no light-time
    t = Time(bjd, format="jd", scale="tdb", location=geo)
    ltt = t.light_travel_time(target, kind="barycentric").to(u.min).value
    out[n] = dict(bjd_tdb=bjd, utc_no_lighttime=naive, barycentric_ltt_min=ltt)
    # audit ephemeris (T0 4000.8658739613, P 11.132745380039134) without light-time, for C093 comparison
    bjd_a = 2461000.8658739613 + n * 11.132745380039134
    out[n]["audit_eph_utc_no_lighttime"] = Time(bjd_a, format="jd", scale="tdb").utc.iso
    out[n]["pubprep_exact_eph_utc_no_lighttime"] = Time(2457000 + 4000.865 + n * 11.132738735768573, format="jd", scale="tdb").utc.iso
json.dump(out, open("t5_ltt.json", "w"), indent=1); print(json.dumps(out, indent=1))
