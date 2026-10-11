"""Independent check of S110 observed times: Roemer delay from JPL ephemeris (astropy get_body_barycentric,
builtin) r_earth.n/c, no Time.light_travel_time. Target: Gaia DR3 2016 position with proper motion ignored
(effect < 1e-3 s). Writes s110_indep.json."""
import json, numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, SkyCoord
import astropy.units as u, astropy.constants as c
n_hat = SkyCoord(80.85983*u.deg, -1.32696*u.deg).cartesian.xyz.value  # ICRS unit vector
def obs(bjd):
    t = Time(bjd, format='jd', scale='tdb')
    for _ in range(5):
        r = get_body_barycentric('earth', t).xyz.to(u.m).value
        d = float(np.dot(r, n_hat)/c.c.value)   # seconds; BJD = t_obs + d
        t = Time(bjd, format='jd', scale='tdb') - d*u.s
    return t.utc.iso[:19], d/60
E = {'A_U':(3.182785,2460998.736,range(109,118)), 'A_pub':(3.1827856,2460998.7375,range(109,118)),
     'A_audit':(3.182785293657263,2457000+3998.7367930333035,range(109,118)),
     'B_U':(11.13274,2461000.865,range(31,34)), 'B_audit':(11.132745380039134,2457000+4000.8658739613074,range(31,34)),
     'B_pub':(11.132738735768573,2457000+4000.8650000483804,range(31,34))}
gap=(Time('2026-11-13 19:55',scale='utc'),Time('2026-11-14 00:55',scale='utc'))
out={}
for k,(P,T0,ns) in E.items():
    rows=[]
    for n in ns:
        s,dm=obs(T0+n*P); rows.append((n,s,round(dm,2)))
    out[k]=rows; print(k, rows)
tA=Time(obs(2460998.736+113*3.182785)[0],scale='utc'); h=1.87/24/2
print('gap margins A n=113 (U): ingress-gapstart [min]', ((tA-h*u.day)-gap[0]).to(u.min).value, ' gapend-egress [min]', (gap[1]-(tA+h*u.day)).to(u.min).value)
out['gap_margin_A113_U_min']=[((tA-h*u.day)-gap[0]).to(u.min).value,(gap[1]-(tA+h*u.day)).to(u.min).value]
json.dump(out,open('s110_indep.json','w'),indent=1)
