"""Independent check of the S110 barycentric light-time and UTC times for B (verifier T5).
Roemer delay = (r_earth_bary . n_hat)/c computed from astropy get_body_barycentric (not Time.light_travel_time)."""
import json, numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric
import astropy.units as u, astropy.constants as c
OUT = "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T5_signal_B/verify/v_ltt_check.json"
tgt = SkyCoord(80.85983*u.deg, -1.32696*u.deg)
n = tgt.cartesian.xyz.value  # ICRS unit vector
P, T0 = 11.13274, 2461000.865
res = {}
for E in (31, 32, 33):
    bjd = T0 + E*P
    t = Time(bjd, format="jd", scale="tdb")
    for _ in range(4):
        r = get_body_barycentric("earth", t).xyz.to(u.m).value
        dt = np.dot(r, n)/c.c.value  # s ; BJD = t_obs + dt
        t = Time(bjd, format="jd", scale="tdb") - dt*u.s
    res[E] = dict(bjd_tdb=bjd, roemer_min=dt/60, utc_obs=t.utc.iso, utc_no_ltt=Time(bjd, format="jd", scale="tdb").utc.iso)
# sanity: time system of an actual TESS S98 file is irrelevant here; also 4-event fit and 7-event fit times
b = dict(fit4=(4000.8655179850807, 11.132739168890842), fit7=(4000.866048942678, 11.132750610389188))
for k,(t0,p) in b.items():
    res[k] = {E: Time(2457000+t0+E*p, format="jd", scale="tdb").utc.iso for E in (31,32,33)}
json.dump(res, open(OUT, "w"), indent=1); print(json.dumps(res, indent=1))
