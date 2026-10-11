import numpy as np, json
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
c=SkyCoord(80.85983,-1.32696,unit='deg')
def bjd_to_utc(bjd):
    # iterate: geocentric UTC such that UTC+ltt_bary = BJD_TDB
    t=Time(bjd,format='jd',scale='tdb')
    for _ in range(3):
        ltt=t.light_travel_time(c,kind='barycentric',location=EarthLocation.from_geocentric(0,0,0,unit='m'))
        t=Time(bjd,format='jd',scale='tdb')-ltt
    return t.utc
start=Time('2026-10-31T00:00:00',scale='utc'); end=Time('2026-11-27T00:00:00',scale='utc')
res={}
for name,T0,sT0,P,sP in [('A',2460998.7375,0.0023,3.1827856,0.0000058),('B',2461000.865,0.003,11.13274,0.00003)]:
    lo=int(np.floor((start.tdb.jd-T0)/P))-1; hi=int(np.ceil((end.tdb.jd-T0)/P))+1
    rows=[]
    for n in range(lo,hi+1):
        tc=T0+n*P; s=np.hypot(sT0,n*sP)
        u_=bjd_to_utc(tc)
        if start<=u_<=end: rows.append((n,u_.iso[:16],round(s*1440,1)))
    res[name]=rows; print(name,len(rows)); [print(' ',r) for r in rows]
json.dump(res,open('v02_s110.json','w'),indent=1)
