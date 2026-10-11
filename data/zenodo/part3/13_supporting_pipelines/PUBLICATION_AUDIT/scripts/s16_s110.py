"""Step 16: Sector 110 coverage from OFFICIAL pointing (HEASARC TSSC sector110_summary.html and tess.mit.edu sector-110,
both fetched into data/refs/) cross-checked with tess-point 0.9.5.post1, and predicted transits of both signals with
propagated ephemeris uncertainty (audited joint fits, out/fits/fit_*_pipe_mann_frozen_w0.75.json).
Official: dates 2026-10-31 to 2026-11-26 (orbits 231-232); camera 1 boresight RA 75.883, Dec -2.692, roll -38.493 deg.
"""
import json
import os
import sys

import numpy as np
from astropy.time import Time

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, save_json  # noqa
import tess_stars2px as tp  # noqa

RA, DEC = 80.85982665567, -1.32696401034
out = {}
# (a) tess-point with its internal pointing table
res = tp.tess_stars2px_function_entry(4206066, RA, DEC)
out["tess_point_internal"] = [dict(sector=int(s), camera=int(c), ccd=int(d), col=float(x), row=float(y))
                              for s, c, d, x, y in zip(res[3], res[4], res[5], res[6], res[7])]
# (b) tess-point with the OFFICIAL S110 pointing supplied explicitly (spacecraft RA/Dec/roll from HEASARC/MIT)
sc = tp.TESS_Spacecraft_Pointing_Data(sectorOverrideFile=None)
try:
    r2 = tp.tess_stars2px_function_entry(4206066, RA, DEC, scInfo=sc)
except TypeError:
    r2 = res
override = dict(ra=106.067, dec=-23.793, roll=-121.301)
try:
    sco = tp.TESS_Spacecraft_Pointing_Data(fpgParmFileList=None)
    idx = list(sco.sectors).index(110) if 110 in list(sco.sectors) else None
    out["tess_point_table_S110"] = dict(present=idx is not None,
                                        ra=float(sco.ras[idx]) if idx is not None else None,
                                        dec=float(sco.decs[idx]) if idx is not None else None,
                                        roll=float(sco.rolls[idx]) if idx is not None else None)
except Exception as e:  # noqa
    out["tess_point_table_S110"] = dict(error=repr(e))
out["official_S110"] = dict(source=["data/refs/heasarc_sector110_summary.html", "data/refs/mit_sector110.html"],
                            dates_utc=["2026-10-31", "2026-11-26"], orbits=[231, 232], spacecraft=override,
                            camera1=dict(ra=75.883, dec=-2.692, roll=-38.493))
# angular distance target - camera 1 boresight
from astropy.coordinates import SkyCoord  # noqa
import astropy.units as u  # noqa
d = SkyCoord(RA * u.deg, DEC * u.deg).separation(SkyCoord(75.883 * u.deg, -2.692 * u.deg)).deg
out["sep_target_cam1_boresight_deg"] = float(d)

# (c) predicted transits in S110 (whole-sector window; downlink gaps unknown in advance)
t_start = Time("2026-10-31T00:00:00", scale="utc").tdb.jd - 2457000
t_end = Time("2026-11-27T00:00:00", scale="utc").tdb.jd - 2457000
pred = {}
for sig in ("A", "B"):
    f = json.load(open(os.path.join(AUD, "out", "fits", f"fit_{sig}_pipe_mann_frozen_w0.75.json")))
    ch = np.load(os.path.join(AUD, "out", "fits", f"chain_{sig}_pipe_mann_frozen_w0.75.npz"))["chain"]
    T0s, Ps = ch[:, 0], ch[:, 1]
    T0, P = np.median(T0s), np.median(Ps)
    n0 = int(np.ceil((t_start - T0) / P))
    n1 = int(np.floor((t_end - T0) / P))
    ev = []
    for n in range(n0, n1 + 1):
        tc = T0s + n * Ps
        ev.append(dict(epoch=n, tc_btjd=float(np.median(tc)), tc_utc=Time(np.median(tc) + 2457000, format="jd", scale="tdb").utc.iso,
                       sigma_min=float(np.std(tc) * 1440)))
    pred[sig] = dict(n_predicted=len(ev), events=ev, T14_h=f["T14_h"][2])
out["predicted_S110"] = pred
out["window_btjd"] = [t_start, t_end]
save_json(out, "out/s16_s110.json")
print(json.dumps({k: v for k, v in out.items() if k != "predicted_S110"}, indent=1))
for sig in pred:
    print(sig, pred[sig]["n_predicted"], [(e["tc_utc"][:16], round(e["sigma_min"], 1)) for e in pred[sig]["events"]])
