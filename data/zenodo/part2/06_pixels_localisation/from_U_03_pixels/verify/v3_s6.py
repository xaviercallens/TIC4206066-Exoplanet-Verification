"""Verifier: S6 3.18 d raw and dump-excised positions with own code (v1_loc functions)."""
import json, numpy as np, v1_loc as V
from astropy.io import fits
out = {}
r = V.run("S6_ffi", "A", False, nboot=50); out["S6_A_raw"] = r; print("S6 A raw", r["n"], np.round(r["off"],1), np.round(r["boot_sd"],1), [tuple(np.round(j,1)) for j in r["jackknife"]])
# dump-excised +-0.5 h
d = fits.getdata(V.F["S6_ffi"], 1); t = d["TIME"]; q = d["QUALITY"]
dumps = t[(q & 32) > 0]
ds = V.load("S6_ffi"); keep = np.ones(len(ds["t"]), bool)
for td in dumps: keep &= np.abs(ds["t"] - td) > 0.5/24
for k in ("t","F","pc1","pc2"): ds[k] = ds[k][keep]
f = V.Fitter(ds); evs = V.events(ds, "A"); Dm, Vm = V.combine(ds, evs, "A", False); x, y, a = f.fit(Dm, Vm)
out["S6_A_dumpx0.5h"] = dict(n=len(evs), off=f.to_arcsec(x, y)); print("S6 A dumpx0.5h", len(evs), np.round(f.to_arcsec(x, y),1))
# distance of each A transit to nearest dump
P, T0 = 3.1827856, 3998.7375
n = np.arange(np.ceil((np.nanmin(t)-T0)/P), np.floor((np.nanmax(t)-T0)/P)+1); tc = T0+n*P
dd = [float(np.min(np.abs(dumps - x))*24) for x in tc]; out["dist_h"] = list(zip(tc.tolist(), dd)); print("dump distances h", np.round(dd,2))
json.dump(out, open("v3_s6.json","w"), indent=1, default=float)
