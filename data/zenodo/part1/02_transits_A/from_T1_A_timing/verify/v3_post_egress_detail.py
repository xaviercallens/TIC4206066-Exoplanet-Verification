"""Verifier: post-egress data layout (minutes after egress, QUALITY==0 cadences) for S98 E-1 and E1 (read-only raw QLP LC)."""
import numpy as np, json
from astropy.io import fits
f = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"
d = fits.open(f)[1].data; t = np.array(d["TIME"], float); q = np.array(d["QUALITY"]).astype(np.int64)
PA, T0A, T14 = 3.1827850, 3998.73599, 1.887/24
out = {}
for n in (-1, 1):
    eg = T0A + n*PA + T14/2
    m = (q == 0) & (t > eg) & (t < eg + 0.3)
    x = np.round((t[m] - eg)*1440, 1)
    # segments separated by > 10 min
    seg = np.split(x, np.where(np.diff(x) > 10)[0] + 1)
    out[f"E{n}"] = [[float(s[0]), float(s[-1]), int(len(s))] for s in seg if len(s)]
    print(f"E{n}: post-egress QUALITY==0 segments [start_min, end_min, n]:", out[f"E{n}"])
json.dump(out, open(__file__.replace(".py", ".json"), "w"), indent=1)
