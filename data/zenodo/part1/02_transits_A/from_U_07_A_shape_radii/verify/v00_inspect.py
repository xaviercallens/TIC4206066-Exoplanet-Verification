"""Verifier: inspect raw files (headers, quality bits, S98 segments around Signal-A epochs). Own code."""
import os, hashlib, json
import numpy as np
from astropy.io import fits
D = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "00_data", "mast")
F = {6: "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits",
     32: "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits",
     98: "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"}
PA, T0A = 3.1827856, 2460998.7375 - 2457000
PB, T0B = 11.13274, 2461000.865 - 2457000
for s, fn in F.items():
    p = os.path.join(D, fn)
    h = fits.open(p)
    print(s, fn, hashlib.sha256(open(p, 'rb').read()).hexdigest()[:16])
    h0, h1 = h[0].header, h[1].header
    for kw in ("TIMESYS", "BJDREFI", "BJDREFF", "TIMEDEL", "CROWDSAP", "FLFRCSAP", "DILUTION", "TESSMAG", "TEFF", "RADIUS"):
        for hh in (h0, h1):
            if kw in hh:
                print("   ", kw, hh[kw]); break
    print("   cols", h[1].columns.names)
    q = np.asarray(h[1].data["QUALITY"]).astype(np.int64)
    bits = {b: int(((q >> b) & 1).sum()) for b in range(32) if ((q >> b) & 1).any()}
    print("   n", len(q), "q==0", int((q == 0).sum()), "bits", bits)
