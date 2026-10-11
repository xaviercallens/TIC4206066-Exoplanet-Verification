# Read-only dump of full QLP FITS headers (S98, S32) for T7 reconciliation
from astropy.io import fits
import sys
U="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/"
for f in ["hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits","hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits"]:
    with fits.open(U+f, memmap=True) as h:
        for i,hd in enumerate(h):
            print("=====",f,"HDU",i)
            print(repr(hd.header))
