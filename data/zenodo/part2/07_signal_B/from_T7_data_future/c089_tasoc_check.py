# C089 minimal check: does verif T1 drop TASOC_S06 at the 'too few points after detrend' branch (which writes no result)?
import sys, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, "/Users/paulsmith/space <3/tess-jev/addenda/independent_verification_tic4206066_SEALED/scripts")
from lcio import load, fast_detrend
V = "/Users/paulsmith/space <3/tess-jev/addenda/independent_verification_tic4206066_SEALED/downloads/"
L = load("TASOC_S06", path=V + "lc_hlsp_tasoc_tess_ffi_tic00004206066-s0006-cam1-ccd1-c1800_tess_v05_ens-lc.fits")
print("TASOC_S06 good points:", len(L["t"]))
if len(L["t"]) >= 50:
    fd, keep = fast_detrend(L["t"], L["f"], window=1.0)
    print("kept after fast_detrend:", int(keep.sum()), "finite fd:", int(np.isfinite(fd).sum()))
