"""Record inputs: product files, sha256, FITS product identifiers, ephemeris, software versions."""
import json, os, sys, platform
import numpy, astropy, scipy, matplotlib
from astropy.io import fits
from common import ROOT, OUT, PKG, DATA, EPHEM_FILE, EPH, sha256
rec = {"ephemeris_file": EPHEM_FILE, "ephemeris_sha256": sha256(os.path.join(ROOT, EPHEM_FILE)),
       "ephemeris_BTJD": {s: {"P_d": v[0], "T0_BTJD": v[1], "T14_d": v[2]} for s, v in EPH.items()},
       "inputs": {}}
for k, p in DATA.items():
    h = fits.open(os.path.join(ROOT, p)); h0, h1 = h[0].header, h[1].header
    rec["inputs"][k] = {"path": p, "sha256": sha256(os.path.join(ROOT, p)),
        "sector": h0.get("SECTOR"), "camera": h0.get("CAMERA"), "ccd": h0.get("CCD"),
        "object": h0.get("OBJECT"), "ticid": h0.get("TICID"), "creator": h0.get("CREATOR"),
        "procver": h0.get("PROCVER"), "data_rel": h0.get("DATA_REL"), "timedel_d": h1.get("TIMEDEL"),
        "n_cadences": int(h1.get("NAXIS2")), "tstart": h1.get("TSTART"), "tstop": h1.get("TSTOP")}
rec["momentum_dumps_source"] = "addenda/independent_verification_tic4206066_SEALED/out/t7_momentum_dumps.json"
rec["momentum_dumps_sha256"] = sha256(os.path.join(ROOT, rec["momentum_dumps_source"]))
rec["old_T9_source"] = "addenda/independent_verification_tic4206066_SEALED/out/t9_aperture_diff.json"
rec["old_T9_sha256"] = sha256(os.path.join(ROOT, rec["old_T9_source"]))
rec["t9_cutouts_11x11"] = {k: {"path": v, "sha256": sha256(os.path.join(ROOT, v))} for k, v in {
    "S98": "addenda/independent_verification_tic4206066_SEALED/downloads/tesscut/tess-s0098-1-2_80.859827_-1.326964_11x11_astrocut.fits",
    "S06": "addenda/independent_verification_tic4206066_SEALED/downloads/tesscut/tess-s0006-1-1_80.859827_-1.326964_11x11_astrocut.fits"}.items()}
from common import QLP_S98
rec["qlp_s98_for_scattered_light_flags"] = {"path": QLP_S98, "sha256": sha256(os.path.join(ROOT, QLP_S98))}
rec["software"] = {"python": sys.version.split()[0], "platform": platform.platform(), "numpy": numpy.__version__,
                   "astropy": astropy.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
                   "interpreter": sys.executable}
json.dump(rec, open(os.path.join(OUT, "inputs.json"), "w"), indent=1)
print(json.dumps(rec["software"]))
