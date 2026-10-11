"""Step 2: fresh TESScut sector list + 21x21 FFI cutouts (S6, S32, S98), compared with the earlier cutouts
array-by-array (TESScut files carry creation timestamps, so byte identity is not expected)."""
import glob
import io
import json
import os
import sys
import zipfile

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, OLD_DEEP, fetch, save_json, utcnow  # noqa

RA, DEC = 80.85982665567, -1.32696401034
rec = fetch(f"https://mast.stsci.edu/tesscut/api/v0.1/sector?ra={RA}&dec={DEC}&radius=0", "data/tesscut/sectors.json",
            note="TESScut sector list")
sectors = json.load(open(os.path.join(AUD, "data/tesscut/sectors.json")))
print(json.dumps(sectors, indent=1))

out = {"sector_list": sectors, "cutouts": []}
for s in [6, 32, 98]:
    url = f"https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={RA}&dec={DEC}&y=21&x=21&sector={s}"
    z = f"data/tesscut/tesscut_s{s}.zip"
    if not os.path.exists(os.path.join(AUD, z)):
        fetch(url, z, note=f"TESScut 21x21 sector {s}", timeout=900)
    with zipfile.ZipFile(os.path.join(AUD, z)) as zf:
        names = zf.namelist()
        zf.extractall(os.path.join(AUD, "data/tesscut"))
    fn = os.path.join(AUD, "data/tesscut", names[0])
    old = glob.glob(os.path.join(OLD_DEEP, "pixels", f"tess-s{s:04d}-*_astrocut.fits"))
    with fits.open(fn) as h:
        new_flux = h[1].data["FLUX"]
        new_time = h[1].data["TIME"]
        hdr0, hdr1, hdr2 = h[0].header, h[1].header, h[2].header
        info = dict(sector=s, file=os.path.relpath(fn, AUD), shape=list(new_flux.shape), camera=hdr0.get("CAMERA"),
                    ccd=hdr0.get("CCD"), procver=hdr0.get("PROCVER"), date=hdr0.get("DATE"),
                    ffi_type=hdr0.get("FFI_TYPE"), n_cad=len(new_time),
                    tstart=float(np.nanmin(new_time)), tstop=float(np.nanmax(new_time)),
                    timedel=hdr1.get("TIMEDEL"), exposure_d=hdr1.get("EXPOSURE"), int_time=hdr1.get("INT_TIME"),
                    frametim=hdr1.get("FRAMETIM"), num_frm=hdr1.get("NUM_FRM"), timesys=hdr1.get("TIMESYS"),
                    bjdrefi=hdr1.get("BJDREFI"), bjdreff=hdr1.get("BJDREFF"), tunit=hdr1.get("TIMEUNIT"),
                    backapp=hdr1.get("BACKAPP"), wcs_1crv=hdr2.get("1CRV4P", None))
    if old:
        with fits.open(old[0]) as h:
            of = h[1].data["FLUX"]
            ot = h[1].data["TIME"]
        info["old_file"] = os.path.relpath(old[0], AUD)
        info["flux_identical"] = bool(np.array_equal(np.nan_to_num(of, nan=-9e9), np.nan_to_num(new_flux, nan=-9e9)))
        info["time_identical"] = bool(np.array_equal(np.nan_to_num(ot, nan=-9e9), np.nan_to_num(new_time, nan=-9e9)))
    out["cutouts"].append(info)
    print(json.dumps(info))
out["utc"] = utcnow()
save_json(out, "data/tesscut/tesscut_check.json")
