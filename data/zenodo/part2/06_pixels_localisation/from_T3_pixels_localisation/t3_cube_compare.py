"""Check whether the audit's and U's S98 TESScut cubes contain identical pixel data (read-only, memmap)."""
import numpy as np, json
from astropy.io import fits
R = "/Users/paulsmith/space <3/tess-jev/addenda/"
a = R + "tic4206066_PUBLICATION_AUDIT/data/tesscut/tess-s0098-1-2_80.859827_-1.326964_21x21_astrocut.fits"
u = R + "opus55_final_audit_tic4206066/00_data/tesscut/21x21/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits"
out = {}
with fits.open(a, memmap=True) as ha, fits.open(u, memmap=True) as hu:
    for k in ("TIME", "QUALITY"):
        x, y = np.asarray(ha[1].data[k]), np.asarray(hu[1].data[k])
        out[k + "_equal"] = bool(x.shape == y.shape and np.array_equal(x, y, equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x, y))
    fa, fu = ha[1].data["FLUX"], hu[1].data["FLUX"]
    out["flux_shape"] = [list(fa.shape), list(fu.shape)]
    n = fa.shape[0]; eq = True; maxdiff = 0.0
    for i in range(0, n, 2000):
        A = np.asarray(fa[i:i + 2000]); B = np.asarray(fu[i:i + 2000])
        m = np.isfinite(A) & np.isfinite(B)
        if not np.array_equal(np.isfinite(A), np.isfinite(B)): eq = False
        if m.any(): maxdiff = max(maxdiff, float(np.max(np.abs(A[m] - B[m]))))
    out["flux_max_abs_diff"] = maxdiff; out["flux_nan_pattern_equal"] = eq
    for k in ("1CRV4P", "2CRV4P"):
        out[k] = [ha[1].header.get(k), hu[1].header.get(k)]
    out["DATE"] = [ha[0].header.get("DATE"), hu[0].header.get("DATE")]
print(json.dumps(out, indent=1))
json.dump(out, open("/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation/t3_cube_compare.json", "w"), indent=1)
