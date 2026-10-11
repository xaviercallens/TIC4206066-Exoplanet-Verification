"""G1: verify the reading of the official MAST PRF files and of the TESS_PRF package.

1. Each PRF FITS is 117x117 = 9x9 interleaved 13x13 sub-PRFs (MAST archive manual; 00README).
   Check: total = 81, every sub-PRF sums to ~1.
2. Own un-interleaving (derived from CRPIX1P = 59, CDELT1P = 1/9): sub-PRF for a source at intra-pixel
   offset s/9 px (s=-4..4) is prf[(4-s)::9] along each axis.  Check the centroid moves by +s/9.
3. TESS_PRF.locate(): centroid of the model at (x0+dx, y0+dy) must move by (dx, dy) (to ~0.02 px).
4. Interpolation across the 5x5 grid of PRF files: TESS_PRF's bilinear result == own bilinear result.
5. Sensitivity to the 44-column collateral offset quoted in the README.
"""
import glob, os
import numpy as np
from astropy.io import fits
from common import *

res = {}
files = sorted(glob.glob(os.path.join(PRFDIR, "cam1_ccd2", "*.fits")))
raw = fits.getdata(files[12]); hdr = fits.getheader(files[12])
res["file"] = os.path.basename(files[12])
res["shape"] = raw.shape
res["total"] = float(raw.sum())
res["hdr"] = {k: hdr.get(k) for k in ["CCD_RREF", "CCD_CREF", "NSAMP", "CRPIX1P", "CRVAL1P", "CDELT1P", "VERSION"]}


def own_sub(prf, sx, sy):
    return prf[(4 - sy)::9, (4 - sx)::9]


def centroid(img):
    yy, xx = np.indices(img.shape)
    return (img * xx).sum() / img.sum(), (img * yy).sum() / img.sum()


sums, cents = [], []
for sy in range(-4, 5):
    for sx in range(-4, 5):
        s = own_sub(raw, sx, sy)
        assert s.shape == (13, 13)
        sums.append(s.sum())
        cx, cy = centroid(s)
        cents.append((sx / 9, sy / 9, cx, cy))
cents = np.array(cents)
res["subprf_sum_min_max"] = [float(np.min(sums)), float(np.max(sums))]
# centroid slope vs intended offset
px = np.polyfit(cents[:, 0], cents[:, 2], 1)
py = np.polyfit(cents[:, 1], cents[:, 3], 1)
res["own_unpack_centroid_slope_x"] = px.tolist()
res["own_unpack_centroid_slope_y"] = py.tolist()
res["own_unpack_centroid_resid_max_px"] = float(max(np.abs(cents[:, 2] - np.polyval(px, cents[:, 0])).max(),
                                                   np.abs(cents[:, 3] - np.polyval(py, cents[:, 1])).max()))

# 3. TESS_PRF.locate
out3 = {}
for key, (col, row) in {"S6_cam1ccd1": (76.0, 889.0), "S32_cam1ccd2": (1372.0, 943.0), "S98_cam1ccd2": (676.0, 982.0)}.items():
    cam, ccd = 1, int(key[-1])
    m = PRFModel(cam, ccd, 10, col, row)
    shape = (25, 25)
    c0 = centroid(m.image(12.0, 12.0, shape))
    devs = []
    for dx in np.linspace(-0.5, 0.5, 11):
        for dy in np.linspace(-0.5, 0.5, 11):
            c = centroid(m.image(12.0 + dx, 12.0 + dy, shape))
            devs.append((dx, dy, c[0] - c0[0] - dx, c[1] - c0[1] - dy))
    devs = np.array(devs)
    # also large shifts (integer + fraction)
    devs2 = []
    for dx, dy in [(1.3, -2.7), (-3.45, 0.2), (2.0, 2.0), (-0.51, 0.49)]:
        c = centroid(m.image(12.0 + dx, 12.0 + dy, shape))
        devs2.append((dx, dy, c[0] - c0[0] - dx, c[1] - c0[1] - dy))
    # 4. own bilinear interpolation of the 4 surrounding files
    fl = sorted(glob.glob(os.path.join(PRFDIR, f"cam{cam}_ccd{ccd}", "*.fits")))
    rows = np.array([int(os.path.basename(f).split("-row")[1][:4]) for f in fl])
    cols = np.array([int(os.path.basename(f).split("-col")[1][:4]) for f in fl])
    ur, uc = np.unique(rows), np.unique(cols)
    r1, r2 = ur[ur < row].max(), ur[ur > row].min()
    c1, c2 = uc[uc < col].max(), uc[uc > col].min()
    get = lambda r, c: fits.getdata(fl[np.where((rows == r) & (cols == c))[0][0]])
    wr, wc = (row - r1) / (r2 - r1), (col - c1) / (c2 - c1)
    mine = (get(r1, c1) * (1 - wr) * (1 - wc) + get(r1, c2) * (1 - wr) * wc + get(r2, c1) * wr * (1 - wc)
            + get(r2, c2) * wr * wc)
    grid_diff = float(np.abs(mine - m.p.prf).max() / mine.max())
    # PRF variation between the grid corners (how much the interpolation matters)
    corner_var = float(np.abs(get(r1, c1) - get(r2, c2)).max() / mine.max())
    # 5. 44-column sensitivity: centroid of the central sub-PRF at col vs col+44
    try:
        m44 = PRFModel(cam, ccd, 10, col + 44.0, row)
        c44 = centroid(m44.image(12.0, 12.0, shape))
        d44 = [c44[0] - c0[0], c44[1] - c0[1]]
        img_a = m.image(12.0, 12.0, shape); img_b = m44.image(12.0, 12.0, shape)
        d44_img = float(np.abs(img_a - img_b).max() / img_a.max())
    except Exception as e:
        d44, d44_img = str(e), None
    out3[key] = dict(colnum=col, rownum=row, grid_corners=[int(r1), int(r2), int(c1), int(c2)],
                     max_abs_centroid_dev_px=float(np.abs(devs[:, 2:]).max()),
                     rms_centroid_dev_px=float(np.sqrt((devs[:, 2:] ** 2).mean())),
                     large_shift_devs=np.array(devs2).tolist(),
                     tessprf_vs_own_bilinear_maxreldiff=grid_diff,
                     corner_prf_maxreldiff=corner_var,
                     centroid_change_if_col_plus44_px=d44, image_maxreldiff_col_plus44=d44_img)
res["locate_tests"] = out3
res["which_set"] = ("start_s0004 for Sectors >= 4 (S6, S32, S98): MAST TESS Archive Manual 'Data Product Overview' "
                    "https://outerspace.stsci.edu/spaces/TESS/pages/14563420 and start_s0004/00README.txt; "
                    "VERSION='UPDATED_2.0' in the FITS header")
dump(res, "g1_prf_check.json")
print(json.dumps(jsonable(res), indent=1))
