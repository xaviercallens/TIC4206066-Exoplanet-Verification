"""Shared helpers for the independent pixel-localisation audit of TIC 4206066 (audit 03_pixels).

Everything here is written from scratch for this audit.  Inputs are only the raw public files in
../00_data (TESScut FFI cutouts, SPOC/TESS-SPOC TPFs, QLP light curves, MAST PRF FITS files, Gaia DR3 and
TIC 8.2 cone searches).  No code or derived table of earlier analyses is imported or read here.
"""
import os, json, glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import PRF as TESSPRF   # TESS_PRF 0.1.3 (Bell & Higgins 2022, ascl:2207.008)

ROOT = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066"
DATA = os.path.join(ROOT, "00_data")
OUT = os.path.join(ROOT, "03_pixels")
PRFDIR = os.path.join(DATA, "prf", "start_s0004")      # Sector >= 4 PRF set (MAST archive manual)

TARGET_GAIA = 3220388198192519424

# Ephemerides quoted in the manuscript (BJD_TDB - 2457000 = BTJD)
EPH = {
    "A": dict(P=3.1827856, T0=2460998.7375 - 2457000.0, T14=1.87 / 24.0, depth=500e-6),
    "B": dict(P=11.13274, T0=2461000.865 - 2457000.0, T14=1.65 / 24.0, depth=1000e-6),
}

# SPOC quality bits recommended as "lesser quality" by the MAST archive manual
# (bits 1,2,3,4,5,6,8,10,13,15 -> values 1,2,4,8,16,32,128,512,4096,16384)
DEFAULT_BITMASK = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 512 | 4096 | 16384
STRAY_BIT = 2048

FILES = {
    "S6_ffi": os.path.join(DATA, "tesscut/21x21/tess-s0006-1-1_80.859830_-1.326960_21x21_astrocut.fits"),
    "S32_ffi": os.path.join(DATA, "tesscut/21x21/tess-s0032-1-2_80.859830_-1.326960_21x21_astrocut.fits"),
    "S98_ffi": os.path.join(DATA, "tesscut/21x21/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits"),
    "S32_tpf": os.path.join(DATA, "mast/tess2020324010417-s0032-0000000004206066-0200-s_tp.fits"),
    "S6_tspoc_tp": os.path.join(DATA, "mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits"),
    "S32_tspoc_tp": os.path.join(DATA, "mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_tp.fits"),
    "S6_tspoc_lc": os.path.join(DATA, "mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits"),
    "S32_spoc_lc": os.path.join(DATA, "mast/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits"),
    "S32_tspoc_lc": os.path.join(DATA, "mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits"),
    "QLP98": os.path.join(DATA, "mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"),
    "QLP32": os.path.join(DATA, "mast/hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits"),
    "QLP6": os.path.join(DATA, "mast/hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits"),
}


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float):
        return None if not np.isfinite(o) else o
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def dump(obj, name):
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(jsonable(obj), f, indent=1)


# ----------------------------------------------------------------------------------------------------------
# pixel data
# ----------------------------------------------------------------------------------------------------------
class Cube:
    """Pixel time series (TESScut cutout or SPOC TPF)."""

    def __init__(self, key):
        self.key = key
        path = FILES[key]
        self.path = path
        with fits.open(path) as h:
            h0, h1, h2 = h[0].header, h[1].header, h[2].header
            d = h[1].data
            self.sector = int(h0["SECTOR"])
            self.cam = int(h0["CAMERA"])
            self.ccd = int(h0["CCD"])
            self.time = np.array(d["TIME"], float)
            self.flux = np.array(d["FLUX"], np.float64)
            self.ferr = np.array(d["FLUX_ERR"], np.float64)
            self.quality = np.array(d["QUALITY"], np.int64)
            self.cadenceno = np.array(d["CADENCENO"], np.int64)
            self.poscorr1 = np.array(d["POS_CORR1"], float)
            self.poscorr2 = np.array(d["POS_CORR2"], float)
            self.aperture = np.array(h[2].data, np.int64)
            self.col0 = int(h1["1CRV4P"])   # CCD column of array column 0
            self.row0 = int(h1["2CRV4P"])   # CCD row of array row 0
            self.wcs = WCS(h2)
            self.timedel = float(h1["TIMEDEL"])
            self.hdr1 = h1
        self.ny, self.nx = self.flux.shape[1:]
        self.is_tpf = "tp.fits" in path

    def good(self, bitmask=DEFAULT_BITMASK):
        ok = np.isfinite(self.time) & ((self.quality & bitmask) == 0)
        ok &= np.all(np.isfinite(self.flux.reshape(len(self.time), -1)), axis=1)
        return ok

    def sky_to_pix(self, ra, dec):
        x, y = self.wcs.all_world2pix(np.atleast_1d(ra), np.atleast_1d(dec), 0)
        return x, y   # x = column index, y = row index (0-based, integer = pixel centre)

    def pix_to_sky(self, x, y):
        ra, dec = self.wcs.all_pix2world(np.atleast_1d(x), np.atleast_1d(y), 0)
        return ra, dec

    def mid_epoch_jyear(self):
        t = np.nanmedian(self.time) + 2457000.0
        return 2000.0 + (t - 2451545.0) / 365.25


# ----------------------------------------------------------------------------------------------------------
# PRF
# ----------------------------------------------------------------------------------------------------------
class PRFModel:
    """Thin wrapper around TESS_PRF (official MAST start_s0004 PRF files, 81 interleaved 13x13 sub-PRFs)."""

    def __init__(self, cam, ccd, sector, colnum, rownum):
        assert sector >= 4
        self.p = TESSPRF.TESS_PRF(cam, ccd, sector, colnum, rownum, localdatadir=PRFDIR)

    def image(self, x, y, shape):
        """Normalised PRF of a point source at 0-based (col=x, row=y) inside an array of `shape`."""
        # TESS_PRF places the 13x13 stamp; sources further than 6 px outside the stamp give zeros.
        if x < -7 or y < -7 or x > shape[1] + 6 or y > shape[0] + 6:
            return np.zeros(shape)
        return self.p.locate(x, y, shape)


def prf_for_cube(cube, xc=None, yc=None):
    if xc is None:
        xc, yc = (cube.nx - 1) / 2.0, (cube.ny - 1) / 2.0
    return PRFModel(cube.cam, cube.ccd, cube.sector, cube.col0 + xc, cube.row0 + yc)


# ----------------------------------------------------------------------------------------------------------
# catalogue
# ----------------------------------------------------------------------------------------------------------
def stassun_T(G, bprp):
    """Stassun et al. 2019 (AJ 158, 138) eqs. 1-2."""
    G = np.asarray(G, float); c = np.asarray(bprp, float)
    T = G - 0.00522555 * c**3 + 0.0891337 * c**2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(c), T, G - 0.430)


def load_gaia(radius_arcsec=240.0):
    import pandas as pd
    g = pd.read_csv(os.path.join(DATA, "catalogs/gaiadr3_gaia_source_cone240as.csv"))
    # GAIA ids must be read as strings: as float64 they lose precision (19-digit ids)
    t = pd.read_csv(os.path.join(DATA, "catalogs/tic82_cone300as.csv"), dtype={"GAIA": str})
    t = t[t["GAIA"].notna()].copy()
    t["GAIA"] = np.array([int(s) for s in t["GAIA"]], dtype=np.int64)
    g = g[g["dist_arcsec_j2016pos"] <= radius_arcsec].copy()
    # some Gaia ids are attached to >1 TIC entry (TIC splits / duplicates): keep, per Gaia id, the non-DUPLICATE
    # TIC entry closest on the sky to the Gaia DR3 position (TIC J2000 vs Gaia 2016.0; faint neighbours have
    # small proper motions).
    t = t[t["disposition"].fillna("") != "DUPLICATE"]
    gpos = g.set_index("source_id")[["ra", "dec"]]
    t = t[t["GAIA"].isin(gpos.index)].copy()
    gra = gpos.loc[t["GAIA"], "ra"].values; gde = gpos.loc[t["GAIA"], "dec"].values
    t["sep_gaia"] = np.hypot((t["ra"].values - gra) * np.cos(np.radians(gde)), t["dec"].values - gde) * 3600
    t = t.sort_values("sep_gaia").drop_duplicates("GAIA", keep="first")
    g = g.merge(t[["GAIA", "ID", "Tmag", "e_Tmag", "disposition"]], left_on="source_id", right_on="GAIA", how="left")
    g["T_stassun"] = stassun_T(g["phot_g_mean_mag"], g["bp_rp"])
    g["Tmag_used"] = np.where(np.isfinite(g["Tmag"]), g["Tmag"], g["T_stassun"])
    g["Tmag_src"] = np.where(np.isfinite(g["Tmag"]), "TIC8.2", "Stassun19")
    g = g[np.isfinite(g["Tmag_used"])]
    return g.reset_index(drop=True)


def propagate(g, jyear):
    """Linear proper-motion propagation from the Gaia DR3 epoch (2016.0) to `jyear`.
    Sources without a 5/6-parameter solution keep their 2016.0 position."""
    dt = jyear - 2016.0
    pmra = np.nan_to_num(np.asarray(g["pmra"], float))      # mas/yr, includes cos(dec)
    pmdec = np.nan_to_num(np.asarray(g["pmdec"], float))
    dec = np.asarray(g["dec"], float) + pmdec * dt / 3.6e6
    ra = np.asarray(g["ra"], float) + pmra * dt / 3.6e6 / np.cos(np.radians(dec))
    return ra, dec


def offsets_arcsec(ra0, dec0, ra, dec):
    """(East, North) offsets in arcsec of (ra,dec) from (ra0,dec0) (gnomonic, small angles)."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c0 = SkyCoord(ra0 * u.deg, dec0 * u.deg)
    c = SkyCoord(np.atleast_1d(ra) * u.deg, np.atleast_1d(dec) * u.deg)
    de, dn = c0.spherical_offsets_to(c)
    return de.to(u.arcsec).value, dn.to(u.arcsec).value


def pix_jacobian_arcsec(cube, x, y, h=0.5):
    """2x2 matrix J with (dE, dN)[arcsec] = J @ (dx, dy)[pix] at pixel (x, y)."""
    ra0, dec0 = cube.pix_to_sky(x, y)
    J = np.zeros((2, 2))
    for k, (dx, dy) in enumerate([(h, 0.0), (0.0, h)]):
        ra1, dec1 = cube.pix_to_sky(x + dx, y + dy)
        e, n = offsets_arcsec(ra0[0], dec0[0], ra1, dec1)
        J[0, k] = e[0] / h
        J[1, k] = n[0] / h
    return J


# ----------------------------------------------------------------------------------------------------------
# transit windows
# ----------------------------------------------------------------------------------------------------------
def transit_centres(eph, tmin, tmax):
    P, T0 = eph["P"], eph["T0"]
    n0 = np.ceil((tmin - T0) / P)
    n1 = np.floor((tmax - T0) / P)
    n = np.arange(n0, n1 + 1)
    return n.astype(int), T0 + n * P


def box_fraction(t, dt, tc, T14):
    """Fraction of each exposure [t-dt/2, t+dt/2] inside the box [tc-T14/2, tc+T14/2]."""
    lo = np.maximum(t - dt / 2, tc - T14 / 2)
    hi = np.minimum(t + dt / 2, tc + T14 / 2)
    return np.clip(hi - lo, 0, None) / dt


def in_any_transit(t, dt, eph, pad):
    """True for cadences overlapping any transit (+/- pad days) of ephemeris `eph`."""
    ph = ((t - eph["T0"]) / eph["P"] + 0.5) % 1.0 - 0.5
    return np.abs(ph * eph["P"]) < eph["T14"] / 2 + pad + dt / 2
