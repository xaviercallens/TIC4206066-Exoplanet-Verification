"""Verifier helpers (written from scratch; no code imported from the auditor or earlier analyses).
Inputs: raw public FITS in ../../00_data/mast (sha256 checked against DATA_MANIFEST.json).
Time: BTJD = BJD_TDB - 2457000 in all products used."""
import os, json, hashlib
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "00_data"))
MAST = os.path.join(DATA, "mast")

F = dict(
    qlp98="hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits",
    qlp6="hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits",
    qlp32="hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits",
    ts6="hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits",
    s32="tess2020324010417-s0032-0000000004206066-0200-s_lc.fits",
    tglc6="hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits",
)

# ephemerides
P_A, T0_A, T14_A = 3.1827856, 3998.7375, 1.87 / 24
P_MS, T0_MS = 11.13274, 4000.865          # manuscript B ephemeris (BTJD)
DUR = 1.6 / 24
# SPOC quality bits treated as bad (SDPDD: 1,2,4,8,16,32,128,1024(2-min: impulsive outlier),2048? no, 4096,8192,32768)
SPOC_BAD = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 1024 | 4096 | 8192 | 32768


def sha(fn):
    h = hashlib.sha256()
    with open(fn, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def check_manifest():
    m = json.load(open(os.path.join(DATA, "DATA_MANIFEST.json")))
    byname = {os.path.basename(e.get("relpath", "")): e for e in m["files"]}
    out = {}
    for k, fn in F.items():
        s = sha(os.path.join(MAST, fn))
        e = byname.get(fn, {})
        out[fn] = dict(sha256=s, manifest=e.get("sha256"), match=(s == e.get("sha256")), url=e.get("url") or e.get("final_url"))
    return out


def inA(t, half=T14_A):
    ph = ((t - T0_A) / P_A + 0.5) % 1 - 0.5
    return np.abs(ph * P_A) < half


def biweight(t, f, mask, window=0.75):
    from wotan import flatten
    use = np.isfinite(f) & ~mask
    _, tr = flatten(t[use], f[use], method="biweight", window_length=window, return_trend=True, break_tolerance=0.4)
    ok = np.isfinite(tr)
    trend = np.interp(t, t[use][ok], tr[ok])
    return f / trend


def rstd(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def load(kind):
    d = fits.open(os.path.join(MAST, F[kind]))[1].data
    t = np.asarray(d["TIME"], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    return d, t, q
