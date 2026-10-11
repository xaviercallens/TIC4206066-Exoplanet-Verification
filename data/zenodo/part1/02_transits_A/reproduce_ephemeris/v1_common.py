"""Verifier's independent timing code (written from scratch; no import from auditor or package code).
Method deliberately differs from the auditor: detrending = quadratic polynomial in time fitted jointly with the
template inside a +-0.3 d window (no wotan); depth constrained >= 0 (negative-depth fits are replaced by the
no-transit chi2); template via batman with supersampling."""
import numpy as np, batman, hashlib, json
from astropy.io import fits
D = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/"
FILES = {6: ("mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX"),
         32: ("mast/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX"),
         98: ("mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX")}
PB, T0B, T14B = 11.13274, 4000.865, 1.66/24
PA, T0A = 3.1827856, 3998.7375
W = 0.30

def load(s):
    rel, col = FILES[s]
    raw = open(D+rel, "rb").read(); sh = hashlib.sha256(raw).hexdigest()
    man = {e["relpath"]: e["sha256"] for e in json.load(open(D+"DATA_MANIFEST.json"))["files"]}
    assert man[rel] == sh
    d = fits.open(D+rel)[1].data
    t = np.array(d["TIME"], float); f = np.array(d[col], float); q = np.array(d["QUALITY"]).astype(np.int64)
    good = np.isfinite(t) & np.isfinite(f) & (((q & 17087) == 0) if s != 98 else (q == 0))
    t, f = t[good], f[good]; f = f/np.nanmedian(f)
    return t, f, float(np.median(np.diff(t))), sh

def model_fn(cad):
    pm = batman.TransitParams(); pm.t0=0; pm.per=PA; pm.rp=0.0212
    rho=2.51*1.40996; G=6.674e-8; pm.a=(rho*G*(PA*86400)**2/(3*np.pi))**(1/3)
    pm.inc=np.degrees(np.arccos(0.35/pm.a)); pm.ecc=0; pm.w=90; pm.u=[0.37,0.28]; pm.limb_dark="quadratic"
    tau = np.arange(-0.35, 0.35, 10/86400)
    ss = max(1, int(round(cad*86400/20)))
    m = batman.TransitModel(pm, tau, supersample_factor=ss, exp_time=cad) if ss > 1 else batman.TransitModel(pm, tau)
    fl = m.light_curve(pm); s = (1-fl)/0.0212**2   # unit ~ depth k^2 scale
    s = s/ s.max()
    return lambda x: np.interp(x, tau, s, left=0, right=0)

T14A = 1.887/24
def bmask(t, pad=0.4/24):
    n = np.round((t-T0B)/PB); return np.abs(t-(T0B+n*PB)) < T14B/2+pad

def scan(t, y, tc, tpl, R=144, step=0.25, dpos=True):
    x = t-tc; dts = np.arange(-R, R+step/2, step)/1440
    P0 = np.vstack([np.ones_like(x), x, x*x]).T
    c0, *_ = np.linalg.lstsq(P0, y, rcond=None); chi_null = np.sum((y-P0@c0)**2)
    chis = np.empty(len(dts)); Ds = np.empty(len(dts))
    for i, dt in enumerate(dts):
        s = tpl(x-dt)
        if s.max() == 0: chis[i]=chi_null; Ds[i]=0; continue
        X = np.column_stack([P0, -s]); c, *_ = np.linalg.lstsq(X, y, rcond=None)
        if dpos and c[3] < 0: chis[i]=chi_null; Ds[i]=0; continue
        chis[i]=np.sum((y-X@c)**2); Ds[i]=c[3]
    return dts*1440, chis, Ds, chi_null

def best(dts, chis, Ds, R):
    m = np.abs(dts) <= R+1e-9; idx = np.where(m)[0]; j = idx[np.argmin(chis[idx])]
    return dts[j], Ds[j], j, (j == idx[0] or j == idx[-1])

def window(t, f, tc):
    w = (np.abs(t-tc) < W) & ~bmask(t)
    return t[w], f[w]

def baseline(t, tc):
    pre = np.sum((t > tc-T14A/2-0.5*T14A) & (t < tc-T14A/2)); post = np.sum((t > tc+T14A/2) & (t < tc+T14A/2+0.5*T14A))
    cov = np.sum(np.abs(t-tc) < T14A/2)
    return pre, post, cov
