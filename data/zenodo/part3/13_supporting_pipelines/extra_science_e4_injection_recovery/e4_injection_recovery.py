"""Injection-recovery completeness of a BLS search in the TESS Sector 98 QLP light curve of TIC 4206066 (the discovery sector,
57-d double sector). Data: U/00_data/mast QLP s0098 llc, SYS_RM_FLUX, QUALITY == 0 only (all QLP-flagged cadences removed, as in
the Note). Signals A and B masked (+-1.5 T14) with the reconciled ephemerides. Injection: batman, circular, quadratic LD
(Claret 2017 PHOENIX, from U/07 common.claret_ld(): 0.316, 0.313), R* 0.623, M* 0.606, b ~ U(0, 0.9), random phase,
exposure 200 s (supersample 5). Detrend: wotan biweight 0.75 d (as the discovery search). Search: astropy BoxLeastSquares,
periods 0.5-28 d (uniform in frequency, 0.3 x duration/baseline spacing), durations 0.5-4 h.
Recovered = one of the top-5 iterative peaks within 0.5 % of P_inj (or P/2, 2P) and in phase (1.5 T14), at S/N >= 8 (the real
signal A reaches S/N 8.2 in this simple search, 4th peak) or S/N >= 12 (above the strongest systematic peak of the null, 15.56 d, S/N 12.0).
S/N = astropy depth_snr (depth / depth_err) with a constant per-point error = 1.4826*MAD of the flattened flux (white noise)."""
import os, sys, json, numpy as np
from multiprocessing import Pool
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
import batman
from wotan import flatten
F = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"
d = fits.open(F)[1].data
m = (d["QUALITY"] == 0) & np.isfinite(d["SYS_RM_FLUX"]) & np.isfinite(d["TIME"])
t0 = np.array(d["TIME"][m], float); f0 = np.array(d["SYS_RM_FLUX"][m], float); f0 /= np.nanmedian(f0)
BT = 2457000.0
EPH = dict(A=(3.182785, 2460998.736 - BT, 1.87/24), B=(11.13274, 2461000.865 - BT, 1.66/24))
keep = np.ones_like(t0, bool)
for P, T0, T14 in EPH.values():
    ph = ((t0 - T0 + 0.5*P) % P) - 0.5*P; keep &= np.abs(ph) > 1.5*T14
t, f = t0[keep], f0[keep]
Rs, Ms = 0.623, 0.606; RE_RS = 6.371e6/6.957e8
G, Msun = 6.674e-11, 1.989e30
base = t.max() - t.min()
def search(tt, ff):
    fl = flatten(tt, ff, method="biweight", window_length=0.75)
    ok = np.isfinite(fl); r = fl[ok] - np.median(fl[ok]); sig = 1.4826*np.median(np.abs(r))
    bls = BoxLeastSquares(tt[ok], fl[ok], dy=np.full(ok.sum(), sig))
    fmin, fmax = 1/28., 1/0.5; df = 0.3*(0.5/24)/base**2*base  # ~0.3*q/base with q~0.5h/base: conservative step
    freqs = np.arange(fmin, fmax, 0.25*(1.0/24)/base/1.0)
    res = bls.power(1/freqs, [0.5/24, 1/24, 1.5/24, 2.5/24, 4/24], objective="snr")
    i = np.nanargmax(res.power); return res.period[i], res.transit_time[i], res.depth_snr[i], res.power[i]
def one(args):
    try:
        return _one(args)
    except Exception as ex:
        return dict(seed=args[0], P=args[1], rp=args[2], error=repr(ex))
def _one(args):
    seed, P, rp = args
    rng = np.random.default_rng(seed)
    b = rng.uniform(0, 0.9); T0 = t.min() + rng.uniform(0, P)
    a = (G*Ms*Msun*(P*86400)**2/(4*np.pi**2))**(1/3)/(Rs*6.957e8)
    if b >= a: return None
    pm = batman.TransitParams(); pm.t0, pm.per, pm.rp, pm.a = T0, P, rp*RE_RS/Rs, a
    pm.inc = np.degrees(np.arccos(b/a)); pm.ecc, pm.w, pm.u, pm.limb_dark = 0., 90., [0.316, 0.313], "quadratic"
    mod = batman.TransitModel(pm, t, supersample_factor=5, exp_time=200/86400.).light_curve(pm)
    ntr = len(np.unique(np.round((t[mod < 1] - T0)/P)))
    T14 = P/np.pi*np.arcsin(min(1.0, np.sqrt(max((1+pm.rp)**2 - b**2, 0))/a))
    ff = f*mod; keep = np.ones_like(t, bool); peaks = []; hit8 = hit12 = False
    for it in range(5):                      # iterative search: top-5 peaks, each masked (+-0.15 d) before the next
        if keep.sum() < 2000: break             # short trial periods can mask most of the data
        Pb, Tb, snr, pw = search(t[keep], ff[keep]); peaks.append((float(Pb), float(snr)))
        okP = any(abs(Pb/(P*r) - 1) < 0.005 for r in (1, 0.5, 2))
        dph = ((Tb - T0 + 0.5*P) % P) - 0.5*P
        okT = abs(dph) < max(1.5*T14, 0.05) or abs(abs(dph) - 0.5*P) < max(1.5*T14, 0.05)
        if okP and okT: hit8 |= snr >= 8.0; hit12 |= snr >= 12.0
        ph = ((t - Tb + 0.5*Pb) % Pb) - 0.5*Pb; keep &= np.abs(ph) > 0.15
    return dict(seed=seed, P=P, rp=rp, b=b, n_tr=int(ntr), peaks=peaks, rec8=bool(hit8), rec12=bool(hit12))
if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    rng = np.random.default_rng(2026)
    jobs = [(i, float(np.exp(rng.uniform(np.log(0.5), np.log(28)))), float(rng.uniform(0.5, 4.0))) for i in range(n)]
    # null (no injection) check
    print("null search:", search(t, f), "N points", len(t), flush=True)
    with Pool(int(os.environ.get("NPROC", "4"))) as pool:
        out = []
        for k, r in enumerate(pool.imap_unordered(one, jobs, chunksize=4)):
            if r: out.append(r)
            if k % 200 == 0: print(k, flush=True); json.dump(out, open("e4_injections.json", "w"))
    json.dump(out, open("e4_injections.json", "w")); print("done", len(out))
