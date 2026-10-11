"""Independent verification of Audit 5 (signal B depth). Written from scratch; reads only raw FITS in ../../00_data/mast.
Method deliberately different from the auditor's: no wotan; per-event local polynomial baseline fitted jointly with a
batman template (prior shape) and, where A overlaps, an A template; errors from (a) white OOT scatter, (b) 'null fits'
of the same template at many random out-of-transit times of the same light curve (empirical, includes red noise)."""
import os, json, sys
import numpy as np
from astropy.io import fits
import batman

H = os.path.dirname(os.path.abspath(__file__))
M = os.path.abspath(os.path.join(H, "..", "..", "00_data", "mast"))
rng = np.random.default_rng(2026)

PA, TA, T14A = 3.1827856, 3998.7375, 1.87 / 24
PB, TB = 11.13274, 4000.865
RHO_SUN = 1.98841e33 / (4 / 3 * np.pi * 6.957e10 ** 3); G = 6.6743e-8
def ars(rho, P): return (rho * RHO_SUN * G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)

def tmodel(k, b, rho, P, u=(0.32, 0.31)):
    pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = k; pm.a = ars(rho, P)
    pm.inc = np.degrees(np.arccos(b / pm.a)); pm.ecc = 0; pm.w = 90; pm.u = list(u); pm.limb_dark = "quadratic"
    return pm

def unit(pm, trel, cad):
    ss = max(1, int(np.ceil(cad * 1440 / 1.0)));  ss += (ss % 2 == 0)
    m = batman.TransitModel(pm, np.asarray(trel, float), supersample_factor=ss, exp_time=cad) if ss > 1 else batman.TransitModel(pm, np.asarray(trel, float))
    m0 = batman.TransitModel(pm, np.array([0.0]))
    return (1 - m.light_curve(pm)) / (1 - m0.light_curve(pm)[0])

PMB = tmodel(0.0330, 0.875, 2.51, PB)
PMA = tmodel(0.0212, 0.35, 2.51, PA)

def dtph(t, P, T0):
    return ((t - T0) / P + 0.5) % 1 - 0.5

def load(sec, which="base", keep512=False):
    if sec == 6:
        f = "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits"; col = "PDCSAP_FLUX"
    elif sec == 32:
        f = "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits"; col = "PDCSAP_FLUX"
    else:
        f = "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"; col = "SYS_RM_FLUX"
    h = fits.open(os.path.join(M, f)); d = h[1].data
    t = np.array(d["TIME"], float); fl = np.array(d[col], float); q = np.array(d["QUALITY"]).astype(np.int64)
    cad = float(h[1].header["TIMEDEL"])
    if sec == 98:
        bm = 1 | 4 | 8 | 16 | 32 | 128 | 2048 | 4096 | (1 << 29) | (1 << 30)
    else:
        bm = 17087 & ~512 if keep512 else 17087
    g = np.isfinite(t) & np.isfinite(fl) & ((q & bm) == 0)
    return t[g], fl[g] / np.nanmedian(fl[g]), q[g], cad

def design(t, tcs, cad, poly, Amode):
    """Amode: 'fixed' subtract 522 ppm A; 'fit' free A amp; 'mask' drop |dtA|<T14A (manuscript), 'maskhalf' drop |dtA|<0.6 T14A"""
    cols = []; rows = np.zeros(len(t), bool); blocks = []
    for j, tc in enumerate(tcs):
        w = np.abs(t - tc) < 0.3
        blocks.append(w)
    return blocks

def fit(t, f, cad, tcs, poly=1, Amode="fixed", perevent=False, A_ppm=522.0, hw=0.3):
    dA = dtph(t, PA, TA) * PA
    keep = np.ones(len(t), bool)
    if Amode == "mask": keep &= np.abs(dA) >= T14A
    if Amode == "maskhalf": keep &= np.abs(dA) >= 0.6 * T14A
    ff = f.copy()
    if Amode == "fixed": ff = ff + A_ppm * 1e-6 * unit(PMA, dA, cad)
    ws = [(np.abs(t - tc) < hw) & keep for tc in tcs]
    ws = [(w, tc) for w, tc in zip(ws, tcs) if w.sum() > 10 and np.sum(np.abs(t[w] - tc) < 0.035) >= 1]
    N = sum(w.sum() for w, _ in ws); ne = len(ws); nb = (poly + 1) * ne; nd = ne if perevent else 1
    needA = Amode == "fit" and any(np.any(np.abs(dA[w]) < 0.6 * T14A) for w, _ in ws)
    X = np.zeros((N, nb + nd + needA)); y = np.zeros(N); i = 0
    for j, (w, tc) in enumerate(ws):
        n = w.sum(); x = (t[w] - tc) / hw
        for p in range(poly + 1): X[i:i + n, j * (poly + 1) + p] = x ** p
        X[i:i + n, nb + (j if perevent else 0)] = -unit(PMB, t[w] - tc, cad)
        if needA: X[i:i + n, -1] = -unit(PMA, dA[w], cad)
        y[i:i + n] = ff[w]; i += n
    c, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ c
    s = 1.4826 * np.median(np.abs(r - np.median(r)))
    cov = np.linalg.pinv(X.T @ X) * s ** 2
    D = c[nb:nb + nd] * 1e6; e = np.sqrt(np.diag(cov)[nb:nb + nd]) * 1e6
    return D, e, ne

def btimes(t):
    n0 = int(np.floor((t.min() - TB) / PB)) - 1; n1 = int(np.ceil((t.max() - TB) / PB)) + 1
    out = []
    for n in range(n0, n1 + 1):
        tc = TB + n * PB
        if np.sum(np.abs(t - tc) < 0.035) >= 1: out.append((n, tc))
    return out

def null_err(t, f, cad, ntc, poly, Amode, nnull=400):
    """same fit at random fake mid-times (away from real A/B), grouped in sets of ntc events; std of fitted depths"""
    dA = np.abs(dtph(t, PA, TA) * PA); dB = np.abs(dtph(t, PB, TB) * PB)
    Ds = []; tries = 0
    while len(Ds) < nnull and tries < nnull * 30:
        tries += 1
        tcs = rng.choice(t, ntc)
        ok = True
        for tc in tcs:
            w = np.abs(t - tc) < 0.3
            if np.any(dB[w] < 0.1) or np.sum(np.abs(t - tc) < 0.035) * cad < 0.035 or w.sum() < 0.3 * 0.6 / cad:
                ok = False; break
            # require full in-transit coverage
            if np.sum(np.abs(t - tc) < 0.0345) * cad < 0.85 * 0.069: ok = False; break
        if not ok: continue
        D, e, ne = fit(t, f, cad, tcs, poly, Amode)
        Ds.append(D[0])
    Ds = np.array(Ds)
    return float(np.std(Ds)), float(np.mean(Ds)), len(Ds)

def detrend_median(t, f, win=0.75):
    """simple own detrend: running median of out-of-transit points (A,B masked), interpolated"""
    m = (np.abs(dtph(t, PA, TA) * PA) > T14A) & (np.abs(dtph(t, PB, TB) * PB) > 0.07)
    tt, ff = t[m], f[m]; tr = np.empty(len(tt))
    for i in range(len(tt)):
        w = np.abs(tt - tt[i]) < win / 2
        tr[i] = np.median(ff[w])
    return f / np.interp(t, tt, tr)

if __name__ == "__main__":
    res = {}
    info = {}
    for sec in (6, 32, 98):
        t, f, q, cad = load(sec)
        bt = btimes(t)
        info[sec] = dict(cad_d=cad, B_epochs=[n for n, _ in bt], B_tc=[tc for _, tc in bt])
        if sec == 32:
            info[sec]["A_minus_B_h"] = []
            for n, tc in bt:
                ph = dtph(tc, PA, TA) * PA * 24
                inB = np.abs(t - tc) < 1.62 / 48
                maskd = inB & (np.abs(dtph(t, PA, TA) * PA) < T14A)
                info[sec]["A_minus_B_h"].append(dict(epoch=n, tB_minus_tA_h=ph, frac_B_intransit_masked_by_pmT14A=float(maskd.sum() / inB.sum())))
    print(json.dumps(info, indent=1))
    res["info"] = info
    # S6 bit 512
    h = fits.open(os.path.join(M, "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits")); d = h[1].data
    q = np.array(d["QUALITY"]).astype(int); t6 = np.array(d["TIME"]); fl = np.array(d["PDCSAP_FLUX"])
    b512 = np.where((q & 512) != 0)[0]
    res["S6_bit512"] = [dict(t=float(t6[i]), flux_rel_ppm=float((fl[i] / np.nanmedian(fl) - 1) * 1e6), q=int(q[i]),
                             dt_from_B_h=float(dtph(t6[i], PB, TB) * PB * 24)) for i in b512]
    print("S6 bit512:", res["S6_bit512"])
    # S98 flagged B epochs
    h = fits.open(os.path.join(M, "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")); d = h[1].data
    t98 = np.array(d["TIME"]); q98 = np.array(d["QUALITY"]).astype(np.int64)
    ev98 = []
    for n in range(-2, 6):
        tc = TB + n * PB; w = np.abs(t98 - tc) < 1.62 / 48
        if w.sum() == 0: continue
        qq = q98[w]
        ev98.append(dict(epoch=n, n_in=int(w.sum()), frac_2048=float(np.mean(qq & 2048 > 0)), frac_b29=float(np.mean(qq & (1 << 29) > 0)),
                         frac_b30=float(np.mean(qq & (1 << 30) > 0)), frac_clean=float(np.mean((qq & (1 | 4 | 8 | 16 | 32 | 128 | 2048 | 4096 | (1 << 29) | (1 << 30))) == 0))))
    res["S98_flags"] = ev98
    print(json.dumps(ev98, indent=1))
    json.dump(res, open(os.path.join(H, "v01_info.json"), "w"), indent=1)
