"""Configurable per-sector depth measurement for signal B (one consistent method, many switches).

Model in each +-HW window around a predicted B mid-time (fixed linear ephemeris unless dt mode says otherwise):
    fd(t) = poly_event(t) - D_sector * s_B(t - tc - dt) - D_A,sector * s_A(t - tcA)
s_B, s_A: unit-mid-depth batman templates (fixed shapes), exposure-integrated (supersampled to <=1 min).
The A term is included only when A in-transit cadences fall in the window (only Sector 32 in practice).
"""
import numpy as np
from b_lib import *

SHAPE_A = dict(k=0.0212, b=0.35, rho=2.51, u1=LD_NOM[0], u2=LD_NOM[1])
_TPL_CACHE = {}


def tplA():
    if "A" not in _TPL_CACHE:
        _TPL_CACHE["A"] = Template(SHAPE_A, P=EPH_A["P"])
    return _TPL_CACHE["A"]


_LC_CACHE = {}


def get_lc(key, s98_mask=MASK_S98_BASE, spoc_mask=16575):
    """product loaders. spoc_mask default 16575 = lightkurve default 17087 without bit 10 (512, PDC impulsive outlier)."""
    ck = (key, s98_mask, spoc_mask)
    if ck in _LC_CACHE:
        return _LC_CACHE[ck]
    if key == "S6_tessspoc_pdc":
        lc = load_spoc(F_S6_SPOC, key, 6, "PDCSAP_FLUX", spoc_mask)
    elif key == "S6_tessspoc_sap":
        lc = load_spoc(F_S6_SPOC, key, 6, "SAP_FLUX", spoc_mask, crowd_correct_sap=True)
    elif key == "S6_qlp_ksp":
        lc = load_qlp(F_S6_QLP, key, 6, "KSPSAP_FLUX", QLP_SPACECRAFT | 4096)
    elif key == "S6_qlp_sap":
        lc = load_qlp(F_S6_QLP, key, 6, "SAP_FLUX", QLP_SPACECRAFT | 4096)
    elif key == "S6_tglc_aper":
        h = fits.open(F_S6_TGLC); d = h[1].data
        t = np.asarray(d["time"], float); f = np.asarray(d["cal_aper_flux"], float)
        q = np.asarray(d["TESS_flags"]).astype(np.int64); tq = np.asarray(d["TGLC_flags"]).astype(np.int64)
        m = ((q & spoc_mask) == 0) & (tq == 0)
        lc = dict(name=key, sector=6, cad=1800 / 86400.0, meta=dict(file=os.path.basename(F_S6_TGLC), col="cal_aper_flux"))
        good = m & np.isfinite(f) & np.isfinite(t) & (f > 0)
        lc.update(t=t[good], f=f[good] / np.median(f[good]), e=None, q=q[good])
    elif key == "S6_eleanor_corr":
        h = fits.open(F_S6_ELEANOR); d = h[1].data
        t = np.asarray(d["TIME"], float); f = np.asarray(d["CORR_FLUX"], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
        good = ((q & spoc_mask) == 0) & np.isfinite(f) & np.isfinite(t) & (f > 0)
        lc = dict(name=key, sector=6, cad=1800 / 86400.0, t=t[good], f=f[good] / np.median(f[good]), e=None, q=q[good],
                  meta=dict(file=os.path.basename(F_S6_ELEANOR), col="CORR_FLUX"))
    elif key == "S32_spoc2m_pdc":
        lc = load_spoc(F_S32_2M, key, 32, "PDCSAP_FLUX", spoc_mask)
    elif key == "S32_spoc2m_sap":
        lc = load_spoc(F_S32_2M, key, 32, "SAP_FLUX", spoc_mask, crowd_correct_sap=True)
    elif key == "S32_tessspoc10_pdc":
        lc = load_spoc(F_S32_SPOC10, key, 32, "PDCSAP_FLUX", spoc_mask)
    elif key == "S32_qlp_ksp":
        lc = load_qlp(F_S32_QLP, key, 32, "KSPSAP_FLUX", (QLP_SPACECRAFT & ~32) | 4096)   # QLP S32 bit 32 misplaced (00_data README)
    elif key == "S32_qlp_sap":
        lc = load_qlp(F_S32_QLP, key, 32, "SAP_FLUX", (QLP_SPACECRAFT & ~32) | 4096)
    elif key == "S98_qlp_sysrm":
        lc = load_qlp(F_S98_QLP, key, 98, "SYS_RM_FLUX", s98_mask)
    elif key == "S98_qlp_det":
        lc = load_qlp(F_S98_QLP, key, 98, "DET_FLUX", s98_mask)
    elif key == "S98_qlp_sap":
        lc = load_qlp(F_S98_QLP, key, 98, "SAP_FLUX", s98_mask)
    elif key.startswith("own"):
        # own{sector}_r{radius}
        sec = int(key[3:].split("_")[0]); rad = float(key.split("_r")[1])
        extra = qlp_s98_bad_time_fn(s98_mask) if sec == 98 else None
        lc = own_phot(sec, rad, bitmask=spoc_mask | 2048, extra_bad_times=extra)
        lc["name"] = key
    else:
        raise KeyError(key)
    _LC_CACHE[ck] = lc
    return lc


def build_windows(lc, evs, tplB, tA, poly, use, dts, supersample, per_event_depth, fitA):
    """returns X (N x p), y, index of B column(s), index of A column (or None)"""
    tt_all, y_all, cols = [], [], []
    ne = len(evs); nb = (poly + 1) * ne
    nd = ne if per_event_depth else 1
    blocks = []
    for j, ev in enumerate(evs):
        tt = lc["t"][ev["w"]]; blocks.append((j, ev, tt)); tt_all.append(tt); y_all.append(lc[use][ev["w"]])
    y = np.concatenate(y_all); N = len(y)
    tcat = np.concatenate(tt_all)
    dA = phase_dt(tcat, EPH_A)
    needA = fitA and np.any(np.abs(dA) < 0.6 * EPH_A["T14"])
    X = np.zeros((N, nb + nd + (1 if needA else 0)))
    i0 = 0
    for j, ev, tt in blocks:
        n = len(tt); x = (tt - ev["tc"]) / 0.3
        for p in range(poly + 1):
            X[i0:i0 + n, j * (poly + 1) + p] = x ** p
        dt = dts[j] if np.ndim(dts) else dts
        X[i0:i0 + n, nb + (j if per_event_depth else 0)] = -tplB.unit(tt - ev["tc"] - dt, lc["cad"], supersample)
        if needA:
            X[i0:i0 + n, -1] = -tA.unit(phase_dt(tt, EPH_A), lc["cad"], supersample)
        i0 += n
    return X, y, nb, needA


def lsq(X, y, sig):
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c
    cov = np.linalg.pinv(X.T @ X) * sig ** 2
    return c, cov, r


def measure(lc_d, tplB, cfg, eph=EPH_B):
    """lc_d: detrended light curve (or raw with use='f'); returns per-sector + per-event results."""
    use = cfg.get("use", "fd"); poly = cfg.get("poly", 1); hw = cfg.get("hw", 0.3); ss = cfg.get("supersample", True)
    A_mode = cfg.get("A_mode", "model")
    dropA = {"model": 0.0, "fixed": 0.0, "mask0.6": 0.6 * EPH_A["T14"], "mask1.0": 1.0 * EPH_A["T14"]}[A_mode]
    if A_mode == "fixed":
        # subtract the A template with a fixed mid-depth (default 522 ppm: k=0.0212, b=0.35, LD_NOM)
        lc_d = dict(lc_d)
        lc_d[use] = lc_d[use] + cfg.get("A_depth_ppm", 522.0) * 1e-6 * tplA().unit(phase_dt(lc_d["t"], EPH_A), lc_d["cad"], ss)
    evs = events_in(lc_d, eph, T14=tplB.t14_d, half_window=hw, min_cov=cfg.get("min_cov", 0.25), drop_A_hw=dropA if dropA > 0 else 1e-9)
    if cfg.get("drop_partial"):
        evs = [e for e in evs if e["cov"] >= 0.9]
    if cfg.get("epochs") is not None:
        evs = [e for e in evs if e["epoch"] in cfg["epochs"]]
    if not evs:
        return None
    fitA = A_mode == "model"
    sig = lc_d["sig"]
    dtmode = cfg.get("dt", "fixed")
    grid = np.arange(-60, 60.01, 1.0) / 1440.0
    if dtmode == "fixed":
        dts = 0.0
    elif dtmode == "sector":
        best = None
        for g in grid:
            X, y, nb, needA = build_windows(lc_d, evs, tplB, tplA(), poly, use, g, ss, False, fitA)
            c, cov, r = lsq(X, y, sig); ch = np.sum(r ** 2)
            if best is None or ch < best[0]:
                best = (ch, g)
        dts = best[1]
    elif dtmode == "event":
        dts = np.zeros(len(evs))
        for j, ev in enumerate(evs):
            best = None
            for g in grid:
                X, y, nb, needA = build_windows(lc_d, [ev], tplB, tplA(), poly, use, g, ss, False, fitA)
                c, cov, r = lsq(X, y, sig); ch = np.sum(r ** 2)
                if best is None or ch < best[0]:
                    best = (ch, g)
            dts[j] = best[1]
    X, y, nb, needA = build_windows(lc_d, evs, tplB, tplA(), poly, use, dts, ss, False, fitA)
    c, cov, r = lsq(X, y, sig)
    D = c[nb] * 1e6; eD = np.sqrt(cov[nb, nb]) * 1e6
    chi2 = float(np.sum((r / sig) ** 2)); dof = len(y) - X.shape[1]
    Xe, ye, nbe, _ = build_windows(lc_d, evs, tplB, tplA(), poly, use, dts, ss, True, fitA)
    ce, cove, re_ = lsq(Xe, ye, sig)
    De = (ce[nbe:nbe + len(evs)] * 1e6).tolist(); eDe = (np.sqrt(np.diag(cove))[nbe:nbe + len(evs)] * 1e6).tolist()
    res = dict(D=float(D), eD_white=float(eD), chi2=chi2, dof=int(dof), chi2_red=float(chi2 / max(dof, 1)),
               eD_white_scaled=float(eD * np.sqrt(max(chi2 / max(dof, 1), 1.0))),
               sig_ppm=float(sig * 1e6), n=int(len(y)), epochs=[e["epoch"] for e in evs], cov=[round(e["cov"], 3) for e in evs],
               D_event=De, eD_event_white=eDe, dt_min=(np.atleast_1d(dts) * 1440).tolist(),
               A_in_window=bool(needA), D_A=float(c[-1] * 1e6) if needA else None)
    return res


def run_config(cfg, tplB):
    """cfg: products {6:key,32:key,98:key}, detrend dict, masks, fit options. Returns per-sector results + summary."""
    out = {}
    for sec in (6, 32, 98):
        key = cfg["products"][sec]
        lc = get_lc(key, cfg.get("s98_mask", MASK_S98_BASE), cfg.get("spoc_mask", 16575))
        det = cfg.get("detrend", dict(method="biweight", window=0.75))
        if det["method"] == "none":
            lcd = dict(lc); lcd["fd"] = lc["f"].copy()
            oot = ~(in_transit(lc["t"], EPH_B, EPH_B["T14"]) | in_transit(lc["t"], EPH_A, EPH_A["T14"]))
            # local noise from first differences (insensitive to slow trends)
            lcd["sig"] = float(1.4826 * np.median(np.abs(np.diff(lc["f"][oot]))) / np.sqrt(2))
        else:
            kw = {k: v for k, v in det.items() if k not in ("method", "window")}
            lcd = detrend(lc, method=det["method"], window=det["window"], **kw)
        r = measure(lcd, tplB, cfg)
        r["product"] = key
        out[sec] = r
    return out


def summarize(out, errkey="eD_white", kappa=None):
    """kappa: optional {sector: inflation factor} applied to errkey"""
    dep = {}
    for s in (6, 32, 98):
        e = out[s][errkey] * (kappa[s] if kappa else 1.0)
        dep[s] = (out[s]["D"], e)
    return s98_vs_archival(dep), dep


def prep_detrended(lc, cfg, extra_mask=None):
    det = cfg.get("detrend", dict(method="biweight", window=0.75))
    if det["method"] == "none":
        lcd = dict(lc); lcd["fd"] = lc["f"].copy()
        oot = ~(in_transit(lc["t"], EPH_B, EPH_B["T14"]) | in_transit(lc["t"], EPH_A, EPH_A["T14"]))
        if extra_mask is not None:
            oot &= ~extra_mask
        lcd["sig"] = float(1.4826 * np.median(np.abs(np.diff(lc["f"][oot]))) / np.sqrt(2))
        return lcd
    kw = {k: v for k, v in det.items() if k not in ("method", "window")}
    return detrend(lc, method=det["method"], window=det["window"], extra_mask=extra_mask, **kw)


def inject_measure(lc, tplB, cfg, eph_inj, D_inj_ppm, hw_mask=None):
    """multiply the raw normalised flux by an injected B-template transit train, re-detrend (masking real A, real B and
    the injected events), and measure the injected depth with the same configuration."""
    lci = dict(lc)
    dt = phase_dt(lc["t"], eph_inj)
    lci["f"] = lc["f"] * (1.0 - D_inj_ppm * 1e-6 * tplB.unit(dt, lc["cad"], cfg.get("supersample", True)))
    hw_mask = hw_mask or EPH_B["T14"]
    lcd = prep_detrended(lci, cfg, extra_mask=np.abs(dt) < hw_mask)
    return measure(lcd, tplB, cfg, eph=eph_inj)
