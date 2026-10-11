"""R-S6 (c) and (d): is the Sector 6 photometric signal consistent with Sectors 32+98, is it affected by the S6 momentum
dumps, does S6 matter for the odd/even, phase-0.5 and V-shape tests, and does cadence integration (supersampling) of the
30-min S6 data matter.

All fits are linear in depth with a fixed-shape unit template (k=0.0212, b=0.35, rho=2.51, Claret LD; D = unbinned
mid-transit depth), cadence-integrated per sector (boxcar over TIMEDEL), one linear baseline per +-0.3 d window, raw flux.
Transit centres from the S32+S98-only timing ephemeris of s1 (so nothing in the S6 test uses S6).
Empirical errors: the same fit repeated with the whole ephemeris shifted by random phases that avoid the real transits
and phase 0.5 ('fake ephemerides', 400 per test).
"""
import json, time
import numpy as np
import batman
from rs6_common import *

rng = np.random.default_rng(606)
U = claret_u()
S1 = json.load(open(os.path.join(OUT, "s1_timing_R60.json")))
E3298 = S1["ephemerides"]["S32_S98_sel"]
PE, TE = E3298["P"], E3298["T0"]
NFAKE = 400


def make_tpl(k=0.0212, b=0.35, rho=2.51):
    tau = np.arange(-0.40, 0.40 + 1e-9, 2.0 / 86400.0)
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp = 0.0, P_A, k
    pm.a = ars_from_rho(rho, P_A); pm.inc = np.degrees(np.arccos(b / pm.a)); pm.ecc, pm.w = 0.0, 90.0
    pm.u, pm.limb_dark = list(U), "quadratic"
    fl = batman.TransitModel(pm, tau).light_curve(pm)
    s = (1 - fl) / (1 - fl[np.argmin(np.abs(tau))])
    return tau, s


TAU, S_UNB = make_tpl()


def binned(s_unb, sec, integrate=True):
    if not integrate:
        return s_unb
    from scipy.signal import fftconvolve
    n = int(round(CAD_S[sec] / 2.0))
    return fftconvolve(s_unb, np.ones(n) / n, mode="same")


TPL = {sec: binned(S_UNB, sec) for sec in CAD_S}
TPL_NOINT = {sec: S_UNB for sec in CAD_S}


def box_unb(width_d):
    return (np.abs(TAU) < width_d / 2).astype(float)


def trap_unb(T14d, T12d):
    a = np.abs(TAU)
    s = np.clip((T14d / 2 - a) / max(T12d, 1e-6), 0, 1)
    return s


LC = {}
for sec in (6, 32, 98):
    t, f, meta = load(sec)
    dB, _ = phase_dt(t, P_B, T0_B)
    keep = np.abs(dB) > T14_B
    LC[sec] = dict(t=t[keep], f=f[keep], meta=meta)
DUMPS = load_dumps_s6()


def collect(sec, phase=0.0, dt_min=0.0, parity=None, excl_real=False, extra_mask=None, half=0.30):
    """Windows around T = TE + (E + phase) PE + dt.  Returns t, f, x (t - centre), window id, epochs."""
    t, f = LC[sec]["t"], LC[sec]["f"]
    c = TE + phase * PE + dt_min / 1440.0
    dd, ep = phase_dt(t, PE, c)
    m = np.abs(dd) < half
    if excl_real:
        dA, _ = phase_dt(t, PE, TE)
        m &= np.abs(dA) > T14_A
        m &= np.abs(phase_dt(t, PE, TE + 0.5 * PE)[0]) > T14_A
    if parity is not None:
        m &= (ep % 2) == parity
    if extra_mask is not None:
        m &= ~extra_mask(t)
    t, f, dd, ep = t[m], f[m], dd[m], ep[m]
    # keep windows with enough points and in-window coverage
    keep = np.zeros(len(t), bool)
    for e in np.unique(ep):
        w = ep == e
        if w.sum() >= 8 and np.sum(np.abs(dd[w]) < T14_A / 2) >= 1:
            keep |= w
    return t[keep], f[keep], dd[keep], ep[keep]


def design(x, ep, s):
    """Columns: per-window offset and slope, then -s (depth)."""
    ue, wid = np.unique(ep, return_inverse=True)
    nw = len(ue)
    X = np.zeros((len(x), 2 * nw + 1))
    X[np.arange(len(x)), wid] = 1.0
    X[np.arange(len(x)), nw + wid] = x
    X[:, -1] = -s
    return X


def fit_depth(sets, shape="tpl", integrate=True, width=None, T12=None, steps=None, sig_fixed=None):
    """sets: list of (sec, t, f, x, ep).  Common depth over all sets; per-sector white sigma from residuals
    (two passes).  steps: list of dump times; adds one free step per dump inside a window."""
    Xs, ys, ws, secs = [], [], [], []
    blocks = []
    for (sec, t, f, x, ep) in sets:
        if shape == "tpl":
            s = np.interp(x, TAU, TPL[sec] if integrate else TPL_NOINT[sec], left=0, right=0)
        elif shape == "box":
            s = np.interp(x, TAU, binned(box_unb(width), sec), left=0, right=0)
        elif shape == "trap":
            s = np.interp(x, TAU, binned(trap_unb(width, T12), sec), left=0, right=0)
        X = design(x, ep, s)
        if steps is not None and sec == 6:
            cols = []
            for td in steps:
                w = np.abs(t - td) < 0.30
                if w.any() and (t[w] < td).any() and (t[w] > td).any():
                    ep_d = ep[np.argmin(np.abs(t - td))]
                    col = ((t > td) & (ep == ep_d)).astype(float)
                    if col.sum() > 0 and (ep == ep_d).sum() > col.sum():
                        cols.append(col)
            if cols:
                X = np.hstack([X[:, :-1], np.array(cols).T, X[:, -1:]])
        blocks.append((sec, X, f))
    # block-diagonal solve: depth common, nuisance per block
    sig = {sec: 1.0 for sec, _, _ in blocks} if sig_fixed is None else dict(sig_fixed)
    for it in range(2 if sig_fixed is None else 1):
        A = 0.0; B = 0.0
        parts = []
        for sec, X, f in blocks:
            w = 1.0 / sig[sec] ** 2
            N = X[:, :-1]; s = X[:, -1]
            # profile the nuisance columns: project them out
            Q, _ = np.linalg.qr(N)
            s_p = s - Q @ (Q.T @ s); f_p = f - Q @ (Q.T @ f)
            A += w * s_p @ s_p; B += w * s_p @ f_p
            parts.append((sec, s_p, f_p, N.shape[1], len(f)))
        D = B / A
        if sig_fixed is not None:
            break
        for sec, s_p, f_p, npar, n in parts:
            r = f_p - D * s_p
            sig[sec] = float(np.sqrt(np.sum(r * r) / max(n - npar - 1, 1)))
        if it == 0:
            continue
        # recompute A, B with the final sigma
        A = 0.0; B = 0.0
        for sec, s_p, f_p, npar, n in parts:
            A += s_p @ s_p / sig[sec] ** 2; B += s_p @ f_p / sig[sec] ** 2
        D = B / A
    eD = float(np.sqrt(1.0 / A))
    chi2 = float(sum(np.sum((f_p - D * s_p) ** 2) / sig[sec] ** 2 for sec, s_p, f_p, _, _ in parts))
    out = dict(D=float(D), eD=eD, sig=sig, chi2=chi2, n=int(sum(p[4] for p in parts)))
    return out


def sets_for(secs, **kw):
    out = []
    for sec in secs:
        t, f, x, ep = collect(sec, **kw)
        out.append((sec, t, f, x, ep))
    return out


def shifted(S, dt_min):
    """Same data points (window membership fixed at dt=0), template centre moved by dt_min."""
    return [(sec, t, f, x - dt_min / 1440.0, ep) for (sec, t, f, x, ep) in S]


def empirical(secs, nfake=NFAKE, **kw):
    vals = []
    ph = rng.uniform(0.08, 0.92, nfake * 3)
    ph = ph[np.abs(ph - 0.5) > 0.08][:nfake]
    for p in ph:
        kw2 = dict(kw); kw2["phase"] = p; kw2["excl_real"] = True
        r = fit_depth(sets_for(secs, **{k: v for k, v in kw2.items() if k in ("phase", "dt_min", "parity", "excl_real", "extra_mask")}),
                      **{k: v for k, v in kw.items() if k in ("shape", "integrate", "width", "T12", "steps")})
        vals.append(r["D"])
    v = np.array(vals)
    return dict(sd=float(np.std(v)), mad_sigma=float(1.4826 * np.median(np.abs(v - np.median(v)))), mean=float(np.mean(v)), n=int(len(v)))


def depth_test(secs, emp=True, **kw):
    ckw = {k: v for k, v in kw.items() if k in ("phase", "dt_min", "parity", "excl_real", "extra_mask")}
    fkw = {k: v for k, v in kw.items() if k in ("shape", "integrate", "width", "T12", "steps")}
    r = fit_depth(sets_for(secs, **ckw), **fkw)
    if emp:
        e = empirical(secs, **kw)
        r["emp"] = e
        r["eD_used"] = max(r["eD"], e["sd"])
    return r


def timing_scan(secs, rng_min=90, step=0.5, extra_mask=None, **kw):
    dts = np.arange(-rng_min, rng_min + 1e-9, step)
    chis, Ds = [], []
    S0 = sets_for(secs, extra_mask=extra_mask)
    sref = fit_depth(S0, **kw)["sig"]
    for dt in dts:
        r = fit_depth(shifted(S0, dt), sig_fixed=sref, **kw)
        chis.append(r["chi2"]); Ds.append(r["D"])
    chis = np.array(chis); j = int(np.argmin(chis))
    ok = np.where(chis - chis[j] <= 1.0)[0]
    return dict(dt_best=float(dts[j]), lo=float(dts[ok[0]]), hi=float(dts[ok[-1]]), e_half=float(0.5 * (dts[ok[-1]] - dts[ok[0]])),
                D_best=float(Ds[j]), dchi2_at_0=float(chis[np.argmin(np.abs(dts))] - chis[j]),
                chi2_curve=[[float(a), float(b)] for a, b in zip(dts[::4], chis[::4] - chis[j])])


def main():
    t0 = time.time()
    R = dict(ephemeris_used=dict(source="s1 S32_S98_sel", P=PE, T0=TE, eP=E3298["eP"], eT0=E3298["eT0"]),
             dumps_S6=DUMPS.tolist())

    # ---------------- (c) depths per sector and S6 vs S32+S98
    dep = {}
    for key, secs in (("S6", (6,)), ("S32", (32,)), ("S98", (98,)), ("S32_S98", (32, 98)), ("all3", (6, 32, 98))):
        dep[key] = depth_test(secs)
        print("depth", key, {k: v for k, v in dep[key].items() if k != "sig"}, f"{time.time()-t0:.0f}s", flush=True)
    d6, d3 = dep["S6"], dep["S32_S98"]
    dep["S6_minus_S32S98_ppm"] = (d6["D"] - d3["D"]) * 1e6
    dep["z_S6_vs_S32S98"] = (d6["D"] - d3["D"]) / np.hypot(d6["eD_used"], d3["eD_used"])
    dep["chi2_three_sectors"] = float(sum(((dep[k]["D"] - dep["all3"]["D"]) / dep[k]["eD_used"]) ** 2 for k in ("S6", "S32", "S98")))
    R["depths"] = dep

    # ---------------- (c) S6 timing relative to the S32+S98 ephemeris (sector-level, depth free)
    R["S6_timing"] = timing_scan((6,))
    # injection-based sector-level error: inject the template (S6 depth) at the fake ephemeris and recover the common dt
    errs = []
    Dinj = d6["D"]
    base_t, base_f = LC[6]["t"].copy(), LC[6]["f"].copy()
    for p in rng.uniform(0.1, 0.9, 600):
        if abs(p - 0.5) < 0.12 or len(errs) >= 200:
            continue
        u = rng.uniform(-10, 10)
        c = TE + p * PE + u / 1440.0
        dd, _ = phase_dt(base_t, PE, c)
        dA, _ = phase_dt(base_t, PE, TE)
        LC[6]["f"] = base_f * (1 - Dinj * np.interp(dd, TAU, TPL[6], left=0, right=0))
        dts = np.arange(-60, 60.01, 1.0)
        ch = []
        S0 = sets_for((6,), phase=p, excl_real=True)
        sref = fit_depth(S0)["sig"]
        for dt in dts:
            ch.append(fit_depth(shifted(S0, dt), sig_fixed=sref)["chi2"])
        errs.append(dts[int(np.argmin(ch))] - u)
    LC[6]["f"] = base_f
    errs = np.array(errs)
    R["S6_timing"]["inj_sector_level"] = dict(n=int(len(errs)), sd_min=float(np.std(errs)),
                                              mad_sigma_min=float(1.4826 * np.median(np.abs(errs - np.median(errs)))),
                                              hw68_min=float(0.5 * (np.percentile(errs, 84.135) - np.percentile(errs, 15.865))))
    E6 = np.array([-795, -794, -793, -791, -790, -789], float)
    sp = np.sqrt(E3298["eT0"] ** 2 + E6 ** 2 * E3298["eP"] ** 2 + 2 * E6 * E3298["cov"]) * 1440
    R["S6_timing"]["pred_sigma_min"] = float(np.mean(sp))
    e_meas = max(R["S6_timing"]["e_half"], R["S6_timing"]["inj_sector_level"]["sd_min"])
    R["S6_timing"]["z"] = float(R["S6_timing"]["dt_best"] / np.hypot(e_meas, np.mean(sp)))
    print("S6 timing", {k: v for k, v in R["S6_timing"].items() if k != "chi2_curve"}, f"{time.time()-t0:.0f}s", flush=True)

    # ---------------- (c) momentum-dump robustness of the S6 photometry
    dm = {}
    for hw in (0.5, 1.0, 1.5):
        msk = lambda t, hw=hw: np.any(np.abs(t[:, None] - DUMPS[None, :]) < hw / 24.0, axis=1)
        r = fit_depth(sets_for((6,), extra_mask=msk))
        # timing with the mask
        dts = np.arange(-60, 60.01, 1.0); ch = []
        S0 = sets_for((6,), extra_mask=msk)
        for dt in dts:
            ch.append(fit_depth(shifted(S0, dt), sig_fixed=r["sig"])["chi2"])
        ch = np.array(ch); j = int(np.argmin(ch)); ok = np.where(ch - ch[j] <= 1)[0]
        dm[f"mask_pm{hw}h"] = dict(D=r["D"], eD=r["eD"], n=r["n"], dt_best=float(dts[j]), e_half=float(0.5 * (dts[ok[-1]] - dts[ok[0]])))
    r = fit_depth(sets_for((6,)), steps=DUMPS)
    dts = np.arange(-60, 60.01, 1.0); ch = []
    S0 = sets_for((6,))
    for dt in dts:
        ch.append(fit_depth(shifted(S0, dt), steps=DUMPS, sig_fixed=r["sig"])["chi2"])
    ch = np.array(ch); j = int(np.argmin(ch)); ok = np.where(ch - ch[j] <= 1)[0]
    dm["free_step_at_each_dump"] = dict(D=r["D"], eD=r["eD"], dt_best=float(dts[j]), e_half=float(0.5 * (dts[ok[-1]] - dts[ok[0]])))
    # out-of-transit flux steps at the dumps (PDCSAP): linear baseline + step, +-0.5 d, A and B in-transit removed
    steps = []
    t6, f6 = LC[6]["t"], LC[6]["f"]
    dA, _ = phase_dt(t6, PE, TE)
    for td in DUMPS:
        w = (np.abs(t6 - td) < 0.5) & (np.abs(dA) > T14_A)
        if (t6[w] < td).sum() < 5 or (t6[w] > td).sum() < 5:
            steps.append(dict(t=float(td), step_ppm=None)); continue
        x = t6[w] - td
        X = np.vstack([np.ones_like(x), x, (x > 0).astype(float)]).T
        c, *_ = np.linalg.lstsq(X, f6[w], rcond=None)
        r_ = f6[w] - X @ c; s2 = np.sum(r_ ** 2) / (len(r_) - 3)
        cov = np.linalg.inv(X.T @ X) * s2
        steps.append(dict(t=float(td), step_ppm=float(c[2] * 1e6), e_step_ppm=float(np.sqrt(cov[2, 2]) * 1e6),
                          nearest_A_dt_h=float((td - (TE + np.round((td - TE) / PE) * PE)) * 24)))
    dm["oot_steps_at_dumps"] = steps
    # dump-locked alternative: would a signal locked to the dumps (constant dump-to-dip lag) fit S6 as well?
    dmid = DUMPS
    lags = np.arange(-4.0, 4.001, 0.05)
    chil = []
    for lag in lags:
        cs = dmid + lag / 24.0
        x = t6 - cs[np.argmin(np.abs(t6[:, None] - cs[None, :]), axis=1)]
        m = np.abs(x) < 0.30
        ep = np.argmin(np.abs(t6[:, None] - cs[None, :]), axis=1)
        r_ = fit_depth([(6, t6[m], f6[m], x[m], ep[m])])
        chil.append((r_["chi2"], r_["D"]))
    chil = np.array(chil)
    jl = int(np.argmin(chil[:, 0]))
    r_lin = fit_depth(sets_for((6,)))
    # same sigma normalisation: recompute linear-ephemeris chi2 with the dump-locked sigma
    snr_l = []
    for lag in lags:
        cs = dmid + lag / 24.0
        idx = np.argmin(np.abs(t6[:, None] - cs[None, :]), axis=1)
        x = t6 - cs[idx]; m = np.abs(x) < 0.30
        r_ = fit_depth([(6, t6[m], f6[m], x[m], idx[m])])
        snr_l.append(r_["D"] / r_["eD"])
    snr_l = np.array(snr_l); jl = int(np.argmax(snr_l))
    dm["dump_locked_model"] = dict(best_lag_h=float(lags[jl]), best_snr=float(snr_l[jl]), D_ppm_at_best=float(chil[jl, 1] * 1e6),
                                   snr_linear_ephemeris=float(r_lin["D"] / r_lin["eD"]), D_ppm_linear=float(r_lin["D"] * 1e6),
                                   dump_spacing_d=float(np.median(np.diff(DUMPS))),
                                   note="best S/N of a dip at a fixed lag after each dump (lags -4..+4 h) vs the S/N of the S32+S98 linear ephemeris")
    R["S6_dumps"] = dm
    print("dumps", json.dumps(dm, default=float)[:1500], f"{time.time()-t0:.0f}s", flush=True)

    # ---------------- (d) cadence integration of the 30-min S6 data
    r_int = fit_depth(sets_for((6,)), integrate=True)
    r_noint = fit_depth(sets_for((6,)), integrate=False)
    R["S6_integration"] = dict(D_integrated_ppm=r_int["D"] * 1e6, D_unintegrated_ppm=r_noint["D"] * 1e6,
                               ratio=r_noint["D"] / r_int["D"], chi2_integrated=r_int["chi2"], chi2_unintegrated=r_noint["chi2"],
                               template_peak_binned=float(TPL[6].max()), template_peak_unbinned=float(S_UNB.max()))
    print("integration", R["S6_integration"], flush=True)

    # ---------------- odd/even and phase 0.5, with and without S6
    oe = {}
    for key, secs in (("all3", (6, 32, 98)), ("S32_S98", (32, 98))):
        od = depth_test(secs, parity=1); ev = depth_test(secs, parity=0)
        oe[key] = dict(odd_ppm=od["D"] * 1e6, e_odd_ppm=od["eD_used"] * 1e6, even_ppm=ev["D"] * 1e6, e_even_ppm=ev["eD_used"] * 1e6,
                       z=float((od["D"] - ev["D"]) / np.hypot(od["eD_used"], ev["eD_used"])))
        print("odd/even", key, oe[key], flush=True)
    R["odd_even"] = oe
    sec05 = {}
    for key, secs in (("all3", (6, 32, 98)), ("S32_S98", (32, 98))):
        r = depth_test(secs, phase=0.5, shape="box", width=T14_A)
        sec05[key] = dict(depth_ppm=r["D"] * 1e6, e_ppm=r["eD_used"] * 1e6, flux_ppm=-r["D"] * 1e6)
        print("phase0.5", key, sec05[key], flush=True)
    R["phase05"] = sec05

    # ---------------- V shape: trapezoid (free T14, T12) vs triangle (T12 = T14/2), with and without S6
    vs = {}
    for key, secs in (("all3", (6, 32, 98)), ("S32_S98", (32, 98))):
        S = sets_for(secs)
        best_tr, best_v = (np.inf, None), (np.inf, None)
        sig_ref = fit_depth(S)["sig"]
        for T14h in np.arange(1.5, 2.6, 0.02):
            T14d = T14h / 24
            c = _chi_fixed(S, "trap", T14d, T14d / 2, sig_ref)
            if c < best_v[0]:
                best_v = (c, T14h)
            for f12 in np.arange(0.02, 0.5, 0.02):
                c = _chi_fixed(S, "trap", T14d, f12 * T14d, sig_ref)
                if c < best_tr[0]:
                    best_tr = (c, (T14h, f12))
        vs[key] = dict(chi2_trap=best_tr[0], trap_T14_h=best_tr[1][0], trap_T12_over_T14=best_tr[1][1], chi2_V=best_v[0],
                       V_T14_h=best_v[1], dchi2=best_v[0] - best_tr[0], sigma_ppm={str(k): v * 1e6 for k, v in sig_ref.items()})
        print("vshape", key, vs[key], f"{time.time()-t0:.0f}s", flush=True)
    R["vshape"] = vs
    R["runtime_s"] = time.time() - t0
    json.dump(R, open(os.path.join(OUT, "s3_s6_checks.json"), "w"), indent=1, default=float)


_QCACHE = {}


def _chi_fixed(S, shape, width, T12, sig):
    """chi2 with fixed per-sector sigma (from the reference template fit), depth and baselines profiled."""
    A = 0.0; B = 0.0; C = 0.0
    for (sec, t, f, x, ep) in S:
        s = np.interp(x, TAU, binned(trap_unb(width, T12), sec), left=0, right=0)
        key = (sec, len(t), float(t[0]))
        if key not in _QCACHE:
            X = design(x, ep, s)
            _QCACHE[key] = np.linalg.qr(X[:, :-1])[0]
        Q = _QCACHE[key]
        s_p = s - Q @ (Q.T @ s); f_p = f - Q @ (Q.T @ f)
        w = 1 / sig[sec] ** 2
        A += w * s_p @ s_p; B += w * s_p @ f_p; C += w * f_p @ f_p
    return float(C - B * B / A)


if __name__ == "__main__":
    main()
