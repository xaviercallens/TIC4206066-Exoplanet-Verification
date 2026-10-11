"""R-S6 (a) and (c)-timing: per-transit times of Signal A with a fixed-shape template, injection-calibrated errors,
linear ephemerides with and without Sector 6, the S32<->S98 cycle count, and the S6 O-C against the S32+S98 ephemeris.

Per transit: raw normalised flux in +-0.3 d windows (B removed), model f = c0 + c1 (t - tc) - D s(t - tc - dt), where s is
the unit-depth template (k=0.0212, b=0.35, rho=2.51, Claret LD) integrated over the cadence (boxcar of TIMEDEL; S6 1800 s).
dt on a 0.25-min grid (+-60 min, as in the manuscript, and +-144 min as a check); c0, c1, D linear at each dt.
Selection (manuscript rule): in-transit coverage >= 0.7, S/N >= 2, Delta-chi2=1 interval not touching the scan edge.
Errors: per-sector SD of (recovered - injected) for 600 injections of the template at random out-of-transit times of the
same light curve (same window, coverage and S/N selection), depth = sector mean depth.
usage: s1_timing.py [RANGE_MIN=60] [NINJ=600]
"""
import sys, json, time
import numpy as np
import batman
from rs6_common import *

RANGE = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
NINJ = int(sys.argv[2]) if len(sys.argv) > 2 else 600
rng = np.random.default_rng(20260925)
U = claret_u()
SHAPE = dict(k=0.0212, b=0.35, rho=2.51)
DTS = np.arange(-RANGE, RANGE + 1e-9, 0.25) / 1440.0


def unit_template():
    tau = np.arange(-0.40, 0.40 + 1e-9, 2.0 / 86400.0)
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp = 0.0, P_A, SHAPE["k"]
    pm.a = ars_from_rho(SHAPE["rho"], P_A)
    pm.inc = np.degrees(np.arccos(SHAPE["b"] / pm.a)); pm.ecc, pm.w = 0.0, 90.0
    pm.u, pm.limb_dark = list(U), "quadratic"
    fl = batman.TransitModel(pm, tau).light_curve(pm)
    s = (1 - fl) / (1 - fl[np.argmin(np.abs(tau))])      # unit depth at mid-transit (unbinned)
    out = {}
    for sec, cad in CAD_S.items():
        n = int(round(cad / 2.0))
        ker = np.ones(n) / n
        out[sec] = (tau, np.convolve(s, ker, mode="same"))   # boxcar over the cadence (centred: TIMEPIXR=0.5)
    return out


TPL = unit_template()


def tmpl(sec, x):
    tau, s = TPL[sec]
    return np.interp(x, tau, s, left=0.0, right=0.0)


def scan(t, f, tc, sec, dts=DTS):
    """Vectorised linear fits over the dt grid.  Returns dict with best dt, D, sigma_D, S/N, edge flag, chi2 curve."""
    x = t - tc
    S = np.array([tmpl(sec, x - dt) for dt in dts])            # (ndt, n)
    n = len(t)
    one = np.ones(n)
    # normal equations for [c0, c1, -D]
    A = np.zeros((len(dts), 3, 3)); B = np.zeros((len(dts), 3))
    A[:, 0, 0] = n; A[:, 0, 1] = A[:, 1, 0] = x.sum(); A[:, 1, 1] = (x * x).sum()
    A[:, 0, 2] = A[:, 2, 0] = -S.sum(1); A[:, 1, 2] = A[:, 2, 1] = -(S @ x); A[:, 2, 2] = (S * S).sum(1)
    B[:, 0] = f.sum(); B[:, 1] = (x * f).sum(); B[:, 2] = -(S @ f)
    try:
        c = np.linalg.solve(A, B[..., None])[..., 0]
    except np.linalg.LinAlgError:
        return None
    model = c[:, 0:1] + c[:, 1:2] * x[None, :] - c[:, 2:3] * S
    chi = ((f[None, :] - model) ** 2).sum(1)
    j = int(np.argmin(chi))
    r = f - model[j]
    s2 = np.sum(r * r) / max(n - 4, 1)
    chir = chi / s2
    ok = np.where(chir - chir[j] <= 1.0)[0]
    edge = bool(ok[0] == 0 or ok[-1] == len(dts) - 1)
    cov = np.linalg.inv(A[j]) * s2
    D, eD = c[j, 2], float(np.sqrt(cov[2, 2]))
    return dict(dt=float(dts[j]), D=float(D), eD=eD, snr=float(D / eD), edge=edge,
                e_dt_white=float(0.5 * (dts[ok[-1]] - dts[ok[0]])), rms=float(np.sqrt(s2)), n=n)


def real_events():
    rows = []
    for sec in (6, 32, 98):
        t, f, ep, meta = windows(sec)
        for e in np.unique(ep):
            m = ep == e
            tc = T0_A + e * P_A
            cov, npre, npost = coverage(t[m], tc, CAD_S[sec])
            if m.sum() < 10:
                continue
            r = scan(t[m], f[m], tc, sec)
            if r is None:
                continue
            r.update(sector=sec, epoch=int(e), tc_pred=float(tc), coverage=cov, n_pre=npre, n_post=npost)
            r["sel"] = bool(cov >= 0.7 and r["snr"] >= 2 and not r["edge"])
            rows.append(r)
    return rows


def injections(sec, depth):
    t, f, meta = load(sec)
    dA, _ = phase_dt(t, P_A, T0_A)
    dB, _ = phase_dt(t, P_B, T0_B)
    keep = (np.abs(dB) > T14_B)
    t, f, dA = t[keep], f[keep], dA[keep]
    res = []
    tries = 0
    while len(res) < NINJ and tries < 50 * NINJ:
        tries += 1
        ti = rng.uniform(t.min() + 0.1, t.max() - 0.1)
        # window must not contain any real A in-transit cadence (A transits +-1 T14 excluded)
        w = np.abs(t - ti) < 0.30
        if np.any(np.abs(dA[w]) < T14_A):
            continue
        tw, fw = t[w], f[w]
        if len(tw) < 10:
            continue
        cov, npre, npost = coverage(tw, ti, CAD_S[sec])
        if cov < 0.7:
            continue
        u = rng.uniform(-10, 10) / 1440.0
        fi = fw * (1.0 - depth * tmpl(sec, tw - ti - u))
        r = scan(tw, fi, ti, sec)
        if r is None:
            continue
        sel = r["snr"] >= 2 and not r["edge"]
        res.append(dict(err_min=(r["dt"] - u) * 1440.0, sel=bool(sel), D=r["D"], eD=r["eD"], snr=r["snr"],
                        n_pre=npre, n_post=npost))
    return res


def wls(ev, sig_by_sec, ref=0):
    E = np.array([e["epoch"] for e in ev], float) - ref
    tobs = np.array([e["tc_pred"] + e["dt"] for e in ev])
    s = np.array([sig_by_sec[e["sector"]] for e in ev]) / 1440.0
    X = np.vstack([np.ones_like(E), E]).T
    W = 1 / s ** 2
    C = np.linalg.inv(X.T @ (X * W[:, None]))
    p = C @ (X.T @ (W * tobs))
    r = tobs - X @ p
    chi2 = float(np.sum(r ** 2 * W))
    return dict(T0=float(p[0]), eT0=float(np.sqrt(C[0, 0])), P=float(p[1]), eP=float(np.sqrt(C[1, 1])),
                cov=float(C[0, 1]), chi2=chi2, dof=int(len(ev) - 2), n=int(len(ev)), ref_epoch=ref,
                rms_min=float(np.std(r) * 1440), oc_min=(r * 1440).tolist(),
                epochs=[int(e["epoch"]) for e in ev], sectors=[int(e["sector"]) for e in ev])


def main():
    t0 = time.time()
    ev = real_events()
    for r in ev:
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    sel = [r for r in ev if r["sel"]]
    print("selected", len(sel), {s: sum(r["sector"] == s for r in sel) for s in (6, 32, 98)})
    # injections at each sector's mean depth of the selected events
    inj = {}
    sig = {}
    for sec in (6, 32, 98):
        ds = np.array([r["D"] for r in sel if r["sector"] == sec]); es = np.array([r["eD"] for r in sel if r["sector"] == sec])
        dep = float(np.sum(ds / es ** 2) / np.sum(1 / es ** 2))
        res = injections(sec, dep)
        e = np.array([x["err_min"] for x in res if x["sel"]])
        e_bb = np.array([x["err_min"] for x in res if x["sel"] and x["n_pre"] > 0 and x["n_post"] > 0])
        inj[sec] = dict(depth_injected=dep, n_inj=len(res), n_sel=int(len(e)), sd=float(np.std(e)),
                        mad_sigma=float(1.4826 * np.median(np.abs(e - np.median(e)))),
                        hw68=float(0.5 * (np.percentile(e, 84.135) - np.percentile(e, 15.865))),
                        frac_gt30=float(np.mean(np.abs(e) > 30)), mean=float(np.mean(e)),
                        sd_two_sided_baseline=float(np.std(e_bb)), sel_frac=float(len(e) / len(res)))
        sig[sec] = inj[sec]["sd"]
        print("inj", sec, inj[sec], f"{time.time()-t0:.0f}s", flush=True)
        np.savetxt(os.path.join(OUT, f"s1_inj_S{sec}_R{int(RANGE)}.txt"),
                   np.array([[x["err_min"], x["sel"], x["snr"], x["n_pre"], x["n_post"]] for x in res]),
                   header="err_min sel snr n_pre n_post")

    R = dict(range_min=RANGE, template=dict(SHAPE, u=U.tolist()), events=ev, injections=inj, sigma_used_min=sig)
    by = lambda secs, extra=lambda r: True: [r for r in sel if r["sector"] in secs and extra(r)]
    fits_ = {}
    fits_["all3_sel"] = wls(by((6, 32, 98)), sig)
    fits_["S32_S98_sel"] = wls(by((32, 98)), sig)
    fits_["S98_sel"] = wls(by((98,)), sig, ref=7)
    fits_["S32_sel"] = wls(by((32,)), sig, ref=-570)
    fits_["S6_sel"] = wls(by((6,)), sig, ref=-792)
    two = lambda r: r["n_pre"] > 0 and r["n_post"] > 0
    fits_["all3_sel_twosided"] = wls(by((6, 32, 98), two), sig)
    fits_["S32_S98_sel_twosided"] = wls(by((32, 98), two), sig)
    # S32+S98 with robust (MAD) sigmas
    sig_mad = {s: inj[s]["mad_sigma"] for s in inj}
    fits_["S32_S98_sel_madsig"] = wls(by((32, 98)), sig_mad)
    fits_["all3_sel_madsig"] = wls(by((6, 32, 98)), sig_mad)
    # decorrelated reference epoch for S32+S98 (minimises cov(T0,P))
    f = fits_["S32_S98_sel"]
    Eref = int(round(-f["cov"] / f["eP"] ** 2))
    fits_["S32_S98_sel_decorr"] = wls(by((32, 98)), sig, ref=Eref)
    R["ephemerides"] = fits_

    # ---- comparison with the manuscript and nested-difference significance
    Pm, ePm = 3.1827856, 0.0000058
    a, s = fits_["all3_sel"], fits_["S32_S98_sel"]
    R["compare"] = dict(
        manuscript=dict(P=Pm, eP=ePm, T0_BJD=2460998.7375, eT0=0.0023),
        S32_S98_minus_manuscript_P=s["P"] - Pm,
        z_vs_manuscript_naive=(s["P"] - Pm) / np.hypot(s["eP"], ePm),
        S32_S98_minus_all3_P=s["P"] - a["P"],
        z_nested=(s["P"] - a["P"]) / np.sqrt(max(s["eP"] ** 2 - a["eP"] ** 2, 1e-30)),
        S32_S98_minus_all3_T0_min=(s["T0"] - a["T0"]) * 1440,
        z_nested_T0=(s["T0"] - a["T0"]) / np.sqrt(max(s["eT0"] ** 2 - a["eT0"] ** 2, 1e-30)),
        eP_ratio=s["eP"] / a["eP"], eT0_ratio=s["eT0"] / a["eT0"],
        S32_S98_T0_BJD=s["T0"] + 2457000.0)

    # ---- cycle count S32 <-> S98 using only in-sector information
    f98, f32 = fits_["S98_sel"], fits_["S32_sel"]
    T98 = f98["T0"] + 0.0; E98 = f98["ref_epoch"]
    T32 = f32["T0"]; E32 = f32["ref_epoch"]
    Pw = (f98["P"] / f98["eP"] ** 2 + f32["P"] / f32["eP"] ** 2) / (1 / f98["eP"] ** 2 + 1 / f32["eP"] ** 2)
    ePw = 1 / np.sqrt(1 / f98["eP"] ** 2 + 1 / f32["eP"] ** 2)
    N0 = int(round((T98 - T32) / Pw))
    cyc = []
    for j in range(-3, 4):
        Pj = (T98 - T32) / (N0 + j)
        ePj = np.hypot(f98["eT0"], f32["eT0"]) / (N0 + j)
        cyc.append(dict(j=j, N=N0 + j, P=Pj, eP=ePj, z_vs_insector=(Pj - Pw) / np.hypot(ePw, ePj),
                        z_vs_S98only=(Pj - f98["P"]) / f98["eP"]))
    R["cycle_count"] = dict(T32_BTJD=T32, E32=E32, eT32_min=f32["eT0"] * 1440, T98_BTJD=T98, E98=E98,
                            eT98_min=f98["eT0"] * 1440, N_expected=E98 - E32, N0=N0,
                            P_insector_weighted=Pw, eP_insector=ePw, P_S98only=f98["P"], eP_S98only=f98["eP"],
                            P_S32only=f32["P"], eP_S32only=f32["eP"], alias_spacing_d=Pw / N0,
                            alias_spacing_over_sigma=(Pw / N0) / ePw,
                            S98only_extrap_to_S32_sigma_d=(E98 - E32) * f98["eP"],
                            aliases=cyc)

    # ---- S6 O-C against the S32+S98 ephemeris (independent of S6)
    s6 = by((6,))
    E = np.array([r["epoch"] for r in s6], float)
    tobs = np.array([r["tc_pred"] + r["dt"] for r in s6])
    tpred = s["T0"] + E * s["P"]
    sp = np.sqrt(s["eT0"] ** 2 + E ** 2 * s["eP"] ** 2 + 2 * E * s["cov"])     # prediction sigma (correlated)
    oc = (tobs - tpred) * 1440
    w = 1 / (sig[6] ** 2)
    mean_oc = float(np.mean(oc)); e_mean = sig[6] / np.sqrt(len(oc))
    sp_mean = float(np.mean(sp) * 1440)
    R["S6_vs_S32S98"] = dict(n=int(len(oc)), epochs=E.astype(int).tolist(), oc_min=oc.tolist(),
                             sigma_event_min=sig[6], pred_sigma_min=(sp * 1440).tolist(),
                             mean_oc_min=mean_oc, e_mean_oc_min=float(e_mean), pred_sigma_mean_min=sp_mean,
                             z_mean=float(mean_oc / np.hypot(e_mean, sp_mean)),
                             chi2_events=float(np.sum(oc ** 2 / (sig[6] ** 2 + (sp * 1440) ** 2))),
                             dof=int(len(oc)))
    # all S6 events with data (not only selected)
    s6all = [r for r in ev if r["sector"] == 6]
    R["S6_events_all"] = [dict(epoch=r["epoch"], coverage=r["coverage"], n_pre=r["n_pre"], n_post=r["n_post"], snr=r["snr"],
                               dt_min=r["dt"] * 1440, oc_vs_S32S98_min=(r["tc_pred"] + r["dt"] - (s["T0"] + r["epoch"] * s["P"])) * 1440,
                               sel=r["sel"], edge=r["edge"]) for r in s6all]
    R["runtime_s"] = time.time() - t0
    json.dump(R, open(os.path.join(OUT, f"s1_timing_R{int(RANGE)}.json"), "w"), indent=1, default=float)
    for k_, v in fits_.items():
        print(k_, {kk: v[kk] for kk in ("n", "P", "eP", "T0", "eT0", "chi2", "dof")})
    print("compare", R["compare"])
    print("cycle", {k_: v for k_, v in R["cycle_count"].items() if k_ != "aliases"})
    for c in cyc:
        print("  alias", c)
    print("S6 vs S32+S98", {k_: v for k_, v in R["S6_vs_S32S98"].items()})


if __name__ == "__main__":
    main()
