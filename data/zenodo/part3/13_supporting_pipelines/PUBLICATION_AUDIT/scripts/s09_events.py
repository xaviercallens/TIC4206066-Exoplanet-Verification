"""Step 9: event-by-event audit for both signals, in every data set (pipeline products and own photometry).

For each predicted transit n (tc = T0 + n P):
  coverage   : (in-transit cadences x cadence) / T14
  depth D    : linear least squares  f = c0 + c1 (t - tc) - D s(t - tc), s = transit template with unit central depth
               (template: batman model from the joint fit with the Mann density prior); errors from the covariance
               scaled by sqrt(chi2_red) and by the red-noise factor beta (time-averaging method, out-of-transit part)
  timing     : scan dt in [-60, 60] min (0.5-min grid); chi2(dt) with D free; error from Delta chi2 = 1 after
               rescaling chi2 by chi2_red*beta^2; events with coverage < 0.7 or S/N < 2 get no timing
  local noise: robust rms of out-of-transit residuals and binned rms on the T14 time scale
Aggregates: weighted mean depth per data set/sector, chi2 of a constant depth, share of the total (S/N)^2 carried by
the one and two strongest events, leave-one-out S/N, linear-ephemeris fit to the event times (O-C) with chi2.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import transit_tools as tt  # noqa
from audit_common import AUD, save_json  # noqa

EPH = {k: dict(v) for k, v in lcio.EPH.items()}
ephf = os.path.join(AUD, "out", "ephem_current.json")
if os.path.exists(ephf):
    EPH.update(json.load(open(ephf)))


def template_params(sig):
    f = os.path.join(AUD, "out", "fits", f"fit_{sig}_pipe_mann_frozen_w0.75.json")
    if os.path.exists(f):
        r = json.load(open(f))
        return dict(k=r["k"][2], b=r["b"][2], rho=r["rho"][2], u1=r["u1"][2], u2=r["u2"][2], src=os.path.basename(f))
    return dict(k=0.0212 if sig == "A" else 0.028, b=0.36 if sig == "A" else 0.8, rho=2.5, u1=0.33, u2=0.30, src="prelim")


def unit_template(sig, t_rel, cad):
    p = template_params(sig)
    e = EPH[sig]
    m = tt.Model(t_rel, cad)
    ars = tt.ars_from_rho(p["rho"], e["P"])
    lc = m.lc(0.0, e["P"], p["k"], p["b"], ars, p["u1"], p["u2"])
    mc = tt.Model(np.array([0.0]), 1e-4).lc(0.0, e["P"], p["k"], p["b"], ars, p["u1"], p["u2"])[0]
    return (1 - lc) / (1 - mc)   # unit depth at mid-transit (for an infinitely short exposure)


def fit_event(t, f, tc, sig, cad, T14):
    trel = t - tc
    s0 = unit_template(sig, trel, cad)
    oot = np.abs(trel) > 0.75 * T14
    X = np.vstack([np.ones_like(trel), trel, -s0]).T
    c, *_ = np.linalg.lstsq(X, f, rcond=None)
    r = f - X @ c
    dof = max(len(f) - 3, 1)
    chi2red = np.sum(r ** 2) / dof / np.var(r[oot]) if oot.sum() > 5 else 1.0
    s2 = np.var(r)
    cov = np.linalg.pinv(X.T @ X) * s2
    beta = tt.beta_red(r[oot], cad) if oot.sum() > 20 else 1.0
    D, eD = c[2], np.sqrt(cov[2, 2]) * beta
    # timing scan
    dts = np.arange(-60, 60.01, 0.5) / 1440.0
    chis = []
    for dt in dts:
        s = unit_template(sig, trel - dt, cad)
        Xd = np.vstack([np.ones_like(trel), trel, -s]).T
        cc, *_ = np.linalg.lstsq(Xd, f, rcond=None)
        chis.append(np.sum((f - Xd @ cc) ** 2))
    chis = np.array(chis) / (s2 * beta ** 2)
    j = int(np.argmin(chis))
    dchi = chis - chis[j]
    ok = np.where(dchi <= 1.0)[0]
    lo, hi = dts[ok[0]], dts[ok[-1]]
    edge = (j == 0) or (j == len(dts) - 1) or ok[0] == 0 or ok[-1] == len(dts) - 1
    return dict(depth_ppm=float(D * 1e6), depth_err_ppm=float(eD * 1e6), snr=float(D / eD), beta=float(beta),
                rms_oot_ppm=float(np.std(r[oot]) * 1e6) if oot.sum() > 2 else None, n=int(len(f)),
                dt_min=float(dts[j] * 1440), dt_lo_min=float(lo * 1440), dt_hi_min=float(hi * 1440),
                e_dt_min=float(0.5 * (hi - lo) * 1440), timing_at_scan_edge=bool(edge), chi2red=float(chi2red))


def analyse(sig, lcs, label):
    e = EPH[sig]
    other = "B" if sig == "A" else "A"
    T14 = e["T14"]
    rows = []
    for lc in lcs:
        t, f = lc["t"], lc["fd"]
        n0 = int(np.floor((t.min() - e["T0"]) / e["P"])) - 1
        n1 = int(np.ceil((t.max() - e["T0"]) / e["P"])) + 1
        for n in range(n0, n1 + 1):
            tc = e["T0"] + n * e["P"]
            w = (np.abs(t - tc) < 0.35) & ~lcio.in_transit(t, other, 1.0, EPH)
            inn = w & (np.abs(t - tc) < T14 / 2)
            cov = inn.sum() * lc["cad"] / T14
            if tc < t.min() - 0.2 or tc > t.max() + 0.2:
                continue
            row = dict(dataset=lc["key"], sector=lc["sector"], epoch=n, tc_pred=float(tc), coverage=float(min(cov, 1.5)),
                       n_window=int(w.sum()),
                       other_signal_within_12h=bool(np.min(np.abs(((tc - EPH[other]["T0"]) / EPH[other]["P"] + 0.5) % 1 - 0.5) * EPH[other]["P"]) < 0.5))
            if cov >= 0.5 and w.sum() > 20:
                row.update(fit_event(t[w], f[w], tc, sig, lc["cad"], T14))
            rows.append(row)
    return rows


def summarise(rows, key):
    R = [r for r in rows if r["dataset"] == key and "depth_ppm" in r]
    if not R:
        return None
    d = np.array([r["depth_ppm"] for r in R])
    e = np.array([r["depth_err_ppm"] for r in R])
    w = 1 / e ** 2
    mu = float(np.sum(w * d) / w.sum())
    emu = float(1 / np.sqrt(w.sum()))
    chi2 = float(np.sum(((d - mu) / e) ** 2))
    snr2 = (d / e) ** 2
    order = np.argsort(snr2)[::-1]
    tot = float(np.sqrt(np.sum(w * d) ** 2 / np.sum(w)))  # S/N of weighted mean
    loo = []
    for i in range(len(R)):
        m = np.ones(len(R), bool)
        m[i] = False
        loo.append(float(np.sum(w[m] * d[m]) / np.sqrt(np.sum(w[m]))))
    return dict(n_events=len(R), mean_depth_ppm=mu, err_ppm=emu, snr_total=tot, chi2_const=chi2, dof=len(R) - 1,
                share_top1=float(snr2[order[0]] / snr2.sum()), share_top2=float(snr2[order[:2]].sum() / snr2.sum()),
                loo_snr_min=float(min(loo)), loo_snr_max=float(max(loo)), epochs=[r["epoch"] for r in R])


def linear_ephem(rows, keys, min_cov=0.7, min_snr=2.0):
    R = [r for r in rows if r["dataset"] in keys and "dt_min" in r and r["coverage"] >= min_cov and r["snr"] >= min_snr
         and not r["timing_at_scan_edge"]]
    if len(R) < 3:
        return None
    n = np.array([r["epoch"] for r in R], float)
    tc = np.array([r["tc_pred"] + r["dt_min"] / 1440 for r in R])
    et = np.array([max(r["e_dt_min"], 0.5) / 1440 for r in R])
    A = np.vstack([np.ones_like(n), n]).T
    Wm = np.diag(1 / et ** 2)
    cov = np.linalg.inv(A.T @ Wm @ A)
    c = cov @ A.T @ Wm @ tc
    oc = (tc - A @ c) * 1440
    chi2 = float(np.sum((oc / (et * 1440)) ** 2))
    dof = len(R) - 2
    scale = np.sqrt(max(chi2 / dof, 1.0))
    return dict(n=len(R), T0=float(c[0]), P=float(c[1]), eT0_formal=float(np.sqrt(cov[0, 0])), eP_formal=float(np.sqrt(cov[1, 1])),
                eT0_scaled=float(np.sqrt(cov[0, 0]) * scale), eP_scaled=float(np.sqrt(cov[1, 1]) * scale), chi2=chi2, dof=dof,
                oc_min=[float(x) for x in oc], err_min=[float(x * 1440) for x in et], epochs=[int(x) for x in n],
                datasets=[r["dataset"] for r in R], rms_oc_min=float(np.std(oc)))


if __name__ == "__main__":
    SETS = {
        "pipe": [lcio.load_product("TESSSPOC_S6_PDC"), lcio.load_product("SPOC2min_S32_PDC"), lcio.load_product("QLP_S98_SYS", "frozen")],
        "pipe_s98default": [lcio.load_product("QLP_S98_SYS", "default")],
        "own": [lcio.load_own("tesscut_s6", 2.5), lcio.load_own("tesscut_s32", 2.0), lcio.load_own("tesscut_s98", 2.0)],
        "own_nostray": [lcio.load_own("tesscut_s98", 2.0, drop_stray=True)],
        "own2min": [lcio.load_own("spoc_tpf_s32", 2.0)],
        "others": [lcio.load_product("QLP_S6"), lcio.load_product("QLP_S32"), lcio.load_product("TESSSPOC_S32_PDC"),
                   lcio.load_product("eleanor_S6"), lcio.load_product("TGLC_S6"), lcio.load_product("SPOC2min_S32_SAP"),
                   lcio.load_product("TESSSPOC_S6_SAP")],
    }
    out = {"ephemeris_used": EPH, "template": {s: template_params(s) for s in ("A", "B")}}
    for sig in ("A", "B"):
        out[sig] = {}
        for sname, lcs in SETS.items():
            d = [lcio.detrend(lc, window=0.75, eph=EPH) for lc in lcs]
            rows = analyse(sig, d, sname)
            summ = {lc["key"]: summarise(rows, lc["key"]) for lc in d}
            out[sig][sname] = dict(rows=rows, summary=summ)
            for k, v in summ.items():
                if v:
                    print(f"{sig} {sname:16s} {k:26s} n={v['n_events']:2d} depth={v['mean_depth_ppm']:6.0f}±{v['err_ppm']:4.0f} "
                          f"S/N={v['snr_total']:5.1f} chi2={v['chi2_const']:5.1f}/{v['dof']} top1={v['share_top1']:.2f} "
                          f"top2={v['share_top2']:.2f} LOO S/N {v['loo_snr_min']:.1f}-{v['loo_snr_max']:.1f}")
        # timing: pipeline set (one light curve per sector)
        keys = [lc["key"] for lc in SETS["pipe"]]
        allrows = out[sig]["pipe"]["rows"]
        out[sig]["ephem_pipe"] = linear_ephem(allrows, keys)
        keys_own = [lc["key"] for lc in SETS["own"]]
        out[sig]["ephem_own"] = linear_ephem(out[sig]["own"]["rows"], keys_own)
        for lab in ("ephem_pipe", "ephem_own"):
            v = out[sig][lab]
            if v:
                print(f"{sig} {lab}: n={v['n']} P={v['P']:.7f}±{v['eP_scaled']:.7f} T0={v['T0']:.5f}±{v['eT0_scaled']:.5f} "
                      f"chi2={v['chi2']:.1f}/{v['dof']} rms O-C={v['rms_oc_min']:.1f} min")
    save_json(out, "out/s09_events.json")
