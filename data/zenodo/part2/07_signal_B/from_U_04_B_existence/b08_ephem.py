"""Step 8: per-event mid-times of B and linear ephemerides for several event sets, with injection-calibrated timing errors
(same fitter as for the real events: trapezoid T14 = 1.6 h, ingress 0.15 T14, local linear baseline, +-90 min scan)."""
import numpy as np, json
import common_b as c
import archdata as ad
import s98data as sd

rng = np.random.default_rng(8)
PM, TM = c.P_B_MS, c.T0_B_MS


def supersample_trap(t, tc, T14, tau, cad):
    if cad < 0.005:
        return c.trapezoid(t - tc, T14, tau)
    k = 15
    off = (np.arange(k) + 0.5) / k - 0.5
    return np.mean([c.trapezoid(t + o * cad - tc, T14, tau) for o in off], axis=0)


def fit_event_ss(t, f, tc0, cad, scan_min=90.0, step_min=0.5, T14=c.T14_B):
    """fit_event with the model integrated over the cadence (needed for 30-min S6 data)."""
    win = 3 * T14
    m = np.isfinite(f) & (np.abs(t - tc0) < win + scan_min / 1440)
    tt, ff = t[m], f[m]
    if len(tt) < 8: return None
    x0 = tt - tc0
    best = None
    for dtc in np.arange(-scan_min, scan_min + 1e-9, step_min) / 1440:
        s = supersample_trap(tt, tc0 + dtc, T14, 0.15 * T14, cad)
        if np.sum(s > 0.3) < 1: continue
        A = np.vstack([np.ones_like(x0), x0, -s]).T
        coef, *_ = np.linalg.lstsq(A, ff, rcond=None)
        chi = np.sum((ff - A @ coef) ** 2)
        if coef[2] > 0 and (best is None or chi < best[0]):
            best = (chi, dtc, coef[2])
    if best is None: return None
    return dict(tc=tc0 + best[1], dt_min=best[1] * 1440, D=best[2], at_edge=bool(abs(best[1]) * 1440 > scan_min - 1))


def inj_err(t, f, cad, depth, avoid, n=300):
    dts = []
    k = 0
    while len(dts) < n and k < 30 * n:
        k += 1
        tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
        if np.min(np.abs(tc - np.asarray(avoid))) < 4 * c.T14_B: continue
        w = np.abs(t - tc) < c.T14_B / 2
        if w.sum() < 0.8 * c.T14_B / cad: continue
        fi = f - depth * supersample_trap(t, tc, c.T14_B, 0.15 * c.T14_B, cad)
        r = fit_event_ss(t, fi, tc, cad)
        if r is not None: dts.append(r["dt_min"])
    dts = np.array(dts)
    return float(np.std(dts)), float(c.robust_std(dts)), float(np.mean(np.abs(dts) > 30)), len(dts)


data = {6: ad.get(6, "pipe"), 32: ad.get(32, "pipe"), 98: sd.get("QLP_DET_bw075")}
# S98 scattered-light epoch 2 needs the relaxed selection (2^29 kept): build it here
q = c.load_qlp(98)
t98, qq = q["t"], q["q"]
rel = np.isfinite(q["f"]) & (((qq >> 11) & 1) == 0) & (((qq >> 5) & 1) == 0) & (((qq >> 30) & 1) == 0)
f98r = np.full_like(t98, np.nan)
f98r[rel] = c.detrend(t98[rel], q["f"][rel], "biweight", 0.75, mask=c.in_A(t98[rel]))
keep = rel & ~c.in_A(t98)
data["98rel"] = (t98[keep], f98r[keep])
cads = {6: 1800 / 86400, 32: 120 / 86400, 98: 200 / 86400, "98rel": 200 / 86400}
events = []
for sec in [6, 32, 98, "98rel"]:
    t, f = data[sec]
    ns, tcs = c.pred_times(PM, TM, t.min(), t.max())
    for n, tc in zip(ns, tcs):
        if sec == "98rel" and n != 2: continue
        if sec == 98 and n == 2: continue
        w = np.abs(t - tc) < c.T14_B / 2
        cov = w.sum() / (c.T14_B / cads[sec])
        if cov < 0.3: 
            events.append(dict(sector=str(sec), epoch=int(n), tc_pred=float(tc), coverage=float(cov), used=False)); continue
        r = fit_event_ss(t, f, tc, cads[sec])
        events.append(dict(sector=str(sec), epoch=int(n), tc_pred=float(tc), coverage=float(cov), used=r is not None,
                           tc=r and float(r["tc"]), dt_min=r and float(r["dt_min"]), D_ppm=r and float(r["D"] * 1e6), at_edge=r and r["at_edge"]))
errs = {}
avoid_all = list(c.pred_times(PM, TM, 1400, 4100)[1])
for sec, dep in [(6, 850e-6), (32, 1000e-6), (98, 1100e-6), ("98rel", 1100e-6)]:
    t, f = data[sec]
    if sec == "98rel":
        m = (t < 4024.4) & (t > 4017.0)   # the 2^29 half-orbit containing epoch 2
        t, f = t[m], f[m]
    errs[str(sec)] = inj_err(t, f, cads[sec], dep, avoid_all)
    print("timing error", sec, errs[str(sec)], flush=True)
# partial S98 epoch 4: 55 % of cadences flagged -> injection with random 55 % dropout (from b01: 7.8 min); recompute here
t, f = data[98]
dts = []
while len(dts) < 300:
    tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
    if np.min(np.abs(tc - np.asarray(avoid_all))) < 4 * c.T14_B: continue
    ww = np.abs(t - tc) < c.T14_B / 2
    if ww.sum() < 0.8 * c.T14_B / cads[98]: continue
    fi = f - 1100e-6 * c.trapezoid(t - tc, c.T14_B, 0.15 * c.T14_B)
    fi = np.where((np.abs(t - tc) < 0.5) & (rng.random(len(t)) < 0.55), np.nan, fi)
    r = fit_event_ss(t, fi, tc, cads[98])
    if r: dts.append(r["dt_min"])
errs["98_partial"] = (float(np.std(dts)), float(c.robust_std(np.array(dts))), float(np.mean(np.abs(np.array(dts)) > 30)), len(dts))
print("timing error 98 partial", errs["98_partial"])
for e in events:
    if not e["used"]: continue
    key = e["sector"]
    if e["sector"] == "98" and e["epoch"] == 4: key = "98_partial"
    e["err_min"] = errs[key][0]
    print(e)


def linfit(ep, tt, ee):
    ep = np.asarray(ep, float); tt = np.asarray(tt); ee = np.asarray(ee)
    A = np.vstack([np.ones_like(ep), ep]).T
    W = np.diag(1 / ee ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ tt
    res = tt - A @ p
    chi2 = float(np.sum((res / ee) ** 2)); dof = len(tt) - 2
    sc = np.sqrt(max(1.0, chi2 / dof)) if dof > 0 else 1.0
    return dict(P=float(p[1]), eP=float(np.sqrt(cov[1, 1])), eP_scaled=float(np.sqrt(cov[1, 1]) * sc), T0=float(p[0]), eT0=float(np.sqrt(cov[0, 0])),
                eT0_scaled=float(np.sqrt(cov[0, 0]) * sc), chi2=chi2, dof=dof, oc_min=(res * 1440).round(2).tolist(), n=len(tt))


def sel(cond):
    E = [e for e in events if e["used"] and not e.get("at_edge") and cond(e)]
    return [e["epoch"] for e in E], [e["tc"] for e in E], [e["err_min"] / 1440 for e in E], [(e["sector"], e["epoch"]) for e in E]


sets = {
    "S6x2 + S98(0,1) [manuscript-like, 4 events]": lambda e: (e["sector"] == "6") or (e["sector"] == "98" and e["epoch"] in (0, 1)),
    "S6x2 + S32x2 + S98(0,1,4)": lambda e: e["sector"] in ("6", "32", "98"),
    "S6x2 + S32x2 + S98(0,1,4) + S98 SL(2)": lambda e: True,
    "S6x2 + S98(0,1,4)": lambda e: e["sector"] in ("6", "98"),
    "S98 only (0,1,4)": lambda e: e["sector"] == "98",
}
eph = {}
for k, fn in sets.items():
    ep, tt, ee, ids = sel(fn)
    eph[k] = dict(linfit(ep, tt, ee), events=ids)
    print("%-45s P = %.6f +- %.6f (scaled %.6f)  T0 = %.4f +- %.4f (scaled %.4f)  chi2 %.1f/%d  O-C %s" % (
        k, eph[k]["P"], eph[k]["eP"], eph[k]["eP_scaled"], eph[k]["T0"], eph[k]["eT0"], eph[k]["eT0_scaled"], eph[k]["chi2"], eph[k]["dof"], eph[k]["oc_min"]))
json.dump(dict(events=events, timing_errors_min=errs, ephemerides=eph), open("out/b08_ephem.json", "w"), indent=1, default=float)
