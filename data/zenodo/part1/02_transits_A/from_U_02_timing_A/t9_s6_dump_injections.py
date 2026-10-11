"""Audit 2: S6 injections that mirror the real events' proximity to reaction-wheel desaturations.
Every real S6 transit of A lies 0.4-3.5 h from a momentum dump, but the generic injections (t4) must avoid the real
transits and therefore never sample the 0-4 h vicinity of a dump. Here the real A signal is first divided out
(template at the measured per-event time and depth), then the template is injected at dump + delta,
|delta| uniform in [0.3, 4] h (either side), and measured exactly as in t4. A control set uses the same A-removed
light curve with uniform random times (checks that the removal itself does not change the statistics)."""
import json
import numpy as np, pandas as pd
import t2_common as C

rng = np.random.default_rng(606)
tpl = C.Template(); T14 = tpl.t14
ev = pd.read_csv(f"{C.OUT}/t3_events_base.csv")
meta = json.load(open(f"{C.OUT}/t3_meta_base.json"))
P, T0 = meta["ephem_window"]["P"], meta["ephem_window"]["T0"]
lc = C.load_lc(6)
t, f = lc["t"], lc["f"].copy()
e6 = ev[ev.sector == 6]
for _, r in e6.iterrows():
    f = f / (1 - r.D * tpl(t - r.tc_obs, lc["cad"]))
Dinj = float(np.average(e6[e6.sel].D, weights=1 / e6[e6.sel].eD_white ** 2))
dumps = np.array([1468.3912, 1471.5161, 1474.6411, 1481.2868, 1484.4117, 1487.5366])
_, tcsA = C.a_epochs(t, P, T0)
rows = []
for mode, N in (("near_dump", 600), ("random_control", 600)):
    k = 0; tries = 0
    while k < N and tries < 50 * N:
        tries += 1
        if mode == "near_dump":
            ti = rng.choice(dumps) + rng.choice([-1, 1]) * rng.uniform(0.3, 4.0) / 24
        else:
            ti = rng.uniform(t.min(), t.max())
        ch = np.abs(t - ti) < C.CFG["chunk"] + 0.05
        fi = f.copy(); fi[ch] = f[ch] * (1 - Dinj * tpl(t[ch] - ti, lc["cad"]))
        # A already removed: the detrending masks the real A epochs as for real events (harmless)
        prep = C.prepare_window(dict(lc, f=fi), ti, tcsA, T14, extra_mask=ti)
        if prep is None:
            continue
        tt, yy, cov = prep
        if cov < 0.5 or len(tt) < C.n_min(lc):
            continue
        r = C.measure(tt, yy, ti, lc["cad"], tpl, rng, n_boot=50, pb_max=40); r.pop("chi_curve")
        r.update(mode=mode, t_inj=ti, coverage=cov, ddump_h=float(np.min(np.abs(ti - dumps)) * 24))
        rows.append(r); k += 1
d = pd.DataFrame(rows)
d["sel"] = (d.coverage >= 0.7) & (d.snr >= 2) & (~d.edge)
d.to_csv(f"{C.OUT}/t9_s6_dump_inj.csv", index=False)
rs = lambda x: float(1.4826 * np.median(np.abs(x - np.median(x))))
res = {}
for m in ("near_dump", "random_control"):
    q = d[(d["mode"] == m) & d.sel]
    res[m] = dict(n=int(len(q)), frac_sel=float(d[d["mode"] == m].sel.mean()), sd=float(q.dt.std()), robust=rs(q.dt.values),
                  sd_R60=float(d[(d["mode"] == m) & (d.coverage >= 0.7)].dt_alt.std()), mean=float(q.dt.mean()),
                  frac_gt_30=float(np.mean(np.abs(q.dt) > 30)), median_e_white=float(q.e_white.median()))
res["D_inj"] = Dinj
json.dump(res, open(f"{C.OUT}/t9_results.json", "w"), indent=1)
print(json.dumps(res, indent=1))
