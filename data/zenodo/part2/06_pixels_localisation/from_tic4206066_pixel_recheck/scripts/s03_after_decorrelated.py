"""AFTER: per-event local correction for pointing drift (centroid decorrelation, cf. Vanderburg & Johnson 2014).
For each transit: window +-0.35 d; OOT fit of relative flux with [1, dt, dt^2, cx, cy, cx^2, cy^2, cx*cy, bkg];
divide; depth = 1 - mean(in-transit). Combine events (inverse variance); aperture differences per event,
error from scatter across events (empirical)."""
import json, os, sys
import numpy as np
from common import ROOT, OUT, EPH, load, radii_grid, centroid
T7 = json.load(open(os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/out/t7_momentum_dumps.json")))
DUMPS = {"S6_ffi": T7["S06"]["dump_times"], "S32_ffi": T7["S32"]["dump_times"], "S32_tpf2min": T7["S32"]["dump_times"], "S98_ffi": T7["S98_tesscut"]["dump_times"]}

RAD = (0.5, 1.0, 2.0, 3.0)
DIFFS = ((1.0, 2.0), (0.5, 2.0), (1.0, 3.0), (2.0, 3.0))   # (1.0, 2.0) is the quantity of the old T9 test

def run(name, sig, in_frac=0.25, min_in=3, dump_rule=True):
    t, F, bk, c, med = load(name)
    rr = radii_grid(F.shape[1:], c); cx, cy = centroid(F, c)
    ap = {r: F[:, rr <= r].sum(1) for r in RAD}
    P, T0, T14 = EPH[sig]; Po, T0o, T14o = EPH["B" if sig == "A" else "A"]
    ev = []
    for n in np.unique(np.round((t - T0) / P)):
        tc = T0 + n * P; w = np.abs(t - tc) < 0.35
        if w.sum() < 20: continue
        dt = t[w] - tc; inn = np.abs(dt) < in_frac * T14; oot = np.abs(dt) > 0.75 * T14
        no = np.round((t[w] - T0o) / Po); oot &= np.abs(t[w] - (T0o + no * Po)) > 0.75 * T14o
        cad = np.median(np.diff(t)); expect = 2 * in_frac * T14 / cad
        if inn.sum() < max(min_in, 0.7 * expect): continue          # drop events with gaps inside the transit window
        if np.max(np.diff(t[w])) > 0.1: continue                     # drop windows with a data gap > 2.4 h
        if dump_rule and any(abs(x - tc) < 0.35 for x in DUMPS[name]): continue   # drop windows containing a momentum dump
        if inn.sum() < min_in or oot.sum() < 20 or (dt[oot] < 0).sum() < 5 or (dt[oot] > 0).sum() < 5: continue
        x = cx[w] - cx[w][oot].mean(); y = cy[w] - cy[w][oot].mean(); b = bk[w] - bk[w][oot].mean()
        X = np.vstack([np.ones_like(dt), dt, dt**2, x, y, x*x, y*y, x*y, b]).T
        row = {"epoch": int(n), "t_mid_BTJD": float(tc), "n_in": int(inn.sum()), "raw": {}, "corrected": {}}
        for r in RAD:
            fo = ap[r][w] / np.median(ap[r][w][oot])
            # raw = only a quadratic in time (no pointing terms); corrected = full model
            for key, cols in (("raw", slice(0, 3)), ("corrected", slice(None))):
                coef, *_ = np.linalg.lstsq(X[oot][:, cols], fo[oot], rcond=None)
                res = fo / (X[:, cols] @ coef)
                row[key][str(r)] = [float(1e6 * (1 - res[inn].mean())), float(1e6 * np.std(res[oot]) / np.sqrt(inn.sum()))]
        ev.append(row)
    summ = {"n_events": len(ev), "in_window_frac_T14": in_frac, "dump_rule": dump_rule}
    if len(ev) < 2:
        return {"summary": summ, "events": ev}
    for key in ("raw", "corrected"):
        s = {}
        for r in RAD:
            v = np.array([e[key][str(r)][0] for e in ev]); e_ = np.array([e[key][str(r)][1] for e in ev]); w_ = 1 / e_**2
            s[str(r)] = {"depth_ppm": float((v*w_).sum()/w_.sum()), "err_formal": float(1/np.sqrt(w_.sum())),
                         "err_event_scatter": float(np.std(v, ddof=1)/np.sqrt(len(v))) if len(v) > 1 else None}
        for a, b in DIFFS:
            d = np.array([e[key][str(a)][0] - e[key][str(b)][0] for e in ev])
            s[f"diff_{a}-{b}"] = {"ppm": float(d.mean()), "err_empirical": float(d.std(ddof=1)/np.sqrt(len(d))) if len(d) > 1 else None}
        summ[key] = s
    return {"summary": summ, "events": ev}

res = {"A": {}, "B": {}, "A_sensitivity_no_dump_rule": {}}
for name in ("S32_tpf2min", "S32_ffi", "S98_ffi"):
    res["A"][name] = run(name, "A")
    res["A_sensitivity_no_dump_rule"][name] = run(name, "A", dump_rule=False)
# S6 FFI has 30-min cadence: the +-T14/4 window holds < 3 points, so a wider in-transit window is used (documented).
# Every S6 A transit lies 0.4-3.5 h from a momentum dump, so the dump rule removes all of them.
res["A"]["S6_ffi"] = run("S6_ffi", "A", in_frac=0.45, min_in=2)
res["A_sensitivity_no_dump_rule"]["S6_ffi"] = run("S6_ffi", "A", in_frac=0.45, min_in=2, dump_rule=False)
res["B"]["S32_tpf2min"] = run("S32_tpf2min", "B")
# Combined single-pixel excess (r=0.5 minus r=2.0) over independent photon sets (S32 FFI shares photons with the
# S32 2-min TPF and is left out), inverse-variance weights from the empirical per-sector errors.
comb = {}
for key in ("raw", "corrected"):
    for dk in ("diff_1.0-2.0", "diff_0.5-2.0"):
        sets = [n for n in ("S6_ffi", "S32_tpf2min", "S98_ffi") if key in res["A"][n]["summary"]]
        v = [(res["A"][n]["summary"][key][dk]["ppm"], res["A"][n]["summary"][key][dk]["err_empirical"]) for n in sets]
        w = [1 / e**2 for _, e in v]
        comb.setdefault(key, {})[dk] = {"ppm": float(sum(x*wi for (x, _), wi in zip(v, w)) / sum(w)), "err": float(1 / np.sqrt(sum(w))), "sets": sets}
res["A_combined"] = comb
print("combined single-pixel excess:", comb)
json.dump(res, open(os.path.join(OUT, "after_decorrelated.json"), "w"), indent=1)
for s in ("A", "B", "A_sensitivity_no_dump_rule"):
    for name, v in res[s].items():
        if "corrected" not in v["summary"]:
            print(s, name, v["summary"]["n_events"], "events - too few"); continue
        c = v["summary"]["corrected"]
        print(s, name, v["summary"]["n_events"], {k: (round(x["ppm"]), round(x["err_empirical"] or 0)) for k, x in c.items() if k.startswith("diff")})
