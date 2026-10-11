"""Step 10c: what produces the largest phase-scan features that s10 found outside the transits?
  B, own photometry: |z| = 7.7 at B-phase 0.925 (+440 ppm)   (pipeline photometry: max |z| = 3.6, at phase 0.24)
  A, own photometry: |z| = 3.8 at A-phase 0.27 (+147 ppm)    (pipeline photometry: max |z| = 2.2)
For each feature: the box depth per cycle (which transit-window-sized stretch of data produces it), per sector, and the
same box in the pipeline photometry. Also: the s10 random boxes were drawn only at phases 0.1-0.9, so the s10
'fraction of random boxes' does not describe phases < 0.1 or > 0.9; here the empirical sigma is re-derived over 0.05-0.95.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import s10_eb_tests as s10  # noqa
from audit_common import save_json  # noqa

rng = np.random.default_rng(7)


def fold_data(sig, setname):
    lcs = s10.load(setname)
    other = "B" if sig == "A" else "A"
    rows = []
    for lc, sec in zip(lcs, (6, 32, 98)):
        t, f, w = lc["t"], lc["fd"], np.full(len(lc["t"]), 1 / lc["sig"] ** 2)
        keep = ~lcio.in_transit(t, other, 1.0, s10.EPH)
        rows.append((sec, t[keep], f[keep], w[keep]))
    return rows


def per_cycle(rows, sig, phase):
    e = s10.EPH[sig]
    P, T0, T14 = e["P"], e["T0"], e["T14"]
    out = []
    for sec, t, f, w in rows:
        ph = ((t - T0) / P - phase + 0.5) % 1 - 0.5
        m = np.abs(ph * P) < T14 / 2
        cyc = np.round((t - T0) / P - phase).astype(int)
        for c in np.unique(cyc[m]):
            mm = m & (cyc == c)
            if mm.sum() < 3:
                continue
            d = 1 - np.sum(w[mm] * f[mm]) / np.sum(w[mm])
            err = 1 / np.sqrt(np.sum(w[mm]))
            out.append(dict(sector=sec, cycle=int(c), t_mid=float(T0 + (c + phase) * P), n=int(mm.sum()),
                            depth_ppm=float(d * 1e6), err_white_ppm=float(err * 1e6), z_white=float(d / err)))
    return out


def emp_sigma_full(rows, sig):
    e = s10.EPH[sig]
    P, T0, T14 = e["P"], e["T0"], e["T14"]
    t = np.concatenate([r[1] for r in rows]); f = np.concatenate([r[2] for r in rows]); w = np.concatenate([r[3] for r in rows])
    vals = []
    for ph in rng.uniform(0.05, 0.95, 3000):
        if abs(ph) < 1.5 * T14 / P or abs(ph - 1) < 1.5 * T14 / P:
            continue
        d, n = s10.box(t, f, w, P, T0, ph, T14)
        if n >= 3:
            vals.append(d)
    vals = np.array(vals)
    return float(1.4826 * np.median(np.abs(vals - np.median(vals)))), vals


res = {}
for sig, phase in (("B", 0.925), ("A", 0.27)):
    for setname in ("own", "pipe"):
        rows = fold_data(sig, setname)
        pc = per_cycle(rows, sig, phase)
        s_rob, vals = emp_sigma_full(rows, sig)
        t = np.concatenate([r[1] for r in rows]); f = np.concatenate([r[2] for r in rows]); w = np.concatenate([r[3] for r in rows])
        e = s10.EPH[sig]
        d, n = s10.box(t, f, w, e["P"], e["T0"], phase, e["T14"])
        # leave-one-cycle-out
        loo = []
        for c in pc:
            cyc = np.round((t - e["T0"]) / e["P"] - phase).astype(int)
            m = cyc != c["cycle"]
            dl, _ = s10.box(t[m], f[m], w[m], e["P"], e["T0"], phase, e["T14"])
            loo.append(dict(sector=c["sector"], cycle=c["cycle"], depth_without_ppm=float(dl * 1e6)))
        top = sorted(pc, key=lambda r: -r["z_white"])[:4]
        res[f"{sig}_{setname}_phase{phase}"] = dict(fold_depth_ppm=float(d * 1e6), emp_sigma_ppm=s_rob * 1e6, z_emp=float(d / s_rob),
                                                    frac_random_ge=float(np.mean(vals >= d)), per_cycle=pc, top_cycles=top, leave_one_out=loo)
        print(f"{sig} {setname:4s} phase {phase}: fold {d*1e6:6.0f} ppm, emp sigma {s_rob*1e6:.0f}, z {d/s_rob:.1f}; top cycles:",
              "; ".join(f"S{r['sector']} c{r['cycle']} t={r['t_mid']:.2f} {r['depth_ppm']:.0f}±{r['err_white_ppm']:.0f}" for r in top),
              "| min leave-one-out", f"{min(l['depth_without_ppm'] for l in loo):.0f}", flush=True)
save_json(res, "out/s10c_phase_features.json")
