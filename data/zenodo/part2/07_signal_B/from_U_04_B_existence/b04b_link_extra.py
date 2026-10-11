"""Step 4 extras: (a) number of distinct winning box configurations in the link-test null (resolution of the FAP);
(b) exponential-tail extrapolation of the null; (c) look-elsewhere check: link test of the other S98 search peaks
(9.97 d, 1.06 d) with the same machinery."""
import numpy as np, json, sys, time
import common_b as c
import archdata as ad

rng = np.random.default_rng(44)


def scan(grids, Pc, Tref, w, nstep):
    Pv = Pc + np.linspace(-w, w, nstep)
    s, nev = ad.combined_snr(grids, Tref, Pv)
    i = int(np.argmax(s))
    return s[i], Pv[i]


def config_key(grids, Tref, P):
    key = []
    for g in grids:
        for n in range(-2000, 1):
            tc = Tref + n * P
            if g.tmin - 0.05 < tc < g.tmax + 0.05:
                d, s = g.lookup(np.array([tc]))
                if np.isfinite(d[0]):
                    key.append(int(np.floor(tc / c.T14_B)))
    return tuple(key)


def run(phot, sectors, Pc, Tref, s6_halfwin=0.908, prange=(9.0, 13.5), ntr=50000, excl=0.1, keys=True):
    n6 = (Tref - 1479.0) / Pc
    w = s6_halfwin / n6
    nstep = int(np.ceil(2 * s6_halfwin / (0.5 / 1440))) + 1
    grids = [ad.BoxGrid(*ad.get(s, phot)) for s in sectors]
    obs, Pm = scan(grids, Pc, Tref, w, nstep)
    nulls = np.empty(ntr); ks = set()
    for k in range(ntr):
        while True:
            P0 = rng.uniform(*prange)
            if abs(P0 - Pc) > excl: break
        T0 = Tref + rng.uniform(-0.5, 0.5) * P0
        w0 = w * P0 / Pc
        v, Pbest = scan(grids, P0, T0, w0, nstep)
        nulls[k] = v
        if keys:
            ks.add(config_key(grids, T0, Pbest))
    top = np.sort(nulls)[-int(0.02 * ntr):]            # top 2 %
    thr = top[0]
    lam = 1.0 / np.mean(top - thr)                       # exponential tail fit above the 98th percentile
    p_tail = 0.02 * np.exp(-lam * (obs - thr)) if obs > thr else float(np.mean(nulls >= obs))
    nge = int(np.sum(nulls >= obs))
    return dict(phot=phot, sectors=sectors, Pc=Pc, Tref=Tref, w=w, obs=float(obs), P_at_max=float(Pm), ntr=ntr, n_ge=nge,
                fap=(nge + 1) / (ntr + 1), n_distinct_configs=len(ks) if keys else None,
                null_q=[float(x) for x in np.percentile(nulls, [50, 90, 99, 99.9, 100])], tail_exp_fit=dict(thr98=float(thr), lam=float(lam), p_extrap=float(p_tail)))


if __name__ == "__main__":
    res = []
    for cf in [dict(phot="pipe", sectors=[6, 32], Pc=c.P_B98, Tref=c.T0_B98),
               dict(phot="pipe", sectors=[6], Pc=c.P_B98, Tref=c.T0_B98),
               dict(phot="own", sectors=[6, 32], Pc=c.P_B98, Tref=c.T0_B98),
               # other S98 search peaks (own-photometry and QLP BLS ephemerides from b04b_other_peaks_s98_ephem.json)
               dict(phot="pipe", sectors=[6, 32], Pc=9.969684460454253, Tref=4002.0337894659297, prange=(8.0, 12.0), ntr=20000),
               dict(phot="pipe", sectors=[6, 32], Pc=9.968488235884825, Tref=4002.0299393705286, prange=(8.0, 12.0), ntr=20000),
               dict(phot="pipe", sectors=[6, 32], Pc=1.0556499747408545, Tref=3996.1361893705284, prange=(0.9, 1.2), excl=0.005, ntr=20000),
               dict(phot="pipe", sectors=[6, 32], Pc=1.0552258187737855, Tref=3996.642122799263, prange=(0.9, 1.2), excl=0.005, ntr=20000)]:
        t0 = time.time()
        r = run(**cf)
        res.append(r)
        print(json.dumps(r), "(%.0fs)" % (time.time() - t0), flush=True)
        json.dump(res, open("out/b04b_link_extra.json", "w"), indent=1)
