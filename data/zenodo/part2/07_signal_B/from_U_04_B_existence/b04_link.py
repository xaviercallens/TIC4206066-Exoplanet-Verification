"""Step 4: archival link test.  Statistic = max over a period scan (centred on an S98-only period) of the combined 1.6-h
box-depth S/N at the predicted Sector 6 / 32 times.  Null: random period in 9-13.5 d (|P - Pc| > 0.1 d) and random epoch,
scan width rescaled so that the scanned time window in the archive is the same; data either kept as observed (exchangeable
permutation null) or with B's archival windows removed."""
import numpy as np, json, sys, time
import common_b as c
import archdata as ad

NTR = int(sys.argv[1]) if len(sys.argv) > 1 else 50000
which = sys.argv[2] if len(sys.argv) > 2 else "main"
rng = np.random.default_rng(4)
WINDOWS = {
    "W1_manuscript": dict(Pc=c.P_B98, Tref=c.T0_B98, w=0.004),
    "W2_S98timing_3sig": dict(Pc=11.13508, Tref=4000.8646, w=0.0045),
    "W3_S98_2event_3sig": dict(Pc=11.13508, Tref=4000.8646, w=0.0113),
}
BARCH = [1473.7330, 1484.8657, 2175.0963, 2186.2290, 2197.3617]   # manuscript ephemeris propagated (for 'removed')


def remove_B(t, f):
    m = np.ones_like(t, bool)
    for tc in BARCH:
        m &= np.abs(t - tc) > (c.T14_B / 2 + c.T14_B + 0.02)
    return t[m], f[m]


def scan(grids, Pc, Tref, w, nstep):
    Pv = Pc + np.linspace(-w, w, nstep)
    s, nev = ad.combined_snr(grids, Tref, Pv)
    i = int(np.argmax(s))
    return s[i], Pv[i], s, nev


def run(phot, sectors, window, null_mode, ntr, detw=0.75):
    W = WINDOWS[window]
    step_time_S6 = 0.5 / 1440      # scan resolution: 0.5 min at Sector 6
    n6 = (W["Tref"] - 1479.0) / W["Pc"]
    nstep = int(np.ceil(2 * W["w"] * n6 / step_time_S6)) + 1
    data = {s: ad.get(s, phot if not (s == 32 and phot == "eleanor") else "pipe", window=detw) for s in sectors}
    grids_obs = [ad.BoxGrid(*data[s]) for s in sectors]
    obs, Pmax, curve, nev = scan(grids_obs, W["Pc"], W["Tref"], W["w"], nstep)
    # per-event values at the maximum
    per = []
    for g, s in zip(grids_obs, sectors):
        for n in range(-300, 1):
            tc = W["Tref"] + n * Pmax
            if g.tmin - 0.05 < tc < g.tmax + 0.05:
                d, e = g.lookup(np.array([tc]))
                per.append(dict(sector=s, n=n, tc=float(tc), depth_ppm=float(d[0] * 1e6) if np.isfinite(d[0]) else None,
                                err_ppm=float(e[0] * 1e6) if np.isfinite(e[0]) else None))
    if null_mode == "removed":
        grids_null = [ad.BoxGrid(*remove_B(*data[s])) for s in sectors]
    else:
        grids_null = grids_obs
    nulls = np.empty(ntr); wins = set()
    for k in range(ntr):
        while True:
            P0 = rng.uniform(9.0, 13.5)
            if abs(P0 - W["Pc"]) > 0.1: break
        T0 = W["Tref"] + rng.uniform(-0.5, 0.5) * P0
        w0 = W["w"] * P0 / W["Pc"]
        nulls[k] = scan(grids_null, P0, T0, w0, nstep)[0]
        # distinct window id: predicted S6 time of the scan centre, in units of the S6 scan width
        t6 = T0 + np.round((1479.0 - T0) / P0) * P0
        wins.add(int(np.floor(t6 / (2 * W["w"] * n6))))
    nge = int(np.sum(nulls >= obs))
    return dict(phot=phot, sectors=sectors, window=window, null=null_mode, detrend_window=detw, obs=float(obs), P_at_max=float(Pmax),
                Pc=W["Pc"], w=W["w"], nstep=nstep, per_event=per, beta=[float(g.beta) for g in grids_obs],
                ntr=ntr, n_ge=nge, fap=(nge + 1) / (ntr + 1), fap_naive=nge / ntr,
                null_q=dict(zip(["50", "90", "99", "99.9", "max"], [float(x) for x in np.percentile(nulls, [50, 90, 99, 99.9, 100])])),
                n_distinct_s6_windows=len(wins),
                curve=dict(P=(W["Pc"] + np.linspace(-W["w"], W["w"], nstep))[::20].tolist(), snr=curve[::20].tolist()))


if __name__ == "__main__":
    configs = []
    if which == "main":
        for null in ["kept", "removed"]:
            for secs in [[6, 32], [6], [32]]:
                configs.append(("pipe", secs, "W1_manuscript", null, NTR))
        configs.append(("own", [6, 32], "W1_manuscript", "kept", NTR))
        configs.append(("own", [6], "W1_manuscript", "kept", NTR))
    elif which == "windows":
        for win in ["W2_S98timing_3sig", "W3_S98_2event_3sig"]:
            for secs in [[6, 32], [6]]:
                configs.append(("pipe", secs, win, "kept", NTR))
    elif which == "phot":
        for phot in ["qlp", "sap", "own05"]:
            configs.append((phot, [6, 32], "W1_manuscript", "kept", NTR))
            configs.append((phot, [6], "W1_manuscript", "kept", NTR))
        configs.append(("eleanor", [6], "W1_manuscript", "kept", NTR))
        configs.append(("tglc_aper", [6], "W1_manuscript", "kept", NTR))
    res = []
    for cf in configs:
        t0 = time.time()
        r = run(*cf)
        res.append(r)
        print("%-9s %-8s %-20s %-8s obs %.2f at P %.6f  null med %.2f p99 %.2f max %.2f  n>=obs %d/%d  FAP %.2e  windows %d  beta %s  (%.0fs)" % (
            r["phot"], "+".join(map(str, r["sectors"])), r["window"], r["null"], r["obs"], r["P_at_max"], r["null_q"]["50"], r["null_q"]["99"],
            r["null_q"]["max"], r["n_ge"], r["ntr"], r["fap"], r["n_distinct_s6_windows"], np.round(r["beta"], 2), time.time() - t0), flush=True)
        print("    events:", [(e["sector"], round(e["tc"], 3), e["depth_ppm"] and round(e["depth_ppm"]), e["err_ppm"] and round(e["err_ppm"])) for e in r["per_event"]], flush=True)
        json.dump(res, open("out/b04_link_%s.json" % which, "w"), indent=1)
