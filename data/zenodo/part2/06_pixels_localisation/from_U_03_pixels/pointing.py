"""(5) Pointing / jitter near the transits.

For a target-dominated scene, a pointing shift delta (px) that is correlated with the transit window adds a dipole
-grad(M).delta to the difference image; a PRF fit then moves by about -delta/depth.  With depth ~400 ppm a
transit-correlated 1e-3 px shift moves the fitted position by ~2.5 px (~50").  Here we measure, per event, the
in-transit minus (linear) baseline value of
  * POS_CORR1/2 of the SPOC S32 2-min TPF (spacecraft/local image motion from SPOC motion polynomials), and
  * the flux-weighted centroid of the 5x5 pixels around the target (all products; for an on-target transit the
    centroid does not move, so this measures pointing plus any off-target contribution),
and the same statistic at random out-of-transit epochs.  Also: centroid steps at the S6 momentum dumps.
"""
from run_real import *


def centroid_series(ds, half=2):
    xc, yc = int(round(ds.xt_wcs)), int(round(ds.yt_wcs))
    sub = ds.F[:, yc - half:yc + half + 1, xc - half:xc + half + 1]
    sub = sub - np.median(ds.F[:, [0, 0, -1, -1], [0, -1, 0, -1]], axis=1)[:, None, None] * (not ds.cube.is_tpf)
    yy, xx = np.mgrid[yc - half:yc + half + 1, xc - half:xc + half + 1]
    tot = sub.sum((1, 2))
    return (sub * xx).sum((1, 2)) / tot, (sub * yy).sum((1, 2)) / tot


def event_shift(ds, series, idx, tc, sig):
    t = ds.t[idx]
    f = box_fraction(t, ds.dt, tc, EPH[sig]["T14"])
    X = np.vstack([np.ones_like(t), t - tc, f]).T
    out = []
    for s in series:
        y = s[idx]
        ok = np.isfinite(y)
        b, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
        r = y[ok] - X[ok] @ b
        cov = np.linalg.pinv(X[ok].T @ X[ok]) * (r @ r) / max(ok.sum() - 3, 1)
        out.append((b[2], np.sqrt(cov[2, 2])))
    return out


res = {}
rng = np.random.default_rng(5)
for name, key, var in [("S32_tpf", "S32_tpf", "default"), ("S32_ffi", "S32_ffi", "default"),
                       ("S98_clean29", "S98_ffi", "clean29"), ("S6_ffi", "S6_ffi", "default"),
                       ("S6_dumpx0.5h", "S6_ffi", "default")]:
    ds = Dataset(key, var, **({"dump_exclude_h": 0.5} if "dumpx" in name else {}))
    sf = ds.scene_fit()
    J = pix_jacobian_arcsec(ds.cube, sf["x_t"], sf["y_t"])
    cx, cy = centroid_series(ds)
    series = [cx, cy]
    labels = ["centroid_col", "centroid_row"]
    if ds.cube.is_tpf:
        pc1 = ds.cube.poscorr1[ds.good]; pc2 = ds.cube.poscorr2[ds.good]
        series += [pc1, pc2]; labels += ["POS_CORR1", "POS_CORR2"]
    res[name] = {}
    for sig in "AB":
        evs = [e for e in ds.event_windows(sig, min_side_frac=MIN_SIDE) if e["ok"]]
        if not evs:
            continue
        per = np.array([[v for pair in event_shift(ds, series, e["idx"], e["tc"], sig) for v in pair] for e in evs])
        vals, errs = per[:, 0::2], per[:, 1::2]
        w = 1 / errs ** 2
        mean = (vals * w).sum(0) / w.sum(0)
        emean = 1 / np.sqrt(w.sum(0))
        # random-epoch distribution of the same weighted mean (same number of events)
        same = in_any_transit(ds.t, ds.dt, EPH[sig], PAD_OTHER)
        rand = []
        for k in range(300):
            fake = []
            while len(fake) < len(evs):
                tc = rng.uniform(ds.t.min(), ds.t.max())
                ev = ds.event_windows(sig, centres=[tc], min_side_frac=MIN_SIDE)
                if ev and ev[0]["ok"] and not same[ev[0]["idx"]].any():
                    fake.append(ev[0])
            pf = np.array([[v for pair in event_shift(ds, series, e["idx"], e["tc"], sig) for v in pair] for e in fake])
            wv = 1 / pf[:, 1::2] ** 2
            rand.append((pf[:, 0::2] * wv).sum(0) / wv.sum(0))
        rand = np.array(rand)
        depth = EPH[sig]["depth"] * 0.8     # box-average depth, order of magnitude
        # implied PRF-fit displacement for a pointing shift equal to the measured centroid shift
        implied_px = -mean[:2] / depth
        implied_arcsec = J @ implied_px
        res[name][sig] = dict(n_events=len(evs), labels=labels, per_event=per, weighted_mean=mean, err=emean,
                              random_sd=rand.std(0), random_p16_p84=np.percentile(rand, [16, 84], axis=0),
                              z_vs_random=(mean - rand.mean(0)) / rand.std(0),
                              implied_prf_offset_arcsec_if_pointing=implied_arcsec)
        print(name, sig, "n", len(evs), " ".join(f"{l}={m*1e3:+.3f}+-{e*1e3:.3f}mpx(rand sd {s*1e3:.3f})"
                                                for l, m, e, s in zip(labels, mean, emean, rand.std(0))),
              "implied(E,N)\"", np.round(implied_arcsec, 1), flush=True)
    # centroid steps at momentum dumps (median 2 h after minus 2 h before)
    steps = []
    for td in ds.dump_times:
        a = (ds.t > td - 2 / 24) & (ds.t < td - 0.02)
        b = (ds.t > td + 0.02) & (ds.t < td + 2 / 24)
        if a.sum() > 2 and b.sum() > 2:
            steps.append([td, float(np.median(cx[b]) - np.median(cx[a])), float(np.median(cy[b]) - np.median(cy[a]))])
    res[name]["dump_centroid_steps_px"] = steps
    if steps:
        st = np.array(steps)
        print(name, "dump centroid steps (mpx):", np.round(st[:, 1:] * 1e3, 2).tolist())
dump(res, "pointing.json")
