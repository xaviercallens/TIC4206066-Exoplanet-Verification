"""S32 2-min TPF: remove pointing-induced pixel variations with SPOC POS_CORR1/2 before difference imaging.

Per pixel, the high-passed flux (minus a 0.5-d running median) is regressed on the high-passed POS_CORR1/2
(out-of-transit cadences of both signals only); the fitted pixel gradients g1, g2 are then used to subtract
g1*POS_CORR1 + g2*POS_CORR2 from the pixel series.  The localisation and its on-target injection calibration
(model-free, same as inject.py) are then redone.  POS_CORR comes from SPOC motion polynomials fitted to many
stars on the CCD, so it is independent of the target's own transit.
"""
import inject
from inject import *
from scipy.ndimage import median_filter


def running_median(t, y, width_d=0.5):
    # cadence-based running median (2-min cadence -> 360 cadences per 0.5 d)
    n = max(3, int(round(width_d / np.median(np.diff(t)))) | 1)
    return median_filter(y, size=n, mode="nearest")


def decorrelate(ds):
    good = ds.good
    pc1 = ds.cube.poscorr1[good].copy(); pc2 = ds.cube.poscorr2[good].copy()
    ok = np.isfinite(pc1) & np.isfinite(pc2)
    oot = ok & ~in_any_transit(ds.t, ds.dt, EPH["A"], 1 / 24) & ~in_any_transit(ds.t, ds.dt, EPH["B"], 1 / 24)
    p1 = pc1 - running_median(ds.t, np.nan_to_num(pc1)); p2 = pc2 - running_median(ds.t, np.nan_to_num(pc2))
    X = np.vstack([p1, p2]).T
    Fn = ds.F.reshape(len(ds.t), -1)
    G = np.zeros((2, Fn.shape[1]))
    for j in range(Fn.shape[1]):
        r = Fn[:, j] - running_median(ds.t, Fn[:, j])
        G[:, j], *_ = np.linalg.lstsq(X[oot], r[oot], rcond=None)
    corr = np.vstack([np.nan_to_num(pc1), np.nan_to_num(pc2)]).T @ G
    ds.F = (Fn - corr).reshape(ds.F.shape)
    return G.reshape((2,) + ds.F.shape[1:])


if __name__ == "__main__":
    import sys
    CFGS = sys.argv[1:] or ["S32_tpf"]
    out = {}
    for cfg in CFGS:
        for sig in "AB":
            job = f"{cfg}_poscorr_{sig}"
            # build a dataset, decorrelate, then run the same injection code on it
            ds = make_dataset(cfg)
            G = decorrelate(ds)
            sf = ds.scene_fit()
            evs = [e for e in ds.event_windows(sig, min_side_frac=MIN_SIDE) if e["ok"]]
            real, Dm, Ve, _ = analyse(ds, sig, sf, events=evs, want_nb=False)
            # jackknife
            jk = []
            for k in range(len(evs)):
                o2, *_ = analyse(ds, sig, sf, events=[e for i, e in enumerate(evs) if i != k], want_nb=False)
                jk.append((evs[k]["tc"], o2["off_E"], o2["off_N"]))
            # injections on the decorrelated cube
            inject.make_dataset = lambda name, _ds=ds: _ds
            inject.JOBS[job] = (cfg, sig, 300, (0, 0), 101 + 2 * CFGS.index(cfg) + (sig == "B"))
            r = inject.run_job(job, os.path.join("inj", f"inj_{job}.json"))
            X = np.array([[t["off_E"], t["off_N"]] for t in r["trials"]])
            out[job] = dict(real_off=[real["off_E"], real["off_N"]], snr=real["snr_target"], dchi2_target=real["dchi2_target"],
                            jackknife=jk, inj_mean=X.mean(0), inj_sd=X.std(0, ddof=1), inj_cov=np.cov(X.T),
                            gradient_rms=[float(np.sqrt((G[0] ** 2).mean())), float(np.sqrt((G[1] ** 2).mean()))])
            print(job, "real", np.round([real["off_E"], real["off_N"]], 1), "S/N", round(real["snr_target"], 1),
                  "inj mean", np.round(X.mean(0), 1), "sd", np.round(X.std(0, ddof=1), 1),
                  "jackknife", [(round(a, 3), round(b, 1), round(c, 1)) for a, b, c in jk], flush=True)
    dump(out, "poscorr_decorr_" + "_".join(CFGS) + ".json")
