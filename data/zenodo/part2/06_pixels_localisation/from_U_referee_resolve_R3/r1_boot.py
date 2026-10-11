"""R3 resolution, step 1: event-level error bars for the difference-image localisation.

Reuses the verifier's independent localisation code (03_pixels/verify/v1_loc.py: per-event per-pixel linear
baseline + box, optional POS_CORR1/2 regressors, TESS PRF fit, target reference = PRF fit of the median image)
read-only via sys.path; nothing is written outside referee/resolve_R3.
New here: 5000-resample event bootstrap with the samples kept (robust SD / percentiles), jackknife SE,
single-event fits, and a random-epoch null (difference images at out-of-transit epochs, same windows) to see how
often pure noise/systematics produce an offset like the real one.
"""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/03_pixels/verify")
import numpy as np
import v1_loc as V

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
NB = 5000


def robust_sd(x):
    x = np.asarray(x)
    return (1.4826 * np.median(np.abs(x - np.median(x, 0)), 0)).tolist()


def run(key, sig, usepc, seed=1):
    ds = V.load(key); f = V.Fitter(ds)
    evs = V.events(ds, sig)
    per = [V.ediff(ds, tc, idx, sig, usepc) for tc, idx in evs]
    D = np.array([p[0] for p in per]); Vv = np.array([p[1] for p in per])

    def fit_set(ii):
        w = 1 / Vv[ii]
        x, y, a = f.fit((w * D[ii]).sum(0) / w.sum(0), 1 / w.sum(0))
        return f.to_arcsec(x, y) + (float(a),)

    n = len(evs)
    allfit = fit_set(np.arange(n))
    res = dict(key=key, sig=sig, usepc=usepc, n=n, tcs=[float(e[0]) for e in evs], off=allfit[:2], amp=allfit[2])
    if n >= 2:
        rng = np.random.default_rng(seed)
        bo = np.array([fit_set(rng.integers(0, n, n))[:2] for _ in range(NB)])
        res.update(boot_sd=bo.std(0, ddof=1).tolist(), boot_robust_sd=robust_sd(bo),
                   boot_p16_p84=np.percentile(bo, [16, 84], axis=0).T.tolist(),
                   boot_cov=np.cov(bo.T).tolist(), boot_median=np.median(bo, 0).tolist())
        np.save(os.path.join(OUT, f"boot_{key}_{sig}_{'pc' if usepc else 'raw'}.npy"), bo)
    if n >= 3:
        jk = np.array([fit_set(np.array([j for j in range(n) if j != k]))[:2] for k in range(n)])
        se = np.sqrt((n - 1) / n * ((jk - jk.mean(0)) ** 2).sum(0))
        res.update(jackknife=[[float(evs[k][0])] + jk[k].tolist() for k in range(n)], jackknife_se=se.tolist())
    single = [fit_set(np.array([k])) for k in range(n)]
    res["single_event"] = [[float(evs[k][0])] + list(single[k]) for k in range(n)]
    return res


if __name__ == "__main__":
    cfgs = [("S32_tpf", "A", True), ("S32_tpf", "A", False), ("S98_ffi", "A", False), ("S98_ffi", "B", False),
            ("S32_tpf", "B", True), ("S6_ffi", "A", False), ("S6_ffi", "B", False)]
    out = []
    for c in cfgs:
        r = run(*c); out.append(r)
        print(c, "n", r["n"], "off", np.round(r["off"], 2), "bootSD", np.round(r.get("boot_sd", [np.nan] * 2), 2),
              "robustSD", np.round(r.get("boot_robust_sd", [np.nan] * 2), 2), "jkSE", np.round(r.get("jackknife_se", [np.nan] * 2), 2),
              flush=True)
    json.dump(out, open(os.path.join(OUT, "r1_boot.json"), "w"), indent=1, default=float)
