import json, os, sys
import numpy as np
from v01_independent import *

def nsig(dep):
    d6, e6 = dep[6]; d32, e32 = dep[32]; d98, e98 = dep[98]
    w = np.array([1 / e6 ** 2, 1 / e32 ** 2]); a = (w[0] * d6 + w[1] * d32) / w.sum(); ea = 1 / np.sqrt(w.sum())
    D = np.array([d6, d32, d98]); E = np.array([e6, e32, e98]); W = 1 / E ** 2; mu = np.sum(W * D) / W.sum()
    chi2 = float(np.sum(((D - mu) / E) ** 2))
    return dict(arch=[a, ea], nsig=float((d98 - a) / np.hypot(e98, ea)), common=[float(mu), float(1 / np.sqrt(W.sum()))], dchi2=chi2)

if __name__ == "__main__":
    out = {}
    cfgs = []
    for det in ("local", "median"):
        for poly in (1, 2):
            for Amode in ("fixed", "fit", "mask", "maskhalf"):
                for keep512 in (False, True):
                    cfgs.append((det, poly, Amode, keep512))
    data = {}
    for keep512 in (False, True):
        for sec in (6, 32, 98):
            t, f, q, cad = load(sec, keep512=keep512)
            data[(sec, keep512, "local")] = (t, f, cad)
            data[(sec, keep512, "median")] = (t, detrend_median(t, f), cad)
    nullcache = {}
    for det, poly, Amode, keep512 in cfgs:
        if keep512 and not (Amode == "fixed"):
            continue
        dep_w = {}; dep_n = {}; per = {}
        for sec in (6, 32, 98):
            t, f, cad = data[(sec, keep512, det)]
            tcs = [tc for n, tc in btimes(t)]
            D, e, ne = fit(t, f, cad, tcs, poly, Amode)
            De, ee, _ = fit(t, f, cad, tcs, poly, Amode, perevent=True)
            key = (sec, det, poly, "fixed" if Amode != "fit" else "fit", ne)
            if key not in nullcache:
                nullcache[key] = null_err(t, f, cad, ne, poly, "fixed" if Amode != "fit" else "fit", nnull=300)
            sn, mn, nn = nullcache[key]
            dep_w[sec] = (float(D[0]), float(e[0])); dep_n[sec] = (float(D[0]), max(sn, float(e[0])))
            per[sec] = dict(D=float(D[0]), e_white=float(e[0]), e_null=sn, null_mean=mn, n_null=nn, D_event=De.tolist(), e_event=ee.tolist())
        name = f"{det}_p{poly}_{Amode}{'_keep512' if keep512 else ''}"
        out[name] = dict(per=per, white=nsig(dep_w), null=nsig(dep_n))
        print(name, " ".join(f"S{s}={per[s]['D']:.0f}+-{per[s]['e_white']:.0f}/{per[s]['e_null']:.0f}" for s in (6, 32, 98)),
              f"nsig_w={out[name]['white']['nsig']:.2f} nsig_null={out[name]['null']['nsig']:.2f} common={out[name]['null']['common'][0]:.0f}+-{out[name]['null']['common'][1]:.0f} dchi2={out[name]['null']['dchi2']:.1f}", flush=True)
    json.dump(out, open(os.path.join(H, "v02_depths.json"), "w"), indent=1)
