"""Step 8b: INDEPENDENT fit configuration (different sampler, LD parameterisation, priors and detrending):
dynesty static nested sampling; Kipping (2013) q1, q2 ~ U(0,1) (no Claret prior); rho ~ log-uniform(0.1, 50);
b ~ U(0, 1+k); k ~ log-uniform(0.005, 0.2); per-data-set offset and jitter; biweight window 1.0 d instead of 0.75 d.
Same likelihood function as s08 (batman), so this checks sampler/prior/detrending sensitivity, not the model.
usage: s08b_dynesty_fit.py SIGNAL SET
"""
import json
import os
import sys
import time

import numpy as np
import dynesty

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import transit_tools as tt  # noqa
from audit_common import AUD, save_json  # noqa

SIGNAL, SET = sys.argv[1], sys.argv[2]
EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
if SET == "pipe":
    raw = [lcio.load_product("TESSSPOC_S6_PDC"), lcio.load_product("SPOC2min_S32_PDC"), lcio.load_product("QLP_S98_SYS", "frozen")]
else:
    raw = [lcio.load_own("tesscut_s6", 2.5), lcio.load_own("tesscut_s32", 2.0), lcio.load_own("tesscut_s98", 2.0, drop_stray=True)]
LCS = [lcio.detrend(lc, window=1.0, eph=EPH) for lc in raw]
W = tt.windows(LCS, SIGNAL, EPH, half=0.30)
MODELS = [tt.Model(w["t"], w["cad"]) for w in W]
NW = len(W)
E = EPH[SIGNAL]
ND = 7 + 2 * NW


def ptform(u):
    x = np.empty(ND)
    x[0] = E["T0"] + (u[0] - 0.5) * 0.02
    x[1] = E["P"] + (u[1] - 0.5) * (2e-4 if SIGNAL == "A" else 6e-4)
    x[2] = np.exp(np.log(0.005) + u[2] * (np.log(0.2) - np.log(0.005)))      # k
    x[3] = u[3] * (1 + x[2])                                                  # b
    x[4] = np.exp(np.log(0.1) + u[4] * (np.log(50) - np.log(0.1)))            # rho
    x[5], x[6] = u[5], u[6]                                                   # q1, q2
    x[7:7 + NW] = 1 + (u[7:7 + NW] - 0.5) * 0.004
    x[7 + NW:] = np.log(1e-6) + u[7 + NW:] * (np.log(1e-2) - np.log(1e-6))
    return x


def loglike(x):
    T0, P, k, b, rho, q1, q2 = x[:7]
    c = x[7:7 + NW]
    lj = x[7 + NW:]
    u1, u2 = 2 * np.sqrt(q1) * q2, np.sqrt(q1) * (1 - 2 * q2)
    ars = tt.ars_from_rho(rho, P)
    if b >= ars:
        return -1e300
    ll = 0.0
    for w, m, ci, lji in zip(W, MODELS, c, lj):
        mod = ci * m.lc(T0, P, k, b, ars, u1, u2)
        s2 = w["sig"] ** 2 + np.exp(2 * lji)
        ll += -0.5 * np.sum((w["f"] - mod) ** 2 / s2 + np.log(2 * np.pi * s2))
    return ll


t0 = time.time()
rs = np.random.default_rng(42)
sampler = dynesty.NestedSampler(loglike, ptform, ND, nlive=600, bound="multi", sample="rwalk", rstate=rs)
sampler.run_nested(dlogz=0.1, print_progress=False)
res = sampler.results
w = np.exp(res.logwt - res.logz[-1])
w /= w.sum()
samp = dynesty.utils.resample_equal(res.samples, w, rstate=rs)
T0, P, k, b, rho, q1, q2 = [samp[:, i] for i in range(7)]
u1, u2 = 2 * np.sqrt(q1) * q2, np.sqrt(q1) * (1 - 2 * q2)
ars = tt.ars_from_rho(rho, P)
T14 = tt.t14_hours(P, k, b, ars)
ST = json.load(open(os.path.join(AUD, "out", "s04_stellar.json")))
Rs = np.random.default_rng(1).normal(ST["radius_mann2015_MK"]["median"], ST["radius_mann2015_MK"]["std"], len(k))
pct = lambda v: [float(x) for x in np.percentile(v, [2.275, 15.865, 50, 84.135, 97.725])]
out = dict(signal=SIGNAL, set=SET, sampler="dynesty static, nlive 600, rwalk, dlogz 0.1", ld="Kipping q1,q2 ~ U(0,1)", detrend="biweight 1.0 d",
           logz=float(res.logz[-1]), logz_err=float(res.logzerr[-1]), niter=int(res.niter), runtime_s=time.time() - t0, n_samples=len(samp),
           T0=pct(T0), P=pct(P), k=pct(k), depth_k2_ppm=pct(k ** 2 * 1e6), b=pct(b), rho=pct(rho), T14_h=pct(T14),
           Rp_Rearth=pct(k * Rs * 109.076), u1=pct(u1), u2=pct(u2),
           datasets=[dict(key=x["key"], n=len(x["t"])) for x in W])
save_json(out, os.path.join(AUD, "out", "fits", f"dynesty_{SIGNAL}_{SET}.json"))
print(json.dumps({k_: out[k_] for k_ in ("signal", "set", "P", "k", "b", "rho", "T14_h", "Rp_Rearth", "u1", "u2", "logz", "runtime_s")}, indent=0))
