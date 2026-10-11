"""Step 3: S98-only null.  Null light curves: B removed (+-T14 around the three clean S98 events), then the flux is cyclically
shifted by an independent random offset within each continuous segment (gaps > 0.5 d); timestamps stay fixed.
Family 'shift' = as above; family 'invshift' = flux inverted about 1 (dips <-> bumps) then shifted; family 'revshift' =
time-reversed within each segment then shifted.  For each null light curve: BLS 1-20 d (same grid, durations and dy as
the observed search), record the global maximum SDE, the maximum depth S/N, and the period of the maximum.
FAP = (1 + #null >= observed) / (N + 1), where observed = SDE of the B peak in the observed spectrum."""
import numpy as np, json, sys, time, os
import common_b as c
import s98data as sd

var = sys.argv[1]
family = sys.argv[2]
N = int(sys.argv[3])
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
rng = np.random.default_rng(seed)
t, f = sd.get(var)
T = t.max() - t.min()
P = sd.period_grid(T)
r, pw, snr = sd.bls_spectrum(t, f, P)
S = sd.sde(pw)
near = np.abs(P - 11.133) < 0.03
iB = np.where(near)[0][np.argmax(pw[near])]
obs = dict(SDE_B=float(S[iB]), P_B=float(P[iB]), snr_B=float(snr[iB]), SDE_globalmax=float(S.max()), P_globalmax=float(P[np.argmax(pw)]),
           snr_max_any=float(np.nanmax(snr)))
dy = c.robust_std(f)
# remove B
keep = ~sd.inB(t, half=c.T14_B)
for tc in [4045.40, 4023.134]:
    keep &= np.abs(t - tc) > c.T14_B
tb, fb = t[keep], f[keep]
gi = np.where(np.diff(tb) > 0.5)[0]
edges = np.r_[0, gi + 1, len(tb)]
segs = [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)]
out = dict(var=var, family=family, obs=obs, n_points=int(len(t)), n_points_null=int(len(tb)), segments=[(float(tb[a]), float(tb[b - 1])) for a, b in segs],
           sde_max=[], snr_max=[], P_max=[])
fn = "out/b03_null_%s_%s_s%d.json" % (var, family, seed)
t0 = time.time()
for k in range(N):
    fn_ = fb.copy()
    if family == "invshift":
        fn_ = 2.0 - fn_
    for a, b in segs:
        x = fn_[a:b]
        if family == "revshift":
            x = x[::-1]
        fn_[a:b] = np.roll(x, rng.integers(0, b - a))
    rr, pwn, snrn = sd.bls_spectrum(tb, fn_, P, dy=dy)
    Sn = sd.sde(pwn)
    j = int(np.argmax(pwn))
    out["sde_max"].append(float(Sn[j])); out["snr_max"].append(float(np.nanmax(snrn))); out["P_max"].append(float(P[j]))
    if (k + 1) % 250 == 0 or k == N - 1:
        a = np.array(out["sde_max"]); s_ = np.array(out["snr_max"])
        out["n_done"] = k + 1
        out["n_ge_SDE"] = int(np.sum(a >= obs["SDE_B"])); out["fap_SDE"] = (out["n_ge_SDE"] + 1) / (k + 2)
        out["n_ge_snr"] = int(np.sum(s_ >= obs["snr_B"])); out["fap_snr"] = (out["n_ge_snr"] + 1) / (k + 2)
        out["null_q_SDE"] = [float(x) for x in np.percentile(a, [50, 90, 95, 99, 100])]
        json.dump(out, open(fn, "w"))
        print(k + 1, "SDE obs %.2f n>= %d fap %.4f | snr obs %.2f n>= %d fap %.4f | %.0fs" % (obs["SDE_B"], out["n_ge_SDE"], out["fap_SDE"],
              obs["snr_B"], out["n_ge_snr"], out["fap_snr"], time.time() - t0), flush=True)
