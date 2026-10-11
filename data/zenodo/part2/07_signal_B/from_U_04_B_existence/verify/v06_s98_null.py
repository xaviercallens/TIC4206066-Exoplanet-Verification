"""V06: S98-only BLS significance of B in QLP S98 (q==0), own detrending (biweight 0.75 d, A masked, A removed).
Flux column: DET_FLUX or SYS_RM_FLUX (argv[1]).  Own log-frequency grid 1-20 d, durations 1-3 h, astropy BLS (likelihood).
Null: B removed (+-1 T14 around the manuscript ephemeris), flux cyclically shifted within each segment (gaps > 0.5 d).
Records SDE and depth S/N of the global maximum per null; FAP = (1 + n_ge)/(N + 1) for both statistics."""
import sys, json, time
import numpy as np
from astropy.timeseries import BoxLeastSquares
import vcommon as v

col = sys.argv[1]; N = int(sys.argv[2]); seed = int(sys.argv[3])
rng = np.random.default_rng(seed)
d, t, q = v.load("qlp98")
f = np.asarray(d[col], float)
ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
t, f = t[ok], f[ok]
# per-segment normalisation (DET_FLUX has per-segment offsets), then biweight with A masked
gi = np.where(np.diff(t) > 0.5)[0]; ed = np.r_[0, gi + 1, len(t)]
for a, b in zip(ed[:-1], ed[1:]):
    f[a:b] /= np.median(f[a:b])
A = v.inA(t)
fd = v.biweight(t, f, A, 0.75)
keep = np.isfinite(fd) & ~A
t, fd = t[keep], fd[keep]
s = v.rstd(fd); keep = fd < 1 + 5 * s; t, fd = t[keep], fd[keep]
T = t.max() - t.min()
durs = np.array([1.0, 1.3, 1.6, 2.0, 2.5, 3.0]) / 24
dlnf = durs.min() / (T * 2.0)
freqs = np.exp(np.arange(np.log(1 / 20.0), np.log(1.0), dlnf))
P = np.sort(1 / freqs)
dy = np.full_like(fd, v.rstd(fd))


def spec(tt, ff):
    r = BoxLeastSquares(tt, ff, dy=dy[:len(tt)]).power(P, durs, objective="likelihood")
    pw = np.asarray(r.power)
    sde = (pw - pw.mean()) / pw.std()
    snr = np.asarray(r.depth) / np.asarray(r.depth_err)
    return sde, snr, r


t0 = time.time()
sde, snr, r = spec(t, fd)
nearB = np.abs(P / 11.1327 - 1) < 0.003
iB = np.flatnonzero(nearB)[np.argmax(sde[nearB])]
order = np.argsort(-sde)
peaks = []
for i in order:
    if all(abs(P[i] / p - 1) > 0.01 for p, _ in peaks):
        peaks.append((float(P[i]), float(sde[i])))
    if len(peaks) >= 6:
        break
obs = dict(col=col, nP=len(P), SDE_B=float(sde[iB]), P_B=float(P[iB]), snrB=float(snr[iB]), depthB_ppm=float(r.depth[iB] * 1e6),
           SDE_max=float(sde.max()), P_max=float(P[np.argmax(sde)]), top_peaks=peaks, rank_B=int(1 + sum(1 for p, sv in peaks if sv > sde[iB])))
print(json.dumps(obs), "%.1fs" % (time.time() - t0), flush=True)
# null
mB = np.zeros_like(t, bool)
for n in range(-1, 5):
    mB |= np.abs(t - (v.T0_MS + n * v.P_MS)) < v.DUR
tn, fn = t[~mB], fd[~mB]
gi = np.where(np.diff(tn) > 0.5)[0]; ed = np.r_[0, gi + 1, len(tn)]
nulls_sde, nulls_snr = [], []
for k in range(N):
    fs = fn.copy()
    for a, b in zip(ed[:-1], ed[1:]):
        fs[a:b] = np.roll(fn[a:b], rng.integers(1, b - a))
    sd_, sn_, _ = spec(tn, fs)
    j = int(np.argmax(sd_))
    nulls_sde.append(float(sd_[j])); nulls_snr.append(float(np.nanmax(sn_)))
    if (k + 1) % 25 == 0 or k == N - 1:
        a1 = np.array(nulls_sde); a2 = np.array(nulls_snr)
        res = dict(obs=obs, n=k + 1, n_ge_sde=int(np.sum(a1 >= obs["SDE_B"])), n_ge_snr=int(np.sum(a2 >= obs["snrB"])),
                   fap_sde=(1 + int(np.sum(a1 >= obs["SDE_B"]))) / (k + 2), fap_snr=(1 + int(np.sum(a2 >= obs["snrB"]))) / (k + 2),
                   null_sde_q=[float(x) for x in np.percentile(a1, [50, 90, 99, 100])], null_snr_q=[float(x) for x in np.percentile(a2, [50, 90, 99, 100])])
        json.dump(dict(res, sde=nulls_sde, snr=nulls_snr), open(v.HERE + "/out/v06_null_%s.json" % col, "w"))
        print(k + 1, {kk: res[kk] for kk in ("n_ge_sde", "fap_sde", "n_ge_snr", "fap_snr", "null_sde_q", "null_snr_q")}, "%.0fs" % (time.time() - t0), flush=True)
