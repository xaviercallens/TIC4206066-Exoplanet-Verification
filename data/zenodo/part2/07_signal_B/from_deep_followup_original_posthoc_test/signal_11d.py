"""Second signal near 11.13 d (seen in S98 by the frozen search). Test it on all data with the 3.18-d transits masked. Exploratory addendum."""
import json, numpy as np, os
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from wotan import flatten
import batman
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A = json.load(open("out/lc_analysis.json")); P1 = A["fit_rho_prior"]["P"][1]; T1 = A["fit_rho_prior"]["t0"][1]; D1 = A["fit_rho_prior"]["T14_h"][1] / 24
src = {"S6 TESS-SPOC 30min": "lc/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "S32 SPOC 2min": "lc/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits",
       "S98 QLP 200s": "lc/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"}
cj = json.load(open("../../candidates/tic4206066_p11.132976/candidate.json"))
P11, T11, D11 = cj["signal"]["period_d"], cj["signal"]["epoch_btjd"], cj["signal"]["duration_hr"] / 24
def load(p):
    d = fits.open(p)[1].data; t = np.asarray(d["TIME"], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
    if "PDCSAP_FLUX" in d.columns.names: f = np.asarray(d["PDCSAP_FLUX"], float); g = q == 0
    else: f = np.asarray(d["SYS_RM_FLUX"], float); g = (q & 1610619071) == 0
    m = np.isfinite(t) & np.isfinite(f) & g; t, f = t[m], f[m] / np.nanmedian(f[m])
    msk = (np.abs((((t - T1) / P1 + 0.5) % 1 - 0.5) * P1) < 1.0 * D1)  # remove 3.18-d transits
    t, f = t[~msk], f[~msk]
    _, tr = flatten(t, f, method="biweight", window_length=0.75, return_trend=True, break_tolerance=0.3)
    ok = np.isfinite(tr); t, f = t[ok], f[ok] / tr[ok]
    s = 1.4826 * np.median(np.abs(f - np.median(f))); k = f < 1 + 5 * s
    return t[k], f[k], float(s)
L = {k: load(v) for k, v in src.items()}
t = np.concatenate([v[0] for v in L.values()]); f = np.concatenate([v[1] for v in L.values()]); e = np.concatenate([np.full(len(v[0]), v[2]) for v in L.values()])
# 1) independent search on S6+S32 only (no S98), long-period range
m_old = t < 3000
bls_old = BoxLeastSquares(t[m_old], f[m_old], e[m_old]); per = np.exp(np.linspace(np.log(1), np.log(25), 60000))
r_old = bls_old.power(per, [0.05, 0.067, 0.09], objective="snr")
top_old = []
for i in np.argsort(r_old.power)[::-1]:
    if all(abs(per[i] - p) / p > 0.01 for p, _ in top_old): top_old.append((float(per[i]), float(r_old.power[i])))
    if len(top_old) >= 6: break
snr_old_at_11 = float(r_old.power[np.argmin(np.abs(per - P11))])
# 2) cycle-count candidates linking S98 to S6/S32: for each integer N near (T11 - t) / P, compute box depth of S6+S32 at those phases
res = []
for dP in np.linspace(-0.004, 0.004, 1601):
    Pp = P11 + dP; ph = ((t[m_old] - T11) / Pp + 0.5) % 1 - 0.5; inn = np.abs(ph * Pp) < D11 / 2
    if inn.sum() < 10: continue
    w = 1 / e[m_old] ** 2; dd = (1 - np.sum(w[inn] * f[m_old][inn]) / np.sum(w[inn])) * 1e6; ee = 1e6 / np.sqrt(np.sum(w[inn])); res.append((Pp, dd, ee, dd / ee, int(inn.sum())))
res = np.array(res); j = int(np.argmax(res[:, 3]))
best = dict(P=float(res[j, 0]), depth=float(res[j, 1]), err=float(res[j, 2]), snr=float(res[j, 3]), n_in=int(res[j, 4]))
noise = dict(median_snr=float(np.median(res[:, 3])), std_snr=float(np.std(res[:, 3])), p99=float(np.percentile(res[:, 3], 99)))
# 3) predicted events with best P, individual depths
Pb, Tb = best["P"], T11
ev = []
for name, (tt, ff, s) in L.items():
    for n in range(int(np.floor((tt.min() - Tb) / Pb)) - 1, int(np.ceil((tt.max() - Tb) / Pb)) + 2):
        tc = Tb + n * Pb; w_ = np.abs(tt - tc) < 0.35
        if w_.sum() < 5: continue
        inn = np.abs(tt[w_] - tc) < D11 / 2; out_ = ~inn
        cad = np.median(np.diff(tt)); cov = inn.sum() * cad / D11
        if cov < 0.5: ev.append(dict(set=name, n=n, tc=float(tc), coverage=float(cov))); continue
        d = (np.mean(ff[w_][out_]) - np.mean(ff[w_][inn])) * 1e6; de = np.std(ff[w_][out_]) / np.sqrt(inn.sum()) * 1e6
        ev.append(dict(set=name, n=n, tc=float(tc), coverage=float(cov), depth=float(d), err=float(de), snr=float(d / de)))
# plot
fig = plt.figure(figsize=(18, 10)); gs = fig.add_gridspec(2, 4)
a0 = fig.add_subplot(gs[0, :2]); a0.plot(per, r_old.power, lw=0.5); a0.axvline(P11, color="C3", ls="--", label="11.13 d (from S98)"); a0.set_xscale("log"); a0.legend(); a0.set_title("S6 + S32 only (no S98): BLS after masking 3.18-d transits"); a0.set_xlabel("period (d)")
a1 = fig.add_subplot(gs[0, 2:]); a1.plot(res[:, 0], res[:, 3], lw=0.7); a1.axhline(noise["p99"], color="k", ls=":", label="99th percentile"); a1.set_xlabel("period (d)"); a1.set_ylabel("S6+S32 box S/N at S98 phase"); a1.legend(); a1.set_title("which cycle count links S98 to S6/S32?")
evv = [x for x in ev if "depth" in x]
for i, x in enumerate(evv[:8]):
    a = fig.add_subplot(gs[1, i % 4]) if i < 4 else None
    if a is None: break
    tt, ff, s = L[x["set"]]; w_ = np.abs(tt - x["tc"]) < 0.3
    a.plot((tt[w_] - x["tc"]) * 24, (ff[w_] - 1) * 1e6, ".", ms=2); a.axvspan(-D11 * 12, D11 * 12, color="C3", alpha=0.15); a.set_ylim(-3000, 2000)
    a.set_title(f"{x['set']} n={x['n']}: {x['depth']:.0f}±{x['err']:.0f} ppm", fontsize=9)
plt.tight_layout(); plt.savefig("plots/signal_11d.png", dpi=100); plt.close()
out = dict(P_s98=P11, T0_s98=T11, dur_h=D11 * 24, bls_S6_S32_top=top_old, bls_S6_S32_snr_at_11d=snr_old_at_11, best_link=best, link_noise=noise, events=ev)
json.dump(out, open("out/signal_11d.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("bls_S6_S32_top", "bls_S6_S32_snr_at_11d", "best_link", "link_noise")}, indent=1))
for x in ev: print(x)
