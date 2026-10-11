"""Null test for the 11.13-d link: repeat the S98->S6/S32 cycle-count scan with the S98 epoch shifted by random phase offsets.
FAP = fraction of null scans whose max S/N >= observed. Exploratory addendum."""
import json, numpy as np, importlib.util
spec = importlib.util.spec_from_file_location("s11", "signal_11d.py")
# reuse data loading without re-running plots: copy the minimal pieces
from astropy.io import fits
from wotan import flatten
A = json.load(open("out/lc_analysis.json")); P1 = A["fit_rho_prior"]["P"][1]; T1 = A["fit_rho_prior"]["t0"][1]; D1 = A["fit_rho_prior"]["T14_h"][1] / 24
S = json.load(open("out/signal_11d.json")); P11, T11, D11 = S["P_s98"], S["T0_s98"], S["dur_h"] / 24
src = ["lc/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "lc/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits"]
T, F, E = [], [], []
for p in src:
    d = fits.open(p)[1].data; t = np.asarray(d["TIME"], float); f = np.asarray(d["PDCSAP_FLUX"], float); q = np.asarray(d["QUALITY"])
    m = np.isfinite(t) & np.isfinite(f) & (q == 0); t, f = t[m], f[m] / np.nanmedian(f[m])
    msk = np.abs((((t - T1) / P1 + 0.5) % 1 - 0.5) * P1) < D1; t, f = t[~msk], f[~msk]
    _, tr = flatten(t, f, method="biweight", window_length=0.75, return_trend=True, break_tolerance=0.3); ok = np.isfinite(tr); t, f = t[ok], f[ok] / tr[ok]
    s = 1.4826 * np.median(np.abs(f - np.median(f))); k = f < 1 + 5 * s; T.append(t[k]); F.append(f[k]); E.append(np.full(k.sum(), s))
t = np.concatenate(T); f = np.concatenate(F); w = 1 / np.concatenate(E) ** 2
dPs = np.linspace(-0.004, 0.004, 1601)
def scan(T0):
    best = -99
    for dP in dPs:
        Pp = P11 + dP; inn = np.abs(((t - T0) / Pp + 0.5) % 1 - 0.5) * Pp < D11 / 2
        if inn.sum() < 10: continue
        snr = (1 - np.sum(w[inn] * f[inn]) / np.sum(w[inn])) * np.sqrt(np.sum(w[inn]))
        best = max(best, snr)
    return best
obs = scan(T11)
rng = np.random.default_rng(7); null = []
for u in rng.uniform(0.1, 0.9, 300):
    null.append(scan(T11 + u * P11))
null = np.array(null); fap = float(np.mean(null >= obs))
out = dict(observed_max_snr=float(obs), null_n=len(null), null_median=float(np.median(null)), null_p95=float(np.percentile(null, 95)), null_p99=float(np.percentile(null, 99)), null_max=float(null.max()), FAP=fap,
           note="null: S98 epoch shifted by U(0.1,0.9) x P; same 1601-period scan over S6+S32 residual data (3.18-d transits masked)")
json.dump(out, open("out/signal_11d_null.json", "w"), indent=1); print(out)
